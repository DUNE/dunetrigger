////////////////////////////////////////////////////////////////////////
// Class:       TriggerROPActivityMakerTPC
// Plugin Type: producer
// File:        TriggerROPActivityMakerTPC_module.cc
//
// Applies a configurable TriggerActivityROPMaker algorithm to each
// Readout Plane (ROP) of every detector element (APA/CRP).
// A fresh algorithm instance is built for each ROP in every event,
// so that each plane has independent internal state.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityROPData.hpp"
#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerPrimitive.hpp"

#include "dunetrigger/experimental/TPCTriggerSim/Verbosity.hh"
#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityROP.hpp"
#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityROPFactory.hpp"

#include <nlohmann/json.hpp>
#include <numeric>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace dunetrigger::experimental {

  using tpctriggeralgs::TriggerActivityROPData;
  using tpctriggeralgs::TriggerPrimitive;

  class TriggerROPActivityMakerTPC : public art::EDProducer {
  public:
    explicit TriggerROPActivityMakerTPC(fhicl::ParameterSet const &p);

    TriggerROPActivityMakerTPC(TriggerROPActivityMakerTPC const &) = delete;
    TriggerROPActivityMakerTPC(TriggerROPActivityMakerTPC &&) = delete;
    TriggerROPActivityMakerTPC &operator=(TriggerROPActivityMakerTPC const &) = delete;
    TriggerROPActivityMakerTPC &operator=(TriggerROPActivityMakerTPC &&) = delete;

    void beginJob() override;
    void produce(art::Event &e) override;

  private:
    const std::string algname_;
    const art::InputTag tp_tag_;
    const std::vector<raw::ChannelID_t> channel_mask_;
    const int verbosity_;
    const bool flush_;

    // Geometry service handle, initialised once in the constructor.
    geo::WireReadoutGeom const *geom_;

    // Per-ROP algorithm config index.
    // Built in the constructor using the geometry service and FHiCL input.
    const std::map<unsigned int, nlohmann::json> algconfigs_;

    // Algorithm factory instance
    std::shared_ptr<tpctriggeralgs::AbstractFactory<tpctriggeralgs::TriggerActivityROPMaker>> factory_ = tpctriggeralgs::TriggerActivityROPFactory::get_instance();

    // Try to read `key` from `pset` as type T; on success store it in `out` and return true.
    template <typename T>
    static bool tryGet(const fhicl::ParameterSet &pset, const std::string &key, nlohmann::json &out) {
      try {
        out = pset.get<T>(key);
        return true;
      } catch (const fhicl::exception &) {
        return false;
      }
    }

    // Convert a FHiCL atom or sequence of atoms to JSON, trying the narrowest type first
    // (signed int, unsigned int for values above INT64_MAX, double, bool, then string).
    static nlohmann::json valueToJson(const fhicl::ParameterSet &pset, const std::string &key) {
      nlohmann::json v;
      if (pset.is_key_to_sequence(key)) {
        if (tryGet<std::vector<int64_t>>(pset, key, v) || tryGet<std::vector<uint64_t>>(pset, key, v) ||
            tryGet<std::vector<double>>(pset, key, v) || tryGet<std::vector<bool>>(pset, key, v) ||
            tryGet<std::vector<std::string>>(pset, key, v))
          return v;
        throw cet::exception("TriggerROPActivityMakerTPC")
          << "FHiCL sequence '" << key << "' must be a homogeneous sequence of atoms.\n";
      }
      if (tryGet<int64_t>(pset, key, v) || tryGet<uint64_t>(pset, key, v) || tryGet<double>(pset, key, v) ||
          tryGet<bool>(pset, key, v) || tryGet<std::string>(pset, key, v))
        return v;
      throw cet::exception("TriggerROPActivityMakerTPC") << "Cannot convert FHiCL parameter '" << key << "' to JSON.\n";
    }

    // Convert a FHiCL parameter set to a JSON object understood by tpctriggeralgs.
    // Nested tables become nested JSON objects. Unknown keys are rejected by the algorithms themselves.
    static nlohmann::json toJson(const fhicl::ParameterSet &pset) {
      nlohmann::json cfg = nlohmann::json::object();
      for (const auto &key : pset.get_names()) {
        if (pset.is_key_to_table(key))
          cfg[key] = toJson(pset.get<fhicl::ParameterSet>(key));
        else
          cfg[key] = valueToJson(pset, key);
      }
      return cfg;
    }

    // Sort and deduplicate the channel mask so it can be searched with std::binary_search.
    static std::vector<raw::ChannelID_t> sortedUnique(std::vector<raw::ChannelID_t> v) {
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      return v;
    }

    // Build the per-ROP config map from FHiCL and the geometry service. Called once during construction. 
    // Only requested ROP IDs specified within "active_rops" list in config fcl are used in TAMaking. 
    static std::map<unsigned int, nlohmann::json> buildAlgConfigs(fhicl::ParameterSet const &p, geo::WireReadoutGeom const &geom) {
      std::vector<unsigned int> active_rops;

      // If "active_rops" is missing from the FHiCL config, default to all ROPs present in the geo
      if (!p.get_if_present("active_rops", active_rops)) {
        active_rops.resize(geom.MaxROPs());
        std::iota(active_rops.begin(), active_rops.end(), 0);
      }

      std::map<unsigned int, nlohmann::json> configs;
      for (unsigned int r : active_rops) {
        if (r >= geom.MaxROPs()) {
          throw cet::exception("TriggerROPActivityMakerTPC")
            << "active_rops includes ROP " << r << " but the geometry has only " << geom.MaxROPs() << " ROPs per TPC set.\n";
        }
        const std::string key = "algconfig_rop" + std::to_string(r);
        if (!p.has_key(key)) {
          throw cet::exception("TriggerROPActivityMakerTPC") 
            << "active_rops includes ROP " << r << " but no matching '" << key << "' was provided.\n";
        }
        configs[r] = toJson(p.get<fhicl::ParameterSet>(key));
      }
      return configs;
    }

    // Return true if two TriggerPrimitives represent the same hit.
    static bool sameTP(const TriggerPrimitive &a, const TriggerPrimitive &b) {
      return a.channel == b.channel && a.time_start == b.time_start && a.adc_integral == b.adc_integral;
    }

    // Comparator: sort art::Ptr<TriggerPrimitive>s by (time_start, channel).
    // the old code sorted by channel - not sure if this is necessary
    static bool earlierTP(const art::Ptr<TriggerPrimitive> &a, const art::Ptr<TriggerPrimitive> &b) {
      return std::tie(a->time_start, a->channel) < std::tie(b->time_start, b->channel);
    }

    // Build art::Ptr associations between one output TA and its input TPs.
    void buildAssociations(
                           const tpctriggeralgs::TriggerActivityROP &ta,
                           const std::vector<art::Ptr<TriggerPrimitive>> &plane_tps,
                           art::Ptr<TriggerActivityROPData> taPtr,
                           art::Assns<TriggerActivityROPData, TriggerPrimitive> &assns) const;
  };

  // Constructor

  TriggerROPActivityMakerTPC::TriggerROPActivityMakerTPC(fhicl::ParameterSet const &p)
    : EDProducer{p}
    , algname_(p.get<std::string>("algorithm"))
    , tp_tag_(p.get<art::InputTag>("tp_tag"))
    , channel_mask_(sortedUnique(p.get<std::vector<raw::ChannelID_t>>("channel_mask", {})))
    , verbosity_(p.get<int>("verbosity", 1))
    , flush_(p.get<bool>("flush", true))
    , geom_(&art::ServiceHandle<geo::WireReadout>()->Get())
    , algconfigs_(buildAlgConfigs(p, *geom_)) {

      // Validate the algorithm name and every per-ROP config at construction time,
      // so configuration errors surface at job start rather than on the first event.
      for (const auto &[rop, cfg] : algconfigs_) {
        auto alg = factory_->build_maker(algname_);
        if (!alg)
          throw cet::exception("TriggerROPActivityMakerTPC") << "Unknown algorithm: '" << algname_ << "'\n";
        try {
          alg->configure(cfg);
        } catch (const std::exception &ex) {
          throw cet::exception("TriggerROPActivityMakerTPC")
            << "Invalid configuration 'algconfig_rop" << rop << "' for " << algname_ << ": " << ex.what() << "\n";
        }
      }

      produces<std::vector<TriggerActivityROPData>>();
      produces<art::Assns<TriggerActivityROPData, TriggerPrimitive>>();
      consumes<std::vector<TriggerPrimitive>>(tp_tag_);
    }

  //  beginJob
  void TriggerROPActivityMakerTPC::beginJob() {

    std::cout << "Producing TAs for ROP planes: ";
    for (const auto &[rop, _] : algconfigs_) std::cout << rop << " ";
    std::cout << "(" << algconfigs_.size() << "/" << geom_->MaxROPs() << ")\n";

    if (verbosity_ >= Verbosity::kInfo) {
      std::cout << "Masked channels:";
      for (auto ch : channel_mask_) std::cout << " " << ch;
      std::cout << "\n";
    }
  }

  // producer
  void TriggerROPActivityMakerTPC::produce(art::Event &e) {

    auto ta_vec = std::make_unique<std::vector<TriggerActivityROPData>>();
    auto assns  = std::make_unique<art::Assns<TriggerActivityROPData, TriggerPrimitive>>();

    art::PtrMaker<TriggerActivityROPData> taPtrMaker{e};
    auto tpHandle = e.getValidHandle<std::vector<TriggerPrimitive>>(tp_tag_);

    // Split TPs by ROPID, such that each TAMaker instance accepts TP data from single APA plane.
    std::map<readout::ROPID, std::vector<art::Ptr<TriggerPrimitive>>> tps_by_plane;
    for (size_t i = 0; i < tpHandle->size(); ++i) {
      tps_by_plane[geom_->ChannelToROP((*tpHandle)[i].channel)].push_back(art::Ptr<TriggerPrimitive>(tpHandle, i));
    }

    // Run the TAMaker alg. on each plane
    for (auto &[plane, plane_tps] : tps_by_plane) {
      const unsigned int rop_idx = plane.ROP;

      // Skip ROPs not listed in active_rops
      if (algconfigs_.find(rop_idx) == algconfigs_.end()) continue;

      std::sort(plane_tps.begin(), plane_tps.end(), earlierTP);

      auto alg = factory_->build_maker(algname_);
      if (!alg)
        throw cet::exception("TriggerROPActivityMakerTPC") << "Unknown algorithm: '" << algname_ << "'\n";

      alg->configure(algconfigs_.at(rop_idx));

      if (verbosity_ >= Verbosity::kInfo) std::cout << "Running " << algname_ << " on plane " << plane << "\n";

      // Feed TPs into the algorithm, skipping any masked channels.
      std::vector<tpctriggeralgs::TriggerActivityROP> created_tas;
      for (const auto &itp : plane_tps) {
        if (std::binary_search(channel_mask_.begin(), channel_mask_.end(), itp->channel)) {
          if (verbosity_ >= Verbosity::kDebug)
            std::cout << "Skipping masked channel " << itp->channel << "\n";
          continue;
        }
        (*alg)(*itp, created_tas);
      }

      // Optionally flush the algorithm window to ensure final in-progress window is evaluated before the event closes.
      if (flush_) alg->flush(std::numeric_limits<uint64_t>::max(), created_tas);

      if (verbosity_ >= Verbosity::kInfo && !created_tas.empty())
        std::cout << "Created " << created_tas.size() << " TAs on plane " << plane << "\n";

      // Store output TAs and build TP associations
      for (const auto &ta : created_tas) {
        auto taPtr = taPtrMaker(ta_vec->size());
        ta_vec->emplace_back(ta);
        buildAssociations(ta, plane_tps, taPtr, *assns);
      }
    }

    e.put(std::move(ta_vec));
    e.put(std::move(assns));
  }

  // buildAssociations 

  void TriggerROPActivityMakerTPC::buildAssociations(
                                                  const tpctriggeralgs::TriggerActivityROP &ta,
                                                  const std::vector<art::Ptr<TriggerPrimitive>> &plane_tps,
                                                  art::Ptr<TriggerActivityROPData> taPtr,
                                                  art::Assns<TriggerActivityROPData, TriggerPrimitive> &assns) const {

    art::PtrVector<TriggerPrimitive> ptrs;
    for (const auto &in_tp : ta.inputs) {
      auto it = std::find_if(plane_tps.begin(), plane_tps.end(), [&](const art::Ptr<TriggerPrimitive> &t) { return sameTP(*t, in_tp); });

      if (it == plane_tps.end()) {
        mf::LogWarning("TriggerROPActivityMakerTPC") << "TP recorded in TA not found in input list -- skipping association.";
        continue;
      }

      size_t matches = 0;
      while (it != plane_tps.end()) {
        ptrs.push_back(*it);
        ++matches;
        it = std::find_if(std::next(it), plane_tps.end(), [&](const art::Ptr<TriggerPrimitive> &t) { return sameTP(*t, in_tp); });
      }
      if (matches > 1)
        mf::LogWarning("TriggerROPActivityMakerTPC") << matches << " matching TPs found -- possible duplicate in input list.";
    }
    assns.addMany(taPtr, ptrs);
  }

} // namespace dunetrigger::experimental

DEFINE_ART_MODULE(dunetrigger::experimental::TriggerROPActivityMakerTPC)
