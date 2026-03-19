////////////////////////////////////////////////////////////////////////
// Class:       TriggerActivityMakerTPC
// Plugin Type: producer
// File:        TriggerActivityMakerTPC_module.cc
//
// Applies a configurable TriggerActivityMaker algorithm to each
// Readout Plane (ROP) of every detector element (APA/CRP).
// One algorithm instance per ROP is built in beginJob() 
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

#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include "dunetrigger/TriggerSim/Verbosity.hh"
#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivity.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivityFactory.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace dunetrigger {

  using dunedaq::trgdataformats::TriggerActivityData;
  using dunedaq::trgdataformats::TriggerPrimitive;

  // Bundles a TriggerPrimitive with its index in the original input collection,
  // This is needed later to build art::Ptr associations without re-searching.
  struct IndexedTP {
    size_t index;
    TriggerPrimitive tp;
  }; 

  class TriggerActivityMakerTPC : public art::EDProducer {
  public:
    explicit TriggerActivityMakerTPC(fhicl::ParameterSet const &p);

    TriggerActivityMakerTPC(TriggerActivityMakerTPC const &) = delete;
    TriggerActivityMakerTPC(TriggerActivityMakerTPC &&) = delete;
    TriggerActivityMakerTPC &operator=(TriggerActivityMakerTPC const &) = delete;
    TriggerActivityMakerTPC &operator=(TriggerActivityMakerTPC &&) = delete;

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
    std::shared_ptr<triggeralgs::AbstractFactory<triggeralgs::TriggerActivityMaker>> factory_ = triggeralgs::TriggerActivityFactory::get_instance();

    // Convert a FHiCL parameter set to a JSON object understood by triggeralgs.
    static nlohmann::json toJson(const fhicl::ParameterSet &pset) {
      nlohmann::json cfg;
      for (const auto &key : pset.get_all_keys()) {
        try {
          cfg[key] = pset.get<uint64_t>(key);
          continue;
        } catch (const fhicl::exception &) {}
        try {
          cfg[key] = pset.get<bool>(key);
          continue;
        } catch (const fhicl::exception &) {
          throw cet::exception("TriggerActivityMakerTPC") << "FHiCL parameter '" << key << "' is neither a uint64_t nor a bool.\n";
        }
      }
      return cfg;
    }

    // Build the per-ROP config map from FHiCL and the geometry service. Called once during construction. 
    // Only requested ROP IDs specified within "active_rops" list in config fcl are used in TAMaking. 
    static std::map<unsigned int, nlohmann::json> buildAlgConfigs(fhicl::ParameterSet const &p, geo::WireReadoutGeom const &geom) {
      const auto active_rops = p.get<std::vector<unsigned int>>("active_rops");

      std::map<unsigned int, nlohmann::json> configs;
      for (unsigned int r : active_rops) {
        const std::string key = "algconfig_rop" + std::to_string(r);
        if (!p.has_key(key))
          throw cet::exception("TriggerActivityMakerTPC") << "active_rops includes ROP " << r << " but no matching '" << key << "' was provided.\n";
        configs[r] = toJson(p.get<fhicl::ParameterSet>(key));
      }
      return configs;
    }

    // Return true if two TriggerPrimitives represent the same hit.
    static bool sameTP(const TriggerPrimitive &a, const TriggerPrimitive &b) {
      return a.channel == b.channel && a.time_start == b.time_start && a.adc_integral == b.adc_integral;
    }

    // Comparator: sort IndexedTPs by (time_start, channel).
    static bool earlierTP(const IndexedTP &a, const IndexedTP &b) {
      return std::tie(a.tp.time_start, a.tp.channel) < std::tie(b.tp.time_start, b.tp.channel);
    }

    // Build art::Ptr associations between one output TA and its input TPs.
    void buildAssociations(
                           const triggeralgs::TriggerActivity &ta,
                           const std::vector<IndexedTP> &scope_tps,
                           art::Ptr<TriggerActivityData> taPtr,
                           const art::ValidHandle<std::vector<TriggerPrimitive>> &tpHandle,
                           art::Assns<TriggerActivityData, TriggerPrimitive> &assns) const;
  };

  // Constructor

  TriggerActivityMakerTPC::TriggerActivityMakerTPC(fhicl::ParameterSet const &p)
    : EDProducer{p}
    , algname_(p.get<std::string>("algorithm"))
    , tp_tag_(p.get<art::InputTag>("tp_tag"))
    , channel_mask_(p.get<std::vector<raw::ChannelID_t>>("channel_mask", {}))
    , verbosity_(p.get<int>("verbosity", 1))
    , flush_(p.get<bool>("flush", false))
    , geom_(&art::ServiceHandle<geo::WireReadout>()->Get())
    , algconfigs_(buildAlgConfigs(p, *geom_)) {

      produces<std::vector<TriggerActivityData>>();
      produces<art::Assns<TriggerActivityData, TriggerPrimitive>>();
      consumes<std::vector<TriggerPrimitive>>(tp_tag_);
    }

  //  beginJob
  void TriggerActivityMakerTPC::beginJob() {

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
  void TriggerActivityMakerTPC::produce(art::Event &e) {

    auto ta_vec = std::make_unique<std::vector<TriggerActivityData>>();
    auto assns  = std::make_unique<art::Assns<TriggerActivityData, TriggerPrimitive>>();

    art::PtrMaker<TriggerActivityData> taPtrMaker{e};
    auto tpHandle = e.getValidHandle<std::vector<TriggerPrimitive>>(tp_tag_);

    // Split TPs by ROPID, such that each TAMaker instance accepts TP data from single APA plane. 
    // Each TP is stored alongside its original vector index so that art::Ptr associations can be built without re-searching the full TP collection.

    std::map<readout::ROPID, std::vector<IndexedTP>> tps_by_scope;
    for (size_t i = 0; i < tpHandle->size(); ++i) {
      const TriggerPrimitive &tp = (*tpHandle)[i];
      tps_by_scope[geom_->ChannelToROP(tp.channel)].push_back({i, tp});
    }

    // Run the TAMaker alg. on each scope 
    for (auto &[scope, indexed_tps] : tps_by_scope) {
      const unsigned int rop_idx = scope.ROP;

      // Skip ROPs not listed in active_rops
      if (algconfigs_.find(rop_idx) == algconfigs_.end()) continue;

      std::sort(indexed_tps.begin(), indexed_tps.end(), earlierTP);

      auto alg = factory_->build_maker(algname_);
      if (!alg)
        throw cet::exception("TriggerActivityMakerTPC") << "Unknown algorithm: '" << algname_ << "'\n";

      alg->configure(algconfigs_.at(rop_idx));

      if (verbosity_ >= Verbosity::kInfo) std::cout << "Running " << algname_ << " on scope " << scope << "\n";

      // Feed TPs into the algorithm, skipping any masked channels.
      std::vector<triggeralgs::TriggerActivity> created_tas;
      for (const auto &itp : indexed_tps) {
        if (std::binary_search(channel_mask_.begin(), channel_mask_.end(), itp.tp.channel)) {
          if (verbosity_ >= Verbosity::kDebug)
            std::cout << "Skipping masked channel " << itp.tp.channel << "\n";
          continue;
        }
        (*alg)(itp.tp, created_tas);
      }

      // Optionally flush the algorithm window with a dummy TP, ensuring the
      // final in-progress window is evaluated before the event closes.
      if (flush_) (*alg)(TriggerPrimitive{}, created_tas);

      if (verbosity_ >= Verbosity::kInfo && !created_tas.empty())
        std::cout << "Created " << created_tas.size() << " TAs on scope " << scope << "\n";

      // Store output TAs and build TP associations
      for (const auto &ta : created_tas) {
        auto taPtr = taPtrMaker(ta_vec->size());
        ta_vec->emplace_back(ta);
        buildAssociations(ta, indexed_tps, taPtr, tpHandle, *assns);
      }
    }

    e.put(std::move(ta_vec));
    e.put(std::move(assns));
  }

  // buildAssociations 

  void TriggerActivityMakerTPC::buildAssociations(
                                                  const triggeralgs::TriggerActivity &ta,
                                                  const std::vector<IndexedTP> &scope_tps,
                                                  art::Ptr<TriggerActivityData> taPtr,
                                                  const art::ValidHandle<std::vector<TriggerPrimitive>> &tpHandle,
                                                  art::Assns<TriggerActivityData, TriggerPrimitive> &assns) const {

    art::PtrVector<TriggerPrimitive> ptrs;
    for (const auto &in_tp : ta.inputs) {
      auto it = std::find_if(scope_tps.begin(), scope_tps.end(), [&](const IndexedTP &t) { return sameTP(t.tp, in_tp); });

      if (it == scope_tps.end()) {
        mf::LogWarning("TriggerActivityMakerTPC") << "TP recorded in TA not found in input list -- skipping association.";
        continue;
      }

      size_t matches = 0;
      while (it != scope_tps.end()) {
        ptrs.push_back(art::Ptr<TriggerPrimitive>(tpHandle, it->index));
        ++matches;
        it = std::find_if(std::next(it), scope_tps.end(), [&](const IndexedTP &t) { return sameTP(t.tp, in_tp); });
      }
      if (matches > 1)
        mf::LogWarning("TriggerActivityMakerTPC") << matches << " matching TPs found -- possible duplicate in input list.";
    }
    assns.addMany(taPtr, ptrs);
  }

} // namespace dunetrigger

DEFINE_ART_MODULE(dunetrigger::TriggerActivityMakerTPC)
