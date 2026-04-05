// RawDigitsReadTest_module.cc
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Selector.h"

#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/WireReadout.h"

#include <regex>    
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace duneana {

nlohmann::json parse_to_json(const std::string &str) {
  auto j = nlohmann::json::parse(str, nullptr, false);
  if (j.is_discarded()) {
    return str;
  }
  return j;
}

nlohmann::json fhicl_to_json(const fhicl::ParameterSet &pset) {
  nlohmann::json j;
  for (const std::string &key : pset.get_names
    ()) {
//   for (const std::string &key : pset.get_all_keys()) {
    // std::cout << " >>> " << key << std::endl;
    if (pset.is_key_to_atom(key)) {
      std::string field_str = pset.get<std::string>(key);
      j[key] = parse_to_json(field_str);
    } else if (pset.is_key_to_sequence(key)) {
      std::vector<std::string> field_vec =
          pset.get<std::vector<std::string>>(key);
      std::vector<json> json_vec;
      for (const std::string &val : field_vec) {
        json_vec.push_back(parse_to_json(val));
      }
      j[key] = json_vec;
    } else if (pset.is_key_to_table(key)) {
      fhicl::ParameterSet sub_pset = pset.get<fhicl::ParameterSet>(key);
      j[key] = fhicl_to_json(sub_pset);
    }
  }
  return j;
}

class RawDigitsReadTest : public art::EDAnalyzer {
public:
    explicit RawDigitsReadTest(fhicl::ParameterSet const& pset);

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;

private:
};

// Constructor — read FHiCL parameters
RawDigitsReadTest::RawDigitsReadTest(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    // , fHitLabel(pset.get<art::InputTag>("HitLabel"))
{}

void RawDigitsReadTest::beginJob() {
    // art::ServiceHandle<art::TFileService> tfs;
    // fHitHist = tfs->make<TH1F>("hNHits", "Number of hits;N hits;Events", 100, 0, 500);
}

void RawDigitsReadTest::analyze(art::Event const& e) {

    geo::WireReadoutGeom const *geom =
        &art::ServiceHandle<geo::WireReadout>()->Get();

    std::regex rawdigi_regex("daq.*");

    // art::SelectorByFunction s([](art::BranchDescription const& p){ return true;}, "pippo");
    art::SelectorByFunction s(
        [rawdigi_regex](art::BranchDescription const& p){
            return std::regex_match(p.inputTag().instance(), rawdigi_regex);
        },
        "pippo"
    );

    std::cout << "---RawDigits Input tags-------------------------------------------------" << std::endl;

    auto input_tags = e.getInputTags<std::vector<raw::RawDigit>>(s);
    for ( auto i : input_tags) {
        std::cout << i.label() << "   " << i.instance() << "   " << i.process() << std::endl;
    }

    std::cout << "---RawDigits Handles----------------------------------------------------" << std::endl;

    auto vec_h = e.getMany<std::vector<raw::RawDigit>>(s);

    for( auto h : vec_h) {
        auto i = h.provenance()->inputTag();
        

        auto ropid = geom->ChannelToROP(h->front().Channel());

        std::cout << i.label() << "   " << i.instance() << "   " << i.process() << " : size=" << h->size() << ", tpcset=" << ropid.asTPCsetID() << std::endl;
        // std::cout << "Provenance: " << *(h.provenance()) << std::endl;
        // std::cout << "PSet: " << fhicl_to_json(h.provenance()->parameterSet()).dump(4) << std::endl;
        std::cout << "wcls_main.structs.process_apa_index : " <<h.provenance()->parameterSet().get<int>("wcls_main.structs.process_apa_index") << std::endl;
    }

    std::cout << "------------------------------------------------------------------------" << std::endl;

}

void RawDigitsReadTest::endJob() {
    // Optional cleanup / summary
}

} // namespace duneana

DEFINE_ART_MODULE(duneana::RawDigitsReadTest)