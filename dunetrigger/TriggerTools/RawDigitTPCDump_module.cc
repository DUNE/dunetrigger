////////////////////////////////////////////////////////////////////////
// Class:       RawDigitTPCDump
// Plugin Type: analyzer
// File:        RawDigitTPCDump_module.cc
//
// Diagnostic: for every std::vector<raw::RawDigit> collection matching
// rawdigit_tag (a regex InputTag, see dunetrigger::vendor::lar::util::getManyByRegexTag), print
// which TPCSets (and cryostat:TPC pairs) its channels map to, and report
// channels present in more than one collection.
//
// Run with:  lar -c rawdigit_tpc_dump_vd10kt.fcl -s <artfile.root> -n 1
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "dunetrigger/vendor/lardata/ArtDataHelper/GetManyByRegexTag.h"
#include "dunetrigger/Utilities/ChannelOverlap.h"
#include "fhiclcpp/types/Atom.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RawData/RawDigit.h"

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace dunetrigger {
  class RawDigitTPCDump;
}

class dunetrigger::RawDigitTPCDump : public art::EDAnalyzer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> rawdigit_tag{
        fhicl::Name("rawdigit_tag"),
        fhicl::Comment("label:instance:process regex patterns; empty process = newest matching process"),
        art::InputTag{"tpcrawdecoder:daq.*"}};
    fhicl::Atom<unsigned> max_events{fhicl::Name("max_events"),
                                     fhicl::Comment("Stop dumping after this many events (0 = all)"), 1u};
    fhicl::Atom<bool> verbose{fhicl::Name("verbose"),
                              fhicl::Comment("Also print per-TPCSet channel counts by view"), false};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit RawDigitTPCDump(Parameters const &p);

  void analyze(art::Event const &e) override;

private:
  art::InputTag rawdigit_tag_;
  unsigned max_events_;
  bool verbose_;
  unsigned n_dumped_ = 0;
};

dunetrigger::RawDigitTPCDump::RawDigitTPCDump(Parameters const &p)
    : EDAnalyzer{p}, rawdigit_tag_(p().rawdigit_tag()), max_events_(p().max_events()),
      verbose_(p().verbose()) {
  consumesMany<std::vector<raw::RawDigit>>();
}

void dunetrigger::RawDigitTPCDump::analyze(art::Event const &e) {
  if (max_events_ && n_dumped_ >= max_events_) return;
  ++n_dumped_;

  geo::WireReadoutGeom const &geom = art::ServiceHandle<geo::WireReadout>()->Get();

  // Only the newest matching process is kept, so RawDigits re-saved by
  // later processes are not counted twice.
  auto handles = dunetrigger::vendor::lar::util::getManyByRegexTag<std::vector<raw::RawDigit>>(e, rawdigit_tag_);

  if (handles.empty()) {
    std::cout << "RawDigitTPCDump: Run " << e.run() << " SubRun " << e.subRun() << " Event " << e.event()
              << ": no std::vector<raw::RawDigit> collections matching rawdigit_tag \""
              << rawdigit_tag_.encode() << "\"" << std::endl;
    return;
  }

  std::ostringstream out;
  out << "=== Run " << e.run() << " SubRun " << e.subRun() << " Event " << e.event() << ": "
      << handles.size() << " RawDigit collection(s) matching \"" << rawdigit_tag_.encode() << "\" ===\n";

  for (auto const &h : handles) {
    std::string const tag = h.provenance()->inputTag().encode();

    std::map<int, std::map<int, unsigned>> tpcset_view_nch;  // tpcset -> view -> nchannels
    std::set<std::pair<unsigned, unsigned>> cryo_tpcs;
    std::set<std::pair<unsigned, unsigned>> tpcset_tpcs;      // (tpcset, tpc) pairs covered
    std::set<ULong64_t> n_samples;                            // distinct readout lengths
    std::set<int> compressions;                               // distinct raw::Compress_t values
    unsigned n_invalid = 0, n_empty = 0, n_with_adcs = 0;
    raw::ChannelID_t ch_min = raw::InvalidChannelID, ch_max = 0;

    for (raw::RawDigit const &rd : *h) {
      raw::ChannelID_t const ch = rd.Channel();
      if (!geom.HasChannel(ch)) {
        ++n_invalid;
        continue;
      }
      if (rd.NADC() == 0)
        ++n_empty;
      else
        ++n_with_adcs;
      n_samples.insert(rd.Samples());
      compressions.insert(rd.Compression());
      if (ch < ch_min) ch_min = ch;
      if (ch > ch_max) ch_max = ch;

      readout::ROPID const rop = geom.ChannelToROP(ch);
      unsigned const tpcset = rop.asTPCsetID().TPCset;
      ++tpcset_view_nch[tpcset][geom.View(rop)];

      // a channel can read wires in more than one TPC of its TPCSet
      for (geo::WireID const &wid : geom.ChannelToWire(ch)) {
        cryo_tpcs.emplace(wid.Cryostat, wid.TPC);
        tpcset_tpcs.emplace(tpcset, wid.TPC);
      }
    }

    out << "\n[" << tag << "] " << h->size() << " channels";
    if (!h->empty()) out << " (ID range " << ch_min << "-" << ch_max << ")";
    out << ", " << n_with_adcs << " with ADCs";
    if (n_empty) out << ", " << n_empty << " with no ADCs";
    if (n_invalid) out << ", " << n_invalid << " NOT IN GEOMETRY";
    out << "\n  Samples:";
    for (ULong64_t n : n_samples) out << " " << n;
    out << "  Compression:";
    for (int c : compressions) out << " " << c;
    out << "\n  TPCSets:";
    for (auto const &[tpcset, views] : tpcset_view_nch) out << " " << tpcset;
    out << "\n  TPCs (cryo:tpc):";
    for (auto const &[cryo, tpc] : cryo_tpcs) out << " " << cryo << ":" << tpc;
    out << "\n  TPCSet:TPC:";
    for (auto const &[tpcset, tpc] : tpcset_tpcs) out << " " << tpcset << ":" << tpc;
    out << "\n";

    if (verbose_) {
      for (auto const &[tpcset, views] : tpcset_view_nch) {
        out << "    TPCSet " << tpcset << ":";
        for (auto const &[view, n] : views) out << " view" << view << "=" << n;
        out << "\n";
      }
    }
  }

  // Same check as TriggerPrimitiveMakerTPC and TriggerAnaTree, reported instead of thrown
  auto overlaps = channel_overlap::find(channel_overlap::from_handles(handles));
  if (overlaps.empty()) {
    out << "\nNo channel is present in more than one collection.\n";
  } else {
    channel_overlap::annotate(overlaps, geom);
    out << "\nOVERLAPS (channels present in more than one collection):\n"
        << channel_overlap::describe(overlaps) << "\n";
  }

  // std::cout rather than mf::LogInfo: DUNE message configs often suppress INFO
  std::cout << out.str() << std::flush;
}

DEFINE_ART_MODULE(dunetrigger::RawDigitTPCDump)
