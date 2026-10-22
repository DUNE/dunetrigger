#include "dunetrigger/Utilities/ChannelOverlap.h"

#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

#include <algorithm>
#include <map>
#include <sstream>
#include <unordered_map>

namespace dunetrigger::channel_overlap {

  std::vector<ChannelRange> to_ranges(std::vector<raw::ChannelID_t> channels)
  {
    std::sort(channels.begin(), channels.end());
    channels.erase(std::unique(channels.begin(), channels.end()), channels.end());

    std::vector<ChannelRange> ranges;
    for (raw::ChannelID_t ch : channels) {
      if (!ranges.empty() && ch == ranges.back().second + 1)
        ranges.back().second = ch;
      else
        ranges.emplace_back(ch, ch);
    }
    return ranges;
  }

  std::vector<Overlap> find(std::vector<ChannelSet> const& sets)
  {
    // Collection (index in sets) in which each channel was first seen
    std::unordered_map<raw::ChannelID_t, std::size_t> first_seen;
    // Channels seen again in a later collection, grouped by (first, later) collection
    std::map<std::pair<std::size_t, std::size_t>, std::vector<raw::ChannelID_t>> shared;

    for (std::size_t i = 0; i < sets.size(); ++i) {
      first_seen.reserve(first_seen.size() + sets[i].channels.size());
      for (raw::ChannelID_t ch : sets[i].channels) {
        auto const [it, inserted] = first_seen.try_emplace(ch, i);
        if (!inserted && it->second != i) shared[{it->second, i}].push_back(ch);
      }
    }

    std::vector<Overlap> overlaps;
    overlaps.reserve(shared.size());
    for (auto& [pair, channels] : shared) {
      Overlap o;
      o.first = sets[pair.first].tag;
      o.second = sets[pair.second].tag;
      o.ranges = to_ranges(std::move(channels));
      for (auto const& [lo, hi] : o.ranges)
        o.n_channels += hi - lo + 1;
      overlaps.push_back(std::move(o));
    }
    return overlaps;
  }

  void annotate(std::vector<Overlap>& overlaps, geo::WireReadoutGeom const& geom)
  {
    for (Overlap& o : overlaps) {
      for (auto const& [lo, hi] : o.ranges) {
        for (raw::ChannelID_t ch = lo; ch <= hi; ++ch) {
          if (!geom.HasChannel(ch)) continue;
          unsigned const tpcset = geom.ChannelToROP(ch).asTPCsetID().TPCset;
          o.tpcsets.insert(tpcset);
          // a channel can read wires in more than one TPC of its TPCSet
          for (geo::WireID const& wid : geom.ChannelToWire(ch))
            o.tpcset_tpcs.emplace(tpcset, wid.TPC);
        }
      }
    }
  }

  std::string describe(std::vector<Overlap> const& overlaps, std::size_t max_ranges)
  {
    std::ostringstream out;
    for (std::size_t n = 0; n < overlaps.size(); ++n) {
      Overlap const& o = overlaps[n];
      if (n) out << '\n';
      out << "  " << o.first << " and " << o.second << ": " << o.n_channels << " channels (";
      for (std::size_t i = 0; i < std::min(o.ranges.size(), max_ranges); ++i) {
        if (i) out << ", ";
        out << o.ranges[i].first;
        if (o.ranges[i].second != o.ranges[i].first) out << "-" << o.ranges[i].second;
      }
      if (o.ranges.size() > max_ranges) out << ", ... " << o.ranges.size() - max_ranges << " more ranges";
      out << ")";
      if (!o.tpcsets.empty()) {
        out << "; TPCSets:";
        for (unsigned s : o.tpcsets) out << " " << s;
        out << "; TPCs (tpcset:tpc):";
        for (auto const& [s, t] : o.tpcset_tpcs) out << " " << s << ":" << t;
      }
    }
    return out.str();
  }

} // namespace dunetrigger::channel_overlap
