////////////////////////////////////////////////////////////////////////
// File:  ChannelOverlap.h
//
// Detection of channels present in more than one collection of a
// segmented data product (one std::vector<raw::RawDigit> or
// std::vector<sim::SimChannel> per TPCSet, TPC, ...).
//
// A channel found in two collections is what breaks the consumers: a TP
// maker emits its TPs twice, a backtracker index keeps one copy and drops
// the other, energy sums count it twice. Detection is therefore keyed on
// the channel and needs no geometry. Geometry comes in afterwards, only
// for the channels that do overlap, to say which TPCSets and TPCs are
// involved.
//
// Typical use in a module:
//
//   auto overlaps = channel_overlap::find(channel_overlap::from_handles(handles));
//   if (!overlaps.empty()) {
//     channel_overlap::annotate(overlaps, geom);
//     throw cet::exception("MyModule") << "tag matches collections sharing channels:\n"
//                                      << channel_overlap::describe(overlaps);
//   }
////////////////////////////////////////////////////////////////////////

#ifndef DUNETRIGGER_UTILITIES_CHANNELOVERLAP_H
#define DUNETRIGGER_UTILITIES_CHANNELOVERLAP_H

#include "art/Framework/Principal/Handle.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include <cstddef>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace geo {
  class WireReadoutGeom;
}

namespace dunetrigger::channel_overlap {

  /// The channels of one collection, with the collection's input tag.
  struct ChannelSet {
    std::string tag;
    std::vector<raw::ChannelID_t> channels;
  };

  /// Inclusive [first, last] range of contiguous channel IDs.
  using ChannelRange = std::pair<raw::ChannelID_t, raw::ChannelID_t>;

  /// Channels shared by two collections.
  struct Overlap {
    std::string first;  ///< tag of the collection where the channels were first seen
    std::string second; ///< tag of the collection where they were seen again
    std::vector<ChannelRange> ranges; ///< the shared channels, sorted and merged
    std::size_t n_channels = 0;

    // Filled by annotate()
    std::set<unsigned> tpcsets;
    std::set<std::pair<unsigned, unsigned>> tpcset_tpcs; ///< (TPCSet, TPC)
  };

  /// Sorts and deduplicates `channels` and merges them into contiguous ranges.
  std::vector<ChannelRange> to_ranges(std::vector<raw::ChannelID_t> channels);

  /// Channels present in more than one of `sets`, one entry per pair of
  /// collections, in the order the collections appear in `sets`. A channel
  /// repeated inside a single collection is not an overlap.
  std::vector<Overlap> find(std::vector<ChannelSet> const& sets);

  /// Adds the TPCSets and (TPCSet, TPC) pairs of the overlapping channels.
  /// Channels unknown to the geometry are skipped.
  void annotate(std::vector<Overlap>& overlaps, geo::WireReadoutGeom const& geom);

  /// One line per overlap: the two tags, the channel count and up to
  /// `max_ranges` channel ranges, then TPCSets and TPCs if annotated.
  std::string describe(std::vector<Overlap> const& overlaps, std::size_t max_ranges = 10);

  /// ChannelSets from art handles. `T` needs a `Channel()` member, as
  /// raw::RawDigit and sim::SimChannel have.
  template <typename T>
  std::vector<ChannelSet> from_handles(std::vector<art::Handle<std::vector<T>>> const& handles)
  {
    std::vector<ChannelSet> sets;
    sets.reserve(handles.size());
    for (auto const& handle : handles) {
      ChannelSet set;
      set.tag = handle.provenance()->inputTag().encode();
      set.channels.reserve(handle->size());
      for (T const& object : *handle)
        set.channels.push_back(object.Channel());
      sets.push_back(std::move(set));
    }
    return sets;
  }

} // namespace dunetrigger::channel_overlap

#endif // DUNETRIGGER_UTILITIES_CHANNELOVERLAP_H
