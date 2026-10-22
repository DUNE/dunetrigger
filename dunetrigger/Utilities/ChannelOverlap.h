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
// the channel and needs no geometry. Geometry is used afterwards, only on
// the channels that do overlap, to say which TPCSets and TPCs are involved.
//
// Typical use in a module: one checker per module, built with the wire
// readout geometry, applied to the handles of every event.
//
//   ChannelOverlapCheck overlap_check_{&art::ServiceHandle<geo::WireReadout>()->Get()};
//   ...
//   if (auto overlaps = overlap_check_(handles); !overlaps.empty()) {
//     throw cet::exception("MyModule") << "tag matches collections sharing channels:\n"
//                                      << overlaps.describe();
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

namespace dunetrigger {

  /// Inclusive [first, last] range of contiguous channel IDs.
  using ChannelRange = std::pair<raw::ChannelID_t, raw::ChannelID_t>;

  /// Result of a ChannelOverlapCheck: the channels shared by pairs of
  /// collections, one entry per pair, in the order the collections were given.
  class ChannelOverlaps {
  public:
    struct Overlap {
      std::string first;  ///< tag of the collection where the channels were first seen
      std::string second; ///< tag of the collection where they were seen again
      std::vector<ChannelRange> ranges; ///< the shared channels, sorted and merged
      std::size_t n_channels = 0;

      // Set by the check when it has a geometry
      std::set<unsigned> tpcsets;
      std::set<std::pair<unsigned, unsigned>> tpcset_tpcs; ///< (TPCSet, TPC)
    };

    bool empty() const { return items_.empty(); }
    std::size_t size() const { return items_.size(); }
    Overlap const& operator[](std::size_t i) const { return items_[i]; }
    auto begin() const { return items_.begin(); }
    auto end() const { return items_.end(); }

    /// One line per overlap: the two tags, the channel count and up to
    /// `max_ranges` (from the check) channel ranges, then TPCSets and TPCs
    /// when a geometry was available. Empty string for no overlaps.
    std::string describe() const;

  private:
    friend class ChannelOverlapCheck;
    std::vector<Overlap> items_;
    std::size_t max_ranges_ = 10;
  };

  /// Finds channels present in more than one collection. Holds the geometry
  /// used to annotate the result and the formatting settings; stateless
  /// otherwise, so one instance per module serves every event.
  class ChannelOverlapCheck {
  public:
    /// The channels of one collection, with the collection's input tag.
    struct ChannelSet {
      std::string tag;
      std::vector<raw::ChannelID_t> channels;
    };

    /// @param geom       geometry for the TPCSet/TPC annotation; nullptr for channels only
    /// @param max_ranges channel ranges listed per overlap in describe()
    explicit ChannelOverlapCheck(geo::WireReadoutGeom const* geom = nullptr, std::size_t max_ranges = 10);

    /// Channels present in more than one of `sets`. A channel repeated inside
    /// a single collection is not an overlap.
    ChannelOverlaps operator()(std::vector<ChannelSet> const& sets) const;

    /// Same, from art handles. `T` needs a `Channel()` member, as
    /// raw::RawDigit and sim::SimChannel have.
    template <typename T>
    ChannelOverlaps operator()(std::vector<art::Handle<std::vector<T>>> const& handles) const
    {
      return (*this)(channel_sets(handles));
    }

    /// Sorts and deduplicates `channels` and merges them into contiguous ranges.
    static std::vector<ChannelRange> to_ranges(std::vector<raw::ChannelID_t> channels);

    template <typename T>
    static std::vector<ChannelSet> channel_sets(std::vector<art::Handle<std::vector<T>>> const& handles)
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

  private:
    /// Adds TPCSets and (TPCSet, TPC) pairs of the overlapping channels;
    /// channels unknown to the geometry are skipped.
    void annotate(ChannelOverlaps& overlaps) const;

    geo::WireReadoutGeom const* geom_;
    std::size_t max_ranges_;
  };

} // namespace dunetrigger

#endif // DUNETRIGGER_UTILITIES_CHANNELOVERLAP_H
