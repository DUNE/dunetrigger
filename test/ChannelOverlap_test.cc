// Unit test for dunetrigger::ChannelOverlapCheck: no art, no geometry.

#include "dunetrigger/Utilities/ChannelOverlap.h"

#include <cassert>
#include <iostream>

using dunetrigger::ChannelOverlapCheck;
using dunetrigger::ChannelOverlaps;
using dunetrigger::ChannelRange;
using ChannelSet = ChannelOverlapCheck::ChannelSet;

namespace {
  std::vector<raw::ChannelID_t> range(raw::ChannelID_t lo, raw::ChannelID_t hi)
  {
    std::vector<raw::ChannelID_t> v;
    for (raw::ChannelID_t ch = lo; ch <= hi; ++ch) v.push_back(ch);
    return v;
  }
}

int main()
{
  ChannelOverlapCheck const check; // no geometry: channel-level only

  // to_ranges: sorts, deduplicates, merges contiguous IDs
  {
    auto const ranges = ChannelOverlapCheck::to_ranges({5, 3, 4, 10, 4, 12, 11});
    assert(ranges.size() == 2);
    assert(ranges[0] == ChannelRange(3, 5));
    assert(ranges[1] == ChannelRange(10, 12));
    assert(ChannelOverlapCheck::to_ranges({}).empty());
  }

  // Disjoint collections: no overlap
  {
    std::vector<ChannelSet> sets{{"a", range(0, 99)}, {"b", range(100, 199)}, {"c", range(200, 299)}};
    ChannelOverlaps const overlaps = check(sets);
    assert(overlaps.empty());
    assert(overlaps.describe().empty());
  }

  // A channel repeated inside one collection is not an overlap
  {
    std::vector<ChannelSet> sets{{"a", {1, 2, 2, 3}}, {"b", {4, 5}}};
    assert(check(sets).empty());
  }

  // Two collections sharing a block: one overlap, attributed to the right pair
  {
    std::vector<ChannelSet> sets{{"a", range(0, 99)}, {"b", range(50, 149)}, {"c", range(150, 199)}};
    ChannelOverlaps const overlaps = check(sets);
    assert(overlaps.size() == 1);
    assert(overlaps[0].first == "a");
    assert(overlaps[0].second == "b");
    assert(overlaps[0].n_channels == 50);
    assert(overlaps[0].ranges.size() == 1);
    assert(overlaps[0].ranges[0] == ChannelRange(50, 99));
    assert(overlaps[0].tpcsets.empty()); // no geometry, no annotation

    auto const text = overlaps.describe();
    assert(text.find("a and b: 50 channels (50-99)") != std::string::npos);
    assert(text.find("TPCSets") == std::string::npos);
  }

  // Three collections, same channels in all: every pair against the first owner
  {
    std::vector<ChannelSet> sets{{"a", range(0, 9)}, {"b", range(0, 9)}, {"c", range(0, 9)}};
    ChannelOverlaps const overlaps = check(sets);
    assert(overlaps.size() == 2);
    assert(overlaps[0].first == "a" && overlaps[0].second == "b");
    assert(overlaps[1].first == "a" && overlaps[1].second == "c");
    std::size_t n = 0;
    for (auto const& o : overlaps) n += o.n_channels; // range-for
    assert(n == 20);
  }

  // max_ranges from the check is carried into describe()
  {
    ChannelOverlapCheck const check_two{nullptr, 2};
    std::vector<ChannelSet> sets{{"a", {0, 2, 4, 6}}, {"b", {0, 2, 4, 6}}};
    auto const text = check_two(sets).describe();
    assert(text.find("4 channels (0, 2, ... 2 more ranges)") != std::string::npos);
  }

  std::cout << "ChannelOverlap_test: all checks passed" << std::endl;
  return 0;
}
