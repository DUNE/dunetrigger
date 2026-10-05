# Vendored lardata code (temporary)

`ArtDataHelper/` contains copies of `lardata/ArtDataHelper/GetManyByRegexTag.h`
and `lardata/ArtDataHelper/RegexTagMatcher.{h,cxx}`, which are not yet part of
a LArSoft release.

- Source: `alessandrothea/lardata`, branch `feature/getmanyregex_bugfix`,
  commit `079e218a7d20`. That is lardata `develop` plus the fix making
  `RegexTagMatcher::match()` build its regexes from the configured pattern.
- Changes: a two-line provenance comment at the top of each file, and
  `#include "lardata/ArtDataHelper/..."` rewritten to
  `#include "dunetrigger/vendor/lardata/ArtDataHelper/..."`.
  Code and namespace (`lar::util`) are unchanged.
- `RegexTagMatcher.cxx` is built as the library
  `dunetrigger_vendor_lardata_ArtDataHelper`.

## Removal, once a LArSoft release ships these files (with the fix)

The build stops with an error as soon as the lardata in use provides
`lardata/ArtDataHelper/RegexTagMatcher.h`. Then:

1. In every user, replace `"dunetrigger/vendor/lardata/ArtDataHelper/` with
   `"lardata/ArtDataHelper/` in the includes
   (`grep -rn vendor/lardata dunetrigger` lists them).
2. In their `CMakeLists.txt`, replace `dunetrigger_vendor_lardata_ArtDataHelper`
   with `lardata::ArtDataHelper` (or drop it where that is already linked).
3. Delete `dunetrigger/vendor/` and its `add_subdirectory` line in
   `dunetrigger/CMakeLists.txt`.
