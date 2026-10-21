# Vendored lardata code (temporary)

`ArtDataHelper/` contains copies of `lardata/ArtDataHelper/GetManyByRegexTag.h`
and `lardata/ArtDataHelper/RegexTagSelector.{h,cxx}`, which are not yet part of
a LArSoft release.

- Source: `alessandrothea/lardata`, branch `feature/getmanybyregex_explicitprocess`,
  commit `8ec849769e3d`. It replaces the earlier `RegexTagMatcher` (vendored from
  `feature/getmanyregex_bugfix`, commit `079e218a7d20`) with an `art::SelectorBase`,
  `RegexTagSelector`, passed to `getInputTags()`. Semantics differ from the old
  matcher:
  - an empty label or instance pattern matches only an empty field (use `".*"`
    for "any");
  - an empty process keeps only the newest process with a *valid* matching
    product, falling back to older processes like `e.getHandle<T>()`;
  - an explicit process is a regex and keeps the products of *every* matching
    process (e.g. `".*"` returns duplicates from all processes).

  See the doc comment in `GetManyByRegexTag.h`.
- Changes: a three-line provenance comment at the top of each file,
  `#include "lardata/ArtDataHelper/..."` rewritten to
  `#include "dunetrigger/vendor/lardata/ArtDataHelper/..."`, namespace
  `lar::util` wrapped as `dunetrigger::vendor::lar::util`, and include guards
  prefixed with `DUNETRIGGER_VENDOR_`. The code is otherwise unchanged. Thanks
  to the namespace and guards, the vendored copy coexists with lardata's: both
  can be linked, and both headers included, without clashes.
- `RegexTagSelector.cxx` is built as the library
  `dunetrigger_vendor_lardata_ArtDataHelper`.

## Removal, once a LArSoft release ships these files (with the fix)

CMake prints a warning as soon as the lardata in use provides
`lardata/ArtDataHelper/RegexTagSelector.h` (check in `ArtDataHelper/CMakeLists.txt`);
the build carries on with the vendored copy.

**Check the fix first.** lardata `v10_01_00` already ships `GetManyByRegexTag.h`
with the older `RegexTagMatcher.{h,cxx}`, *without* the fix: its `match()` builds
the regexes from the candidate tag (`tag.instance()`, ...) instead of the
configured pattern, so a pattern such as `"tpcrawdecoder:simpleSC.*"` matches
nothing. That release does not trigger the warning (no `RegexTagSelector.h`), but do
not switch to it. Remove the vendored copy only once the installed lardata
ships `RegexTagSelector`:

```
ls $LARDATA_DIR/source/lardata/ArtDataHelper/RegexTagSelector.cxx
```

(`$LARDATA_DIR` from `ups active`.) Also check that its process semantics match
the vendored commit (see `GetManyByRegexTag.h`), since the jobs here rely on
them. Then:

1. In every user, replace `"dunetrigger/vendor/lardata/ArtDataHelper/` with
   `"lardata/ArtDataHelper/` in the includes, and
   `dunetrigger::vendor::lar::util::getManyByRegexTag` with
   `lar::util::getManyByRegexTag` (`grep -rn vendor/lardata dunetrigger` lists them).
2. In their `CMakeLists.txt`, replace `dunetrigger_vendor_lardata_ArtDataHelper`
   with `lardata::ArtDataHelper` (or drop it where that is already linked).
3. Delete `dunetrigger/vendor/` and its `add_subdirectory` line in
   `dunetrigger/CMakeLists.txt`.
