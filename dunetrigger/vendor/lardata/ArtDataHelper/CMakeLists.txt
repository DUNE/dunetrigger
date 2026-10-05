# TEMPORARY vendored copy of lardata/ArtDataHelper/{GetManyByRegexTag.h,RegexTagMatcher.*},
# see ../README.md.

# Once lardata ships these files, both copies would define the same lar::util
# symbols: stop the build so that this directory gets removed instead.
get_target_property(_lardata_include_dirs lardata::ArtDataHelper INTERFACE_INCLUDE_DIRECTORIES)
foreach(_dir IN LISTS _lardata_include_dirs)
  if(EXISTS "${_dir}/lardata/ArtDataHelper/RegexTagMatcher.h")
    message(FATAL_ERROR "lardata now provides lardata/ArtDataHelper/RegexTagMatcher.h: "
                        "remove dunetrigger/vendor/lardata as described in its README.md")
  endif()
endforeach()

cet_make_library(
  SOURCE RegexTagMatcher.cxx
  LIBRARIES PUBLIC
  canvas::canvas
)

install_headers()
install_source()
