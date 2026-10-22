# TEMPORARY vendored copy of lardata/ArtDataHelper/{GetManyByRegexTag.h,RegexTagSelector.*},
# see ../README.md.

# The vendored code lives in namespace dunetrigger::vendor, so it coexists with
# lardata's copy; once lardata ships these files, remind to switch over.
get_target_property(_lardata_include_dirs lardata::ArtDataHelper INTERFACE_INCLUDE_DIRECTORIES)
foreach(_dir IN LISTS _lardata_include_dirs)
  if(EXISTS "${_dir}/lardata/ArtDataHelper/RegexTagSelector.h")
    message(WARNING "lardata now provides lardata/ArtDataHelper/RegexTagSelector.h: "
                    "consider switching to it and removing dunetrigger/vendor/lardata "
                    "as described in its README.md")
  endif()
endforeach()

cet_make_library(
  SOURCE RegexTagSelector.cxx
  LIBRARIES PUBLIC
  art::Framework_Principal
  canvas::canvas
)

install_headers()
install_source()
