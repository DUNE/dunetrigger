# dunetrigger/experimental/TPCTriggerObjs/CMakeLists.txt


include_directories("${dunedaqdataformats_DIR}/../../../include")
include_directories("${dunedetdataformats_DIR}/../../../include")

cet_make_library(
  LIBRARY_NAME TPCTriggerObjs
  INTERFACE
  SOURCE
  TriggerActivityData.hpp # no .cc files — header-only data product
  # lardataobj::RecoBase
  # ROOT::EG
)


build_dictionary(DICTIONARY_LIBRARIES PRIVATE
  dunetrigger::TPCTriggerObjs
  canvas::canvas
)

install_headers()
install_source()
