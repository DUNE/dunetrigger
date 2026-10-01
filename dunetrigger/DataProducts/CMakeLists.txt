# dunesw/dune/DataProducts/CMakeLists.txt


include_directories("${dunedaqdataformats_DIR}/../../../include")
include_directories("${dunedetdataformats_DIR}/../../../include")

cet_make_library(
  LIBRARY_NAME DataProducts
  INTERFACE
  SOURCE MyProduct.h
  TriggerActivityData.hpp # no .cc files — header-only data product
  # lardataobj::RecoBase
  # ROOT::EG
)


build_dictionary(DICTIONARY_LIBRARIES PRIVATE
  dunetrigger::DataProducts
  canvas::canvas
)

install_headers()
install_source()
