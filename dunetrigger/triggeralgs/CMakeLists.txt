add_compile_options(-O2 -flto -Wno-error=sign-compare -Wno-error=unused-variable)

include(BasicPlugin)
include_directories("${dunedetdataformats_DIR}/../../../include")

basic_plugin(triggeralgs module
  BASENAME_ONLY
  LIBRARIES
  OfflineTPCChannelMap_module
  dunetrigger::DataProducts
  SOURCE
  src/TriggerActivityMakerSWIFT.cpp
  src/TriggerActivityMakerPrescale.cpp
  src/TriggerCandidateMakerPrescale.cpp
)
