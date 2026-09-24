add_compile_options(-O2 -flto -Wno-error=sign-compare -Wno-error=unused-variable)

include(BasicPlugin)
include_directories("${dunedetdataformats_DIR}/../../../include")

basic_plugin(triggeralgs module
  BASENAME_ONLY
  LIBRARIES
  OfflineTPCChannelMap_module
  SOURCE
  src/TriggerActivityMakerSWIFT.cpp
  src/TriggerCandidateMakerPrescale.cpp
)
