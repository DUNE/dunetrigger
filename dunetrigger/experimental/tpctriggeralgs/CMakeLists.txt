add_compile_options(-O2 -flto -Wno-error=sign-compare -Wno-error=unused-variable)

include(BasicPlugin)
include_directories("${dunedetdataformats_DIR}/../../../include")

basic_plugin(tpctriggeralgs module
  BASENAME_ONLY
  LIBRARIES
  dunetrigger::TPCTriggerObjs
  SOURCE
  src/TriggerActivityROPMakerSWIFT.cpp
  src/TriggerActivityROPMakerPassthrough.cpp
)
