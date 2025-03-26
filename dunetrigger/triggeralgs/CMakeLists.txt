add_compile_options(-O3 -flto -march=native -funroll-loops -Wno-error=sign-compare -Wno-error=unused-variable)

include(BasicPlugin)
include_directories("${dunedetdataformats_DIR}/../../../include")
include_directories ( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

basic_plugin(triggeralgs module
  BASENAME_ONLY
  LIBRARIES
  OfflineTPCChannelMap_module
  pthread
  TensorFlow::cc
  TensorFlow::framework
  SOURCE
  src/TAWindow.cpp
  src/TPWindow.cpp
  src/TriggerActivityMakerADCSimpleWindow.cpp
  src/TriggerActivityMakerAEAnomalyWindow.cpp
  src/TriggerActivityMakerBundleN.cpp
  src/TriggerActivityMakerChannelAdjacency.cpp
  src/TriggerActivityMakerChannelDistance.cpp
  src/TriggerActivityMakerDBSCAN.cpp
  src/TriggerActivityMakerHorizontalMuon.cpp
  src/TriggerActivityMakerMichelElectron.cpp
  src/TriggerActivityMakerPlaneCoincidence.cpp
  src/TriggerActivityMakerPrescale.cpp
  src/TriggerActivityMakerSupernova.cpp
  src/TriggerCandidateMakerADCSimpleWindow.cpp
  src/TriggerCandidateMakerAEAnomalyWindow.cpp
  src/TriggerCandidateMakerBundleN.cpp
  src/TriggerCandidateMakerChannelAdjacency.cpp
  src/TriggerCandidateMakerChannelDistance.cpp
  src/TriggerCandidateMakerDBSCAN.cpp
  src/TriggerCandidateMakerHorizontalMuon.cpp
  src/TriggerCandidateMakerMichelElectron.cpp
  src/TriggerCandidateMakerPlaneCoincidence.cpp
  src/TriggerCandidateMakerPrescale.cpp
  src/TriggerCandidateMakerSupernova.cpp
  src/TriggerDecisionMakerSupernova.cpp
  src/dbscan/dbscan.cpp
  src/dbscan/Hit.cpp
)
