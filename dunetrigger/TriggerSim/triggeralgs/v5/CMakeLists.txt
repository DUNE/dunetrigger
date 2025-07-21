add_compile_options(-O2 -flto -Wno-error=sign-compare -Wno-error=unused-variable)

include(BasicPlugin)
include_directories("${dunedetdataformats_DIR}/../../../include")

basic_plugin(triggeralgs_v5 module
  BASENAME_ONLY
  LIBRARIES
  OfflineTPCChannelMap_module
  SOURCE
  src/TAMakerADCSimpleWindowAlgorithm.cpp
  src/TCMakerADCSimpleWindowAlgorithm.cpp
  src/TAMakerHorizontalMuonAlgorithm.cpp
  src/TCMakerHorizontalMuonAlgorithm.cpp
  src/TCMakerPlaneCoincidenceAlgorithm.cpp
  src/TAMakerPlaneCoincidenceAlgorithm.cpp
  src/TAMakerMichelElectronAlgorithm.cpp
  src/TCMakerMichelElectronAlgorithm.cpp
  src/TAMakerPrescaleAlgorithm.cpp
  src/TCMakerPrescaleAlgorithm.cpp
  src/TAMakerSupernovaAlgorithm.cpp
  src/TCMakerSupernovaAlgorithm.cpp
  src/TDMakerSupernovaAlgorithm.cpp
  src/TAMakerDBSCANAlgorithm.cpp
  src/TCMakerDBSCANAlgorithm.cpp
  src/TAMakerChannelDistanceAlgorithm.cpp
  src/TCMakerChannelDistanceAlgorithm.cpp
  src/TAMakerBundleNAlgorithm.cpp
  src/TCMakerBundleNAlgorithm.cpp
  src/TAMakerChannelAdjacencyAlgorithm.cpp
  src/TCMakerChannelAdjacencyAlgorithm.cpp
  src/TAWindow.cpp
  src/TPWindow.cpp
  src/dbscan/dbscan.cpp
  src/dbscan/Hit.cpp
)

