include(BasicPlugin)

basic_plugin(OfflineTPCChannelMap module
            BASENAME_ONLY
            LIBRARIES
            lardataobj::RawData
            dunecore::DuneObj
            lardata::ArtDataHelper
            nusimdata::SimulationBase
            SOURCE
            OfflineTPCChannelMap_module.cc
            LArSoftTPCChannelMap.cc
            # ${ART_FRAMEWORK_CORE}
            # ${ART_FRAMEWORK_PRINCIPAL}
            # ${ART_FRAMEWORK_SERVICES_REGISTRY}
)
