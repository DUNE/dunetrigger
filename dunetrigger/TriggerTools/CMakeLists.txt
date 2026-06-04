include_directories("${dunedetdataformats_DIR}/../../../include") # there's a neater way?

add_subdirectory(fcl)
add_subdirectory(exe)

art_make(BASENAME_ONLY MODULE_LIBRARIES
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  larcorealg::Geometry
  lardataobj::RawData
  dunecore::DuneObj
  lardata::ArtDataHelper
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  canvas::canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::Tree
  triggeralgs_module
  
)

# Example for header-only lib
# cet_make_library(LIBRARY_NAME UtilitiesHeaders INTERFACE
#   SOURCE
#   constexpr_math.h
#   intervals_fhicl.h
#   intervals.h
#   MappedContainer.h
#   MultipleChoiceSelection.h
#   quantities_fhicl.h
#   quantities.h
#   StatCollector.h
#   TrackingTypeUtils.h
#   LIBRARIES INTERFACE
#   larcorealg::headers
#   lardataobj::RecoBase
#   fhiclcpp::fhiclcpp
#   Boost::boost
# )


install_fhicl()
install_headers()
install_source()




