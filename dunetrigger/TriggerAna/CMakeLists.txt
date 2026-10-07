include_directories("${dunedetdataformats_DIR}/../../../include") # there's a neater way?

add_subdirectory(fcl)
add_subdirectory(exe)

# TPBacktracker.cc is the only non-plugin source here: art_make builds it into a
# library (dunetrigger_TriggerAna) that the modules link against.
art_make(BASENAME_ONLY
  LIB_LIBRARIES
  Boost::headers
  canvas::canvas
  lardataobj::Simulation
  lardataobj::RawData
  nusimdata::SimulationBase
  larsim::MCCheater_ParticleInventoryService_service
  nurandom::RandomUtils_NuRandomService_service  # brings CLHEP (pulled in by ParticleInventoryService headers)
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Provenance
  art::Persistency_Common
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
  MODULE_LIBRARIES
  dunetrigger::TriggerAna
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  dunecore::DuneObj
  lardata::ArtDataHelper
  lardata::DetectorPropertiesService
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
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
  larsim::MCCheater_ParticleInventoryService_service
  triggeralgs_module
  dunetrigger_vendor_lardata_ArtDataHelper
)


install_fhicl()
install_headers()
install_source()
