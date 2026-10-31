include_directories("${dunedetdataformats_DIR}/../../../include") # there's a neater way?

add_subdirectory(fcl)
add_subdirectory(exe)

# Library with the analysis helpers shared by the modules (dunetrigger::TriggerAna)
cet_make_library(
  SOURCE
  TPBacktracker.cc
  LIBRARIES
  PUBLIC
  lardataobj::Simulation
  larsim::MCCheater_ParticleInventoryService_service
  nurandom::RandomUtils_NuRandomService_service  # brings CLHEP (pulled in by ParticleInventoryService headers)
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  PRIVATE
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
)

# Libraries used by all the modules of this directory
set(triggerana_module_libraries
  Boost::headers
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::Simulation
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
  triggeralgs_module
  dunetrigger_vendor_lardata_ArtDataHelper
)

cet_build_plugin(TriggerAnaTree art::module BASENAME_ONLY
  LIBRARIES PRIVATE
  dunetrigger::TriggerAna
  ${triggerana_module_libraries}
)

cet_build_plugin(TriggerActivityTxtDump art::module BASENAME_ONLY
  LIBRARIES PRIVATE ${triggerana_module_libraries}
)

cet_build_plugin(TriggerTPCInfoComparator art::module BASENAME_ONLY
  LIBRARIES PRIVATE ${triggerana_module_libraries}
)

cet_build_plugin(TriggerTPCInfoDisplay art::module BASENAME_ONLY
  LIBRARIES PRIVATE ${triggerana_module_libraries}
)

install_fhicl()
install_headers()
install_source()
