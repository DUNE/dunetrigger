include_directories("${dunedetdataformats_DIR}/../../../include") # there's a neater way?

add_subdirectory(TPAlgTools)
add_subdirectory(TAAlgTools)
add_subdirectory(TCAlgTools)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  lardataobj::RawData
  dunecore::DuneObj
  lardata::ArtDataHelper
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
)

install_fhicl()
install_headers()
install_source()
