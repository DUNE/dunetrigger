include_directories("${dunedetdataformats_DIR}/../../../include")

add_subdirectory(fcl)

cet_build_plugin(TPCTriggerAnaTree art::module BASENAME_ONLY
  LIBRARIES PRIVATE
  dunetrigger::TriggerAna
  Boost::headers
  lardata::ArtDataHelper
  larcore::Geometry_Geometry_service
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  canvas::canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::Core
  ROOT::Tree
  dunetrigger::TPCTriggerObjs
)

install_headers()
install_source()
