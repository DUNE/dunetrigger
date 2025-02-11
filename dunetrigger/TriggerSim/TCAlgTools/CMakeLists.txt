cet_build_plugin(TCAlgTPCExample art::tool
  LIBRARIES
  dunecore::DuneObj
  fhiclcpp::fhiclcpp 
  cetlib::cetlib 
  cetlib_except::cetlib_except
)

install_fhicl()
install_headers()
install_source()
