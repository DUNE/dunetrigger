cet_build_plugin(TPAlgTPCExample art::tool
  LIBRARIES
  fhiclcpp::fhiclcpp 
  cetlib::cetlib 
  cetlib_except::cetlib_except
)
cet_build_plugin(TPAlgTPCSimpleThreshold art::tool
  LIBRARIES
  larcore::Geometry_Geometry_service
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)
install_fhicl()
install_headers()
install_source()
