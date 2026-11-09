# Wire-Cell Toolkit components for running TPG inside a WCT graph
# (wclsTPGFilter, wclsTPGSaver). The library name is what goes in the
# WireCellToolkit "plugins" list in FHiCL.
cet_make_library(LIBRARY_NAME WireCellDuneTrigger
  SOURCE
  TPGFilter.cxx
  TPGSaver.cxx
  LIBRARIES
  PRIVATE
  larwirecell::IArtEventVisitor
  WireCell::Aux
  WireCell::Iface
  WireCell::Util
  art::Framework_Core
  art::Framework_Principal
  art_plugin_support::toolMaker
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)

add_subdirectory(cfg)
add_subdirectory(fcl)

install_headers()
install_source()
