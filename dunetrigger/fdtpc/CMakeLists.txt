# basic source code CMakeLists.txt

add_subdirectory(tpctriggeralgs)
add_subdirectory(TPCTriggerSim)
add_subdirectory(TPCTriggerObjs)


install_headers()
install_source()
install_fhicl()
