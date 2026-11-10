# basic source code CMakeLists.txt

add_subdirectory(tpctriggeralgs)
add_subdirectory(TPCTriggerSim)
add_subdirectory(TPCTriggerObjs)
add_subdirectory(TPCTriggerAna)


install_headers()
install_source()
install_fhicl()
