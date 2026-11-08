# basic source code CMakeLists.txt

add_subdirectory(channelmaps)
add_subdirectory(triggeralgs)
add_subdirectory(TriggerSim)
add_subdirectory(TriggerAna)
add_subdirectory(TriggerTools)

install_headers()
install_source()
install_fhicl()
