# Jsonnet for the Wire-Cell TPG graphs. Installed under wire-cell-cfg/dunetrigger,
# which the dunetrigger UPS table prepends to WIRECELL_PATH.
file(GLOB jsonnet_files dune-vd/*.jsonnet)
install_wp(LIST ${jsonnet_files} SUBDIRNAME dunetrigger/dune-vd)
