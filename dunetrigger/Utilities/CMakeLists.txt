# Shared helpers for the TriggerSim / TriggerAna / TriggerTools modules.
# Builds the library dunetrigger_Utilities.

cet_make_library(
  SOURCE
    ChannelOverlap.cxx
  LIBRARIES PUBLIC
    larcorealg::Geometry
    art::Framework_Principal   # art::Handle in ChannelOverlap.h (from_handles)
    canvas::canvas
)

install_headers()
install_source()
