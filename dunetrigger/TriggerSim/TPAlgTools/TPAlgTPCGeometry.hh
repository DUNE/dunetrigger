#ifndef DUNETRIGGER_TRIGGERSIM_TPALGTPCGEOMETRY_hh
#define DUNETRIGGER_TRIGGERSIM_TPALGTPCGEOMETRY_hh

// Geometry lookup shared by the art-facing TPAlgTPC tools. Kept separate so
// the algorithms themselves (and TPAlgTPCTool.hh) do not depend on art
// services.

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "detdataformats/trigger/TriggerPrimitive.hpp"

namespace dunetrigger {

  // Readout plane index (0, 1, 2) of an offline channel.
  inline unsigned int tpc_plane_of_channel(dunedaq::trgdataformats::channel_t const channel)
  {
    geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    return wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(channel)).at(0).Plane;
  }

}

#endif
