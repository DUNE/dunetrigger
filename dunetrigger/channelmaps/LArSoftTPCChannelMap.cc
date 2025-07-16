
#include "dunetrigger/channelmaps/LArSoftTPCChannelMap.hpp"
#include <larcoreobj/SimpleTypesAndConstants/RawTypes.h>
#include "larcore/Geometry/WireReadout.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <memory>

using namespace dunedaq::detchannelmaps;

uint LArSoftTPCChannelMap::get_plane_from_offline_channel(uint offchannel){
    geo::WireReadoutGeom const* geom = &art::ServiceHandle<geo::WireReadout>()->Get();
    return geom->ChannelToROP(static_cast<raw::ChannelID_t>(offchannel)).deepestIndex();
}

uint LArSoftTPCChannelMap::get_element_id_from_offline_channel(uint offchannel){
    geo::WireReadoutGeom const* geom = &art::ServiceHandle<geo::WireReadout>()->Get();
    return geom->ChannelToROP(static_cast<raw::ChannelID_t>(offchannel)).parentID().TPCset;
}

uint LArSoftTPCChannelMap::get_first_channel_on_plane(uint offchannel){
    uint c = 0; // hard coded for protodune - only ever 1 cryostat
    uint s = get_element_id_from_offline_channel(offchannel);
    uint r = get_plane_from_offline_channel(offchannel);

    readout::ROPID rop(c, s, r);

    geo::WireReadoutGeom const* geom = &art::ServiceHandle<geo::WireReadout>()->Get();
    return geom->FirstChannelInROP(rop);
}

uint LArSoftTPCChannelMap::get_nchannels_on_plane(uint offchannel){
    uint c = 0; // hard coded for protodune - only ever 1 cryostat
    uint s = get_element_id_from_offline_channel(offchannel);
    uint r = get_plane_from_offline_channel(offchannel);

    readout::ROPID rop(c, s, r);

    geo::WireReadoutGeom const* geom = &art::ServiceHandle<geo::WireReadout>()->Get();
    return geom->Nchannels(rop);
}
