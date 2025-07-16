#ifndef DUNETRIGGER_CHANNELMAPS_LARSOFTTPCCHANNELMAP_HPP
#define DUNETRIGGER_CHANNELMAPS_LARSOFTTPCCHANNELMAP_HPP

#include "dunetrigger/channelmaps/OfflineTPCChannelMap.hpp"

namespace dunedaq::detchannelmaps {

class LArSoftTPCChannelMap : public TPCChannelMap {
  public:
    LArSoftTPCChannelMap() = default;
    ~LArSoftTPCChannelMap() override = default;
    // overrides of virtual methods like get_first_channel_on_plane(...)
    uint get_plane_from_offline_channel(uint offchannel) override;
    uint get_element_id_from_offline_channel(uint offchannel) override;
    uint get_first_channel_on_plane(uint offchannel) override;
    uint get_nchannels_on_plane(uint offchannel) override;

    // No implementation yet
    uint get_offline_channel_from_crate_slot_fiber_chan(uint crate,
                                                        uint slot,
                                                        uint fiber,
                                                        uint channel) { return 0; }

    uint get_offline_channel_from_crate_slot_stream_chan(uint crate,
                                                        uint slot,
                                                        uint fiber,
                                                        uint channel) { return 0; }

};

}

#endif
