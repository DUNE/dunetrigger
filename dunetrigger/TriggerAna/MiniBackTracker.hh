#ifndef DUNETRIGGER_MINIBACKTRACKER_HH
#define DUNETRIGGER_MINIBACKTRACKER_HH

#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <algorithm>
#include <vector>

namespace dunetrigger {

class MiniBackTracker {
public:
  MiniBackTracker(art::Handle<std::vector<sim::SimChannel>>& simchans_handle)
    : fSimChannelsHandle(simchans_handle) {}

  art::Ptr<sim::SimChannel> findSimChannelPtr(raw::ChannelID_t channel) const;
  std::vector<double> simIDEsToXYZ(std::vector<sim::IDE> const& ides) const;

private:
  art::Handle<std::vector<sim::SimChannel>> fSimChannelsHandle;
  mutable std::vector<art::Ptr<sim::SimChannel>> fSimChannels;

  std::vector<art::Ptr<sim::SimChannel>>& simchannels() const;
};

//-----------------------------------------------------------------------
inline art::Ptr<sim::SimChannel>
MiniBackTracker::findSimChannelPtr(raw::ChannelID_t channel) const
{
  auto ilb = std::lower_bound(simchannels().begin(), simchannels().end(), channel,
    [](art::Ptr<sim::SimChannel> const& a, raw::ChannelID_t ch) { return (a->Channel() < ch); });
  return ((ilb != simchannels().end()) && ((*ilb)->Channel() == channel)) ? *ilb : art::Ptr<sim::SimChannel>{};
}


//------------------------------------------------------------------------------
std::vector<double>
MiniBackTracker::simIDEsToXYZ(std::vector<sim::IDE> const& ides) const
{
  std::vector<double> xyz(3, 0.0);
  double w = 0.0;
  for (auto const& ide : ides) {
    double weight = ide.numElectrons;
    w += weight;
    xyz[0] += (weight * ide.x);
    xyz[1] += (weight * ide.y);
    xyz[2] += (weight * ide.z);
  }
  if (w < 1.e-5)
    throw cet::exception("BackTracker") << "No sim::IDEs providing non-zero number of electrons"
                                        << " can't determine originating location from truth\n";
  xyz[0] = xyz[0] / w;
  xyz[1] = xyz[1] / w;
  xyz[2] = xyz[2] / w;
  return xyz;
}



//-----------------------------------------------------------------------
inline std::vector<art::Ptr<sim::SimChannel>>&
MiniBackTracker::simchannels() const
{
  if (fSimChannels.empty()) {
    art::fill_ptr_vector(fSimChannels, fSimChannelsHandle);

    auto comparesclambda = [](art::Ptr<sim::SimChannel> a, art::Ptr<sim::SimChannel> b) {
      return (a->Channel() < b->Channel());
    };

    if (!std::is_sorted(fSimChannels.begin(), fSimChannels.end(), comparesclambda))
      std::sort(fSimChannels.begin(), fSimChannels.end(), comparesclambda);
  }
  return fSimChannels;
}



} // namespace dunetrigger

#endif // DUNETRIGGER_MINIBACKTRACKER_HH
