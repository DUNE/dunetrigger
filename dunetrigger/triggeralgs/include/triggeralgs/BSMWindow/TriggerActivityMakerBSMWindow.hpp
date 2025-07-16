/**
 * @file TriggerActivityMakerBSMWindow.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRIGGERALGS_BSMWINDOW_TRIGGERACTIVITYMAKERBSMWINDOW_HPP_
#define TRIGGERALGS_BSMWINDOW_TRIGGERACTIVITYMAKERBSMWINDOW_HPP_

#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivityFactory.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/Types.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/BSMWindow/BSMWindow.hpp"
//#include "dunetrigger/triggeralgs/include/triggeralgs/BSMWindow/TreeliteModelInterface.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/BSMWindow/CompiledModelInterface.hpp"
#include "dunetrigger/channelmaps/OfflineTPCChannelMap.hpp"

//#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

#include <fstream>
#include <vector>
#include <algorithm>

#define safe_treelite(call) {  \
  int err = (call); \
  if (err != 0) { \
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
      ": error in " + #call + ":" + TreeliteGetLastError());  \
  } \
}

namespace triggeralgs {
class TriggerActivityMakerBSMWindow : public TriggerActivityMaker
{

public:
  void operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta);
  
  void configure(const nlohmann::json &config);

  ~TriggerActivityMakerBSMWindow() override;

private:

  bool compute_treelite_classification();

  bool compute_dt_filter();

  TriggerActivity construct_ta() const;

  BSMWindow m_current_window;

  timestamp_t m_last_pred_time;
  uint64_t m_primitive_count = 0;

  // Try batching the input - introduce some latency,
  // but maybe speed up the algorithm
  const int nbatch = 1;
  // BDT batching takes a row-major flat array
  std::vector<float> flat_batched_inputs;
  // row-major input for Entry objects used for compiled model
  std::vector<Entry> flat_batched_Entries;
  // Keep track of the number of windows in the current batch
  int nbatch_iterator = 0;

  // Configurable parameters.
  uint32_t m_adc_threshold = 1200000;
  float m_ratio_threshold = 0.65;
  float m_bdt_threshold = 0.95;
  timestamp_t m_window_length = 20000;
  // now the length of the bin
  timestamp_t m_bin_length = 1000;
  int nbins = 20;
  int m_algtype = 0;

  // Geometry information for binning
  geo::WireReadoutGeom const *m_geom;
  std::string m_channel_map_name = "VDColdboxChannelMap";
  std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> channelMap = 
    dunedaq::detchannelmaps::make_map(m_channel_map_name);
  //std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> channelMap; 
  readout::ROPID m_rop;
  //unsigned int m_first_channel;
  //unsigned int m_last_channel;
  channel_t m_first_channel;
  channel_t m_last_channel;
  channel_t m_chan_bin_length = 50;
  int nchanbins = 10;

  // Treelite model
  //std::unique_ptr<TreeliteModelInterface> m_treelite_model_interface;
  // Compiled treelite model interface
  std::unique_ptr<CompiledModelInterface> m_compiled_model_interface;

};
} // namespace triggeralgs

#endif // TRIGGERALGS_BSMWINDOW_TRIGGERACTIVITYMAKERBSMWINDOW_HPP_
