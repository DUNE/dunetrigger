// Out-of-line method implementations for the row structs declared in TriggerAnaRows.hh.
// Built into the dunetrigger::TriggerAna library so that every analyzer can use them
// (a definition inside a module source would only be available to that plugin).

#include "dunetrigger/TriggerAna/TriggerAnaRows.hh"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"

void dunetrigger::TriggerPrimitiveRow::from_tp(const dunedaq::trgdataformats::TriggerPrimitive &tp) {
  version = 2; // temp, since variables below are converted to v2 version while the TP version in TriggerSim is still 1.
               // Go back to "= tp.version" after changing triggeralgs to v5 (and using TriggerPrimitive2.hpp as header)
  flag = 0;
  detid = tp.detid;
  channel = tp.channel;
  samples_over_threshold = tp.time_over_threshold / dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  time_start = tp.time_start;
  samples_to_peak = (tp.time_peak - tp.time_start) / dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  adc_integral = tp.adc_integral;
  adc_peak = tp.adc_peak;
}
