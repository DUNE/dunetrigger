#ifndef TRIGGERALGS_WINDOWBIN_HPP_
#define TRIGGERALGS_WINDOWBIN_HPP_

#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerPrimitive.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/Types.hpp"

#include <ostream>
#include <unordered_map>
#include <vector>

namespace triggeralgs {
  
// Make a class for a mini-window for each bin that goes into the AE input
// Re-used from original Window class in ADCSimpleWindow
// Represents a single bin of TPs in time that goes into a window that is an input to an AE

class WindowBin {
public:
  bool is_empty() const;
      
  void add(TriggerPrimitive const &input_tp);
      
  void clear();
      
  void move(TriggerPrimitive const &input_tp, timestamp_t const &window_length);
      
  void reset(TriggerPrimitive const &input_tp);
      
  friend std::ostream& operator<<(std::ostream& os, const WindowBin& window);

  timestamp_t time_start;
  uint32_t adc_integral;
  std::vector<TriggerPrimitive> tp_list;
};

} // namespace triggeralgs

#endif
