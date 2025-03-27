#include "dunetrigger/triggeralgs/include/triggeralgs/AEAnomalyWindow/BinnedWindow.hpp"

namespace triggeralgs {

void BinnedWindow::addbin(WindowBin const &input_bin) {
  tp_window_bins.push_back(input_bin);
  ae_input.push_back(static_cast<float>(input_bin.adc_integral));
}

void BinnedWindow::resetwindow(WindowBin const &input_bin) {
        window_time_start = input_bin.time_start;
        tp_window_bins.clear();
        ae_input.clear();
}

void BinnedWindow::movebin(WindowBin const &input_bin) {
  // Add the next bin
  // All bins of equal length in time - so just pop out the front bin of TPs
  tp_window_bins.erase(tp_window_bins.begin());
  tp_window_bins.push_back(input_bin);

  ae_input.erase(ae_input.begin());
  ae_input.push_back(input_bin.adc_integral);
}

float BinnedWindow::sumadc() const {
  float sum = std::accumulate(ae_input.begin(), ae_input.end(), 0);
  return sum;
}

int BinnedWindow::bincount() const {
  return (int)ae_input.size();
}

} // namespace triggeralgs
