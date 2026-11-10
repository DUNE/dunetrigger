/**
 * @file TriggerActivityMaker.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERACTIVITYMAKER_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERACTIVITYMAKER_HPP_

#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/Issues.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/Logging.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/TriggerActivity.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/TriggerPrimitive.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/Types.hpp"

#include <atomic>
#include <chrono>
#include <nlohmann/json.hpp>
#include <vector>

namespace tpctriggeralgs {

class TriggerActivityMaker
{
public:
  virtual ~TriggerActivityMaker() = default;
  virtual void operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta) = 0;
  virtual void flush(timestamp_t /* until */, std::vector<TriggerActivity>&) {}
  virtual void configure(const nlohmann::json&) {}
};

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERACTIVITYMAKER_HPP_
