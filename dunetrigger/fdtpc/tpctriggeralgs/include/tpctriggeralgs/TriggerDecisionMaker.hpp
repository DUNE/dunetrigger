/**
 * @file TriggerDecisionMaker.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERDECISIONMAKER_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERDECISIONMAKER_HPP_

#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/Issues.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/Logging.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/TriggerCandidate.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/TriggerDecision.hpp"

#include <atomic>
#include <chrono>
#include <nlohmann/json.hpp>
#include <vector>

namespace tpctriggeralgs {

class TriggerDecisionMaker
{
public:
  virtual ~TriggerDecisionMaker() = default;
  virtual void operator()(const TriggerCandidate& input_tc, std::vector<TriggerDecision>& output_tds) = 0;
  virtual void flush(std::vector<TriggerDecision>&) {}
  virtual void configure(const nlohmann::json&) {}
};

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERDECISIONMAKER_HPP_
