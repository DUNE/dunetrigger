/**
 * @file TriggerCandidate.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERCANDIDATE_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERCANDIDATE_HPP_

#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityData.hpp"
#include "dunetrigger/fdtpc/tpctriggeralgs/include/tpctriggeralgs/TriggerCandidateData.hpp"

#include <vector>

namespace tpctriggeralgs {

struct TriggerCandidate : public TriggerCandidateData
{
  std::vector<TriggerActivityData> inputs;
};

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERCANDIDATE_HPP_
