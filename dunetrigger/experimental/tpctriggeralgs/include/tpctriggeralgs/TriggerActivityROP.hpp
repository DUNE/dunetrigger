/**
 * @file TriggerActivityROP.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERACTIVITYROP_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERACTIVITYROP_HPP_

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityROPData.hpp"
#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerPrimitive.hpp"

#include <vector>

namespace tpctriggeralgs {

struct TriggerActivityROP : public TriggerActivityROPData
{
  std::vector<TriggerPrimitive> inputs;
};

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERACTIVITYROP_HPP_
