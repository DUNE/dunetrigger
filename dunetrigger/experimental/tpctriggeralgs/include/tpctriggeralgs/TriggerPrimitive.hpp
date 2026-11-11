/**
 * @file TriggerPrimitive.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERPRIMITIVE_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERPRIMITIVE_HPP_

#include "detdataformats/trigger/TriggerPrimitive.hpp"

namespace tpctriggeralgs {
// TriggerPrimitive looks a bit difference to TriggerActivity and
// TriggerCandidate: the latter two have distinct overlay and
// non-overlay versions, but TriggerPrimitive does not. So whereas TA
// and TC have (non-overlay) definitions in tpctriggeralgs, TP has a
// definition only in detdataformats. So we just copy it here so that
// references to TP, TA and TC in tpctriggeralgs and downstream code
// "look" the same
using TriggerPrimitive = dunedaq::trgdataformats::TriggerPrimitive;

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_TRIGGERPRIMITIVE_HPP_
