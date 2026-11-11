/**
 * @file TriggerActivityData.hpp
 *
 * Experimental TriggerActivity data product with a minimal set of fields.
 * The per-ROP counterpart is TriggerActivityROPData.
 */

#ifndef DUNETRIGGER_EXPERIMENTAL_TPCTRIGGEROBJS_TRIGGERACTIVITYDATA_HPP_
#define DUNETRIGGER_EXPERIMENTAL_TPCTRIGGEROBJS_TRIGGERACTIVITYDATA_HPP_

#include "detdataformats/trigger/Types.hpp"

namespace dunedaq::trgdataformats::exp {

struct TriggerActivityData
{
  // Update this version number if there are any changes to the in-memory representation of this class!
  static constexpr version_t s_trigger_activity_version = 1; // NOLINT(build/unsigned)

  version_t version = s_trigger_activity_version; // NOLINT(build/unsigned)
  timestamp_t time_start = INVALID_TIMESTAMP;
  timestamp_t time_end = INVALID_TIMESTAMP;
};

} // namespace dunedaq::trgdataformats::exp

#endif // DUNETRIGGER_EXPERIMENTAL_TPCTRIGGEROBJS_TRIGGERACTIVITYDATA_HPP_
