/**
 * @file Types.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_TYPES_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_TYPES_HPP_

#include "detdataformats/trigger/Types.hpp"

namespace tpctriggeralgs {

using timestamp_t = dunedaq::trgdataformats::timestamp_t;
using timestamp_diff_t = dunedaq::trgdataformats::timestamp_diff_t;
using detid_t = dunedaq::trgdataformats::detid_t;
using trigger_number_t = dunedaq::trgdataformats::trigger_number_t;
using channel_t = dunedaq::trgdataformats::channel_t;
using version_t = dunedaq::trgdataformats::version_t;

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_TYPES_HPP_
