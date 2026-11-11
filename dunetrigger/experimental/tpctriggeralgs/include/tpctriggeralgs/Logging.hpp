/**
 * @file Logging.hpp Common logging declarations in tpctriggeralgs
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef tpctriggeralgs_INCLUDE_LOGGING_HPP_
#define tpctriggeralgs_INCLUDE_LOGGING_HPP_

namespace tpctriggeralgs {

/**
* @brief Common name used by TRACE TLOG calls from this package
*/
enum Logging
{
  TLVL_VERY_IMPORTANT = 1,
  TLVL_IMPORTANT      = 2,
  TLVL_GENERAL        = 3,
  TLVL_DEBUG_INFO     = 4,
  TLVL_DEBUG_LOW      = 5,
  TLVL_DEBUG_MEDIUM   = 10,
  TLVL_DEBUG_HIGH     = 15,
  TLVL_DEBUG_ALL      = 20
};

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_LOGGING_HPP_
