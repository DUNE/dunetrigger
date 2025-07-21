/**
 * @file TriggerCandidate.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRIGGERALGS_INCLUDE_TRIGGERALGS_TRIGGERCANDIDATE_HPP_
#define TRIGGERALGS_INCLUDE_TRIGGERALGS_TRIGGERCANDIDATE_HPP_

#include "detdataformats/trigger/TriggerActivityData2.hpp"
#include "detdataformats/trigger/TriggerCandidateData2.hpp"

#include <vector>

namespace triggeralgs {

struct TriggerCandidate : public dunedaq::trgdataformats2::TriggerCandidateData
{
  TriggerCandidate() = default;
  TriggerCandidate(const TriggerCandidate&) = default;
  TriggerCandidate& operator=(const TriggerCandidate&) = default;
  TriggerCandidate& operator=(TriggerCandidate&&) = default;
  ~TriggerCandidate() = default;

  TriggerCandidate(dunedaq::trgdataformats2::TriggerCandidateData&& data)
      : dunedaq::trgdataformats2::TriggerCandidateData(std::move(data)) {}
  TriggerCandidate(const dunedaq::trgdataformats2::TriggerCandidateData &data)
      : dunedaq::trgdataformats2::TriggerCandidateData(data) {}

  std::vector<dunedaq::trgdataformats2::TriggerActivityData> inputs;
};

} // namespace triggeralgs

#endif // TRIGGERALGS_INCLUDE_TRIGGERALGS_TRIGGERCANDIDATE_HPP_
