#include "dunetrigger/TriggerSim/WireCell/TPGSaver.h"

#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"

#include "WireCellUtil/NamedFactory.h"

#include <algorithm>
#include <memory>
#include <tuple>

WIRECELL_FACTORY(wclsTPGSaver,
                 dunetrigger::TPGSaver,
                 wcls::IArtEventVisitor,
                 WireCell::IConfigurable)

using namespace WireCell;

dunetrigger::TPGSaver::TPGSaver() : Aux::Logger("wclsTPGSaver", "tpg") {}

dunetrigger::TPGSaver::~TPGSaver() {}

Configuration dunetrigger::TPGSaver::default_configuration() const
{
  Configuration cfg;
  cfg["instance"] = m_instance;
  cfg["sort"] = m_sort;
  return cfg;
}

void dunetrigger::TPGSaver::configure(const Configuration& cfg)
{
  m_instance = get<std::string>(cfg, "instance", m_instance);
  m_sort = get<bool>(cfg, "sort", m_sort);
}

void dunetrigger::TPGSaver::produces(art::ProducesCollector& collector)
{
  collector.produces<std::vector<tp_t>>(m_instance);
}

void dunetrigger::TPGSaver::add(std::vector<tp_t>&& tps)
{
  std::lock_guard<std::mutex> lock(m_mutex);
  if (m_tps.empty()) {
    m_tps = std::move(tps);
    return;
  }
  m_tps.insert(m_tps.end(),
               std::make_move_iterator(tps.begin()),
               std::make_move_iterator(tps.end()));
}

void dunetrigger::TPGSaver::visit(art::Event& event)
{
  auto out = std::make_unique<std::vector<tp_t>>();
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    out->swap(m_tps);
  }

  if (m_sort) {
    std::sort(out->begin(), out->end(), [](const tp_t& a, const tp_t& b) {
      return std::tie(a.time_start, a.channel, a.time_peak, a.time_over_threshold, a.adc_integral) <
             std::tie(b.time_start, b.channel, b.time_peak, b.time_over_threshold, b.adc_integral);
    });
  }

  log->debug("saving {} TPs to art event {} (instance \"{}\")",
             out->size(), event.event(), m_instance);
  event.put(std::move(out), m_instance);
}
