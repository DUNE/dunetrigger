#include "dunetrigger/TriggerSim/WireCell/TPGFilter.h"
#include "dunetrigger/TriggerSim/WireCell/TPGSaver.h"

#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/IChannel.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <cmath>
#include <unordered_set>

WIRECELL_FACTORY(wclsTPGFilter,
                 dunetrigger::TPGFilter,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

namespace {

  // Convert a JSON object (as written in Jsonnet) to a FHiCL ParameterSet so
  // it can configure an art tool.
  fhicl::ParameterSet json_to_pset(const Json::Value& jobj)
  {
    fhicl::ParameterSet ps;
    for (auto const& key : jobj.getMemberNames()) {
      auto const& val = jobj[key];
      if (val.isObject()) { ps.put(key, json_to_pset(val)); }
      else if (val.isBool()) {
        ps.put(key, val.asBool());
      }
      else if (val.isIntegral()) {
        ps.put(key, static_cast<long long>(val.asInt64()));
      }
      else if (val.isDouble()) {
        ps.put(key, val.asDouble());
      }
      else if (val.isString()) {
        ps.put(key, val.asString());
      }
      else {
        THROW(ValueError() << errmsg{"wclsTPGFilter: unsupported value type for tpalg key \"" +
                                     key + "\""});
      }
    }
    return ps;
  }

}

dunetrigger::TPGFilter::TPGFilter() : Aux::Logger("wclsTPGFilter", "tpg") {}

dunetrigger::TPGFilter::~TPGFilter() {}

Configuration dunetrigger::TPGFilter::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = "AnodePlane";
  cfg["saver"] = "wclsTPGSaver";
  cfg["tpalg"] = Json::objectValue;
  cfg["detid"] = m_detid;
  cfg["default_timestamp"] = Json::UInt64(m_default_timestamp);
  cfg["nticks"] = m_nticks;
  cfg["scale"] = m_scale;
  cfg["pass_through"] = m_pass_through;
  return cfg;
}

void dunetrigger::TPGFilter::configure(const Configuration& cfg)
{
  // anode: channel -> plane index lookup, so the TP tool never needs the
  // LArSoft geometry service.
  m_anode = Factory::find_tn<IAnodePlane>(get<std::string>(cfg, "anode", "AnodePlane"));
  m_plane.clear();
  for (int ch : m_anode->channels()) {
    const int index = m_anode->channel(ch)->planeid().index();
    if (index < 0) {
      log->warn("channel {} of anode {} has no well defined plane, it will be skipped",
                ch, m_anode->ident());
      continue;
    }
    m_plane[ch] = index;
  }

  // the saver is shared by all TPG filters of the graph
  const std::string saver_tn = get<std::string>(cfg, "saver", "wclsTPGSaver");
  auto visitor = Factory::find_tn<wcls::IArtEventVisitor>(saver_tn);
  m_saver = std::dynamic_pointer_cast<TPGSaver>(visitor);
  if (!m_saver) {
    THROW(ValueError() << errmsg{"wclsTPGFilter: \"" + saver_tn + "\" is not a wclsTPGSaver"});
  }

  // the TP algorithm: any TPAlgTPCTool, configured as in FHiCL
  if (!cfg["tpalg"].isObject() || !cfg["tpalg"].isMember("tool_type")) {
    THROW(ValueError() << errmsg{"wclsTPGFilter: \"tpalg\" must be an object with a \"tool_type\""});
  }
  m_tpalg = art::make_tool<TPAlgTPCTool>(json_to_pset(cfg["tpalg"]));

  m_detid = get<int>(cfg, "detid", m_detid);
  m_default_timestamp = cfg["default_timestamp"].isNull() ?
                          m_default_timestamp :
                          cfg["default_timestamp"].asUInt64();
  m_nticks = get<int>(cfg, "nticks", m_nticks);
  m_scale = get<double>(cfg, "scale", m_scale);
  m_pass_through = get<bool>(cfg, "pass_through", m_pass_through);

  log->debug("anode {}: {} channels, tpalg {}, detid {}, pass_through {}",
             m_anode->ident(), m_plane.size(), cfg["tpalg"]["tool_type"].asString(),
             m_detid, m_pass_through);
}

bool dunetrigger::TPGFilter::operator()(const input_pointer& in, output_pointer& out)
{
  out = nullptr;
  if (!in) { // end of stream
    log->debug("EOS at call {}", m_count);
    return true;
  }

  std::vector<TPGSaver::tp_t> tps;
  size_t nchannels = 0;

  auto traces = in->traces();
  if (traces) {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::unordered_set<int> seen;
    std::vector<short> adcs;
    for (auto const& trace : *traces) {
      const int ch = trace->channel();
      // like wclsFrameSaver, use only the first trace of a channel
      if (!seen.insert(ch).second) continue;

      auto pit = m_plane.find(ch);
      if (pit == m_plane.end()) {
        log->warn("channel {} is not in anode {}, skipping", ch, m_anode->ident());
        continue;
      }

      // build the waveform exactly as wclsFrameSaver::save_as_raw does
      const auto& charge = trace->charge();
      const int tbin = trace->tbin();
      size_t ncharge = charge.size();
      int nticks = tbin + ncharge;
      if (m_nticks) {
        if (m_nticks < nticks) { ncharge = m_nticks > tbin ? m_nticks - tbin : 0; }
        nticks = m_nticks;
      }
      adcs.assign(nticks, 0);
      for (size_t ind = 0; ind < ncharge; ++ind) {
        if (std::isfinite(charge[ind])) { adcs[tbin + ind] = m_scale * charge[ind]; }
      }

      m_tpalg->process_waveform(adcs, ch, pit->second, m_detid, m_default_timestamp, tps);
      ++nchannels;
    }
  }

  log->debug("call={} frame={} anode={} channels={} TPs={}",
             m_count, in->ident(), m_anode->ident(), nchannels, tps.size());
  m_saver->add(std::move(tps));

  if (m_pass_through) { out = in; }
  else {
    out = std::make_shared<Aux::SimpleFrame>(in->ident(), in->time(), ITrace::vector{}, in->tick());
  }
  ++m_count;
  return true;
}
