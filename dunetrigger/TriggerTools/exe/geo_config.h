// Shared geometry/wire-readout configuration for the TriggerTools standalone executables.
#ifndef TRIGGERTOOLS_EXE_GEO_CONFIG_H
#define TRIGGERTOOLS_EXE_GEO_CONFIG_H

#include <map>
#include <sstream>
#include <string>

namespace trgtools {

inline std::map<std::string, std::string> const& geo_shortcuts()
{
  static std::map<std::string, std::string> const shortcuts = {
    {"vd_1x8x6", "dunevd10kt_1x8x6_3view_30deg_geo"},
    {"vd_1x8x14", "dunevd10kt_1x8x14_3view_30deg_geo"},
    {"hd_1x2x6", "dune10kt_1x2x6_geo"},
    {"hd_10kt", "dune10kt_geo"},
    {"vd_10kt", "dunevd10kt_geo"},
    {"protodunevd", "protodunevd_v5_geo"}};
  return shortcuts;
}

// Shortcut -> geometry table of geometry_dune.fcl; unknown names are passed through unchanged.
inline std::string resolve_geo(std::string const& name)
{
  auto const& shortcuts = geo_shortcuts();
  auto const it = shortcuts.find(name);
  return it == shortcuts.end() ? name : it->second;
}

// The wire readout class is selected by detector name in DUNEWireReadout_service; the CRP-based
// geometries (protoDUNE-VD, full dunevd10kt) need a table that also configures the PDMapTool.
inline std::string default_wire_readout(std::string const& geo_id)
{
  auto const has = [&](char const* s) { return geo_id.find(s) != std::string::npos; };
  if (has("protodunevd")) return "protodunevd_wire_readout";
  if (has("dunevd10kt") && !has("1x8x") && !has("1x6x")) return "dunevd10kt_wire_readout";
  if (has("dune35t")) return "dune35t_wire_readout";
  return "dune_wire_readout";
}

// Multi-line list of the shortcuts and their targets, for usage messages.
inline std::string shortcut_help()
{
  std::ostringstream os;
  os << "Geometry shortcuts:\n";
  for (auto const& [name, target] : geo_shortcuts()) {
    os << "  " << name << " -> " << target << "\n";
  }
  return os.str();
}

// FHiCL configuration for ArtServiceHelper::load_services.
inline std::string service_config(std::string const& geo_id, std::string const& wire_readout)
{
  std::ostringstream os;
  os << "#include \"geometry_dune.fcl\"\n"
     << "services.Geometry:     @local::" << geo_id << "\n"
     << "services.WireReadout:      @local::" << wire_readout << "\n";
  return os.str();
}

} // namespace trgtools

#endif
