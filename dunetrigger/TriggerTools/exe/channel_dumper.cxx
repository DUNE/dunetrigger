// channel_dumper: dump the channel <-> wire mapping of a TPC (or plane) or of a
// TPC set (or readout plane, ROP) for a given geometry.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dunecore/ArtSupport/ArtServiceHelper.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "geo_config.h"
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

namespace po = boost::program_options;
using json = nlohmann::json;

namespace {

enum class Mode { None, Tpc, TpcSet };

struct Options {
  std::string geo;
  int cryo = 0;
  Mode mode = Mode::None;
  int tpc = -1;
  int tpcset = -1;
  std::optional<int> plane;
  std::optional<int> rop;
  std::string wire_readout;
  bool table = false;
  std::string format = "none";
  std::string prefix;
};

struct ChannelRecord {
  std::vector<geo::WireID> wires;
  std::optional<unsigned int> rop;
};

using ChannelMap = std::map<raw::ChannelID_t, ChannelRecord>;

std::optional<int> optional_int(po::variables_map const& vm, char const* name)
{
  if (vm.count(name) == 0) return std::nullopt;
  int const v = vm[name].as<int>();
  if (v < 0) throw po::error(std::string("--") + name + " must be non-negative");
  return v;
}

// Fills opt from the command line; throws po::error on invalid usage.
// Returns false if help was requested.
bool parse_args(int argc, char** argv, po::options_description const& desc, Options& opt)
{
  po::options_description hidden;
  hidden.add_options()("geo", po::value<std::string>(), "geometry");
  po::positional_options_description positional;
  positional.add("geo", 1);

  po::options_description all;
  all.add(desc).add(hidden);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);

  if (vm.count("help")) return false;
  po::notify(vm);

  if (vm.count("geo") == 0) throw po::error("missing <geo> argument");
  opt.geo = vm["geo"].as<std::string>();
  opt.wire_readout = vm["wire-readout"].as<std::string>();
  opt.table = vm["table"].as<bool>();
  opt.format = vm["format"].as<std::string>();
  opt.prefix = vm["output"].as<std::string>();
  opt.cryo = vm["cryo"].as<int>();
  if (opt.cryo < 0) throw po::error("--cryo must be non-negative");

  auto const tpc = optional_int(vm, "tpc");
  auto const tpcset = optional_int(vm, "tpcset");
  opt.plane = optional_int(vm, "plane");
  opt.rop = optional_int(vm, "rop");

  if (opt.plane && !tpc) throw po::error("--plane requires --tpc");
  if (opt.rop && !tpcset) throw po::error("--rop requires --tpcset");
  if (tpc.has_value() == tpcset.has_value()) {
    throw po::error("exactly one of --tpc / --tpcset is required");
  }

  if (tpc) {
    opt.mode = Mode::Tpc;
    opt.tpc = *tpc;
  }
  else {
    opt.mode = Mode::TpcSet;
    opt.tpcset = *tpcset;
  }
  if (opt.format != "none" && opt.format != "csv" && opt.format != "json" &&
      opt.format != "both") {
    throw po::error("invalid --format '" + opt.format + "' (none|csv|json|both)");
  }
  return true;
}

po::options_description make_options()
{
  po::options_description desc("Options");
  desc.add_options()
    ("cryo", po::value<int>()->default_value(0), "cryostat ID")
    ("tpc", po::value<int>(), "dump the channels of this TPC")
    ("plane", po::value<int>(), "restrict --tpc to this plane")
    ("tpcset", po::value<int>(), "dump the channels of this TPC set")
    ("rop", po::value<int>(), "restrict --tpcset to this readout plane")
    ("format", po::value<std::string>()->default_value("none"),
     "output files: none (stdout summary only), csv, json, both")
    ("table", po::bool_switch(), "print a channel vs z0 table to stdout after the summary")
    ("wire-readout", po::value<std::string>()->default_value(""),
     "wire readout fcl table from geometry_dune.fcl (default: chosen from the geometry)")
    ("output,o", po::value<std::string>()->default_value(""),
     "output file prefix (default <geoid>_tpc<T>[_plane<P>] or <geoid>_tpcset<S>[_rop<R>])")
    ("help,h", "show this help");
  return desc;
}

void print_usage(std::ostream& os, po::options_description const& desc)
{
  os << "Usage: channel_dumper <geo> (--tpc T [--plane P] | --tpcset S [--rop R]) [options]\n"
     << "  <geo>: shortcut or geometry_dune.fcl table\n\n"
     << trgtools::shortcut_help() << "\n"
     << desc << "\n";
}

std::string default_prefix(Options const& opt, std::string const& geo_id)
{
  std::ostringstream os;
  if (opt.mode == Mode::Tpc) {
    os << geo_id << "_tpc" << opt.tpc;
    if (opt.plane) os << "_plane" << *opt.plane;
  }
  else {
    os << geo_id << "_tpcset" << opt.tpcset;
    if (opt.rop) os << "_rop" << *opt.rop;
  }
  return os.str();
}

std::string wire_string(geo::WireID const& w)
{
  std::ostringstream os;
  os << w.Cryostat << ':' << w.TPC << ':' << w.Plane << ':' << w.Wire;
  return os.str();
}

std::string view_name(geo::WireReadoutGeom const& readout, raw::ChannelID_t ch)
{
  return geo::PlaneGeo::ViewName(readout.View(ch));
}

// Validates the selection; returns false (message printed) if it does not exist.
bool validate(Options const& opt, geo::GeometryCore const& geom, geo::WireReadoutGeom const& readout)
{
  if (!geom.HasCryostat(geo::CryostatID{static_cast<geo::CryostatID::CryostatID_t>(opt.cryo)})) {
    std::cerr << "Error: cryostat " << opt.cryo << " does not exist\n";
    return false;
  }
  geo::CryostatID const cryo{static_cast<geo::CryostatID::CryostatID_t>(opt.cryo)};

  if (opt.mode == Mode::Tpc) {
    geo::TPCID const tpcid{cryo, static_cast<geo::TPCID::TPCID_t>(opt.tpc)};
    if (!geom.HasTPC(tpcid)) {
      std::cerr << "Error: TPC " << opt.tpc << " does not exist in cryostat " << opt.cryo << "\n";
      return false;
    }
    if (opt.plane &&
        !readout.HasPlane(geo::PlaneID{tpcid, static_cast<geo::PlaneID::PlaneID_t>(*opt.plane)})) {
      std::cerr << "Error: plane " << *opt.plane << " does not exist in TPC " << opt.tpc << "\n";
      return false;
    }
    return true;
  }

  readout::TPCsetID const setid{cryo, static_cast<readout::TPCsetID::TPCsetID_t>(opt.tpcset)};
  if (!readout.HasTPCset(setid)) {
    std::cerr << "Error: TPC set " << opt.tpcset << " does not exist in cryostat " << opt.cryo
              << "\n";
    return false;
  }
  if (opt.rop &&
      !readout.HasROP(readout::ROPID{setid, static_cast<readout::ROPID::ROPID_t>(*opt.rop)})) {
    std::cerr << "Error: ROP " << *opt.rop << " does not exist in TPC set " << opt.tpcset
              << " (it has " << readout.NROPs(setid) << " ROPs)\n";
    return false;
  }
  return true;
}

ChannelMap collect_tpc(Options const& opt, geo::WireReadoutGeom const& readout)
{
  geo::TPCID const tpcid{geo::CryostatID{static_cast<geo::CryostatID::CryostatID_t>(opt.cryo)},
                         static_cast<geo::TPCID::TPCID_t>(opt.tpc)};
  ChannelMap channels;
  auto add_wires = [&](auto const& id) {
    for (geo::WireID const& wire : readout.Iterate<geo::WireID>(id)) {
      channels[readout.PlaneWireToChannel(wire)].wires.push_back(wire);
    }
  };
  if (opt.plane) {
    add_wires(geo::PlaneID{tpcid, static_cast<geo::PlaneID::PlaneID_t>(*opt.plane)});
  }
  else {
    add_wires(tpcid);
  }
  return channels;
}

void collect_rop(readout::ROPID const& rop, geo::WireReadoutGeom const& readout, ChannelMap& channels)
{
  raw::ChannelID_t const first = readout.FirstChannelInROP(rop);
  for (raw::ChannelID_t ch = first; ch < first + readout.Nchannels(rop); ++ch) {
    ChannelRecord& rec = channels[ch];
    rec.wires = readout.ChannelToWire(ch);
    rec.rop = rop.ROP;
  }
}

ChannelMap collect_tpcset(Options const& opt, geo::WireReadoutGeom const& readout)
{
  readout::TPCsetID const setid{geo::CryostatID{static_cast<geo::CryostatID::CryostatID_t>(opt.cryo)},
                                static_cast<readout::TPCsetID::TPCsetID_t>(opt.tpcset)};
  ChannelMap channels;
  if (opt.rop) {
    collect_rop(readout::ROPID{setid, static_cast<readout::ROPID::ROPID_t>(*opt.rop)}, readout, channels);
  }
  else {
    for (unsigned int r = 0; r < readout.NROPs(setid); ++r) {
      collect_rop(readout::ROPID{setid, static_cast<readout::ROPID::ROPID_t>(r)}, readout, channels);
    }
  }
  return channels;
}

// Channels (and the z0 of their wires) belonging to one plane (TPC mode) or ROP (ROP mode).
struct Group {
  int id = 0;
  std::string view;
  std::set<raw::ChannelID_t> channels;
  double z0_min = std::numeric_limits<double>::infinity();
  double z0_max = -std::numeric_limits<double>::infinity();
  size_t n_non_finite_z0 = 0;
  double dz_min = std::numeric_limits<double>::infinity();
  double dz_max = -std::numeric_limits<double>::infinity();
  std::vector<int> tpcs; // TPCs read out by the ROP (ROP mode only)
  std::map<geo::PlaneID, double> pitches; // wire pitch (cm) of each plane the group's wires use

  double pitch_min() const { return std::min_element(pitches.begin(), pitches.end(), by_pitch)->second; }
  double pitch_max() const { return std::max_element(pitches.begin(), pitches.end(), by_pitch)->second; }
  bool dz_uniform() const { return dz_max - dz_min < 1e-4; }
  bool pitch_uniform() const { return pitch_max() - pitch_min() < 1e-4; }

  static bool by_pitch(std::pair<geo::PlaneID const, double> const& a,
                       std::pair<geo::PlaneID const, double> const& b)
  {
    return a.second < b.second;
  }

  void add_dz(double dz)
  {
    if (!std::isfinite(dz)) return;
    dz_min = std::min(dz_min, dz);
    dz_max = std::max(dz_max, dz);
  }

  void add_z0(double z0)
  {
    if (!std::isfinite(z0)) {
      ++n_non_finite_z0;
      return;
    }
    z0_min = std::min(z0_min, z0);
    z0_max = std::max(z0_max, z0);
  }
};

struct Summary {
  size_t n_channels = 0;
  size_t n_wires = 0;
  size_t n_multi_wire = 0;
  raw::ChannelID_t min_channel = 0;
  raw::ChannelID_t max_channel = 0;
  std::map<std::string, size_t> per_view;
  std::vector<Group> groups;
  std::vector<int> set_tpcs; // TPCs of the selected TPC set (ROP mode only)
};

// z0 (cm) of each wire of the channel, in the order of the channel's wire list.
std::vector<double> channel_z0s(ChannelRecord const& rec, geo::WireReadoutGeom const& readout)
{
  std::vector<double> z0s;
  for (geo::WireID const& wire : rec.wires) z0s.push_back(readout.Wire(wire).ComputeZatY0());
  return z0s;
}

bool wire_in_selected_tpc(Options const& opt, geo::WireID const& wire)
{
  return static_cast<int>(wire.Cryostat) == opt.cryo && static_cast<int>(wire.TPC) == opt.tpc;
}

// Wire spacing along z: pitch of the wire's plane over |sin(ThetaZ)|; nan for wires parallel to z.
double wire_dz(geo::WireReadoutGeom const& readout, geo::WireID const& wire)
{
  double const sin_theta = std::abs(readout.Wire(wire).SinThetaZ());
  if (sin_theta < 1e-6) return std::nan("");
  return readout.Plane(geo::PlaneID(wire)).WirePitch() / sin_theta;
}

int group_id(Options const& opt, ChannelRecord const& rec, geo::WireID const& wire)
{
  return opt.mode == Mode::Tpc ? static_cast<int>(wire.Plane) : static_cast<int>(*rec.rop);
}

// TPC mode groups by the plane of each wire (wrapped wires of other TPCs still count towards the
// channel); ROP mode groups by the ROP recorded for the channel.
std::vector<Group> group_channels(Options const& opt, ChannelMap const& channels,
                                  geo::WireReadoutGeom const& readout)
{
  std::map<int, Group> groups;
  for (auto const& [ch, rec] : channels) {
    auto const z0s = channel_z0s(rec, readout);
    for (size_t i = 0; i < rec.wires.size(); ++i) {
      int const id = group_id(opt, rec, rec.wires[i]);
      Group& g = groups[id];
      g.id = id;
      g.view = view_name(readout, ch);
      g.channels.insert(ch);
      g.add_z0(z0s[i]);
      g.add_dz(wire_dz(readout, rec.wires[i]));
      geo::PlaneID const plane = rec.wires[i];
      bool const in_selected_tpc = wire_in_selected_tpc(opt, rec.wires[i]);
      if (opt.mode == Mode::TpcSet || in_selected_tpc) {
        g.pitches[plane] = readout.Plane(plane).WirePitch();
      }
    }
  }
  std::vector<Group> out;
  for (auto& [id, g] : groups) out.push_back(std::move(g));
  return out;
}

std::string fixed_string(double v, int precision)
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(precision) << v;
  return os.str();
}

std::string range_string(double min, double max, bool uniform)
{
  return uniform ? fixed_string(min, 4) : fixed_string(min, 4) + "-" + fixed_string(max, 4);
}

std::string pitch_string(Group const& g)
{
  return range_string(g.pitch_min(), g.pitch_max(), g.pitch_uniform());
}

std::string dz_string(Group const& g)
{
  return range_string(g.dz_min, g.dz_max, g.dz_uniform());
}

std::string group_label(Options const& opt)
{
  return opt.mode == Mode::Tpc ? "plane" : "rop";
}

std::vector<int> tpc_numbers(std::vector<geo::TPCID> const& ids)
{
  std::vector<int> out;
  for (geo::TPCID const& id : ids) out.push_back(static_cast<int>(id.TPC));
  return out;
}

std::string join_ints(std::vector<int> const& v, std::string const& sep)
{
  std::string out;
  for (size_t i = 0; i < v.size(); ++i) out += (i ? sep : "") + std::to_string(v[i]);
  return out;
}

Summary summarize(Options const& opt, ChannelMap const& channels, geo::WireReadoutGeom const& readout)
{
  Summary s;
  s.groups = group_channels(opt, channels, readout);
  if (opt.mode == Mode::TpcSet) {
    geo::CryostatID const cryo{static_cast<geo::CryostatID::CryostatID_t>(opt.cryo)};
    readout::TPCsetID const setid{cryo, static_cast<readout::TPCsetID::TPCsetID_t>(opt.tpcset)};
    s.set_tpcs = tpc_numbers(readout.TPCsetToTPCs(setid));
    for (Group& g : s.groups) {
      g.tpcs = tpc_numbers(readout.ROPtoTPCs(readout::ROPID{setid, static_cast<readout::ROPID::ROPID_t>(g.id)}));
    }
  }
  s.n_channels = channels.size();
  if (channels.empty()) return s;
  s.min_channel = channels.begin()->first;
  s.max_channel = channels.rbegin()->first;
  for (auto const& [ch, rec] : channels) {
    s.n_wires += rec.wires.size();
    if (rec.wires.size() > 1) ++s.n_multi_wire;
    ++s.per_view[view_name(readout, ch)];
  }
  return s;
}

std::string selection_string(Options const& opt)
{
  std::ostringstream os;
  os << "cryo " << opt.cryo;
  if (opt.mode == Mode::Tpc) {
    os << ", TPC " << opt.tpc;
    os << (opt.plane ? ", plane " + std::to_string(*opt.plane) : ", all planes");
  }
  else {
    os << ", TPC set " << opt.tpcset;
    os << (opt.rop ? ", ROP " + std::to_string(*opt.rop) : ", all ROPs");
  }
  return os.str();
}

// Plain-text table with column widths computed from the header and the cells. Columns are
// right-aligned by default. The last column can instead be left-aligned and unpadded, either sized
// by its content (Left) or with a rule as wide as its header only (Free, for long free-form text).
// A row may carry a note, printed indented below it.
enum class LastColumn { Right, Left, Free };

class TextTable {
public:
  explicit TextTable(std::vector<std::string> header, LastColumn last = LastColumn::Right)
    : header_(std::move(header)), last_(last)
  {}

  void add_row(std::vector<std::string> cells, std::string note = "")
  {
    rows_.push_back({std::move(cells), std::move(note)});
  }

  void print(std::ostream& os, std::string const& indent = "  ") const
  {
    auto const widths = column_widths();
    print_cells(os, indent, header_, widths);
    std::vector<std::string> rule;
    for (size_t i = 0; i < widths.size(); ++i) rule.push_back(std::string(widths[i], '-'));
    print_cells(os, indent, rule, widths);
    for (auto const& row : rows_) {
      print_cells(os, indent, row.cells, widths);
      if (!row.note.empty()) os << indent << "  " << row.note << "\n";
    }
  }

private:
  struct Row {
    std::vector<std::string> cells;
    std::string note;
  };

  std::vector<size_t> column_widths() const
  {
    std::vector<size_t> widths;
    for (auto const& h : header_) widths.push_back(h.size());
    for (auto const& row : rows_) {
      for (size_t i = 0; i < widths.size(); ++i) widths[i] = std::max(widths[i], row.cells[i].size());
    }
    // A free-form last column keeps its header width so the rule does not follow the longest row.
    if (last_ == LastColumn::Free) widths.back() = header_.back().size();
    return widths;
  }

  void print_cells(std::ostream& os, std::string const& indent, std::vector<std::string> const& cells,
                   std::vector<size_t> const& widths) const
  {
    os << indent;
    for (size_t i = 0; i < cells.size(); ++i) {
      bool const last = i + 1 == cells.size();
      if (i) os << "  ";
      if (last && last_ != LastColumn::Right) {
        os << cells[i];
      }
      else {
        os << std::setw(widths[i]) << cells[i];
      }
    }
    os << "\n";
  }

  std::vector<std::string> header_;
  LastColumn last_;
  std::vector<Row> rows_;
};

std::string pitch_note(Group const& g)
{
  std::ostringstream os;
  os << "pitch by plane:";
  for (auto const& [plane, pitch] : g.pitches) {
    os << ' ' << plane.Cryostat << ':' << plane.TPC << ':' << plane.Plane << '=' << fixed_string(pitch, 4);
  }
  return os.str();
}

void print_group_table(Options const& opt, std::vector<Group> const& groups)
{
  bool const any_non_finite =
    std::any_of(groups.begin(), groups.end(), [](Group const& g) { return g.n_non_finite_z0 > 0; });
  std::vector<std::string> header = {group_label(opt), "view", "n_ch", "ch_min", "ch_max",
                                     "pitch [cm]", "dz [cm]", "z0_min [cm]", "z0_max [cm]"};
  if (any_non_finite) header.push_back("n_nonfinite");
  bool const with_tpcs = opt.mode == Mode::TpcSet;
  if (with_tpcs) header.push_back("tpcs");
  TextTable table(header, with_tpcs ? LastColumn::Left : LastColumn::Right);
  for (Group const& g : groups) {
    std::vector<std::string> row = {std::to_string(g.id),
                                    g.view,
                                    std::to_string(g.channels.size()),
                                    std::to_string(*g.channels.begin()),
                                    std::to_string(*g.channels.rbegin()),
                                    pitch_string(g),
                                    dz_string(g),
                                    fixed_string(g.z0_min, 3),
                                    fixed_string(g.z0_max, 3)};
    if (any_non_finite) row.push_back(std::to_string(g.n_non_finite_z0));
    if (with_tpcs) row.push_back(join_ints(g.tpcs, ","));
    table.add_row(std::move(row), g.pitch_uniform() ? "" : pitch_note(g));
  }
  table.print(std::cout);
}

void print_summary(std::string const& detector, Options const& opt, Summary const& s)
{
  std::cout << "Detector        : " << detector << "\n"
            << "Selection       : " << selection_string(opt) << "\n";
  if (opt.mode == Mode::TpcSet) {
    std::cout << "TPCs in set     : " << join_ints(s.set_tpcs, ", ") << "\n";
  }
  std::cout << "Channels        : " << s.n_channels << "\n"
            << "Wires           : " << s.n_wires << "\n"
            << "Channel range   : [" << s.min_channel << ", " << s.max_channel << "]\n"
            << "Channels / view :";
  for (auto const& [view, n] : s.per_view) std::cout << " " << view << "=" << n;
  std::cout << "\nChannels with >1 wire: " << s.n_multi_wire << "\n";
  print_group_table(opt, s.groups);
  std::cout << std::flush;
}

// Keeps, for each channel, only the wires in the given TPC; channels without any are dropped.
ChannelMap restrict_to_tpc(ChannelMap const& channels, int cryo, int tpc)
{
  ChannelMap out;
  for (auto const& [ch, rec] : channels) {
    for (geo::WireID const& wire : rec.wires) {
      if (static_cast<int>(wire.Cryostat) == cryo && static_cast<int>(wire.TPC) == tpc) {
        ChannelRecord& r = out[ch];
        r.rop = rec.rop;
        r.wires.push_back(wire);
      }
    }
  }
  return out;
}

// For a TPC set selection, prints the TPC-mode summary of each TPC of the set, restricted to the
// selected channels.
void print_tpc_summaries(std::string const& detector, Options const& opt, ChannelMap const& channels,
                         Summary const& set_summary, geo::WireReadoutGeom const& readout)
{
  for (int const tpc : set_summary.set_tpcs) {
    ChannelMap const tpc_channels = restrict_to_tpc(channels, opt.cryo, tpc);
    if (tpc_channels.empty()) continue;
    Options tpc_opt = opt;
    tpc_opt.mode = Mode::Tpc;
    tpc_opt.tpc = tpc;
    tpc_opt.plane.reset();
    tpc_opt.rop.reset();
    std::cout << "\n--- TPC " << tpc << " of TPC set " << opt.tpcset << " ---\n";
    print_summary(detector, tpc_opt, summarize(tpc_opt, tpc_channels, readout));
  }
}

json selection_json(Options const& opt, std::vector<int> const& set_tpcs)
{
  json sel = {{"cryo", opt.cryo}};
  if (opt.mode == Mode::Tpc) {
    sel["mode"] = "tpc";
    sel["tpc"] = opt.tpc;
    if (opt.plane) sel["plane"] = *opt.plane;
  }
  else {
    sel["mode"] = "tpcset";
    sel["tpcset"] = opt.tpcset;
    sel["tpcs"] = set_tpcs;
    if (opt.rop) sel["rop"] = *opt.rop;
  }
  return sel;
}

Group const* find_group(std::vector<Group> const& groups, int id)
{
  auto const it = std::find_if(groups.begin(), groups.end(), [&](Group const& g) { return g.id == id; });
  return it == groups.end() ? nullptr : &*it;
}

// i_z = (z0 - group min z0) / dz: the wire index along z, in units of the wire spacing along z.
double wire_iz(Group const& g, double dz, double z0)
{
  return (z0 - g.z0_min) / dz;
}

// i_z rounded to the nearest integer; empty if it is not finite.
std::optional<long> rounded_iz(double iz)
{
  if (!std::isfinite(iz)) return std::nullopt;
  return std::lround(iz);
}

// channel + sign * i_z; empty if i_z is.
std::optional<long> channel_offset(raw::ChannelID_t ch, std::optional<long> iz, int sign)
{
  if (!iz) return std::nullopt;
  return static_cast<long>(ch) + sign * *iz;
}

// Per-wire quantities of a channel, shared by the table and the CSV/JSON files.
struct WireInfo {
  geo::WireID id;
  double z0;                    // z at y = 0 (cm)
  std::optional<long> iz;       // rounded i_z; empty if not finite
  std::optional<long> ch_minus; // channel - i_z
  std::optional<long> ch_plus;  // channel + i_z
  geo::Point_t center;          // wire center (cm)
};

struct ChannelInfo {
  raw::ChannelID_t channel;
  std::string view;
  int group; // plane (TPC mode) or ROP (ROP mode) of the channel's first wire
  std::vector<WireInfo> wires;
};

std::vector<ChannelInfo> channel_infos(Options const& opt, ChannelMap const& channels,
                                       std::vector<Group> const& groups,
                                       geo::WireReadoutGeom const& readout)
{
  std::vector<ChannelInfo> out;
  for (auto const& [ch, rec] : channels) {
    ChannelInfo info{ch, view_name(readout, ch), group_id(opt, rec, rec.wires.front()), {}};
    auto const z0s = channel_z0s(rec, readout);
    for (size_t i = 0; i < rec.wires.size(); ++i) {
      geo::WireID const& wire = rec.wires[i];
      Group const* g = find_group(groups, group_id(opt, rec, wire));
      auto const iz = rounded_iz(g ? wire_iz(*g, wire_dz(readout, wire), z0s[i]) : std::nan(""));
      info.wires.push_back({wire, z0s[i], iz, channel_offset(ch, iz, -1), channel_offset(ch, iz, +1),
                            readout.Wire(wire).GetCenter()});
    }
    out.push_back(std::move(info));
  }
  return out;
}

std::string int_string(std::optional<long> v)
{
  return v ? std::to_string(*v) : "nan";
}

// "(y,z)" of a wire center.
std::string center_string(geo::Point_t const& c)
{
  return "(" + fixed_string(c.Y(), 3) + "," + fixed_string(c.Z(), 3) + ")";
}

// Joins f(wire) over the channel's wires with sep.
template <typename F>
std::string join_wires(ChannelInfo const& info, std::string const& sep, F f)
{
  std::string out;
  for (size_t i = 0; i < info.wires.size(); ++i) out += (i ? sep : "") + f(info.wires[i]);
  return out;
}

std::string z0_of(WireInfo const& w) { return fixed_string(w.z0, 3); }
std::string iz_of(WireInfo const& w) { return int_string(w.iz); }
std::string ch_minus_of(WireInfo const& w) { return int_string(w.ch_minus); }
std::string ch_plus_of(WireInfo const& w) { return int_string(w.ch_plus); }
std::string center_of(WireInfo const& w) { return center_string(w.center); }
std::string center_y_of(WireInfo const& w) { return fixed_string(w.center.Y(), 3); }
std::string center_z_of(WireInfo const& w) { return fixed_string(w.center.Z(), 3); }
std::string id_of(WireInfo const& w) { return wire_string(w.id); }

// Single-wire channels: just the z0. Multi-wire channels:
// z0 list | i_z | ch-i_z | ch+i_z | center (y,z) lists (wire IDs).
std::string z0_cell(ChannelInfo const& info)
{
  std::string out = join_wires(info, ", ", z0_of);
  if (info.wires.size() > 1) {
    out += " | i_z " + join_wires(info, ", ", iz_of) + " | ch-i_z " + join_wires(info, ", ", ch_minus_of) +
           " | ch+i_z " + join_wires(info, ", ", ch_plus_of) + " | center (y,z) " +
           join_wires(info, ", ", center_of) + " (" + join_wires(info, ", ", id_of) + ')';
  }
  return out;
}

void print_table(Options const& opt, std::vector<ChannelInfo> const& infos)
{
  TextTable table({"channel", "view", group_label(opt), "n_wires", "i_z", "ch-i_z", "ch+i_z",
                   "center (y,z) [cm]", "z0 [cm]"},
                  LastColumn::Free);
  for (ChannelInfo const& info : infos) {
    bool const single = info.wires.size() == 1;
    WireInfo const& w = info.wires.front();
    table.add_row({std::to_string(info.channel), info.view, std::to_string(info.group),
                   std::to_string(info.wires.size()), single ? iz_of(w) : "-", single ? ch_minus_of(w) : "-",
                   single ? ch_plus_of(w) : "-", single ? center_of(w) : "-", z0_cell(info)});
  }
  std::cout << "\n";
  table.print(std::cout);
  std::cout << std::flush;
}

json groups_json(Options const& opt, std::vector<Group> const& groups)
{
  json out = json::array();
  for (Group const& g : groups) {
    out.push_back({{group_label(opt), g.id},
                   {"view", g.view},
                   {"n_channels", g.channels.size()},
                   {"channel_min", *g.channels.begin()},
                   {"channel_max", *g.channels.rbegin()},
                   {"z0_min", g.z0_min},
                   {"z0_max", g.z0_max},
                   {"n_non_finite_z0", g.n_non_finite_z0},
                   {"pitch", g.pitch_min()},
                   {"dz", g.dz_min}});
    if (opt.mode == Mode::TpcSet) {
      json& entry = out.back();
      entry["tpcs"] = g.tpcs;
      entry["pitch_min"] = g.pitch_min();
      entry["pitch_max"] = g.pitch_max();
      entry["planes"] = json::array();
      for (auto const& [plane, pitch] : g.pitches) {
        std::ostringstream id;
        id << plane.Cryostat << ':' << plane.TPC << ':' << plane.Plane;
        entry["planes"].push_back({{"plane_id", id.str()}, {"pitch", pitch}});
      }
    }
  }
  return out;
}

json optional_json(std::optional<long> v)
{
  return v ? json(*v) : json(nullptr);
}

// Per channel: the table's columns, with one list entry per wire in the channel's wire order.
json channel_json(Options const& opt, ChannelInfo const& info)
{
  json c = {{"channel", info.channel}, {"view", info.view}, {group_label(opt), info.group}};
  for (char const* key : {"wires", "z0", "i_z", "ch-i_z", "ch+i_z", "center_y", "center_z"}) {
    c[key] = json::array();
  }
  for (WireInfo const& w : info.wires) {
    c["wires"].push_back(wire_string(w.id));
    c["z0"].push_back(w.z0);
    c["i_z"].push_back(optional_json(w.iz));
    c["ch-i_z"].push_back(optional_json(w.ch_minus));
    c["ch+i_z"].push_back(optional_json(w.ch_plus));
    c["center_y"].push_back(w.center.Y());
    c["center_z"].push_back(w.center.Z());
  }
  return c;
}

void write_json(std::string const& path, std::string const& detector, Options const& opt,
                std::vector<ChannelInfo> const& infos, Summary const& s)
{
  json j;
  j["detector"] = detector;
  j["selection"] = selection_json(opt, s.set_tpcs);
  j["summary"] = {{"n_channels", s.n_channels},
                  {"n_wires", s.n_wires},
                  {"min_channel", s.min_channel},
                  {"max_channel", s.max_channel},
                  {"channels_per_view", s.per_view},
                  {"n_channels_multi_wire", s.n_multi_wire},
                  {"groups", groups_json(opt, s.groups)}};
  j["channels"] = json::array();
  for (ChannelInfo const& info : infos) j["channels"].push_back(channel_json(opt, info));
  std::ofstream(path) << std::setw(2) << j << std::endl;
  std::cout << "Wrote " << path << std::endl;
}

// The table's columns; per-wire values are joined by ';' in the channel's wire order.
void write_csv(std::string const& path, Options const& opt, std::vector<ChannelInfo> const& infos)
{
  std::ofstream out(path);
  out << "channel,view,n_wires,wires," << group_label(opt)
      << ",z0,i_z,ch-i_z,ch+i_z,center_y,center_z\n";
  for (ChannelInfo const& info : infos) {
    out << info.channel << ',' << info.view << ',' << info.wires.size() << ','
        << join_wires(info, ";", id_of) << ',' << info.group << ',' << join_wires(info, ";", z0_of) << ','
        << join_wires(info, ";", iz_of) << ',' << join_wires(info, ";", ch_minus_of) << ','
        << join_wires(info, ";", ch_plus_of) << ',' << join_wires(info, ";", center_y_of) << ','
        << join_wires(info, ";", center_z_of) << "\n";
  }
  std::cout << "Wrote " << path << std::endl;
}

} // namespace

int main(int argc, char** argv)
{
  Options opt;
  auto const desc = make_options();
  try {
    if (!parse_args(argc, argv, desc, opt)) {
      print_usage(std::cout, desc);
      return 0;
    }
  }
  catch (po::error const& e) {
    std::cerr << "Error: " << e.what() << "\n\n";
    print_usage(std::cerr, desc);
    return 1;
  }

  std::string const geo_id = trgtools::resolve_geo(opt.geo);

  std::string const readout_id =
    opt.wire_readout.empty() ? trgtools::default_wire_readout(geo_id) : opt.wire_readout;
  std::cout << "Geometry: " << geo_id << ", wire readout: " << readout_id << std::endl;

  std::stringstream config(trgtools::service_config(geo_id, readout_id));
  ArtServiceHelper::load_services(config);

  auto const& geom = *art::ServiceHandle<geo::Geometry>();
  auto const& readout = art::ServiceHandle<geo::WireReadout>()->Get();

  if (!validate(opt, geom, readout)) {
    return 2;
  }

  bool const tpc_mode = opt.mode == Mode::Tpc;
  ChannelMap const channels = tpc_mode ? collect_tpc(opt, readout) : collect_tpcset(opt, readout);
  Summary const summary = summarize(opt, channels, readout);
  print_summary(geom.DetectorName(), opt, summary);
  if (!tpc_mode) print_tpc_summaries(geom.DetectorName(), opt, channels, summary, readout);
  std::vector<ChannelInfo> const infos = channel_infos(opt, channels, summary.groups, readout);
  if (opt.table) print_table(opt, infos);

  std::string const prefix = opt.prefix.empty() ? default_prefix(opt, geo_id) : opt.prefix;
  if (opt.format == "csv" || opt.format == "both") {
    write_csv(prefix + ".csv", opt, infos);
  }
  if (opt.format == "json" || opt.format == "both") {
    write_json(prefix + ".json", geom.DetectorName(), opt, infos, summary);
  }
  return 0;
}
