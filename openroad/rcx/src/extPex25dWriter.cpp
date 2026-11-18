// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "rcx/extPex25dWriter.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "utl/Logger.h"

namespace rcx {

namespace {

// PEX25D grid in um, the same precision as the UF output.
constexpr double kGrid = 0.001;
// Two heights closer than this are treated as the same height.
constexpr double kTol = 1e-6;

// Dielectric bands that reach above the ground plane. A band that ends at or
// below the top of the ground plane lies inside it and has no effect.
std::vector<Pex25dPattern::Dielectric> bandsAboveGround(
    const Pex25dPattern& pattern)
{
  std::vector<Pex25dPattern::Dielectric> bands;
  for (const auto& band : pattern.dielectrics) {
    if (band.z.hi > pattern.ground.hi + kTol) {
      bands.push_back(band);
    }
  }
  return bands;
}

// Simulation domain in PEX25D coordinates, following the UF convention.
struct Domain
{
  double lo[3];
  double hi[3];
};

// x: the window widened by its own width on both sides (SIM_WIN_EXT).
// y: the wire length plus the same length again.
// z: from the top of the ground plane to the top of the dielectric stack.
// The window comes from the UF output and only covers the wires of the main
// layer. In UnderDiag patterns the wires on the layer above can lie outside
// it; findProblem() then refuses the pattern instead of writing a domain that
// cuts those wires off.
// Expects a pattern with a window and at least one band above the ground.
Domain domainOf(const Pex25dPattern& pattern)
{
  const double dx = pattern.window_x_hi - pattern.window_x_lo;
  const double z_top = bandsAboveGround(pattern).back().z.hi;
  return {{pattern.window_x_lo - dx, 0.0, pattern.ground.hi},
          {pattern.window_x_hi + dx, 2 * pattern.window_length, z_top}};
}

// Returns why the pattern cannot be written, or an empty string if it can.
std::string findProblem(const Pex25dPattern& pattern)
{
  if (!pattern.valid) {
    return "wires of one layer have different heights";
  }
  if (pattern.wires.empty()) {
    return "it has no wires";
  }
  if (!pattern.has_ground) {
    return "it has no ground plane";
  }
  if (pattern.has_upper_plane) {
    return "patterns with an upper ground plane are not supported yet";
  }
  const auto bands = bandsAboveGround(pattern);
  if (bands.empty()) {
    return "it has no dielectrics";
  }
  if (bands.front().z.lo > pattern.ground.hi + kTol) {
    return "there is a gap below dielectric " + bands.front().name;
  }
  for (std::size_t i = 0; i < bands.size(); ++i) {
    if (bands[i].z.hi - bands[i].z.lo <= kTol) {
      return "dielectric " + bands[i].name + " has no thickness";
    }
    if (i + 1 < bands.size()
        && std::abs(bands[i].z.hi - bands[i + 1].z.lo) > kTol) {
      return "dielectric " + bands[i + 1].name + " does not start where "
             + bands[i].name + " ends";
    }
  }

  if (!pattern.has_window) {
    return "it has no simulation window";
  }
  const Domain domain = domainOf(pattern);
  for (const auto& w : pattern.wires) {
    if (w.x_lo < domain.lo[0] - kTol || w.x_hi > domain.hi[0] + kTol
        || w.length > domain.hi[1] + kTol) {
      return "wire " + w.name + " is outside the simulation domain";
    }
  }
  for (const auto& m : pattern.metals) {
    if (m.z.hi > domain.hi[2] + kTol) {
      return "layer " + m.name + " reaches above the simulation domain";
    }
  }

  return "";
}

}  // namespace

void Pex25dPattern::addWire(const int met,
                            const uint32_t num,
                            const double x_lo,
                            const double x_hi,
                            const double z_lo,
                            const double z_hi,
                            const double length,
                            const bool focus)
{
  const std::string layer = "M" + std::to_string(met);
  const auto it = std::ranges::find_if(
      metals, [&layer](const Layer& m) { return m.name == layer; });
  if (it == metals.end()) {
    metals.push_back({layer, {z_lo, z_hi}});
  } else if (std::abs(it->z.lo - z_lo) > 1e-9
             || std::abs(it->z.hi - z_hi) > 1e-9) {
    valid = false;  // one layer, two different z ranges
  }
  wires.push_back(
      {layer + "_w" + std::to_string(num), layer, x_lo, x_hi, length, focus});
}

void Pex25dPattern::setGround(const double z_lo, const double z_hi)
{
  has_ground = true;
  ground = {z_lo, z_hi};
}

void Pex25dPattern::setUpperPlane(const int met,
                                  const double z_lo,
                                  const double z_hi)
{
  has_upper_plane = true;
  upper_plane = {"M" + std::to_string(met), {z_lo, z_hi}};
}

void Pex25dPattern::addDielectric(const std::string& name,
                                  const double epsilon,
                                  const double z_lo,
                                  const double z_hi)
{
  dielectrics.push_back({name, epsilon, {z_lo, z_hi}});
}

void Pex25dPattern::setWindow(const double x_lo,
                              const double x_hi,
                              const double length)
{
  has_window = true;
  window_x_lo = x_lo;
  window_x_hi = x_hi;
  window_length = length;
}

Pex25dWriter::Pex25dWriter(utl::Logger* logger) : logger_(logger)
{
}

bool Pex25dWriter::write(const std::string& path,
                         const Pex25dPattern& pattern) const
{
  const std::string problem = findProblem(pattern);
  if (!problem.empty()) {
    logger_->warn(utl::RCX,
                  541,
                  "Pattern {} is not written as PEX25D: {}.",
                  pattern.name,
                  problem);
    return false;
  }

  std::ofstream out(path);
  if (!out) {
    logger_->warn(utl::RCX, 540, "Cannot open PEX25D file {}.", path);
    return false;
  }
  out << "PEX25D 1.0-rc1\n";
  out << "UNITS LENGTH um GRID 0.001\n";
  out << "META technology \"" << pattern.process << "\"\n";
  out << "META source_cell \"" << pattern.name << "\"\n";
  out << "META generator \"OpenRCX gen_solver_patterns\"\n";
  for (const auto& w : pattern.wires) {
    if (w.focus) {
      out << "META focus_conductor " << w.name << "\n";
    }
  }

  writeStack(out, pattern);
  writeConductors(out, pattern);
  writeDomain(out, pattern);
  return true;
}

void Pex25dWriter::writeStack(std::ostream& out,
                              const Pex25dPattern& pattern) const
{
  const auto bands = bandsAboveGround(pattern);

  // PEX25D dielectrics span from the bottom face of one profile to the bottom
  // face of another, so every band bottom needs a profile starting there: the
  // ground plane for the first band, otherwise a wire layer starting at that
  // height, otherwise a shape-less anchor layer.
  std::vector<std::pair<std::string, double>> anchors;  // name, z
  const auto profile_at = [&](const double z) -> std::string {
    for (const auto& m : pattern.metals) {
      if (std::abs(m.z.lo - z) < kTol) {
        return m.name;
      }
    }
    for (const auto& [name, anchor_z] : anchors) {
      if (std::abs(anchor_z - z) < kTol) {
        return name;
      }
    }
    anchors.emplace_back("anchor_" + std::to_string(anchors.size()), z);
    return anchors.back().first;
  };
  std::vector<std::string> bottom(bands.size());
  bottom[0] = "gnd";
  for (std::size_t i = 1; i < bands.size(); ++i) {
    bottom[i] = profile_at(bands[i].z.lo);
  }

  out << "\nGROUND_PLANE gnd Z_OFFSETS " << len(pattern.ground.lo) << " "
      << len(pattern.ground.hi) << "\n";
  for (const auto& m : pattern.metals) {
    out << "METAL " << m.name << " Z_OFFSETS " << len(m.z.lo) << " "
        << len(m.z.hi) << "\n";
  }
  for (const auto& [name, z] : anchors) {
    out << "METAL " << name << " Z_OFFSETS " << len(z) << " " << len(z + kGrid)
        << "\n";
  }

  out << "\n";
  for (std::size_t i = 0; i + 1 < bands.size(); ++i) {
    out << "DIELECTRIC_SIMPLE " << bands[i].name << " WRAPS " << bottom[i]
        << " PERMITTIVITY " << bands[i].epsilon << " BETWEEN " << bottom[i]
        << " " << bottom[i + 1] << "\n";
  }
  out << "DIELECTRIC_BACKGROUND " << bands.back().name << " PERMITTIVITY "
      << bands.back().epsilon << "\n";
}

void Pex25dWriter::writeConductors(std::ostream& out,
                                   const Pex25dPattern& pattern) const
{
  out << "\n";
  for (const auto& w : pattern.wires) {
    out << "CONDUCTOR " << w.name << " " << w.name << "\n";
    out << "BOX CONDUCTOR " << w.name << " LAYER " << w.layer << " LL "
        << len(w.x_lo) << " " << len(0.0) << " UR " << len(w.x_hi) << " "
        << len(w.length) << "\n";
  }
}

void Pex25dWriter::writeDomain(std::ostream& out,
                               const Pex25dPattern& pattern) const
{
  const Domain d = domainOf(pattern);
  out << "\nDOMAIN_BOX LL " << len(d.lo[0]) << " " << len(d.lo[1]) << " "
      << len(d.lo[2]) << " UR " << len(d.hi[0]) << " " << len(d.hi[1]) << " "
      << len(d.hi[2]) << "\n";
}

std::string Pex25dWriter::len(const double v) const
{
  const double snapped = std::llround(v / kGrid) * kGrid;
  if (std::abs(snapped - v) > kTol) {
    logger_->warn(
        utl::RCX, 542, "PEX25D value {} is not on the {} um grid.", v, kGrid);
  }
  std::ostringstream s;
  s << std::fixed << std::setprecision(3) << snapped;
  return s.str();
}

}  // namespace rcx
