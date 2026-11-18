// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "rcx/extPex25dWriter.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "utl/Logger.h"

namespace rcx {

Pex25dWriter::Pex25dWriter(utl::Logger* logger) : logger_(logger)
{
}

bool Pex25dWriter::write(const std::string& path,
                         const Pex25dPattern& pattern) const
{
  if (!pattern.valid || pattern.wires.empty()) {
    logger_->warn(
        utl::RCX, 541, "Pattern {} cannot be written as PEX25D.", pattern.name);
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

  out << "\n";
  for (const auto& m : pattern.metals) {
    out << "METAL " << m.name << " Z_OFFSETS " << len(m.z.lo) << " "
        << len(m.z.hi) << "\n";
  }
  out << "\n";
  for (const auto& w : pattern.wires) {
    out << "CONDUCTOR " << w.name << " " << w.name << "\n";
    out << "BOX CONDUCTOR " << w.name << " LAYER " << w.layer << " LL "
        << len(w.x_lo) << " " << len(0.0) << " UR " << len(w.x_hi) << " "
        << len(w.length) << "\n";
  }

  return true;
}

std::string Pex25dWriter::len(const double v) const
{
  constexpr double kGrid = 0.001;  // um, same precision as the UF
  const double snapped = std::llround(v / kGrid) * kGrid;

  if (std::abs(snapped - v) > 1e-6) {
    logger_->warn(
        utl::RCX, 542, "PEX25D value {} is not on the {} um grid.", v, kGrid);
  }
  std::ostringstream s;
  s << std::fixed << std::setprecision(3) << snapped;
  return s.str();
}

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

}  // namespace rcx
