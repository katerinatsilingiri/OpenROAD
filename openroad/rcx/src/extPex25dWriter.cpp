// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "rcx/extPex25dWriter.h"

#include <fstream>
#include <string>

#include "utl/Logger.h"

namespace rcx {

Pex25dWriter::Pex25dWriter(utl::Logger* logger) : logger_(logger)
{
}

bool Pex25dWriter::write(const std::string& path,
                         const Pex25dPattern& pattern) const
{
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
  return true;
}

}  // namespace rcx
