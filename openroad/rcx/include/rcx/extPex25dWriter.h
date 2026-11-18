// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <string>
#include <vector>

namespace utl {
class Logger;
}

namespace rcx {

// Geometry of one solver pattern, collected while the pattern is generated.
struct Pex25dPattern
{
  struct ZRange
  {
    double lo = 0.0;
    double hi = 0.0;
  };
  struct Layer
  {
    std::string name;
    ZRange z;
  };
  struct Dielectric
  {
    std::string name;
    double epsilon = 0.0;
    ZRange z;
  };
  struct Wire
  {
    std::string name;
    std::string layer;
    double x_lo = 0.0;
    double x_hi = 0.0;
    double length = 0.0;
    bool focus = false;
  };

  std::string name;     // pattern directory
  std::string process;  // process name
  bool has_ground = false;
  ZRange ground;
  bool has_upper_plane = false;
  Layer upper_plane;
  std::vector<Layer> metals;
  std::vector<Dielectric> dielectrics;
  std::vector<Wire> wires;
  bool valid = true;
};

// Writes a Pex25dPattern as a PEX25D text file.
class Pex25dWriter
{
 public:
  explicit Pex25dWriter(utl::Logger* logger);
  // Returns false and writes nothing if the pattern cannot be written.
  bool write(const std::string& path, const Pex25dPattern& pattern) const;

 private:
  utl::Logger* logger_;
};

}  // namespace rcx
