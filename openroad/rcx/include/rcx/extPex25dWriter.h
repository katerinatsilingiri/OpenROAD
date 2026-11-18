// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <iosfwd>
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

  // Adds one wire and, the first time it appears, its layer.
  void addWire(int met,
               uint32_t num,
               double x_lo,
               double x_hi,
               double z_lo,
               double z_hi,
               double length,
               bool focus);
  // Sets the ground plane at the bottom of the pattern.
  void setGround(double z_lo, double z_hi);
  // Records the metal plane above the wires (Under and OverUnder patterns).
  void setUpperPlane(int met, double z_lo, double z_hi);
  // Adds one dielectric band. Bands are added bottom-up.
  void addDielectric(const std::string& name,
                     double epsilon,
                     double z_lo,
                     double z_hi);
};

// Writes a Pex25dPattern as a PEX25D text file.
class Pex25dWriter
{
 public:
  explicit Pex25dWriter(utl::Logger* logger);
  // Returns false and writes nothing if the pattern cannot be written.
  bool write(const std::string& path, const Pex25dPattern& pattern) const;

 private:
  // Writes GROUND_PLANE, METAL (wire layers and anchors) and dielectrics.
  void writeStack(std::ostream& out, const Pex25dPattern& pattern) const;
  // Writes one CONDUCTOR and one BOX per wire.
  void writeConductors(std::ostream& out, const Pex25dPattern& pattern) const;
  // Rounds v to the PEX25D grid and returns it with 3 decimals.
  std::string len(double v) const;

  utl::Logger* logger_;
};

}  // namespace rcx
