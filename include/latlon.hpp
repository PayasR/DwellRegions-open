#pragma once

#include "constants.hpp"

#include <cassert>
#include <cstdint>

namespace DwellRegions {

/**
 * @brief Stores the Latitude of a point at a fixed resolution in degrees
 */
class Latitude
{
public:
  explicit Latitude(int32_t val)
    : _val(val)
  {}

  static Latitude from_degrees(double val)
  {
    assert(val >= -90.0 && val <= 90.0);
    return Latitude(static_cast<int32_t>(val * Constants::resolution));
  }

  static Latitude from_radians(double val)
  {
    const int32_t val_at_resolution = static_cast<int32_t>(val * Constants::resolution);
    return Latitude{static_cast<int32_t>(val_at_resolution * 180 / Constants::PI)};
  }

  int32_t val() { return _val; }
  double in_degrees() const { return static_cast<double>(_val / Constants::resolution); }
  double in_radians() const
  {
    return static_cast<double>((_val * Constants::PI / 180) / Constants::resolution);
  }

  bool operator==(const Latitude& rhs) const { return (_val == rhs._val); }
  bool operator!=(const Latitude& rhs) const { return !(*this == rhs); }

private:
  int32_t _val; // with RESOLUTION set to 10^4, this is stored in units of 10^-4 degrees, so 1 val = 100 micro-degrees
};

// ---------------------------------------------------------------------------------------------------------------------

/**
 * @brief Stores the Longitude of a point at a fixed resolution in degrees
 */
class Longitude
{
public:
  Longitude(int32_t val)
    : _val(val)
  {}

  static Longitude from_degrees(double val)
  {
    assert(val >= -180.0 && val <= 180.0);
    return Longitude(static_cast<int32_t>(val * Constants::resolution));
  }

  static Longitude from_radians(double val)
  {
    const int32_t val_at_resolution = static_cast<int32_t>(val * Constants::resolution);
    return Longitude(val_at_resolution * 180 / Constants::PI);
  }

  int32_t val() { return _val; }
  double in_degrees() const { return static_cast<double>(_val / Constants::resolution); }
  double in_radians() const
  {
    return static_cast<double>((_val * Constants::PI / 180) / Constants::resolution);
  }

  bool operator==(const Longitude& rhs) const { return (_val == rhs._val); }
  bool operator!=(const Longitude& rhs) const { return !(*this == rhs); }

private:
  int32_t _val; // Stored in units of of 10^-4 degrees, so 1 val  = 100 micro-degrees
};
} // namespace DwellRegions
