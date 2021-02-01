#pragma once

#include <cstdint>

#include "constants.hpp"

namespace DwellRegions {

/**
 * @brief Stores an angle at fixed resolution in radians.
 */
class Angle
{
public:
  Angle(Angle const& angle) = default;
  Angle(Angle && angle) = default;

  static Angle from_double_degrees(double val) { return Angle{ static_cast<int32_t>((val * Constants::PI / 180) * Constants::resolution) }; }
  static Angle from_double_radians(double val)
  {
    return Angle{ static_cast<int32_t>(val * Constants::resolution) };
  }

  int32_t val() const { return _val; }

  // Operator Overloads
  // note: Maybe these asserts can be replaced by the -ftrapv flag on Linux?
  // Addition
  Angle operator+=(const Angle& rhs)
  {
    _val += rhs._val;
    return static_cast<Angle>(*this);
  }
  friend Angle operator+(Angle lhs, const Angle& rhs) { return lhs += rhs; }

  // Subtraction
  Angle operator-=(const Angle& rhs)
  {
    _val -= rhs._val;
    return static_cast<Angle>(*this);
  }
  friend Angle operator-(Angle lhs, const Angle& rhs) { return lhs -= rhs; }

  // Multiplication
//  Angle operator*=(const Angle& rhs)
//  {
//    _val *= rhs._val;
//    return static_cast<Angle>(*this);
//  }
//  friend Angle operator*(Angle lhs, const Angle& rhs) { return lhs *= rhs; }
  // Allow multiplication with a constant float
  Angle operator*=(const float rhs)
  {
    _val *= rhs;
    return static_cast<Angle>(*this);
  }
  friend Angle operator*(Angle lhs, const float rhs) { return lhs *= rhs; }

  // Allow multiplication with a constant integer
  Angle operator*=(const int32_t rhs)
  {
    _val *= rhs;
    return static_cast<Angle>(*this);
  }
  friend Angle operator*(Angle lhs, const int32_t rhs) { return lhs *= rhs; }


  // Division
  Angle operator/=(const Angle& rhs)
  {
    int64_t div_result = static_cast<int64_t>(_val) * Constants::resolution / static_cast<int64_t>(rhs.val());
    _val = div_result;
    return *this;
  }
  friend Angle operator/(Angle lhs, const Angle& rhs) { return lhs /= rhs; }

  // Trigonometric functions overload
  friend double sin(const Angle& theta);
  friend double cos(const Angle& theta);
  friend double tan(const Angle& theta);

private:
  explicit Angle(int32_t val)
    : _val(val)
  {}

  double as_double() const { return static_cast<double>(_val) / Constants::resolution; }

  int32_t _val;
};

} // namespace DwellRegions