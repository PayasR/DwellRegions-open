#pragma once

#include <cassert>
#include <cstdint>
#include <numeric>

#include "constants.hpp"

namespace DwellRegions {

/**
 * @brief Stores a Coefficient of the equation of a line segment at a fixed resolution.
 */
class Coefficient
{
public:
  static Coefficient from_val(const int32_t val) { return Coefficient{ val }; }
  static Coefficient from_double(const double val) { return Coefficient{ static_cast<int32_t>(val * Constants::resolution) }; }

  int32_t val() const { return _val; }
  float as_float() const { return static_cast<float>(_val) / Constants::resolution; }

  // Operator Overloads
  // note: Maybe these asserts can be replaced by the -ftrapv flag on Linux?
  // Subtraction
  Coefficient operator-=(const Coefficient& rhs)
  {
    _val -= rhs._val;
    return *this;
  }
  friend Coefficient operator-(Coefficient lhs, const Coefficient& rhs) { return lhs -= rhs; }

  // Multiplication
  Coefficient operator*=(const Coefficient& rhs)
  {
    int64_t mul_result = static_cast<int64_t>(_val) * static_cast<int64_t>(rhs.val());
    _val = mul_result / Constants::resolution;
    return *this;
  }
  friend Coefficient operator*(Coefficient lhs, const Coefficient& rhs) { return lhs *= rhs; }
//  Coefficient operator*=(const int32_t rhs)
//  {
//    _val *= rhs;
//    return *this;
//  }
//  friend Coefficient operator*(Coefficient lhs, const int32_t rhs) { return lhs *= rhs; }
//  Coefficient operator*=(const float rhs)
//  {
//    _val *= rhs;
//    return *this;
//  }
//  friend Coefficient operator*(Coefficient lhs, const float rhs) { return lhs *= rhs; }


  // Division
  Coefficient operator/=(const Coefficient& rhs)
  {
    int64_t div_result = static_cast<int64_t>(_val ) * Constants::resolution / static_cast<int64_t>(rhs.val());
    _val = div_result;
    return *this;
  }
  friend Coefficient operator/(Coefficient lhs, const Coefficient& rhs) { return lhs /= rhs; }

  Coefficient operator/=(const int32_t rhs)
  {
    _val /= rhs;
    return *this;
  }
  friend Coefficient operator/(Coefficient lhs, const int32_t rhs) { return lhs /= rhs; }

  Coefficient operator/=(const float rhs)
  {
    _val /= rhs;
    return *this;
  }
  friend Coefficient operator/(Coefficient lhs, const float rhs) { return lhs /= rhs; }

private:
  explicit Coefficient(const int32_t val)
    : _val(val)
  {}

  int32_t _val;
};
} // namespace DwellRegions