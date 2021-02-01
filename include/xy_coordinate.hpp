#pragma once

#include "constants.hpp"

#include <cassert>
#include <cstdint>

namespace DwellRegions {

/**
 * @brief Represents a coordinate on the X-Y plane. Also provides safe arithmetic operations.
 */
class Coordinate
{
public:
  Coordinate(Coordinate const& rhs) = default;
  Coordinate(Coordinate const&& rhs) noexcept
    : _val(rhs._val)
  {}

  static Coordinate zero() { return Coordinate{ 0 }; }
  static Coordinate from_val(int32_t val) { return Coordinate{ val }; }
  static Coordinate from_double(double val) { return Coordinate{ static_cast<int32_t>(val * Constants::resolution) }; }

  [[nodiscard]] int32_t val() const { return _val; }
  [[nodiscard]] float as_float() const { return static_cast<float>(_val) / Constants::resolution; }

  Coordinate& operator=(Coordinate const& rhs)
  {
    assert(this != &rhs);
    this->_val = rhs._val;
    return *this;
  }

  Coordinate& operator=(Coordinate const&& rhs)
  {
    this->_val = rhs._val;
    return *this;
  }

  bool operator==(const Coordinate& rhs) const { return _val == rhs._val; }
  bool operator!=(const Coordinate& rhs) const { return !(*this == rhs); }

  // Operator Overloads
  // note: Maybe these asserts can be replaced by the -ftrapv flag on Linux?
  // Addition
  Coordinate operator+=(const Coordinate& rhs)
  {
    _val += rhs._val;
    return *this;
  }
  friend Coordinate operator+(Coordinate lhs, const Coordinate& rhs) { return lhs += rhs; }

  // Subtraction
  Coordinate operator-=(const Coordinate& rhs)
  {
    _val -= rhs._val;
    return *this;
  }
  friend Coordinate operator-(Coordinate lhs, const Coordinate& rhs) { return lhs -= rhs; }

  // Multiplication
  Coordinate operator*=(const Coordinate& rhs)
  {
    int64_t mul_result = static_cast<int64_t>(_val) * static_cast<int64_t>(rhs.val());
    _val = mul_result / Constants::resolution;
    return *this;
  }
  friend Coordinate operator*(Coordinate lhs, const Coordinate& rhs) { return lhs *= rhs; }

  // Allow multiplication with a constant float
  Coordinate operator*=(const float rhs)
  {
    int64_t mul_result = static_cast<int64_t>(_val) * static_cast<int64_t>(rhs);
    _val = mul_result;
    return *this;
  }
  friend Coordinate operator*(Coordinate lhs, const float rhs) { return lhs *= rhs; }

  // Allow multiplication with a constant integer
  Coordinate operator*=(const int32_t rhs)
  {
    int64_t mul_result = static_cast<int64_t>(_val) * static_cast<int64_t>(rhs);
    _val = mul_result;
    return *this;
  }
  friend Coordinate operator*(Coordinate lhs, const int32_t rhs) { return lhs *= rhs; }

  // Division
  Coordinate operator/=(const Coordinate& rhs)
  {
    int64_t div_result = static_cast<int64_t>(_val) * Constants::resolution / static_cast<int64_t>(rhs.val());
    _val = div_result;
    return *this;
  }
  friend Coordinate operator/(Coordinate lhs, const Coordinate& rhs) { return lhs /= rhs; }

  // Allow division with a constant float
  Coordinate operator/=(const float rhs)
  {
    _val /= rhs;
    return *this;
  }
  friend Coordinate operator/(Coordinate lhs, const float rhs) { return lhs /= rhs; }

  // Allow multiplication with a constant integer
  Coordinate operator/=(const int32_t rhs)
  {
    _val /= rhs;
    return *this;
  }
  friend Coordinate operator/(Coordinate lhs, const int32_t rhs) { return lhs /= rhs; }

  // Comparison operators
  friend bool operator<(const Coordinate& lhs, const Coordinate& rhs) { return lhs._val < rhs._val; }
  friend bool operator>(const Coordinate& lhs, const Coordinate& rhs) { return rhs < lhs; }
  friend bool operator<=(const Coordinate& lhs, const Coordinate& rhs) { return !(lhs > rhs); }
  friend bool operator>=(const Coordinate& lhs, const Coordinate& rhs) { return !(lhs < rhs); }

private:
  explicit Coordinate(int32_t val)
    : _val(val)
  {}

  int32_t _val; // Stored as 10^-4 units, so 1 val = 100 micro-coordinates
};
} // namespace DwellRegions
