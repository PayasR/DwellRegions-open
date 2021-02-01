#include "../include/angle.hpp"

#include <cmath>

namespace DwellRegions {

double
sin(const Angle& theta)
{
  return std::sin(theta.as_double());
}

double
cos(const Angle& theta)
{
  return std::cos(theta.as_double());
}

double
tan(const Angle& theta)
{
  return std::tan(theta.as_double());
}
} // namespace DwellRegions
