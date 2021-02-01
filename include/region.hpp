#pragma once

#include "point.hpp"
#include <vector>

namespace DwellRegions {
struct Region
{
  Region(std::vector<Point>& p)
    : points(p)
  {}

  const std::vector<Point> points;
};
}
