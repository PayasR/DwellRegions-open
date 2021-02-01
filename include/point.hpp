/// This file contains an immutable Point struct
#pragma once

#include "constants.hpp"
#include <cmath>
#include <tuple>

namespace DwellRegions {
// Let's not pretend we care about higher dimensions
// See https://stackoverflow.com/questions/18636564/lat-long-or-long-lat
struct Point
{
  double Lon;
  double Lat;

  // double x, y, z; // Cartesian coordinates
  double x, y; // Cartesian coordinates (using Equirectangular projection)

  // Point(double lon, double lat): Lon(lon), Lat(lat) {}

  bool operator==(const Point& rhs) const
  {
    return ((std::abs(Lon - rhs.Lon) < Constants::tolerance) && (std::abs(Lat - rhs.Lat) < Constants::tolerance) &&
            (std::abs(x - rhs.x) < Constants::tolerance) && (std::abs(y - rhs.y) < Constants::tolerance));
  }
  bool operator!=(const Point& rhs) const { return !(*this == rhs); }
};

struct PointXYComparator
{
  bool operator()(const Point& p1, const Point& p2)
  {
    if (std::abs(p1.x - p2.x) < Constants::tolerance) {
      return std::isless((p1.y - p2.y), -Constants::tolerance);
    }
    return std::isless((p1.x - p2.x), -Constants::tolerance);
  }
};

struct PointXYEqualityComparator
{
  bool operator()(const Point& p1, const Point& p2)
  {
    return ((std::abs(p1.x - p2.x) < Constants::tolerance) && (std::abs(p1.y - p2.y) < Constants::tolerance));
  }
};

struct PointAnglePairComparator
{
  bool operator()(const std::pair<double, Point> pair1, const std::pair<double, Point> pair2)
  {
    if (pair1.first == pair2.first)
      return pair1.first == pair2.first;
    return pair1.first < pair2.first;
  }
};

Point
Point_from_lon_lat(double lon, double lat);

Point
Point_from_X_Y(double X, double Y);

bool
collinear(const Point& p1, const Point& p2, const Point& p3);

} // namespace DwellRegions
