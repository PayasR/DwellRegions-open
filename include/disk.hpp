/// This file contains the definition of a Disk class
#pragma once

#include "point.hpp"
#include "util.hpp"

#include <cassert>
#include <optional>
#include <vector>

namespace DwellRegions {

class Disk
{
public:
  Disk(const Point& center, const double radius)
    : center(center)
    , radius(radius)
  {}

  // This is the twoPointDisk function in Reaz's code
  Disk(const Point& p1, const Point& p2)
  {
    center = Point_from_X_Y((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    assert(std::islessequal(
      (Util::Geo::distance_in_miles_from_X_Y(p1, center) - Util::Geo::distance_in_miles_from_X_Y(p2, center)),
      Constants::tolerance));
    radius = Util::Geo::distance_in_miles_from_X_Y(p1, center);
  }

  [[nodiscard]] auto get_radius() const { return radius; }
  auto set_radius(double r) { this->radius = r; }
  [[nodiscard]] auto get_center() const { return center; }

  Disk(const Point& p1, const Point& p2, const Point& p3);

  /// Returns true if the point lies on the disk, false otherwise
  [[nodiscard]] bool in_disk(const Point& point) const
  {
    auto distance = Util::Geo::distance_in_miles_from_X_Y(this->center, point);
    return distance <= (radius + Constants::maxdist);
  }

  // TODO: Check the veracity of this function!
  /// Returns a point on disk at theta radians
  [[nodiscard]] Point point_at_angle(double theta) const;

  // Returns the intersection points between disk and disk1
  [[nodiscard]] std::optional<std::vector<Point>> disk_disk_intersection(const Disk& disk1) const;

private:
  Point center;
  double radius; // Stored in miles
};

std::optional<Disk>
welzl(std::vector<Point> points, std::vector<Point> boundary);

std::optional<Disk>
brute_force_SEC(std::vector<Point> points);

bool
check_circle_validity(Disk circle, std::vector<Point> points);

} // namespace DwellRegions
