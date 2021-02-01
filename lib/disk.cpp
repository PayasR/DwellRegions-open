#include "../include/disk.hpp"

#include <cassert>
#include <iostream>
#include <optional>

namespace DwellRegions {
using DwellRegions::Disk;
using DwellRegions::Point;

std::optional<Disk>
brute_force_SEC(std::vector<Point> points)
{
  double radius = 1000000.0;
  double center_x = 0.0, center_y = 0.0;
  bool is_valid_circle = false;

  // all possible 3-points circle
  for (size_t i = 0; i < points.size() - 2; i++) {
    for (size_t j = i + 1; j < points.size() - 1; j++) {
      for (size_t k = j + 1; k < points.size(); k++) {
        Point p1 = Point_from_X_Y(points[i].x, points[i].y);
        Point p2 = Point_from_X_Y(points[j].x, points[j].y);
        Point p3 = Point_from_X_Y(points[k].x, points[k].y);

        if (collinear(p1, p2, p3)) // linearity check
          continue;

        // check if two points are vertical
        if (std::fabs(p1.x - p2.x) < Constants::tolerance) {
          Point temp = p2;
          p2 = p3;
          p3 = temp;
        } else if (std::fabs(p2.x - p3.x) < Constants::tolerance) {
          Point temp = p1;
          p1 = p2;
          p2 = temp;
        }

        // points are not collinear and vertical
        double m_1_2 = (p1.y - p2.y) / (p1.x - p2.x); // Slope of the line joining points 1 and 2
        double m_2_3 = (p2.y - p3.y) / (p2.x - p3.x); // Slope of the line joining points 2 and 3

        // Bisectors of the lines joining points 1_2 and 2_3
        Point bisect_1_2 = Point_from_X_Y((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
        Point bisect_2_3 = Point_from_X_Y((p2.x + p3.x) / 2.0, (p2.y + p3.y) / 2.0);

        double temp_y = (bisect_1_2.x - bisect_2_3.x + m_1_2 * bisect_1_2.y - m_2_3 * bisect_2_3.y) / (m_1_2 - m_2_3);
        double temp_x = bisect_1_2.x - m_1_2 * (temp_y - bisect_1_2.y);
        Point temp_p = Point_from_X_Y(temp_x, temp_y);
        double temp_radius = Util::Geo::distance_in_miles_from_X_Y(p1, temp_p);

        if (std::isgreater((temp_radius - radius), Constants::tolerance))
          continue;
        is_valid_circle = check_circle_validity(Disk(temp_p, temp_radius), points);
        if (is_valid_circle) {
          radius = temp_radius;
          center_x = temp_p.x;
          center_y = temp_p.y;
        }
      }
    }
  }

  // all possible 2-points circle
  for (size_t i = 0; i < points.size() - 1; i++) {
    for (size_t j = i + 1; j < points.size(); j++) {
      Disk circle(points[i], points[j]);
      if (std::isgreater((circle.get_radius() - radius), Constants::tolerance))
        continue;
      is_valid_circle = check_circle_validity(circle, points);
      if (is_valid_circle) {
        radius = circle.get_radius();
        center_x = circle.get_center().x;
        center_y = circle.get_center().y;
      }
    }
  }
  Point center = Point_from_X_Y(center_x, center_y);
  return Disk(center, radius);
}

bool
check_circle_validity(Disk circle, std::vector<Point> points)
{
  for (size_t i = 0; i < points.size(); i++) {
    if (!circle.in_disk(points[i]))
      return false;
  }
  return true;
}

// This implementation of Welzl's algorithm follows the pseudocode from wikipedia
// https://en.wikipedia.org/wiki/Smallest-circle_problem#Welzl's_algorithm
// algorithm welzl is [7]
//     input: Finite sets P and R of points in the plane |R|≤ 3.
//     output: Minimal disk enclosing P with R on the boundary.
//     if P is empty or |R| = 3 then
//         return trivial(R)
//     choose p in P (randomly and uniformly)
//     D := welzl(P - { p }, R)
//     if p is in D then
//         return D
//     return welzl(P - { p }, R ∪ { p })
std::optional<Disk>
welzl(std::vector<Point> points, std::vector<Point> boundary)
{
  // Replace a boundary point with the point that is outside the circle
  // boundary should contain at most 3-points
  // (extension of welzl's algorithm)
  if (boundary.size() == 4) {
    size_t l = 3;
    for (size_t i = 0; i < boundary.size() - 2; i++) {
      for (size_t j = i + 1; j < boundary.size() - 1; j++) {
        for (size_t k = j + 1; k < boundary.size(); k++) {
          Disk circle(boundary[i], boundary[j], boundary[k]);
          if (circle.in_disk(boundary[l])) {
            auto non_boundary_point = boundary[l];
            boundary.erase(boundary.begin() + l);
            points.push_back(non_boundary_point);
            return circle;
          }
          l--;
        }
      }
    }
  }

  // Beginning of original welzl's algorithm
  if (points.empty() || boundary.size() == 3) {
    if (boundary.size() == 0)
      return std::nullopt;
    else if (boundary.size() == 1)
      return std::optional<Disk>(Disk(boundary[0], boundary[0]));
    else if (boundary.size() == 2)
      return std::optional<Disk>(Disk(boundary[0], boundary[1]));
    else if (boundary.size() == 3)
      return std::optional<Disk>(Disk(boundary[0], boundary[1], boundary[2]));
  }

  auto random_pos = Util::rand_size_t(0, points.size() - 1);
  auto random_point = points[random_pos];

  points.erase(points.begin() + random_pos);
  auto disk = welzl(points, boundary);

  if (disk.has_value() && disk->in_disk(random_point)) {
    return disk;
  }

  boundary.emplace_back(random_point);
  return welzl(points, boundary);
}

Disk::Disk(const Point& p1, const Point& p2, const Point& p3)
{
  std::vector<Point> basis_points;
  basis_points.push_back(p1);
  basis_points.push_back(p2);
  basis_points.push_back(p3);
  size_t bsize = 3;

  // Check if the points are collinear
  if (collinear(p1, p2, p3)) {
    size_t k = 2;
    for (size_t i = 0; i < bsize - 1; i++) {
      for (size_t j = i + 1; j < bsize; j++) {
        Disk circle(basis_points[i], basis_points[j]);
        if (circle.in_disk(basis_points[k])) {
          center = circle.get_center();
          radius = circle.get_radius();
          return;
        }
        k--;
      }
    }
  }

  // check if two points are aligned in vertical axis
  if (std::fabs(basis_points[0].x - basis_points[1].x) < Constants::tolerance) {
    Point temp = basis_points[1];
    basis_points[1] = basis_points[2];
    basis_points[2] = temp;
  } else if (std::fabs(basis_points[1].x - basis_points[2].x) < Constants::tolerance) {
    Point temp = basis_points[0];
    basis_points[0] = basis_points[1];
    basis_points[1] = temp;
  }

  // neither points are collinear or vertical
  double m_1_2 = (basis_points[0].y - basis_points[1].y) /
                 (basis_points[0].x - basis_points[1].x); // Slope of the line joining points 1 and 2
  double m_2_3 = (basis_points[1].y - basis_points[2].y) /
                 (basis_points[1].x - basis_points[2].x); // Slope of the line joining points 2 and 3

  // Bisectors of the lines joining points 1_2 and 2_3
  Point bisect_1_2 =
    Point_from_X_Y((basis_points[0].x + basis_points[1].x) / 2.0, (basis_points[0].y + basis_points[1].y) / 2.0);
  Point bisect_2_3 =
    Point_from_X_Y((basis_points[1].x + basis_points[2].x) / 2.0, (basis_points[1].y + basis_points[2].y) / 2.0);

  center.y = (bisect_1_2.x - bisect_2_3.x + m_1_2 * bisect_1_2.y - m_2_3 * bisect_2_3.y) / (m_1_2 - m_2_3);
  center.x = bisect_1_2.x - m_1_2 * (center.y - bisect_1_2.y);

  assert(std::islessequal((Util::Geo::distance_in_miles_from_X_Y(basis_points[0], center) -
                           Util::Geo::distance_in_miles_from_X_Y(basis_points[1], center)),
                          Constants::tolerance));
  assert(std::islessequal((Util::Geo::distance_in_miles_from_X_Y(basis_points[1], center) -
                           Util::Geo::distance_in_miles_from_X_Y(basis_points[2], center)),
                          Constants::tolerance));
  assert(std::islessequal((Util::Geo::distance_in_miles_from_X_Y(basis_points[0], center) -
                           Util::Geo::distance_in_miles_from_X_Y(basis_points[2], center)),
                          Constants::tolerance));

  radius = Util::Geo::distance_in_miles_from_X_Y(basis_points[0], center);
}

[[nodiscard]] Point
Disk::point_at_angle(double theta) const
{
  // Compute the distance from the center in the Y direction, convert to
  // degrees
  const auto distance_from_center_in_y = (radius * std::sin(theta));
  auto y = center.y + distance_from_center_in_y; // final latitude of the point in degrees

  // The same operation applied in the X axis
  const auto distance_from_center_in_x = (radius * std::cos(theta));
  auto x = center.x + distance_from_center_in_x; // final longitude of the
                                                 // point in degrees

  return Point_from_X_Y(x, y);
}

// Returns two circles intersection points
// Source:
// http://www.ambrsoft.com/TrigoCalc/Circles2/circle2intersection/CircleCircleIntersection.htm
[[nodiscard]] std::optional<std::vector<Point>>
Disk::disk_disk_intersection(const Disk& disk1) const
{
  const auto center_x1 = center.x;
  const auto center_y1 = center.y;
  const auto center_x2 = disk1.get_center().x;
  const auto center_y2 = disk1.get_center().y;

  // Distance between two circles' centers
  const auto distance_betweeen_centers =
    std::sqrt((center_x1 - center_x2) * (center_x1 - center_x2) + (center_y1 - center_y2) * (center_y1 - center_y2));

  // Two disks intersect or have tangents (either inner or outer)
  if (std::isgreaterequal(((radius + disk1.get_radius()) - distance_betweeen_centers), Constants::tolerance) &&
      std::isgreaterequal((distance_betweeen_centers - std::fabs(radius - disk1.get_radius())), Constants::tolerance)) {
    // Area according to Heron's formula
    const auto a1 = distance_betweeen_centers + radius + disk1.get_radius();
    const auto a2 = distance_betweeen_centers + radius - disk1.get_radius();
    const auto a3 = distance_betweeen_centers - radius + disk1.get_radius();
    const auto a4 = -distance_betweeen_centers + radius + disk1.get_radius();
    const auto area = std::sqrt(a1 * a2 * a3 * a4) / 4.0;

    // Calculating x axis intersection values
    auto value1 = ((center_x1 + center_x2) / 2.0) +
                  ((center_x1 - center_x2) * (radius * radius - disk1.get_radius() * disk1.get_radius()) /
                   (2.0 * distance_betweeen_centers * distance_betweeen_centers));
    auto value2 = (2.0 * (center_y1 - center_y2) * area) / (distance_betweeen_centers * distance_betweeen_centers);
    const auto x1 = value1 + value2;
    const auto x2 = value1 - value2;

    // Calculating y axis intersection values
    value1 = ((center_y1 + center_y2) / 2.0) +
             ((center_y1 - center_y2) * (radius * radius - disk1.get_radius() * disk1.get_radius()) /
              (2.0 * distance_betweeen_centers * distance_betweeen_centers));
    value2 = (2.0 * (center_x1 - center_x2) * area) / (distance_betweeen_centers * distance_betweeen_centers);
    const auto y1 = value1 + value2;
    const auto y2 = value1 - value2;

    std::vector<Point> intersection_points;
    // Intersection points are (x1, y1) and (x2, y2)
    // Because for every x we have two values of y, and the same thing for y,
    // we have to verify that the intersection points are on the
    // circle otherwise we have to swap between the points
    if (std::fabs((x1 - center_x1) * (x1 - center_x1) + (y1 - center_y1) * (y1 - center_y1) - (radius * radius)) >
        Constants::tolerance) {
      intersection_points.push_back(Point_from_X_Y(x1, y2));
      intersection_points.push_back(Point_from_X_Y(x2, y1));
    } else {
      intersection_points.push_back(Point_from_X_Y(x1, y1));
      intersection_points.push_back(Point_from_X_Y(x2, y2));
    }
    return intersection_points;
  }

  return std::nullopt;
}
} // namespace DwellRegions
