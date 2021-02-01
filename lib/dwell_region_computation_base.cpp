#include "../include/dwell_region_computation_base.hpp"
#include "../include/disk.hpp"
#include "../include/line.hpp"
#include "../include/util.hpp"

#include <algorithm>
#include <iostream>
#include <optional>
#include <set>
#include <unordered_set>

namespace DwellRegions {
double
dot_product_with_unit_vector(Point point, double theta)
{
  // This is the correct version.
  // return (point.Lon * std::cos(theta)) + (point.Lat * std::sin(theta));
  return (point.x * std::cos(theta)) + (point.y * std::sin(theta));
}

// Returns the points in counterclockwise orientation
std::vector<Point>
oriented_counterclockwise(const std::vector<Point>& points)
{
  double center_x = 0.0, center_y = 0.0;
  for (auto point : points) {
    center_x += point.x;
    center_y += point.y;
  }
  center_x = center_x / points.size();
  center_y = center_y / points.size();

  std::vector<DwellRegions::Point> oriented_points;
  std::vector<std::pair<double, DwellRegions::Point>> points_with_angle;
  for (size_t i = 0; i < points.size(); i++) {
    double angle = std::atan2((points[i].y - center_y), (points[i].x - center_x));
    points_with_angle.emplace_back(std::make_pair(angle, points[i]));
  }
  std::sort(points_with_angle.begin(), points_with_angle.end(), DwellRegions::PointAnglePairComparator());
  for (size_t i = 0; i < points_with_angle.size(); i++)
    oriented_points.push_back(points_with_angle[i].second);

  return oriented_points;
}

bool
is_outlier_trajectory(Point start_point_of_traj)
{
  bool outlier = false;
  if (std::isless((start_point_of_traj.Lat - Constants::boundary_min_lat), Constants::tolerance) ||
      std::isgreater((start_point_of_traj.Lat - Constants::boundary_max_lat), Constants::tolerance)) {
    outlier = true;
  }
  if (std::isless((start_point_of_traj.Lon - Constants::boundary_min_lon), Constants::tolerance) ||
      std::isgreater((start_point_of_traj.Lon - Constants::boundary_max_lon), Constants::tolerance)) {
    outlier = true;
  }
  return outlier;
}

// Returns the inner SEC (i.e., the under-circle SEC)
std::optional<Disk>
DwellRegionComputationBase::compute_inner_SEC(const traj_ptr_type& trajectory_ptr,
                                          const std::vector<MaxIDQueue<double>>& directional_heaps)
{
  // Boundary points are the ones that are subtended on unit vectors by points that have
  // maximum dot products for that direction.
  std::vector<Point> boundary_points;
  for (size_t i = 0; i < num_heaps; i++) {
    auto max_dot_product_point = trajectory_ptr->at(directional_heaps[i].peek().value().id);
    boundary_points.push_back(max_dot_product_point);
  }

  // Frontier points are the same as boundary points, but with duplicate points removed
  std::vector<Point> frontier_points(boundary_points.begin(), boundary_points.end());
  std::sort(frontier_points.begin(), frontier_points.end(), PointXYComparator());
  auto last_iter = std::unique(frontier_points.begin(), frontier_points.end(), PointXYEqualityComparator());
  frontier_points.erase(last_iter, frontier_points.end());

  return welzl(frontier_points, {});
  // return brute_force_SEC(frontier_points);
}

// Returns the outer SEC (i.e., the over-circle SEC)
std::optional<Disk>
DwellRegionComputationBase::compute_outer_SEC(const traj_ptr_type& trajectory_ptr,
                                          const std::vector<MaxIDQueue<double>>& directional_heaps)
{
  // Boundary points the furthest points in each directions
  std::vector<Point> boundary_points;
  for (size_t i = 0; i < num_heaps; i++) {
    auto max_dot_product_point = trajectory_ptr->at(directional_heaps[i].peek().value().id);
    boundary_points.push_back(max_dot_product_point);
  }
  // Draw the boundary lines perpendicular to the unit vectors through the frontier points
  std::vector<Line> boundary_lines;
  for (size_t i = 0; i < num_heaps; i++) {
    auto boundary_point = boundary_points[i];
    auto theta = i * angle_between_heaps;
    auto a = std::cos(theta);
    auto b = std::sin(theta);
    auto c = -dot_product_with_unit_vector(boundary_point, theta);

    // equation of line (intercept form) is ax + by + c = 0
    boundary_lines.emplace_back(a, b, c);
    // std::cout << "Boundary point[" << i << "]: " << boundary_point.x << " " << boundary_point.y
    // << std::endl;
    // std::cout << "Boundary line[" << i << "]: " << a << " " << b << " " << c << std::endl;
  }

  // Polygon points are where the boundary lines intersect,
  // i.e. they are the vertices of the minimum bounding polygon with k sides
  std::vector<Point> polygon_points;
  for (size_t i = 0; i < num_heaps; i++) {
    auto next_line_idx = (i + 1) % num_heaps;
    auto intersection_point = boundary_lines[i].intersects_with(boundary_lines[next_line_idx]);
    if (!intersection_point.has_value()) {
      std::cerr << "ERROR: Could not find the intersection between boundary lines " << i << " and " << next_line_idx
                << std::endl;
      exit(1);
    }

    polygon_points.emplace_back(intersection_point.value());
    // std::cout << "Polygon point[" << i << "]: " << intersection_point.value().x << ", "
    // << intersection_point.value().y << std::endl;
  }

  // Remove the duplicate polygon points
  std::sort(polygon_points.begin(), polygon_points.end(), PointXYComparator());
  auto last_iter = std::unique(polygon_points.begin(), polygon_points.end(), PointXYEqualityComparator());
  polygon_points.erase(last_iter, polygon_points.end());

  return welzl(polygon_points, {});
  // return brute_force_SEC(polygon_points);
}

// This is confirmUnconfirmed() in Reaz's code
// Read Section 5.A in the conference paper
std::pair<std::optional<Disk>, std::vector<Point>>
DwellRegionComputationBase::compute_SEC_S(const traj_ptr_type& trajectory_ptr,
                                      std::vector<MaxIDQueue<double>>& directional_heaps,
                                      const Disk& inner_circle)
{
  std::vector<double> threshold_dot_product_at_unit_vectors;
  threshold_dot_product_at_unit_vectors.reserve(num_heaps);
  for (size_t i = 0; i < num_heaps; i++) {
    auto bisector_angle = i * angle_between_heaps + 0.5 * angle_between_heaps;
    auto point = inner_circle.point_at_angle(bisector_angle);
    auto heap_angle = i * angle_between_heaps;
    auto dot_product = dot_product_with_unit_vector(point, heap_angle);
    threshold_dot_product_at_unit_vectors.emplace_back(dot_product);
  }

  // Ring points are the points on the convex hull.
  std::vector<Point> ring_points;
  ring_points.reserve(directional_heaps[0].size());

  std::vector<short> taken_ring_point_indexes(trajectory_ptr->size());

  // Trying to avoid too many malloc syscalls
  std::vector<IDKeyPair<double>> elements_removed_from_heap;
  elements_removed_from_heap.reserve(directional_heaps[0].size());

  for (size_t j = 0; j < num_heaps; j++) {
    size_t elements_removed_count = 0;
    while (!directional_heaps[j].empty()) {
      // Get the point with highest dot product in the direction of the unit vector
      auto furthest_point_id = directional_heaps[j].peek().value().id;
      auto furthest_point_dot_product = directional_heaps[j].peek().value().key;
      directional_heaps[j].pop();

      // If the current point hasn't been already added to the convex hull,
      // and if the dot product of it is greater than the dot product of frontier point and
      // the unit vector
      if (std::isgreaterequal((furthest_point_dot_product - threshold_dot_product_at_unit_vectors[j]),
                              Constants::tolerance) &&
          (taken_ring_point_indexes[furthest_point_id] == 0)) {
        ring_points.emplace_back(trajectory_ptr->at(furthest_point_id));
        taken_ring_point_indexes[furthest_point_id]=1; // Mark as taken
      }

      elements_removed_from_heap[elements_removed_count++] = (IDKeyPair<double>{ furthest_point_id, furthest_point_dot_product });

      if (std::isless((furthest_point_dot_product - threshold_dot_product_at_unit_vectors[j]), Constants::tolerance)) {
        break;
      }
    }

    for (size_t k=0; k < elements_removed_count; k++) {
      directional_heaps[j].push(elements_removed_from_heap[k]);
    }
  }

  std::sort(ring_points.begin(), ring_points.end(), PointXYComparator());
  auto ring_size = std::unique(ring_points.begin(), ring_points.end(), PointXYEqualityComparator());
  ring_points.erase(ring_size, ring_points.end());

  if (ring_points.size() > num_heaps) {
    return std::make_pair(welzl(ring_points, {}), ring_points);
    // return std::make_pair(brute_force_SEC(ring_points), ring_points);
  } else {
    return std::make_pair(inner_circle, ring_points);
  }
}

std::vector<Point>
DwellRegionComputationBase::trajectory_points_on_SEC(const Disk& SEC, const std::vector<Point>& ring_points)
{
  // The subset of ring points those are on the SEC disk,
  // and contribute to form the dwell region.
  std::vector<Point> traj_points_on_disk;
  for (size_t k = 0; k < ring_points.size(); k++) {
    if (std::fabs((Util::Geo::distance_in_miles_from_X_Y(SEC.get_center(), ring_points[k]) - SEC.get_radius())) <
        Constants::tolerance) {
      traj_points_on_disk.push_back(ring_points[k]);
    }
  }

  return traj_points_on_disk;
}

// Returns the points to form the dwell region
// i.e., the points are the vertices of the dwell region.
// The points are in counterclockwise orientation
std::vector<Point>
DwellRegionComputationBase::get_approximate_dwell_region(const double& query_radius, const std::vector<Point>& points)
{
  std::vector<Disk> circles;
  for (auto point : points) {
    Disk circle(point, query_radius);
    circles.push_back(circle);
  }

  // If there are only two circles, there will be only two intersection points
  if (circles.size() == 2) {
    auto result = circles[0].disk_disk_intersection(circles[1]);
    if (!result.has_value()) {
      std::cerr << "ERROR: Could not compute disk-disk intersection. Dwell Region computation aborted." << std::endl;
      exit(1);
    }
    return result.value();
  }

  std::vector<Point> dwell_region_points;
  // nC2 circles intersection points
  for (size_t i = 0; i < circles.size(); i++) {
    for (size_t j = i + 1; j < circles.size(); j++) {
      auto intersection_points = circles[i].disk_disk_intersection(circles[j]).value();
      for (size_t k = 0; k < intersection_points.size(); k++) {
        // The intersection point that contributes to form the dwell region
        // should be inside the other circles
        bool inside_point = true;
        for (size_t l = 0; l < circles.size(); l++) {
          if (l == i || l == j)
            continue;
          if (!circles[l].in_disk(intersection_points[k])) {
            inside_point = false;
            break;
          }
        }
        if (inside_point)
          dwell_region_points.push_back(intersection_points[k]);
      }
    }
  }
  // Dwell region points in a counterclockwise orientation
  return oriented_counterclockwise(dwell_region_points);
}

} // namespace DwellRegions
