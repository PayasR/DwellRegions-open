#pragma once

#include "disk.hpp"
#include "id_queue.hpp"
#include "region.hpp"
#include "trajectory.hpp"
#include "gsl-lite.hpp"

namespace DwellRegions {

// This class contains utility functions for all the online/offline versions of the DwellRegionComputation
// Note 1: LIFETIME: Please mind the lifetime of a DwellRegionComputation object.
// Since it holds a constant pointer to a Trajectory object, it's lifetime is dictated by that of
// the Trajectory object's
//
// Note 2: ASSUMPTION: The heaps representing the frontier points in different directions are all
// equally spaced around the PI radians in a circle. In other words, if num_heaps heaps are used,
// then the nth heap represents the (n-1)*(2*PI/num_heaps) radians from the positive x-axis:
//            0     0        The second heap represents this direction, at (2*PI/num_heaps) degrees
//       0            /  0
//                   /
//   0              /        0
//                 /
// 0              /           0
//               O----------->  The first heap starts from here, 0 degrees
// 0                          0
//
//  0                        0
//
//     0                  0
//          0        0
//
class DwellRegionComputationBase
{
public:
  using traj_ptr_type = gsl::not_null<Trajectory*>;

  DwellRegionComputationBase(size_t num_heaps)
    : num_heaps(num_heaps)
    , angle_between_heaps(2 * Constants::PI / num_heaps)
  {}

protected:
  // Returns the inner SEC (i.e., the under-circle SEC)
  virtual std::optional<Disk> compute_inner_SEC(const traj_ptr_type& trajectory_ptr,
                                        const std::vector<MaxIDQueue<double>>& directional_heaps);

  // Returns the outer SEC (i.e., the over-circle SEC)
  virtual std::optional<Disk> compute_outer_SEC(const traj_ptr_type& trajectory_ptr,
                                        const std::vector<MaxIDQueue<double>>& directional_heaps);

  // Returns the SEC_S that contains all the points of trajectory points set S
  // and the ring points that are on the convex hull of S
  virtual std::pair<std::optional<Disk>, std::vector<Point>> compute_SEC_S(const traj_ptr_type& trajectory_ptr,
                                                                   std::vector<MaxIDQueue<double>>& directional_heaps,
                                                                   const Disk& inner_circle);

  // Returns the trajectory points on the actual SEC_S,
  // the returned vector of trajectory points is a subset of ring points
  virtual std::vector<Point> trajectory_points_on_SEC(const Disk& SEC, const std::vector<Point>& ring_points);

  // Returns the approximate dwell region (R+ in the paper)
  virtual std::vector<Point> get_approximate_dwell_region(const double& query_radius, const std::vector<Point>& points);

  const size_t num_heaps;
  const double angle_between_heaps;
};

double
dot_product_with_unit_vector(Point point, double theta);

std::vector<Point>
oriented_counterclockwise(const std::vector<Point>& points);

bool
is_outlier_trajectory(Point start_point_of_traj);
}