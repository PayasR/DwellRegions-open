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


// The radii bin entries, i.e., the subtrajectory entries in the Rho-index
struct RadiiBinEntry
{
  RadiiBinEntry(size_t traj_id,
                size_t start_idx,
                size_t end_idx,
                size_t duration,
                size_t next_duration,
                size_t next_last_idx)
    : trajectory_id(traj_id)
    , first_idx_subtraj(start_idx)
    , last_idx_subtraj(end_idx)
    , duration_subtraj(duration)
    , duration_of_subtraj_in_next_bin_entry(next_duration)
    , last_idx_of_subtraj_in_next_bin_entry(next_last_idx)
  {}

  size_t trajectory_id;
  size_t first_idx_subtraj;
  size_t last_idx_subtraj;
  // Duration of the subtrajectory in this radii_bin entry
  size_t duration_subtraj;

  // Duration of the subtrajectory in next radii_bin with same start index,
  // this duration is used to extend the subtrajectory in this radii_bin entry
  size_t duration_of_subtraj_in_next_bin_entry;
  // Last index of the subtrajectory in next radii_bin with same start index,
  // this duration is used to extend the subtrajectory in this radii_bin entry
  size_t last_idx_of_subtraj_in_next_bin_entry;
};

// Sort subtrajectories based on duration in increasing order,
// solve ties with start index of subtrajectories (Rho-index)
struct RadiiBinEntryComparator
{
  bool operator()(const RadiiBinEntry& entry1, const RadiiBinEntry& entry2)
  {
    // Compare the duration of the subtraj in radii_bin entries
    if (entry1.duration_subtraj != entry2.duration_subtraj)
      return entry1.duration_subtraj < entry2.duration_subtraj;

    // Compare the trajectory ids if the duration are equal
    if (entry1.trajectory_id != entry2.trajectory_id)
      return entry1.trajectory_id < entry2.trajectory_id;

    // Compare the start indexes if the duration and trajectory ids are equal
    return entry1.first_idx_subtraj < entry2.first_idx_subtraj;
  }
};


// The cell entries, i.e., the subtrajectory entries in the cells in the Tau-Index
struct TauCellEntry
{
  TauCellEntry(size_t traj_id,
               size_t start_idx,
               size_t end_idx,
               size_t duration,
               double radius,
               double center_x,
               double center_y,
               size_t first_time,
               size_t end_time)
    : trajectory_id(traj_id)
    , first_idx_subtraj(start_idx)
    , last_idx_subtraj(end_idx)
    , duration_subtraj(duration)
    , SEC_S_radius(radius)
    , SEC_S_center_x(center_x)
    , SEC_S_center_y(center_y)
    , first_idx_timestamp(first_time)
    , last_idx_timestamp(end_time)
  {}

  size_t trajectory_id;
  size_t first_idx_subtraj;
  size_t last_idx_subtraj;
  // Duration of the subtrajectory in this temporal_bin entry
  size_t duration_subtraj;
  // SEC_S disk radius of the subtrajectory
  double SEC_S_radius;
  // SEC_S disk center (x, y) of the subtrajectory
  double SEC_S_center_x;
  double SEC_S_center_y;

  // First and last timestamps are not stored in the Tau-index,
  // they are required to build the index only
  size_t first_idx_timestamp;
  size_t last_idx_timestamp;
};

// Sort subtrajectories based on duration in increasing order,
// solve ties with start index of subtrajectories (Tau-Index)
struct TauCellEntryComparator
{
  bool operator()(const TauCellEntry& entry1, const TauCellEntry& entry2)
  { 
    // Compare the disk radius of the subtraj in temporal_bin entries
    if (std::isgreater(std::fabs(entry1.SEC_S_radius - entry2.SEC_S_radius), 0.000001))
      return entry1.SEC_S_radius < entry2.SEC_S_radius;

    // Compare the trajectory ids if the radius are equal
    if (entry1.trajectory_id != entry2.trajectory_id)
      return entry1.trajectory_id < entry2.trajectory_id;

    // Compare the start indexes if the radius and ids are equal
    if (entry1.first_idx_subtraj != entry2.first_idx_subtraj)
      return entry1.first_idx_subtraj < entry2.first_idx_subtraj;
    
    // Compare the last indexes if the radius, id, first_index are equal
    if (entry1.last_idx_subtraj != entry2.last_idx_subtraj)
      return entry1.last_idx_subtraj < entry2.last_idx_subtraj;
    
    // Compare the duration of subtrajs (default)
    return entry1.duration_subtraj < entry2.duration_subtraj;
  }
};

}