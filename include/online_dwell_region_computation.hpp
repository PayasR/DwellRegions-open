#pragma once

#include "dwell_region_computation_base.hpp"

namespace DwellRegions {

class OnlineDwellRegionComputation : DwellRegionComputationBase
{
public:
  OnlineDwellRegionComputation(DwellRegionComputationBase::traj_ptr_type T,
                               unsigned num_heaps,
                               unsigned max_window_size,
                               double query_radius)
    : DwellRegionComputationBase(num_heaps)
    , trajectory_ptr(T)
    , max_window_size(max_window_size)
    , query_radius(query_radius)
  {
    for (size_t i = 0; i < num_heaps; i++) {
      directional_heaps.emplace_back(MaxIDQueue<double>(trajectory_ptr->size()));
    }
  }

  struct QueryPerformanceCounters
  {
    unsigned num_SEC_computed = 0;
    unsigned num_query_answered_without_computing_SEC = 0;
    size_t num_filtered_points = 0;
    double heap_update_time_in_seconds = 0.0;

    [[nodiscard]]
    std::string get_header_line() const
    {
      return "num_SEC_computed, num_query_answered_without_computing_SEC, num_filtered_points, "
             "heap_update_time_in_seconds";
    }

    friend std::ostream& operator<<(std::ostream& os, const QueryPerformanceCounters& obj)
    {
      os << obj.num_SEC_computed << ", " << obj.num_query_answered_without_computing_SEC << ", "
         << obj.num_filtered_points << ", " << obj.heap_update_time_in_seconds;
      return os;
    }
  };

  bool has_dwell_region();
  [[nodiscard]] auto num_dwell_regions() const { return sub_trajectories_that_form_dwell_regions.size(); }

  // std::vector<Region> get_dwell_regions() {}
  [[nodiscard]] auto get_subtrajectories_that_form_dwell_regions() const
  {
    return sub_trajectories_that_form_dwell_regions;
  }

  auto get_query_performance_counters() const { return query_perf_counters; }

protected:
  void process_points_and_set_dwell_region_flag();

  const traj_ptr_type trajectory_ptr;
  const size_t max_window_size;
  const double query_radius;
  
  size_t current_time_window = 0;    // Time window starting from the start point to the current point in the trajectory
  size_t current_window_size = 0;    // Current number of points/locations of the object in the trajectory
  size_t current_trajectory_idx = 0; // For all points [0..current_trajectory_idx] in the trajectory we
                                     // have already tried to see if a dwell region exists for the object.
                                     // For all the points from current_trajectory_idx + 1 to the end of
                                     // the trajectory, call process_points_and_set_dwell_region_flag()

  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction

  bool dwell_region_flag = false;
  std::vector<SubTrajectory> sub_trajectories_that_form_dwell_regions;

  QueryPerformanceCounters query_perf_counters;
};
}