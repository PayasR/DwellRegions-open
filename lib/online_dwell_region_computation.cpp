#include "../include/disk.hpp"
#include "../include/online_dwell_region_computation.hpp"
#include "../include/line.hpp"
#include "../include/util.hpp"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <optional>
#include <set>

namespace DwellRegions {

bool
OnlineDwellRegionComputation::has_dwell_region()
{
  if (current_trajectory_idx < trajectory_ptr->size())
    process_points_and_set_dwell_region_flag();

  return dwell_region_flag;
}

// This function sets the dwell_region_flag to true or false depending on whether the trajectory has
// a dwell region or not.
// It processes only the newly added points in the trajectory (points in the range
// [current_trajectory_idx...last_point_index]) Therefore, it acts as both Algorithms 2 and 3 in the
// conference paper. This is queryEvalBF() in Reaz's code
void
OnlineDwellRegionComputation::process_points_and_set_dwell_region_flag()
{
  // At the beginning, there is no dwell region
  // which can be set true, if at least a dwell region exists
  dwell_region_flag = false;
  // Each iteration of the loop processes one point previously unprocessed in the trajectory
  for (; current_trajectory_idx < trajectory_ptr->size(); current_trajectory_idx++) {

    // ------------------ STEP 1: Add point from the trajectory to the heaps -----------------------
    auto curr_point = trajectory_ptr->at(current_trajectory_idx);

    // Iterate over all the heaps and add the dot products of their directions and Points to them
    auto start = clock();
    for (size_t i = 0; i < num_heaps; i++) {
      auto theta = i * angle_between_heaps;
      auto dot_product = dot_product_with_unit_vector(curr_point, theta);
      directional_heaps[i].push(IDKeyPair<double>{ current_trajectory_idx, dot_product });
    }
    auto end = clock();
    query_perf_counters.heap_update_time_in_seconds += (end - start) / (double)CLOCKS_PER_SEC;
    current_window_size++;

    // Add the difference between timestamps of current point and the previous point
    // to compute the current time window
    if (current_window_size > 1) {
      size_t prev_timestamp = trajectory_ptr->timestamp_at(current_trajectory_idx - 1);
      size_t curr_timestamp = trajectory_ptr->timestamp_at(current_trajectory_idx);
      current_time_window += (curr_timestamp - prev_timestamp);
    }

    // Do we have enough points in the window to compute a dwell region?
    // If not, get another point from the trajectory
    if (current_time_window < max_window_size) {
      continue;
    }

    // ------------------ STEP 2: Compute the lower bound of radius -----------------------
    double R_SEC_lower_bound = 0.0, R_SEC_upper_bound = 0.0;

    // Compute the lower bound of SEC_S
    auto SEC_inner = compute_inner_SEC(trajectory_ptr, directional_heaps);
    if (!SEC_inner.has_value()) {
      std::cerr << "ERROR: Welzl's algorithm failed. Could not compute the inner SEC." << std::endl;
      exit(1);
    }
    R_SEC_lower_bound = SEC_inner.value().get_radius();

    // ------------------- STEP 3: Compute the upper bound of radius ----------------------

    // Compute the upper bound of SEC_S
    auto SEC_outer = compute_outer_SEC(trajectory_ptr, directional_heaps);
    if (!SEC_outer.has_value()) {
      std::cerr << "ERROR: Welzl's algorithm failed. Could not compute the outer SEC" << std::endl;
      exit(1);
    }
    R_SEC_upper_bound = SEC_outer.value().get_radius();

    // Sanity Check
    if (std::isgreater((R_SEC_lower_bound - R_SEC_upper_bound), Constants::tolerance)) {
      std::cerr << "ERROR: LOWER BOUND > UPPER BOUND!!" << std::endl;
      exit(1);
    }
    if (R_SEC_upper_bound < 0.0) {
      std::cout << "Upper bound negative. R = " << R_SEC_upper_bound << std::endl;
    }

    // Is the window full?
    // If yes, remove the first point in the window to make space for the next point
    if (std::isgreaterequal((current_time_window - max_window_size), Constants::tolerance)) {
      // ---- CASE 1: If the query radius > upper bound
      // we know that dwell region exists and there's no need to compute the actual SEC
      if (std::isgreater((query_radius - R_SEC_upper_bound), Constants::tolerance)) {

        // Get the first point in the window, create a subtrajectory using the current window and
        // add it to result
        auto first_point_in_window_idx = (current_trajectory_idx - current_window_size + 1);
        sub_trajectories_that_form_dwell_regions.emplace_back(
          *trajectory_ptr, first_point_in_window_idx, current_window_size);

        dwell_region_flag = true;
        query_perf_counters.num_query_answered_without_computing_SEC++;

        // Compute the actual SEC_S, and the approximate dwell region using points on the SEC_S
        auto SEC_ring_points_pair = compute_SEC_S(trajectory_ptr, directional_heaps, SEC_inner.value());
        auto SEC_S = SEC_ring_points_pair.first;
        auto ring_points = SEC_ring_points_pair.second;
        auto dwell_region =
          get_approximate_dwell_region(query_radius, trajectory_points_on_SEC(SEC_S.value(), ring_points));

        // std::cout << "Dwell region points in counterclockwise orientation: " << std::endl;
        // for (size_t t = 0; t < dwell_region.size(); t++)
        //   std::cout << "Dwell region point[" << t << "]: "
        //             << dwell_region[t].x << " " << dwell_region[t].y << std::endl;
      }

      // ---- CASE 2: If the query radius < lower bound
      // The dwell region does not exist
      else if (std::isgreater((R_SEC_lower_bound - query_radius), Constants::tolerance)) {
        // dwell_region_flag = false;
        query_perf_counters.num_query_answered_without_computing_SEC++;
        // std::cout << "NO: Query Radius too small, no SEC computation needed." << std::endl;
      }

      // ---- CASE 3: If the lower bound < query radius < upper bound
      // Compute the actual SEC_S, and make the decision on the basis of it
      else if (std::isgreaterequal((query_radius - R_SEC_lower_bound), Constants::tolerance) &&
               std::islessequal((query_radius - R_SEC_upper_bound), -Constants::tolerance)) {

        query_perf_counters.num_SEC_computed++;
        auto SEC_ring_points_pair = compute_SEC_S(trajectory_ptr, directional_heaps, SEC_inner.value());
        auto SEC_S = SEC_ring_points_pair.first;
        auto ring_points = SEC_ring_points_pair.second;
        query_perf_counters.num_filtered_points += ring_points.size();

        if (!SEC_S.has_value()) {
          std::cerr << "ERROR: Could not compute SEC_S. Exiting..." << std::endl;
          exit(1);
        }

        auto SEC_S_disk = SEC_S.value();

        // Query radius is greater than the radius of the disk or approximately equal to
        // disk, then we have a dwell region, no otherwise.
        if (std::isgreaterequal((query_radius - SEC_S_disk.get_radius()), Constants::tolerance)) {
          // Get the first point in the window, create a subtrajectory using the current window and
          // add it to result
          auto first_point_in_window_idx = (current_trajectory_idx - current_window_size + 1);
          sub_trajectories_that_form_dwell_regions.emplace_back(
            *trajectory_ptr, first_point_in_window_idx, current_window_size);

          dwell_region_flag = true;

          auto dwell_region =
            get_approximate_dwell_region(query_radius, trajectory_points_on_SEC(SEC_S_disk, ring_points));

          // std::cout << "Dwell region points in counterclockwise: " << std::endl;
          // for (size_t t = 0; t < dwell_region.size(); t++)
          // std::cout << "Dwell region point[" << t << "]: "
          // << dwell_region[t].x << " " << dwell_region[t].y << std::endl;
        } else {
          // dwell_region_flag = false;
          // std::cout << "NO: Computed SEC, but no dwell region found." << std::endl;
        }
      }

      // Remove the first element in the window from all the directional heaps (Algorithm 3 in
      // the paper)
      auto first_point_in_window_idx = (current_trajectory_idx - current_window_size + 1);
      auto start = clock();
      for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
        directional_heaps[heap_idx].remove_id(first_point_in_window_idx);
      }
      auto end = clock();
      query_perf_counters.heap_update_time_in_seconds += (end - start) / (double)CLOCKS_PER_SEC;

      current_window_size--;
      size_t first_point_time_window = (trajectory_ptr->timestamp_at(first_point_in_window_idx + 1) -
                                        trajectory_ptr->timestamp_at(first_point_in_window_idx));
      current_time_window -= first_point_time_window;
    }
  }
}
} // namespace DwellRegions
