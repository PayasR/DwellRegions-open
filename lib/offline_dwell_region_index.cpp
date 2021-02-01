#include "../include/disk.hpp"
#include "../include/offline_dwell_region_computation.hpp"
#include "../include/line.hpp"
#include "../include/util.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <optional>
#include <set>
#include <string>

namespace DwellRegions {

// Prepare and setup the bins of different radii
// based on the number of radii bins
void
OfflineDwellRegionIndex::prepare_radii_bins()
{
  double width_of_radii_bins = (max_radius - min_radius) / (num_radii_bins - 1);
  double radius = min_radius;
  while (radius <= max_radius) {
    radii_bins.emplace_back(radius);
    radius += width_of_radii_bins;
    std::vector<RadiiBinEntry> current_bin;
    radii_bins_with_subtrajectories.emplace_back(current_bin);
  }
}

// Process the trajectory and compute the radii bins
// which will be used to compute the dwell regions
// for the dwell region query
void
OfflineDwellRegionIndex::process_trajectory_and_compute_bins(
  size_t trajectory_id,
  size_t radii_idx,
  std::vector<std::vector<RadiiBinEntry>>& current_bins)
{
  double current_radius = radii_bins[radii_idx];
  size_t last_idx_of_last_bin_entry = 0;

  // Move the start point_idx of each subtrajectory over the trajectory points
  for (size_t first_point_idx = 0; first_point_idx < (trajectory_ptr->size() - 1); first_point_idx++) {
    // Initialize variables for new subtrajectory of the trajectory
    size_t current_window_size = 0;
    size_t current_time_window = 0;
    for (auto& heap : directional_heaps) {
      heap.clear();
    }

    // Each iteration of the loop processes one point previously unprocessed in the subtrajectory
    size_t current_point_idx = first_point_idx;
    for (; current_point_idx < trajectory_ptr->size(); current_point_idx++) {
      // ------------------ STEP 1: Add point from the trajectory to the heaps -----------------------
      auto curr_point = trajectory_ptr->at(current_point_idx);
      for (size_t i = 0; i < num_heaps; i++) {
        auto theta = i * angle_between_heaps;
        auto dot_product = dot_product_with_unit_vector(curr_point, theta);
        directional_heaps[i].push(IDKeyPair<double>{ current_point_idx, dot_product });
      }
      current_window_size++;

      // Add the difference between timestamps of current point and the previous point
      // to compute the current time window of the current subtrajectory
      if (current_point_idx > first_point_idx) {
        size_t prev_timestamp = trajectory_ptr->timestamp_at(current_point_idx - 1);
        size_t curr_timestamp = trajectory_ptr->timestamp_at(current_point_idx);
        current_time_window += (curr_timestamp - prev_timestamp);
      }

      // At least three points are required in the window to compute a dwell region,
      // the last added point can be discarded for new radii_bin_entry
      if (current_window_size < 3) {
        continue;
      }

      // For next eligible valid radii_bin_entry, current_point_idx should be
      // greater than the last index of the last radii_bin_entry
      if (current_point_idx <= last_idx_of_last_bin_entry) {
        continue;
      }

      // -----------------STEP 2: Compute the lower and upper bounds of radius ---------------------
      //-------------------   for current subtrajectory with current time window  ------------------
      double R_SEC_lower_bound = 0.0, R_SEC_upper_bound = 0.0;

      // Compute the lower bound of SEC_S
      auto SEC_inner = compute_inner_SEC(trajectory_ptr, directional_heaps);
      if (!SEC_inner.has_value()) {
        std::cerr << "ERROR: Welzl's algorithm failed. Could not compute the inner SEC." << std::endl;
        exit(1);
      }
      R_SEC_lower_bound = SEC_inner.value().get_radius();

      // Compute the upper bound of SEC_S
      auto SEC_outer = compute_outer_SEC(trajectory_ptr, directional_heaps);
      if (!SEC_outer.has_value()) {
        std::cerr << "ERROR: Welzl's algorithm failed. Could not compute the outer SEC" << std::endl;
        exit(1);
      }
      R_SEC_upper_bound = SEC_outer.value().get_radius();

      // Sanity Checks
      if (std::isgreater((R_SEC_lower_bound - R_SEC_upper_bound), Constants::tolerance)) {
        std::cerr << "ERROR: LOWER BOUND > UPPER BOUND!!" << std::endl;
        exit(1);
      }
      if (R_SEC_upper_bound < 0.0) {
        std::cout << "Upper bound negative. R = " << R_SEC_upper_bound << std::endl;
      }

      // -------------- Step 3: Make an entry of a current_radius-maximal subtrajectory--------------
      // Is upper bound > current radius? If yes, compute actual SEC_S;
      // and check SEC_S > current radius Otherwise, continue adding points from trajectory

      // TODO: Argue about using the using the tolerances everywhere for floating point comparisons
      // if (std::isgreater((R_SEC_upper_bound - current_radius), Constants::tolerance))
      if (R_SEC_upper_bound > current_radius) // don't change this condition
      {
        auto SEC_ring_points_pair = compute_SEC_S(trajectory_ptr, directional_heaps, SEC_inner.value());

        auto SEC_S = SEC_ring_points_pair.first;
        if (!SEC_S.has_value()) {
          std::cerr << "ERROR: Could not compute SEC_S. Exiting..." << std::endl;
          exit(1);
        }
        auto SEC_S_disk = SEC_S.value();

        // If the radius of SEC_S is greater than current radius,
        // discard the last added point, and make an entry in radii bin for the subtrajectory,
        // if (std::isgreater((SEC_S_disk.get_radius() - current_radius), Constants::tolerance))
        if (SEC_S_disk.get_radius() > current_radius) // don't change this condition
        {
          // Discard the last added point from the current window for the next radii_bin entry
          current_window_size--;
          size_t last_point_time_window =
            (trajectory_ptr->timestamp_at(current_point_idx) - trajectory_ptr->timestamp_at(current_point_idx - 1));

          for (auto& heap : directional_heaps) {
            heap.remove_id(current_point_idx);
          }

          // Check whether the subtrajectory is current_radius-maximal or not
          // If yes, add it to the current radii_bin
          bool is_subtraj_current_radius_maximal = false;
          auto SEC_outer_subtraj = compute_outer_SEC(trajectory_ptr, directional_heaps);
          double SEC_upper_bound_subtraj = SEC_outer_subtraj.value().get_radius();
          // if (std::islessequal((SEC_upper_bound_subtraj - current_radius), -Constants::tolerance))
          if (SEC_upper_bound_subtraj <= current_radius) // don't change this condition
          {
            is_subtraj_current_radius_maximal = true;
            // std::cout << "Upper bound of subtrajectory: " << SEC_upper_bound_subtraj <<
            // std::endl;
          } else {
            auto SEC_inner_subtraj = compute_inner_SEC(trajectory_ptr, directional_heaps);
            auto SEC_S_subtraj = compute_SEC_S(trajectory_ptr, directional_heaps, SEC_inner_subtraj.value());
            double SEC_S_disk_subtraj = SEC_S_subtraj.first.value().get_radius();
            // if (std::islessequal((SEC_S_disk_subtraj - current_radius), -Constants::tolerance)) {
            if (SEC_S_disk_subtraj <= current_radius) { // don't change this condition
              is_subtraj_current_radius_maximal = true;
              // std::cout << "Lower bound of subtrajectory: " << SEC_lower_bound_subtraj <<
              // std::endl;
            }
          }

          if (is_subtraj_current_radius_maximal && last_idx_of_last_bin_entry < (current_point_idx - 1)) {
            // last index and the time window of the entry subtrajectory
            auto last_point_in_window_idx = current_point_idx - 1;
            auto duration = (current_time_window - last_point_time_window);

            // Make an entry in radii bin for the subtrajectory
            RadiiBinEntry current_bin_entry(trajectory_id,
                                            first_point_idx,
                                            last_point_in_window_idx,
                                            duration,
                                            0, last_point_in_window_idx); // default next duration = 0,
                                                                       // default next last_idx = current last_index
            // Push the new radii_bin entry and update the last index of the last radii_bin entry
            current_bins[radii_idx].emplace_back(current_bin_entry);
            last_idx_of_last_bin_entry = last_point_in_window_idx;
            /*std::cout << "Trajectory id: " << trajectory_id
                      << ", First: " << first_point_idx << ", Last: " << last_point_in_window_idx
                      << ", Duration: " << duration << std::endl;
            std::cout << "Current bin size: " << current_bins[radii_idx].size() << std::endl;*/
          }
          break;
        }
      }
      // If lower bound <= current radius <= upper bound,
      // continue adding points from trajectory
    }
  }
}

// Set the current trajectory, preprocess it, compute the bins for it
// and update the radii_bins that contains subtrajectories from all trajectories
void
OfflineDwellRegionIndex::set_trajectory_and_update_radii_bins(Trajectory* traj, const size_t traj_idx)
{
  // Set the current trajectory
  trajectory_ptr = traj;

  //---------- Step 1: Initialize the directional heaps for the current trajectory-------
  if (traj_idx == 0) {
    for (size_t i = 0; i < num_heaps; i++) {
      directional_heaps.emplace_back(MaxIDQueue<double>(trajectory_ptr->size()));
    }
  } else {
    // Make sure there's enough space in the directional heaps to accommodate all points in the
    // current traj
    if (directional_heaps[0].size() < trajectory_ptr->size()) {
      directional_heaps.clear();
      for (size_t i = 0; i < num_heaps; i++) {
        directional_heaps.emplace_back(MaxIDQueue<double>(trajectory_ptr->size()));
      }
    }
  }

  //---------- Step 2: Preprocess the current trajectory to compute the radii bins--------
  std::vector<std::vector<RadiiBinEntry>> current_bins(radii_bins.size());
  for (size_t idx = 0; idx < radii_bins.size(); idx++) {
    process_trajectory_and_compute_bins(traj_idx, idx, current_bins);

    // Reset all the directional heaps to prepare them form next radius of radii_bin
    for (auto& heap : directional_heaps) {
      heap.clear();
    }
  }

  // Update each bin entries with the duration of
  // valid entries of the next_bin with same start index
  for (size_t i = 0; i < radii_bins.size() - 1; i++) {
    size_t j = 0, k = 0;
    while (j < current_bins[i].size() && k < current_bins[i + 1].size()) {
      auto current_bin_start_idx = current_bins[i][j].first_idx_subtraj;
      auto next_bin_start_idx = current_bins[i + 1][k].first_idx_subtraj;
      if (current_bin_start_idx == next_bin_start_idx) {
        auto next_duration = current_bins[i + 1][k].duration_subtraj;
        current_bins[i][j].duration_of_subtraj_in_next_bin_entry = next_duration;
        auto next_last_idx = current_bins[i + 1][k].last_idx_subtraj;
        current_bins[i][j].last_idx_of_subtraj_in_next_bin_entry = next_last_idx;
        j++;
        k++;
      } else if (current_bin_start_idx < next_bin_start_idx)
        j++;
      else
        k++;
    }
  }

  //----------- Step 3: Insert the computed radii_bin entries for current trajectory-----------
  for (size_t i = 0; i < radii_bins.size(); i++) {
    radii_bins_with_subtrajectories[i].insert(
      radii_bins_with_subtrajectories[i].end(), current_bins[i].begin(), current_bins[i].end());
  }
}

// Preprocess all the trajectories in the full dataset
// and compute the radii_bins for all trajectories
void
OfflineDwellRegionIndex::preprocess_trajectory_dataset(std::filesystem::path filename, std::string data_path)
{
  // Read Points from the file and add_point() to traj here...
  std::ifstream trajFile(filename, std::ios::in);
  std::string line;

  if (!trajFile.is_open()) {
    std::cerr << "ERROR: Could not open file. Exiting." << std::endl;
    exit(1);
  }

  auto start = clock();
  std::getline(trajFile, line);
  size_t num_trajectories = std::stoi(line);

  for (size_t traj_idx = 0; traj_idx < num_trajectories; traj_idx++) {
    // Read out the "Reading FILE.log"
    std::getline(trajFile, line);

    // Read the trajectory in sstream
    std::getline(trajFile, line);
    std::stringstream ss;
    ss << line;

    size_t trajectory_ID, num_points_in_trajectory;
    ss >> trajectory_ID >> num_points_in_trajectory;

    /*std::string filestring = data_path + "GeoLifeDataByID/" + std::to_string(traj_idx) + ".txt";
    auto file = std::filesystem::path(filestring);
    std::ofstream trajByIdFile(file, std::ios::out);
    trajByIdFile << num_points_in_trajectory << std::endl;*/

    // Populate the trajectory object with points from the file
    Trajectory traj;
    while (ss.good()) {
      double lon = 0.0, lat = 0.0;
      size_t timestamp = 0;
      ss >> lat >> lon >> timestamp;
      /*if (timestamp > 0) {
        trajByIdFile << lat << "\t" << lon << "\t" << timestamp << std::endl;
      }*/
      traj.emplace_back(Point_from_lon_lat(lon, lat), timestamp);
    }
    // trajByIdFile.close();

    if (is_outlier_trajectory(traj.at(0))) {
      continue;
    }
    set_trajectory_and_update_radii_bins(&traj, traj_idx);
    if (traj_idx == 8)
      break;
  }

  // Sort the subtrajectories in radii_bins by duration/ time_window
  // (according to paper, sort by tau_i)
  for (size_t idx = 0; idx < radii_bins.size(); idx++) {
    std::sort(radii_bins_with_subtrajectories[idx].begin(),
              radii_bins_with_subtrajectories[idx].end(),
              RadiiBinEntryComparator());

    /*std::cout << "------------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "Current radii bin: " << radii_bins[idx] << ", Size: " << radii_bins_with_subtrajectories[idx].size()
              << std::endl;
    for (size_t j = 0; j < radii_bins_with_subtrajectories[idx].size(); j++) {
      std::cout << "Trajectory id: " << radii_bins_with_subtrajectories[idx][j].trajectory_id
                << ", First: " << radii_bins_with_subtrajectories[idx][j].first_idx_subtraj
                << ", Last: " << radii_bins_with_subtrajectories[idx][j].last_idx_subtraj
                << ", Duration: " << radii_bins_with_subtrajectories[idx][j].duration_subtraj
                << ", Next duration: " << radii_bins_with_subtrajectories[idx][j].duration_of_subtraj_in_next_bin_entry
                << ", Next last index: "
                << radii_bins_with_subtrajectories[idx][j].last_idx_of_subtraj_in_next_bin_entry << std::endl;
    }*/

    int current_radius_bin = (int) (radii_bins[idx] * 10.0 + 0.5);
    std::string index_string = data_path + "Index/" + std::to_string(num_radii_bins) + "_bins/" +
                             "radius_" + std::to_string(current_radius_bin) +
                             "_k_" + std::to_string(num_heaps) + ".txt";
    std::cout << "Writing index file: " << index_string << std::endl;
    auto index_filename = std::filesystem::path(index_string);
    std::ofstream indexFile(index_filename, std::ios::out);
    indexFile << "traj_id" << "\t" << "start_idx" << "\t" << "last_idx" << "\t"
              << "duration" << "\t" << "next_last_idx" << "\t" << "next_duration" << std::endl;
    
    for (auto result : radii_bins_with_subtrajectories[idx]) {
      indexFile << result.trajectory_id << "\t" << result.first_idx_subtraj << "\t"
                << result.last_idx_subtraj << "\t" << result.duration_subtraj << "\t"
                << result.last_idx_of_subtraj_in_next_bin_entry << "\t"
                << result.duration_of_subtraj_in_next_bin_entry << std::endl;
    }
    indexFile.close();
  }
  auto end = clock();
  preprocess_time = (end - start) / (double) CLOCKS_PER_SEC;
  //std::cout << "Preprocessing time: " << preprocess_time << " seconds" << std::endl;
}
} // namespace DwellRegions
