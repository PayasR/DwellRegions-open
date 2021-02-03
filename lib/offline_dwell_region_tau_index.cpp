#include "../include/disk.hpp"
#include "../include/line.hpp"
#include "../include/offline_dwell_region_computation.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

namespace DwellRegions {

// Prepare and setup the bins of different time window
// based on the number of temporal bins
// For EXPONENTIAL partitioning scheme
void
OfflineDwellRegionTauIndex::prepare_temporal_bins()
{
  temporal_bins_with_radius.reserve(num_temporal_bins);

  size_t idx = 0;
  while (idx < num_temporal_bins) {
    std::vector<std::vector<TauCellEntry>> current_bin;
    current_bin.reserve(num_radii); // for radii = 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0
    for (size_t i = 0; i < num_radii; i++) {
      std::vector<TauCellEntry> temp;
      temp.reserve(30000); // default bin size
      current_bin.push_back(temp);
    }
    temporal_bins_with_radius.push_back(current_bin);
    idx++;
  }
}

/* // FOR LINEAR partitioning scheme
void
OfflineDwellRegionIndex::prepare_temporal_bins()
{
  size_t width_of_time_windows = (max_window - min_window) / (num_temporal_bins - 1);
  size_t time_window = min_window;
  time_windows.reserve(num_temporal_bins + 1);
  temporal_bins_with_radius.reserve(num_temporal_bins);

  size_t idx = 0;
  while (idx < num_temporal_bins) {
    time_windows.emplace_back(time_window);
    time_window += width_of_time_windows;
    std::vector<std::vector<TauCellEntry>> current_bin;
    current_bin.reserve(num_radii); // for radii = 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0
    for (size_t i = 0; i < num_radii; i++) {
      std::vector<TauCellEntry> temp;
      temp.reserve(30000); // default bin size
      current_bin.push_back(temp);
    }
    temporal_bins_with_radius.push_back(current_bin);
    idx++;
  }

  if (time_window > max_tolerable_window) {
    time_window = max_tolerable_window;
  }
  time_windows.emplace_back(time_window);
}
*/

// Process the trajectory and compute the temporal bins
// which will be used to compute the dwell regions
// for the dwell region query
void
OfflineDwellRegionTauIndex::process_trajectory_and_compute_bins(size_t trajectory_id, size_t time_idx)
{
  unsigned time_window = time_windows[time_idx];
  size_t current_time_window = 0;
  for (auto& heap : directional_heaps) {
    heap.clear();
  }

  size_t first_point_idx = 0;
  size_t current_point_idx = first_point_idx;
  bool discard_subtraj = false;
  for (; current_point_idx < trajectory_ptr->size(); current_point_idx++) {
    // ------------------ STEP 1: Add point from the trajectory to the heaps -----------------------
    auto curr_point = trajectory_ptr->at(current_point_idx);
    for (size_t i = 0; i < num_heaps; i++) {
      auto theta = i * angle_between_heaps;
      auto dot_product = dot_product_with_unit_vector(curr_point, theta);
      directional_heaps[i].push(IDKeyPair<double>{ current_point_idx, dot_product });
    }

    // Add the difference between timestamps of current point and the previous point
    // to compute the current time window of the current subtrajectory
    bool inconsistent_timestamp = false;
    bool intolerable_time_window = false;
    if (current_point_idx > first_point_idx) {
      size_t prev_timestamp = trajectory_ptr->timestamp_at(current_point_idx - 1);
      size_t curr_timestamp = trajectory_ptr->timestamp_at(current_point_idx);

      if (curr_timestamp < prev_timestamp) {
        inconsistent_timestamp = true;
      } else {
        current_time_window += (curr_timestamp - prev_timestamp);
        if (current_time_window > time_windows[time_idx + 1]) {
          intolerable_time_window = true;
        }
      }
    }

    if (inconsistent_timestamp || intolerable_time_window) {
      // inconsistent timestamp or intolerable time_window,
      // so terminate the current sub-trajectory, i.e.,
      // move the first point of subtrajectory to the current point,
      // and clear the directional heaps (i.e., reset everything)
      first_point_idx = current_point_idx;
      current_point_idx--;
      current_time_window = 0;
      for (auto& heap : directional_heaps) {
        heap.clear();
      }
      discard_subtraj = false;
    }

    if (current_time_window < time_window) {
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

    // -------------Sanity Checks------------
    if (R_SEC_upper_bound < 0.0) {
      std::cout << "Upper bound negative. R = " << R_SEC_upper_bound << std::endl;
      discard_subtraj = true; // discard this subtrajectory
    }
    if (R_SEC_upper_bound < min_tolerable_radius || R_SEC_upper_bound > max_radius) {
      discard_subtraj = true; // discard this subtrajectory
    }
    if (!discard_subtraj && std::isgreater((R_SEC_lower_bound - R_SEC_upper_bound), Constants::tolerance)) {
      std::cout << trajectory_id << " " << first_point_idx << " " << current_point_idx << std::endl;
      std::cout << R_SEC_upper_bound << " " << R_SEC_lower_bound << std::endl;
      std::cerr << "ERROR: LOWER BOUND > UPPER BOUND!!" << std::endl;
      exit(1); // discard this subtrajectory
    }

    if (!discard_subtraj) {
      // -------------------------- Step 3: Compute actual SEC_S ------------------------------------
      auto SEC_ring_points_pair = compute_SEC_S(trajectory_ptr, directional_heaps, SEC_inner.value());
      auto SEC_S = SEC_ring_points_pair.first;
      if (!SEC_S.has_value()) {
        std::cerr << "ERROR: Could not compute SEC_S. Exiting..." << std::endl;
        exit(1);
      }
      auto SEC_S_disk = SEC_S.value();
      double disk_radius = SEC_S_disk.get_radius();

      double temp_min_radius = min_tolerable_radius;
      double current_radius = min_radius;
      size_t radii_idx = 0;
      for (size_t idx = 0; idx < num_radii; idx++) {
        if (disk_radius <= radii_values[idx] && disk_radius > temp_min_radius) {
          current_radius = radii_values[idx];
          radii_idx = idx;
          break;
        }
        temp_min_radius = radii_values[idx];
      }

      // ------------------- Step 4: Create an entry for the subtrajectory in the index --------------
      if (current_radius == min_radius && disk_radius < min_tolerable_radius) {
        discard_subtraj = true; // discard this subtrajectory
      }

      if (!discard_subtraj) {
        Point center = SEC_S_disk.get_center();
        TauCellEntry current_bin_entry(trajectory_id,
                                       first_point_idx,
                                       current_point_idx,
                                       current_time_window,
                                       disk_radius,
                                       center.x,
                                       center.y,
                                       trajectory_ptr->timestamp_at(first_point_idx),
                                       trajectory_ptr->timestamp_at(current_point_idx));
        temporal_bins_with_radius[time_idx][radii_idx].push_back(current_bin_entry);
      }
    }

    // move the time window by one point
    for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
      directional_heaps[heap_idx].remove_id(first_point_idx);
    }
    unsigned first_idx_time =
      (trajectory_ptr->timestamp_at(first_point_idx + 1) - trajectory_ptr->timestamp_at(first_point_idx));
    current_time_window -= first_idx_time;
    first_point_idx++;
    discard_subtraj = false;
  }
}

// Set the current trajectory, preprocess it, compute the bins for it
// and update the temporal_bins that contains subtrajectories from all trajectories
void
OfflineDwellRegionTauIndex::set_trajectory_and_update_temporal_bins(Trajectory* traj, const size_t traj_idx)
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

  //---------- Step 2: Preprocess the current trajectory to compute the temporal bins--------
  for (size_t time_idx = 0; time_idx < num_temporal_bins; time_idx++) {
    process_trajectory_and_compute_bins(traj_idx, time_idx);

    // Reset all the directional heaps to prepare them form next time window
    for (auto& heap : directional_heaps) {
      heap.clear();
    }
  }
}

// Preprocess all the trajectories in the full dataset
// and compute the radii_bins for all trajectories
void
OfflineDwellRegionTauIndex::preprocess_trajectory_dataset(std::filesystem::path filename, std::string data_path)
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

    // Populate the trajectory object with points from the file
    Trajectory traj;
    while (ss.good()) {
      double lon = 0.0, lat = 0.0;
      size_t timestamp = 0;
      ss >> lat >> lon >> timestamp;
      traj.emplace_back(Point_from_lon_lat(lon, lat), timestamp);
    }

    if (is_outlier_trajectory(traj.at(0))) {
      continue;
    }
    set_trajectory_and_update_temporal_bins(&traj, traj_idx);
    // if (traj_idx == 8)
    // break;
  }

  // Sort the subtrajectories in temporal_bins by radius
  for (size_t time_idx = 0; time_idx < num_temporal_bins; time_idx++) {
    for (size_t radii_idx = 0; radii_idx < num_radii; radii_idx++) {
      if (!temporal_bins_with_radius[time_idx][radii_idx].size()) {
        continue;
      }
      std::sort(temporal_bins_with_radius[time_idx][radii_idx].begin(),
                temporal_bins_with_radius[time_idx][radii_idx].end(),
                TauCellEntryComparator());

      remove_duplicate_temporal_entries(time_idx, radii_idx);

      int current_radius_bin = (int)(radii_values[radii_idx] * 10.0);
      std::string index_dir = data_path + "Tau_Index_Exp/" + std::to_string(num_temporal_bins) + "_bins/";
      //std::string index_dir = data_path + "Tau_Index_Linear/" + std::to_string(num_temporal_bins) + "_bins/";
      std::string index_string = index_dir +
                                 std::to_string(time_windows[time_idx]) + "_radius_" +
                                 std::to_string(current_radius_bin) + "_k_" + std::to_string(num_heaps) + ".txt";
      std::cout << "Writing index file: " << index_string << " "
                << temporal_bins_with_radius[time_idx][radii_idx].size() << std::endl;
      auto index_filename = std::filesystem::path(index_string);
      std::ofstream indexFile(index_filename, std::ios::out);
      indexFile << "traj_id"
                << "\t"
                << "start_idx"
                << "\t"
                << "last_idx"
                << "\t"
                << "duration"
                << "\t"
                << "SEC_S_radius"
                << "\t"
                << "SEC_S_center_x"
                << "\t"
                << "SEC_S_center_y" << std::endl;

      for (auto result : temporal_bins_with_radius[time_idx][radii_idx]) {
        indexFile << result.trajectory_id << "\t" << result.first_idx_subtraj << "\t" << result.last_idx_subtraj << "\t"
                  << result.duration_subtraj << "\t" << result.SEC_S_radius << "\t" << result.SEC_S_center_x << "\t"
                  << result.SEC_S_center_y << std::endl;
      }
      indexFile.close();
    }
  }
  auto end = clock();
  preprocess_time = (end - start) / (double)CLOCKS_PER_SEC;
}

void
OfflineDwellRegionTauIndex::remove_duplicate_temporal_entries(size_t time_idx, size_t radii_idx)
{
  size_t comparable_idx = -1;
  bool same_dwell_region = false;
  for (size_t idx = 0; idx < temporal_bins_with_radius[time_idx][radii_idx].size() - 1; idx++) {
    if (temporal_bins_with_radius[time_idx][radii_idx][idx].trajectory_id !=
        temporal_bins_with_radius[time_idx][radii_idx][idx + 1].trajectory_id) {
      if (same_dwell_region) {
        if (temporal_bins_with_radius[time_idx][radii_idx][idx].duration_subtraj <=
            temporal_bins_with_radius[time_idx][radii_idx][comparable_idx].duration_subtraj) {
          temporal_bins_with_radius[time_idx][radii_idx].erase(temporal_bins_with_radius[time_idx][radii_idx].begin() +
                                                               comparable_idx);
        } else {
          temporal_bins_with_radius[time_idx][radii_idx].erase(temporal_bins_with_radius[time_idx][radii_idx].begin() +
                                                               idx);
        }
        idx--;
      }
      same_dwell_region = false;
      continue;
    }
    TauCellEntry entry1 = temporal_bins_with_radius[time_idx][radii_idx][idx];
    TauCellEntry entry2 = temporal_bins_with_radius[time_idx][radii_idx][idx + 1];
    // Both entries are for same SEC_S with same center and radius
    if (std::islessequal(std::fabs(entry1.SEC_S_radius - entry2.SEC_S_radius), 0.00001) &&
        std::islessequal(std::fabs(entry1.SEC_S_center_x - entry2.SEC_S_center_x), 0.00001) &&
        std::islessequal(std::fabs(entry1.SEC_S_center_y - entry2.SEC_S_center_y), 0.00001)) {
      size_t diff_timestamp = (entry2.last_idx_timestamp - entry1.last_idx_timestamp);
      if (diff_timestamp < 50) {
        temporal_bins_with_radius[time_idx][radii_idx][idx].last_idx_timestamp = entry2.last_idx_timestamp;
        temporal_bins_with_radius[time_idx][radii_idx][idx].last_idx_subtraj = entry2.last_idx_subtraj;
        temporal_bins_with_radius[time_idx][radii_idx][idx].duration_subtraj += diff_timestamp;

        // delete the second entry
        temporal_bins_with_radius[time_idx][radii_idx].erase(temporal_bins_with_radius[time_idx][radii_idx].begin() +
                                                             idx + 1);
        idx--;
      } else {
        if (!same_dwell_region) {
          comparable_idx = idx;
          same_dwell_region = true;
        } else {
          if (temporal_bins_with_radius[time_idx][radii_idx][idx].duration_subtraj <=
              temporal_bins_with_radius[time_idx][radii_idx][comparable_idx].duration_subtraj) {
            temporal_bins_with_radius[time_idx][radii_idx].erase(
              temporal_bins_with_radius[time_idx][radii_idx].begin() + comparable_idx);
          } else {
            temporal_bins_with_radius[time_idx][radii_idx].erase(
              temporal_bins_with_radius[time_idx][radii_idx].begin() + idx);
          }
          idx--;
        }
      }
    }
  }
}
} // namespace DwellRegions