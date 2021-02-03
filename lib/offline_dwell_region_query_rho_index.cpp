#include "../include/disk.hpp"
#include "../include/offline_dwell_region_computation.hpp"
#include "../include/line.hpp"
#include "../include/util.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <set>
#include <string>

namespace {
void
find_distinct_subtrajectory_entries(
  std::vector<DwellRegions::OfflineDwellRegionQueryRhoIndex::ResultSubtrajectoryEntry>& subtraj_entries)
{
  for (size_t i = 0; i < (subtraj_entries.size() - 1); i++) {
    if (subtraj_entries[i].last_index == subtraj_entries[i + 1].last_index) {
      subtraj_entries.erase(subtraj_entries.begin() + i + 1);
      i--;
    }
  }
  for (size_t i = 0; i < (subtraj_entries.size() - 1); i++) {
    if (subtraj_entries[i].start_index == subtraj_entries[i + 1].start_index) {
      subtraj_entries.erase(subtraj_entries.begin() + i);
      i--;
    }
  }
}

void
read_trajectory(DwellRegions::Trajectory* traj, std::filesystem::path filename)
{
  std::ifstream trajFile(filename, std::ios::in);
  std::string line;

  if (!trajFile.is_open()) {
    std::cerr << "ERROR: Could not open file. Exiting." << std::endl;
    exit(1);
  }

  std::getline(trajFile, line);
  size_t num_points = std::stoi(line);

  for (size_t point_idx = 0; point_idx < num_points; point_idx++) {
    double lon = 0.0, lat = 0.0;
    size_t timestamp = 0;
    trajFile >> lat >> lon >> timestamp;
    traj->emplace_back(DwellRegions::Point_from_lon_lat(lon, lat), timestamp);
  }
}
}

namespace DwellRegions {

// Prepare and setup the bins of different radii
// based on the number of radii bins
void
OfflineDwellRegionQueryRhoIndex::prepare_radii_bins()
{
  double width_of_radii_bins = (max_radius - min_radius) / (num_radii_bins - 1);
  double radius = min_radius;
  while (radius <= max_radius) {
    radii_bins.emplace_back(radius);
    radius += width_of_radii_bins;
  }
}

void
OfflineDwellRegionQueryRhoIndex::read_indexed_radii_bin(size_t radius_index)
{
  int current_radius_bin = (int) (radii_bins[radius_index] * 10.0 + 0.5);
  std::string index_string = data_path + "Rho_Index/" + std::to_string(num_radii_bins) + "_bins/" +
                            "radius_" + std::to_string(current_radius_bin) +
                            "_k_" + std::to_string(num_heaps) + ".txt";
  auto index_filename = std::filesystem::path(index_string);
  std::ifstream indexFile(index_filename, std::ios::in);
  std::string line;
  
  size_t traj_id, start_idx, last_idx, duration;
  size_t next_last_idx, next_duration;
  std::getline(indexFile, line);

  radii_bin_with_subtrajectories.clear();
  while (!indexFile.eof()) {
    indexFile >> traj_id >> start_idx >> last_idx >> duration
              >> next_last_idx >> next_duration;
    RadiiBinEntry bin_entry(traj_id, start_idx, last_idx, duration, next_duration, next_last_idx);
    radii_bin_with_subtrajectories.emplace_back(bin_entry);
  }
  indexFile.close();
}
//--------------------------------------EXECUTE OFFLINE DWELL REGION QUERY-------------------------------------------

// Shorten subtrajectories to compute the dwell regions
// those satisfy the query parameters to compute dwell regions
std::vector<OfflineDwellRegionQueryRhoIndex::ResultSubtrajectoryEntry>
OfflineDwellRegionQueryRhoIndex::shorten_subtrajs_to_compute_dwell_regions(size_t traj_id, size_t time_idx)
{
  // size_t num_dwell_regions = 0;
  std::vector<ResultSubtrajectoryEntry> subtrajs_for_results;

  // Read the corresponding trajectory with traj_id
  std::string file = data_path + "GeoLifeDataByID/" + std::to_string(traj_id) + ".txt";
  auto filename = std::filesystem::path(file);
  Trajectory traj;
  read_trajectory(&traj, filename);

  // Prepare directional heaps for the trajectory
  directional_heaps.clear();
  for (size_t i = 0; i < num_heaps; i++) {
    directional_heaps.emplace_back(MaxIDQueue<double>(traj.size()));
  }

  // Iterate over all subtrajectories of the trajectory with traj_id
  size_t window_idx = time_idx;
  for (; window_idx < radii_bin_with_subtrajectories.size(); window_idx++) {
    if (radii_bin_with_subtrajectories[window_idx].trajectory_id != traj_id) {
      continue;
    }

    size_t start_idx = radii_bin_with_subtrajectories[window_idx].first_idx_subtraj;
    size_t last_idx = radii_bin_with_subtrajectories[window_idx].last_idx_subtraj;
    for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
      directional_heaps[heap_idx].clear();
    }

    // Update directional heaps with the subtrajctory that needs to be shorten
    // and compute the time window of the subtrajectory
    size_t subtraj_time_window = 0;
    for (size_t point_idx = start_idx; point_idx < (last_idx + 1); point_idx++) {
      auto point = traj.at(point_idx);
      for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
        auto theta = heap_idx * angle_between_heaps;
        auto dot_product = dot_product_with_unit_vector(point, theta);
        directional_heaps[heap_idx].push(IDKeyPair<double>{ point_idx, dot_product });
      }
      if (point_idx > start_idx) {
        subtraj_time_window += (traj.timestamp_at(point_idx) - traj.timestamp_at(point_idx - 1));
      }
    }

    auto inner_circle = compute_inner_SEC(&traj, directional_heaps);
    auto actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
    auto actual_disk_radius = actual_circle.first.value().get_radius();

    // -------------------CASE 1: The full subtrajectory satisfies the query -------------------
    // Check if the subtrajectory satisfies the query, add it in results;
    // Otherwise, shorten the subtrajectory from both side
    if (actual_disk_radius <= query_radius) {
      ResultSubtrajectoryEntry result_entry(traj_id, start_idx, last_idx, actual_disk_radius, subtraj_time_window);
      subtrajs_for_results.emplace_back(result_entry);
    } else {
      std::vector<size_t> removed_traj_indexes; // Track the removed last points for further usage
      auto temp_radius = actual_disk_radius;
      auto temp_window = subtraj_time_window;

      //----------------------CASE 2: Shorten the subtrajctory from right side-------------------
      // Shorten the subtrajectory from right side, delete last points one by one
      size_t current_idx = last_idx;
      while (current_idx > start_idx && actual_disk_radius > query_radius) {
        for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
          directional_heaps[heap_idx].remove_id(current_idx);
        }
        removed_traj_indexes.emplace_back(current_idx);
        inner_circle = compute_inner_SEC(&traj, directional_heaps);
        actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
        actual_disk_radius = actual_circle.first.value().get_radius();
        temp_window -= (traj.timestamp_at(current_idx) - traj.timestamp_at(current_idx - 1));
        current_idx--;
      }
      if (actual_disk_radius <= query_radius && temp_window > max_window_size) {
        ResultSubtrajectoryEntry result_entry(traj_id, start_idx, current_idx, actual_disk_radius, temp_window);
        subtrajs_for_results.emplace_back(result_entry);
      }
      // Add the removed last points back to the subtrajectory
      for (auto idx : removed_traj_indexes) {
        auto point = traj.at(idx);
        for (size_t i = 0; i < num_heaps; i++) {
          auto theta = i * angle_between_heaps;
          auto dot_product = dot_product_with_unit_vector(point, theta);
          directional_heaps[i].push(IDKeyPair<double>{ idx, dot_product });
        }
      }
      //---------------------CASE 3: Shorten the subtrajecory from left side----------------------
      // Shorten the subtrajectory from left side, delete start points one by one
      current_idx = start_idx;
      actual_disk_radius = temp_radius;
      temp_window = subtraj_time_window;
      while (current_idx < last_idx && actual_disk_radius > query_radius) {
        for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
          directional_heaps[heap_idx].remove_id(current_idx);
        }
        inner_circle = compute_inner_SEC(&traj, directional_heaps);
        actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
        actual_disk_radius = actual_circle.first.value().get_radius();
        temp_window -= (traj.timestamp_at(current_idx + 1) - traj.timestamp_at(current_idx));
        current_idx++;
      }
      if (actual_disk_radius <= query_radius && temp_window > max_window_size) {
        ResultSubtrajectoryEntry result_entry(traj_id, current_idx, last_idx, actual_disk_radius, temp_window);
        subtrajs_for_results.emplace_back(result_entry);
      }
    }
  }
  // Sort the resultant subtrajectories to remove the completely overlapped ones
  std::sort(subtrajs_for_results.begin(), subtrajs_for_results.end(), SubtrajectoryEntryComparator());
  find_distinct_subtrajectory_entries(subtrajs_for_results);

  // num_dwell_regions = subtrajs_for_results.size();
  // std::cout << "Trajectory id: " << traj_id << ", Number of dwell regions: " << num_dwell_regions
  // << std::endl;
  return subtrajs_for_results;
}

// Extend subtrajectories to compute the dwell regions
// those satisfy the query parameters to compute dwell regions
std::vector<OfflineDwellRegionQueryRhoIndex::ResultSubtrajectoryEntry>
OfflineDwellRegionQueryRhoIndex::extend_subtrajs_to_compute_dwell_regions(size_t traj_id, size_t time_idx)
{
  // size_t num_dwell_regions = 0;
  std::vector<ResultSubtrajectoryEntry> subtrajs_for_results;

  // Read the corresponding trajectory with traj_id
  std::string file = data_path + "GeoLifeDataByID/" + std::to_string(traj_id) + ".txt";
  auto filename = std::filesystem::path(file);
  Trajectory traj;
  read_trajectory(&traj, filename);

  // Prepare directional heaps for the trajectory
  directional_heaps.clear();
  for (size_t i = 0; i < num_heaps; i++) {
    directional_heaps.emplace_back(MaxIDQueue<double>(traj.size()));
  }

  // Iterate over all subtrajectories of the trajectory with traj_id
  size_t window_idx = time_idx;
  for (; window_idx < radii_bin_with_subtrajectories.size(); window_idx++) {
    if (radii_bin_with_subtrajectories[window_idx].trajectory_id != traj_id) {
      continue;
    }

    size_t start_idx = radii_bin_with_subtrajectories[window_idx].first_idx_subtraj;
    size_t last_idx = radii_bin_with_subtrajectories[window_idx].last_idx_subtraj;
    for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
      directional_heaps[heap_idx].clear();
    }

    // Update directional heaps with the subtrajctory that needs to be extended
    // and compute the time window of the subtrajectory
    size_t subtraj_time_window = 0;
    for (size_t point_idx = start_idx; point_idx < (last_idx + 1); point_idx++) {
      auto point = traj.at(point_idx);
      for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
        auto theta = heap_idx * angle_between_heaps;
        auto dot_product = dot_product_with_unit_vector(point, theta);
        directional_heaps[heap_idx].push(IDKeyPair<double>{ point_idx, dot_product });
      }
      if (point_idx > start_idx) {
        subtraj_time_window += (traj.timestamp_at(point_idx) - traj.timestamp_at(point_idx - 1));
      }
    }
    size_t next_bin_idx = radii_bin_with_subtrajectories[window_idx].last_idx_of_subtraj_in_next_bin_entry;

    auto inner_circle = compute_inner_SEC(&traj, directional_heaps);
    auto actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
    auto actual_disk_radius = actual_circle.first.value().get_radius();

    // If the subtrajectory cannot be extended, add it in results;
    // Otherwise, extend the subtrajectory from both side
    if (last_idx < next_bin_idx) {
      ResultSubtrajectoryEntry result_entry(traj_id, start_idx, last_idx, actual_disk_radius, subtraj_time_window);
      subtrajs_for_results.emplace_back(result_entry);

      std::vector<size_t> removed_traj_indexes;
      auto temp_radius = actual_disk_radius;
      auto temp_window = subtraj_time_window;

      //----------------------CASE 1: Extend the subtrajctory to the right side-------------------
      // Extend the subtrajectory to the right side, add points one by one
      // upto the next_last_index of the next_bin entry of the same subtrajectory
      // or, as long as the extended subtrajectory satisfies the query
      size_t current_idx = last_idx + 1;
      while (current_idx <= next_bin_idx && actual_disk_radius < query_radius) {
        auto point = traj.at(current_idx);
        for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
          auto theta = heap_idx * angle_between_heaps;
          auto dot_product = dot_product_with_unit_vector(point, theta);
          directional_heaps[heap_idx].push(IDKeyPair<double>{ current_idx, dot_product });
        }

        inner_circle = compute_inner_SEC(&traj, directional_heaps);
        actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
        actual_disk_radius = actual_circle.first.value().get_radius();
        temp_window += (traj.timestamp_at(current_idx) - traj.timestamp_at(current_idx - 1));
        current_idx++;
      }
      current_idx--;

      // Adding points can make the actual_disk more than query_radius,
      // If so, delete one point to get the actual_disk that satifies the query
      // Otherwise, add the subtrajectory to the resultant set
      if (actual_disk_radius > query_radius) {
        for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
          directional_heaps[heap_idx].remove_id(current_idx);
        }
        inner_circle = compute_inner_SEC(&traj, directional_heaps);
        actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
        actual_disk_radius = actual_circle.first.value().get_radius();
        temp_window -= (traj.timestamp_at(current_idx) - traj.timestamp_at(current_idx - 1));
        // subtraj_ranges_for_results.emplace_back(std::make_pair(start_idx, current_idx));
      }
      result_entry = ResultSubtrajectoryEntry(traj_id, start_idx, current_idx, actual_disk_radius, temp_window);
      subtrajs_for_results.emplace_back(result_entry);

      //------------------CASE 2: Shorten the extended subtrajctory from the left side-----------------
      // Extend the subtrajectory to the right side upto the next_last_index,
      // then delete start points one by one from left side to satisfy the query
      actual_disk_radius = temp_radius;
      temp_window = subtraj_time_window;
      if (current_idx < next_bin_idx) {
        while (current_idx <= next_bin_idx) {
          auto point = traj.at(current_idx);
          for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
            auto theta = heap_idx * angle_between_heaps;
            auto dot_product = dot_product_with_unit_vector(point, theta);
            directional_heaps[heap_idx].push(IDKeyPair<double>{ current_idx, dot_product });
          }
          temp_window += (traj.timestamp_at(current_idx) - traj.timestamp_at(current_idx - 1));
          current_idx++;
        }

        current_idx = start_idx;
        inner_circle = compute_inner_SEC(&traj, directional_heaps);
        actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
        actual_disk_radius = actual_circle.first.value().get_radius();
        while (current_idx < next_bin_idx && actual_disk_radius > query_radius) {
          for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
            directional_heaps[heap_idx].remove_id(current_idx);
          }
          inner_circle = compute_inner_SEC(&traj, directional_heaps);
          actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
          actual_disk_radius = actual_circle.first.value().get_radius();
          temp_window -= (traj.timestamp_at(current_idx + 1) - traj.timestamp_at(current_idx));
          current_idx++;
        }
        if (actual_disk_radius <= query_radius && temp_window > max_window_size) {
          result_entry = ResultSubtrajectoryEntry(traj_id, current_idx, last_idx, actual_disk_radius, temp_window);
          subtrajs_for_results.emplace_back(result_entry);
        }
      }
    }
    //-------------------CASE 3: The full subtrajectory satisfies the query-------------------
    else {
      ResultSubtrajectoryEntry result_entry(traj_id, start_idx, last_idx, actual_disk_radius, subtraj_time_window);
      subtrajs_for_results.emplace_back(result_entry);
    }
  }
  // Sort the resultant subtrajectories to remove the completely overlapped ones
  std::sort(subtrajs_for_results.begin(), subtrajs_for_results.end(), SubtrajectoryEntryComparator());
  find_distinct_subtrajectory_entries(subtrajs_for_results);

  // num_dwell_regions = subtrajs_for_results.size();
  // std::cout << "Trajectory id: " << traj_id << ", Number of dwell regions: " << num_dwell_regions
  // << std::endl;
  return subtrajs_for_results;
}

// Execute the dwell region query using the preprocessed radii bins
// and get the dwell regions for the dwell region query
void
OfflineDwellRegionQueryRhoIndex::execute_dwell_region_query()
{
  // Note 3: Radii_bins quadrants according to the journal paper (Figure 11)
  //    NW   |   NE
  // --------|--------
  //    SW   |   SE

  size_t radius_index = 0;
  bool closest_to_left_radius = false;
  bool closest_to_right_radius = false;
  bool equal_to_bin_radius = false;

  auto start = clock();
  size_t index = 0;

  // Get the index of radii_bin with radius = query radius
  // or, the index of the closest radius radii_bin, if query radius is in between two bins radii
  while (radii_bins[index] <= query_radius) {
    index++;
    if (index == radii_bins.size())
      break;
  }
  if (index > 0)
    index--;

  // Determine the index of radii_bins to be explored
  if (std::fabs(radii_bins[index] - query_radius) < Constants::tolerance) {
    equal_to_bin_radius = true;
    radius_index = index;
  } else if (index < (radii_bins.size() - 1)) {
    // This is a heuristic. We extend the smaller bin or shorten the larger bin
    // Depending on how close the query radius is to the radii of the two bins.
    // HEURISTIC: Shorten the larger bin if query radius is closer to the larger bin.
    auto diff_with_left_radius = std::fabs(radii_bins[index] - query_radius);
    auto diff_with_right_radius = std::fabs(radii_bins[index + 1] - query_radius);
    if (diff_with_left_radius < diff_with_right_radius && query_radius > radii_bins[index]) {
      closest_to_left_radius = true;
      radius_index = index;
    } else if (diff_with_left_radius < diff_with_right_radius && query_radius < radii_bins[index]) {
      closest_to_right_radius = true;
      radius_index = index;
    } else {
      closest_to_right_radius = true;
      radius_index = index + 1;
    }
  } else {
    closest_to_left_radius = true;
    radius_index = index;
  }

  read_indexed_radii_bin(radius_index);

  // Determine the index of radii_bin[radius_index] with duration
  // greater or equal to max_window_size
  index = 0;
  while (radii_bin_with_subtrajectories[index].duration_subtraj < max_window_size) {
    index++;
    if (index == radii_bin_with_subtrajectories.size())
      break;
  }
  auto duration_index = index;

  size_t num_dwell_regions = 0;
  //-------------------------CASE 1: No subtrajectories satisfy the query-----------------------
  if (duration_index == radii_bin_with_subtrajectories.size()) {
    std::cout << "NO: There is no eligible dwell region!" << std::endl;
    return;
  }

  //------------------CASE 2: Subtrajectories found directly from a preprocessed bin--------------
  if (equal_to_bin_radius) {
    num_dwell_regions = (radii_bin_with_subtrajectories.size() - duration_index);
    std::cout << "YES: Found directly from preprocessed radii_bins!" << std::endl;
    std::cout << "Query radius: " << query_radius << ", Time window: " << max_window_size
              << ", Number of dwell regions: " << num_dwell_regions << std::endl;
    return;
  }

  // Explore the subtrajectories with trajectory ids from preprocessed bins
  std::vector<size_t> traj_ids_list;
  for (size_t i = duration_index; i < radii_bin_with_subtrajectories.size(); i++) {
    traj_ids_list.emplace_back(radii_bin_with_subtrajectories[i].trajectory_id);
  }
  std::set<size_t> ids_set(traj_ids_list.begin(), traj_ids_list.end());
  traj_ids_list.assign(ids_set.begin(), ids_set.end());

  //------------------CASE 3: Explore the subtrajectories of right preprocessed bin---------------
  if (closest_to_right_radius) {
    for (auto traj_id : traj_ids_list) {
      std::vector<ResultSubtrajectoryEntry> temp_results;
      temp_results = shorten_subtrajs_to_compute_dwell_regions(traj_id, duration_index);
      for (auto entry : temp_results) {
        result_subtrajectories_with_dwell_regions.emplace_back(entry);
      }
    }
  }

  //------------------CASE 4: Explore the subtrajectories of left preprocessed bin---------------
  if (closest_to_left_radius) {
    for (auto traj_id : traj_ids_list) {
      std::vector<ResultSubtrajectoryEntry> temp_results;
      temp_results = extend_subtrajs_to_compute_dwell_regions(traj_id, duration_index);
      for (auto entry : temp_results) {
        result_subtrajectories_with_dwell_regions.emplace_back(entry);
      }
    }
  }
  auto end = clock();
  execution_time = (end - start) / (double)CLOCKS_PER_SEC;

  std::cout << "-------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Execution time: " << execution_time << " seconds" << std::endl;
  std::cout << "-------------------------Results of offline query execution--------------------------" << std::endl;
  std::cout << "Total number of dwell regions: " << result_subtrajectories_with_dwell_regions.size() << std::endl;
  for (auto entry : result_subtrajectories_with_dwell_regions) {
    std::cout << "Trajectory id: " << entry.trajectory_id << ", Start index: " << entry.start_index
              << ", Last index: " << entry.last_index << ", actual radius: " << entry.actual_disk_radius
              << ", time window: " << entry.time_window << std::endl;
  }
}

} // namespace DwellRegions
