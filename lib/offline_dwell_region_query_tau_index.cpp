#include "../include/line.hpp"
#include "../include/offline_dwell_region_computation.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>


namespace DwellRegions {

// Prepare and setup the bins of different radii
// based on the number of radii bins
void
OfflineDwellRegionQueryTauIndex::prepare_temporal_bins()
{
  /* // For LINEAR partitioning scheme
  size_t width_of_time_windows = (max_window - min_window) / (num_temporal_bins - 1);
  size_t time_window = min_window;
  time_windows.reserve(num_temporal_bins + 1);
  
  size_t idx = 0;
  while (idx < num_temporal_bins) {
    time_windows.emplace_back(time_window);
    time_window += width_of_time_windows;
    idx++;
  }

  if (time_window > max_tolerable_window) {
    time_window = max_tolerable_window;
  }
  time_windows.emplace_back(time_window);*/
}

void
OfflineDwellRegionQueryTauIndex::read_indexed_temporal_bins(size_t time_index, size_t radius_index)
{
  if (radii_values[radius_index] > query_radius) {
    radius_index--;
  }
  std::cout << time_index << " " << radius_index << std::endl;

  std::string index_dir = data_path + "Tau_Index_Exp/" + std::to_string(num_temporal_bins) + "_bins/";
  //std::string index_dir = data_path + "Tau_Index_Linear/" + std::to_string(num_temporal_bins) + "_bins/";
  for (size_t t_idx = time_index; t_idx < num_temporal_bins; t_idx++) {
    for (size_t r_idx = 0; r_idx < radius_index + 1; r_idx++) {
      int current_radius_bin = (int)(radii_values[r_idx] * 10.0);
      std::string index_string = index_dir + std::to_string(time_windows[t_idx]) + "_radius_" +
                                 std::to_string(current_radius_bin) + "_k_" + std::to_string(num_heaps) + ".txt";
      auto index_filename = std::filesystem::path(index_string);
      std::ifstream indexFile(index_filename, std::ios::in);
      std::string line;

      size_t traj_id, start_idx, last_idx, duration;
      double disk_radius, center_x, center_y;
      std::getline(indexFile, line);

      while (!indexFile.eof()) {
        indexFile >> traj_id >> start_idx >> last_idx >> duration >> disk_radius >> center_x >> center_y;
        TauCellEntry bin_entry(traj_id, start_idx, last_idx, duration, disk_radius, center_x, center_y, 0, 0);
        cell_with_subtrajectories.emplace_back(bin_entry);
      }
      indexFile.close();
    }
  }
}

//--------------------------------------EXECUTE OFFLINE DWELL REGION QUERY-------------------------------------------

// Execute the dwell region query using the preprocessed radii bins
// and get the dwell regions for the dwell region query
void
OfflineDwellRegionQueryTauIndex::execute_dwell_region_query()
{
  // Note 3: Radii_bins quadrants according to the journal paper (Figure 11)
  //    NW   |   NE
  // --------|--------
  //    SW   |   SE

  size_t time_index = 0;
  size_t radius_index = 0;
  bool dwell_region_possible = true;

  auto start = clock();
  if (max_window < max_window_size || query_radius > max_radius) {
    dwell_region_possible = false;
  }

  if (dwell_region_possible) {
    size_t index = 0;
    for (index = 0; index < num_temporal_bins; index++) {
      if (time_windows[index] >= max_window_size) {
        time_index = index;
        break;
      }
    }

    for (index = 0; index < num_radii; index++) {
      if (radii_values[index] >= query_radius) {
        radius_index = index;
        break;
      }
    }

    // read indexed bin_entries from time_index to last and 0 to radius_index
    // discard "radius_index"th bin, if query_radius is less than that bin_radius
    cell_with_subtrajectories.clear();
    read_indexed_temporal_bins(time_index, radius_index);

    if (query_radius < radii_values[radius_index]) {
      // need to explore the last bins with that bin_radius
      std::string index_dir = data_path + "Tau_Index_Exp/" + std::to_string(num_temporal_bins) + "_bins/";
      //std::string index_dir = data_path + "Tau_Index_Linear/" + std::to_string(num_temporal_bins) + "_bins/";
      int current_radius_bin = (int)(radii_values[radius_index] * 10.0);
      for (size_t t_idx = time_index; t_idx < num_temporal_bins; t_idx++) {
        std::string index_string = index_dir + std::to_string(time_windows[t_idx]) + "_radius_" +
                                   std::to_string(current_radius_bin) + "_k_" + std::to_string(num_heaps) + ".txt";
        auto index_filename = std::filesystem::path(index_string);
        std::ifstream indexFile(index_filename, std::ios::in);
        std::string line;

        size_t traj_id, start_idx, last_idx, duration;
        double disk_radius, center_x, center_y;
        std::getline(indexFile, line);

        while (!indexFile.eof()) {
          indexFile >> traj_id >> start_idx >> last_idx >> duration >> disk_radius >> center_x >> center_y;
          if (disk_radius > query_radius) {
            break;
          }
          TauCellEntry bin_entry(traj_id, start_idx, last_idx, duration, disk_radius, center_x, center_y, 0, 0);
          cell_with_subtrajectories.emplace_back(bin_entry);
        }
        indexFile.close();
      }
    }

    if (cell_with_subtrajectories.empty()) {
      dwell_region_possible = false;
    }
  }

  //-------------------------CASE 1: No subtrajectories satisfy the query-----------------------
  if (!dwell_region_possible) {
    std::cout << "NO: There is no eligible dwell region!" << std::endl;
    return;
  }

  //------------------CASE 2: Subtrajectories found directly from a preprocessed bin--------------
  std::sort(cell_with_subtrajectories.begin(), cell_with_subtrajectories.end(), ResultEntryComparator());

  auto end = clock();
  execution_time = (end - start) / (double)CLOCKS_PER_SEC;

  std::cout << "-------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Execution time: " << execution_time << " seconds" << std::endl;
  std::cout << "-------------------------Results of offline query execution--------------------------" << std::endl;
  std::cout << "Total number of dwell regions: " << cell_with_subtrajectories.size() << std::endl;
  for (auto entry : cell_with_subtrajectories) {
    std::cout << "Trajectory id: " << entry.trajectory_id << ", Start index: " << entry.first_idx_subtraj
              << ", Last index: " << entry.last_idx_subtraj << ", actual radius: " << entry.SEC_S_radius
              << ", time window: " << entry.duration_subtraj << ", SEC_S_x: " << entry.SEC_S_center_x
              << ", SEC_S_y: " << entry.SEC_S_center_y << std::endl;
  }
}

} // namespace DwellRegions