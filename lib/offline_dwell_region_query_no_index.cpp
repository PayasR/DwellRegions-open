#include "../include/line.hpp"
#include "../include/offline_dwell_region_computation.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>

namespace {

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


//-----------------------------------------TEST OFFLINE DWELL REGION QUERY---------------------------------------------
void
OfflineDwellRegionQueryNoIndex::test_offline_query_per_traj(size_t traj_id)
{
  DwellRegions::Trajectory traj;
  std::string file = data_path + "GeoLifeDataByID/" + std::to_string(traj_id) + ".txt";
  auto filename = std::filesystem::path(file);

  read_trajectory(&traj, filename);
  directional_heaps.clear();
  for (size_t i = 0; i < num_heaps; i++) {
    directional_heaps.emplace_back(MaxIDQueue<double>(traj.size()));
  }

  size_t counter = 0;
  size_t last_entry_idx = 0;
  // Iterate over the trajcetory starting from start_idx
  for (size_t start_idx = 0; start_idx < (traj.size() - 1); start_idx++) {
    // Initialization for current iteration
    size_t current_time_window = 0;
    for (size_t i = 0; i < num_heaps; i++) {
      directional_heaps[i].clear();
    }

    for (size_t i = start_idx; i < traj.size(); i++) {
      auto point = traj.at(i);
      for (size_t j = 0; j < num_heaps; j++) {
        auto theta = j * angle_between_heaps;
        auto dot_product = dot_product_with_unit_vector(point, theta);
        directional_heaps[j].push(IDKeyPair<double>{ i, dot_product });
      }

      if (i > start_idx) {
        current_time_window += (traj.timestamp_at(i) - traj.timestamp_at(i - 1));
      }
      if (current_time_window < max_window_size) {
        continue;
      }
      if (i <= last_entry_idx) {
        continue;
      }

      auto inner_circle = compute_inner_SEC(&traj, directional_heaps);
      auto outer_circle = compute_outer_SEC(&traj, directional_heaps);
      auto actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
      auto actual_disk_radius = actual_circle.first.value().get_radius();
      
      if (actual_disk_radius <= query_radius) // don't change this condition
      {
        if (i == (traj.size() - 1) && last_entry_idx != i) {
          counter++;
          last_entry_idx = i;
          std::cout << "Start idx: " << start_idx << ", end idx: " << i << ", counter: " << counter
                    << ", actual_disk: " << actual_disk_radius << ", time window: " << current_time_window
                    << ", actual_SEC_S_x: " << actual_circle.first.value().get_center().x
                    << ", actual_SEC_S_y: " << actual_circle.first.value().get_center().y
                    << ", upper: " << outer_circle.value().get_radius() << std::endl;
          break;
        }
      } else {
        for (size_t heap_idx = 0; heap_idx < num_heaps; heap_idx++) {
          directional_heaps[heap_idx].remove_id(i);
        }
        current_time_window -= (traj.timestamp_at(i) - traj.timestamp_at(i - 1));
        inner_circle = compute_inner_SEC(&traj, directional_heaps);
        outer_circle = compute_outer_SEC(&traj, directional_heaps);
        actual_circle = compute_SEC_S(&traj, directional_heaps, inner_circle.value());
        actual_disk_radius = actual_circle.first.value().get_radius();

        // if (std::islessequal((actual_disk_radius - query_radius), Constants::tolerance)
        if (actual_disk_radius <= query_radius // don't change this condition
            && current_time_window >= max_window_size && last_entry_idx != (i - 1)) {
          counter++;
          last_entry_idx = (i - 1);
          std::cout << "Start idx: " << start_idx << ", end idx: " << (i - 1) << ", counter: " << counter
                    << ", actual_disk: " << actual_disk_radius << ", time window: " << current_time_window
                    << ", actual_SEC_S_x: " << actual_circle.first.value().get_center().x
                    << ", actual_SEC_S_y: " << actual_circle.first.value().get_center().y
                    << ", upper: " << outer_circle.value().get_radius() << std::endl;
        }
        break;
      }
    }
  }
  std::cout << "-----------------------In test function---------------------" << std::endl;
  std::cout << "Trajectory id: " << traj_id << ", Total points in traj: " << traj.size()
            << ", Total dwell regions: " << counter << std::endl;
}

} // namespace DwellRegions