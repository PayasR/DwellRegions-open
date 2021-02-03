#pragma once

#include "disk.hpp"
#include "id_queue.hpp"
#include "region.hpp"
#include "trajectory.hpp"

#include "dwell_region_computation_base.hpp"

#include <filesystem>
#include <queue>
#include <type_traits>
#include <utility>

namespace DwellRegions {  

// Preprocessing using the Rho-Index
class OfflineDwellRegionRhoIndex : public DwellRegionComputationBase
{
public:
  OfflineDwellRegionRhoIndex(unsigned num_heaps, unsigned num_bins)
    : DwellRegionComputationBase(num_heaps)
    , num_radii_bins(num_bins)
  {
    prepare_radii_bins();
  }

  // Preprocess the full trajectory dataset
  void preprocess_trajectory_dataset(std::filesystem::path filename, std::string data_path);

  double get_preprocess_time() {
    return preprocess_time;
  }

  const std::vector<double>& get_radii_bins() {
    return radii_bins;
  }

  double get_radii_bins(size_t idx) {
    return radii_bins[idx];
  }

  const std::vector<std::vector<RadiiBinEntry>>& get_radii_bins_with_subtrajectories() {
    return radii_bins_with_subtrajectories;
  }
  
  const std::vector<RadiiBinEntry>& get_radii_bins_with_subtrajectories(size_t radius_idx) {
    return radii_bins_with_subtrajectories[radius_idx];
  }

  const RadiiBinEntry& get_radii_bins_with_subtrajectories(size_t radius_idx, size_t window_idx) {
    return radii_bins_with_subtrajectories[radius_idx][window_idx];
  }

protected:
  void prepare_radii_bins();
  void set_trajectory_and_update_radii_bins(Trajectory* traj, const size_t traj_idx);
  void process_trajectory_and_compute_bins(size_t trajectory_id,
                                           size_t radii_idx,
                                           std::vector<std::vector<RadiiBinEntry>>& current_bins);

  const size_t num_radii_bins;    // number of total radii bins
  const double min_radius = 0.6;  // minimum radius in miles
  const double max_radius = 2.0;  // maximum radius in miles
  std::vector<double> radii_bins; // different radii values

  // According to the journal paper, this is rho_i for r_i,
  // Each bin contains all the subtrajectories from trajectory dataset
  std::vector<std::vector<RadiiBinEntry>> radii_bins_with_subtrajectories;

  Trajectory* trajectory_ptr;
  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction

  double preprocess_time = 0.0;

};

// Query evaluation using the Rho-Index
class OfflineDwellRegionQueryRhoIndex : DwellRegionComputationBase {
  using index_ptr_type = gsl::not_null<OfflineDwellRegionRhoIndex*>;
public:
  OfflineDwellRegionQueryRhoIndex(
    size_t num_heaps,
    size_t max_window_size,
    double query_radius,
    unsigned num_bins,
    std::string data_path
  ) : DwellRegionComputationBase(num_heaps)
    , num_radii_bins(num_bins)
    , max_window_size(max_window_size)
    , query_radius(query_radius)
    , data_path(data_path)
  {
    prepare_radii_bins();
  }

  // Execute the offline dwell region query with query_radius and max_window_size
  void execute_dwell_region_query();

  double get_execution_time() {
    return execution_time;
  }

  size_t get_num_dwell_regions() {
    return result_subtrajectories_with_dwell_regions.size();
  }

  struct ResultSubtrajectoryEntry
  {
    ResultSubtrajectoryEntry(size_t traj_id,
                             size_t start_idx,
                             size_t last_idx,
                             double actual_radius,
                             size_t time_window)
      : trajectory_id(traj_id)
      , start_index(start_idx)
      , last_index(last_idx)
      , actual_disk_radius(actual_radius)
      , time_window(time_window)
    {}

    size_t trajectory_id;
    size_t start_index;
    size_t last_index;
    double actual_disk_radius;
    size_t time_window;
  };

protected:
  void prepare_radii_bins();
  void read_indexed_radii_bin(size_t radius_index);
  
  const size_t num_radii_bins;    // number of total radii bins
  const double min_radius = 0.6;  // minimum radius in miles
  const double max_radius = 2.0;  // maximum radius in miles
  std::vector<double> radii_bins; // different radii values

  // According to the journal paper, this is rho_i for r_i,
  // Each bin contains all the subtrajectories from trajectory dataset
  std::vector<RadiiBinEntry> radii_bin_with_subtrajectories;

  const size_t max_window_size;
  const double query_radius;
  const std::string data_path;

  std::vector<ResultSubtrajectoryEntry> result_subtrajectories_with_dwell_regions;
  double execution_time = 0.0;

  Trajectory* trajectory_ptr;
  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction

  struct SubtrajectoryEntryComparator
  {
    bool operator()(const ResultSubtrajectoryEntry& subtraj1, ResultSubtrajectoryEntry& subtraj2)
    {
      if (subtraj1.last_index != subtraj2.last_index)
        return subtraj1.last_index < subtraj2.last_index;
      return subtraj1.start_index < subtraj2.start_index;
    }
  };
  
private:
  std::vector<ResultSubtrajectoryEntry> shorten_subtrajs_to_compute_dwell_regions(size_t traj_id,
                                                                                  size_t time_idx);
  std::vector<ResultSubtrajectoryEntry> extend_subtrajs_to_compute_dwell_regions(size_t traj_id,
                                                                                 size_t time_idx);
};


// Preprocessing using the Tau-Index
class OfflineDwellRegionTauIndex : public DwellRegionComputationBase
{
public:
  OfflineDwellRegionTauIndex(unsigned num_heaps, unsigned num_bins)
    : DwellRegionComputationBase(num_heaps)
    , num_temporal_bins(num_bins)
  {
    if (num_temporal_bins > max_num_bins) {
      std::cerr << "Number of bins should be within " << max_num_bins << "!" << std::endl;
      exit(1);
    }
    max_window = time_windows[num_temporal_bins - 1];
    max_tolerable_window = time_windows[num_temporal_bins];
    prepare_temporal_bins();
  }

  // Preprocess the full trajectory dataset
  void preprocess_trajectory_dataset(std::filesystem::path filename, std::string data_path);

  double get_preprocess_time() {
    return preprocess_time;
  }

protected:
  void prepare_temporal_bins();
  void set_trajectory_and_update_temporal_bins(Trajectory* traj, const size_t traj_idx);
  void process_trajectory_and_compute_bins(size_t trajectory_id,
                                           size_t time_idx);
  void remove_duplicate_temporal_entries(size_t time_idx, size_t radii_idx);

  // For EXPONENTIAL partitioning scheme
  const size_t num_temporal_bins;   // number of total temporal bins (3, 4, 5, 6, 7)
  const size_t min_window = 3600;   // minimum time window in seconds, 1 hour
  size_t max_window;                // maximum time window in seconds
  size_t max_tolerable_window;
  size_t time_windows[8] = {3600, 7200, 14400, 28800, 57600, 115200, 230400, 259200};
  const size_t max_num_bins = 7;

  /* // For LINEAR partitioning scheme
  const size_t num_temporal_bins;   // number of total temporal bins (3, 4, 5, 6, 7)
  const size_t min_window = 3600;   // minimum time window in seconds, 1 hour
  size_t max_window = 115200;                // maximum time window in seconds
  size_t max_tolerable_window = 144000;
  std::vector<size_t> time_windows;
  const size_t max_num_bins = 7;*/

  const size_t num_radii = 8;  // number of radius values
  const double min_tolerable_radius = 0.4;
  const double min_radius = 0.6;  // minimum radius in miles
  const double max_radius = 2.0;  // maximum radius in miles
  double radii_values[8] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

  std::vector<std::vector<std::vector<TauCellEntry>>> temporal_bins_with_radius;

  Trajectory* trajectory_ptr;
  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction

  double preprocess_time = 0.0;

};

// Query evaulation using the Tau-Index
class OfflineDwellRegionQueryTauIndex : DwellRegionComputationBase {
  using index_ptr_type = gsl::not_null<OfflineDwellRegionTauIndex*>;
public:
  OfflineDwellRegionQueryTauIndex(
    size_t num_heaps,
    size_t max_window_size,
    double query_radius,
    unsigned num_bins,
    std::string data_path
  ) : DwellRegionComputationBase(num_heaps)
    , num_temporal_bins(num_bins)
    , max_window_size(max_window_size)
    , query_radius(query_radius)
    , data_path(data_path)
  {
    if (num_temporal_bins > max_num_bins) {
      std::cerr << "Number of bins should be within " << max_num_bins << "!" << std::endl;
      exit(1);
    }
    max_window = time_windows[num_temporal_bins - 1];
    max_tolerable_window = time_windows[num_temporal_bins];
    prepare_temporal_bins();
  }

  // Execute the offline dwell region query with query_radius and max_window_size
  void execute_dwell_region_query();

  double get_execution_time() {
    return execution_time;
  }

  size_t get_num_dwell_regions() {
    return cell_with_subtrajectories.size();
  }

protected:
  void prepare_temporal_bins();
  void read_indexed_temporal_bins(size_t time_index, size_t radius_index);

  // For EXPONENTIAL partitioning scheme
  const size_t num_temporal_bins;   // number of total temporal bins (3, 4, 5, 6, 7)
  const size_t min_window = 3600;   // minimum time window in seconds, 1 hour
  size_t max_window;                // maximum time window in seconds
  size_t max_tolerable_window;
  size_t time_windows[8] = {3600, 7200, 14400, 28800, 57600, 115200, 230400, 259200};
  const size_t max_num_bins = 7;

  /*//For LINEAR partitioning scheme
  const size_t num_temporal_bins;   // number of total temporal bins (3, 4, 5, 6, 7)
  const size_t min_window = 3600;   // minimum time window in seconds, 1 hour
  size_t max_window = 115200;                // maximum time window in seconds
  size_t max_tolerable_window = 144000;
  std::vector<size_t> time_windows;
  const size_t max_num_bins = 7;*/

  const size_t num_radii = 8;  // number of radius values
  const double min_tolerable_radius = 0.4;
  const double min_radius = 0.6;  // minimum radius in miles
  const double max_radius = 2.0;  // maximum radius in miles
  double radii_values[8] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

  // According to the journal paper, this is lambda_i for r_i,
  // Each bin contains all the subtrajectories from trajectory dataset
  std::vector<TauCellEntry> cell_with_subtrajectories;

  const size_t max_window_size;
  const double query_radius;
  const std::string data_path;

  double execution_time = 0.0;

  Trajectory* trajectory_ptr;
  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction

  struct ResultEntryComparator
  {
    bool operator()(const TauCellEntry& entry1, TauCellEntry& entry2)
    {
      // Compare trajectory ids
      if (entry1.trajectory_id != entry2.trajectory_id)
        return entry1.trajectory_id < entry2.trajectory_id;
      
      // Compare radius of the dwell regions if the ids are same
      if (std::isgreater(std::fabs(entry1.SEC_S_radius - entry2.SEC_S_radius), 0.000001)) {
        return entry1.SEC_S_radius < entry2.SEC_S_radius;
      }

      // Compare the x-coordinate of the dwell regions if id and radius are equal
      if (std::isgreater(std::fabs(entry1.SEC_S_center_x - entry2.SEC_S_center_x), 0.00001)) {
        return entry1.SEC_S_center_x < entry2.SEC_S_center_x;
      }

      // Compare the y-coordinate of the dwell regions if id, radius and x values are equal
      if (std::isgreater(std::fabs(entry1.SEC_S_center_y - entry2.SEC_S_center_y), 0.00001)) {
        return entry1.SEC_S_center_y < entry2.SEC_S_center_y;
      }

      // Compare the duration of the dwell region if all above are equal
      return entry1.duration_subtraj < entry2.duration_subtraj;
    }
  };
};

// Test Query evaulation without any index (NO INDEX case)
class OfflineDwellRegionQueryNoIndex : DwellRegionComputationBase {
  //using index_ptr_type = gsl::not_null<OfflineDwellRegionTauIndex*>;
public:
  OfflineDwellRegionQueryNoIndex(
    size_t num_heaps,
    size_t max_window_size,
    double query_radius,
    std::string data_path
  ) : DwellRegionComputationBase(num_heaps)
    , max_window_size(max_window_size)
    , query_radius(query_radius)
    , data_path(data_path)
  { }

  void test_offline_query_per_traj(size_t traj_id);

protected:
  const size_t max_window_size;
  const double query_radius;
  const std::string data_path;

  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction
};

}