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
    // solve ties with start index of subtrajectories
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

// Preprocessing
class OfflineDwellRegionIndex : public DwellRegionComputationBase
{
public:
  OfflineDwellRegionIndex(unsigned num_heaps, unsigned num_bins)
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

  // According to the journal paper, this is lambda_i for r_i,
  // Each bin contains all the subtrajectories from trajectory dataset
  std::vector<std::vector<RadiiBinEntry>> radii_bins_with_subtrajectories;

  Trajectory* trajectory_ptr;
  std::vector<MaxIDQueue<double>> directional_heaps; // Max heaps storing dot products in each direction

  double preprocess_time = 0.0;

};

// Query
class OfflineDwellRegionQuery : DwellRegionComputationBase {
  using index_ptr_type = gsl::not_null<OfflineDwellRegionIndex*>;
public:
  OfflineDwellRegionQuery(
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
  void execute_dwell_region_query(); // Maybe this should return something?
  void test_offline_query_per_traj(size_t traj_id);

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

  // According to the journal paper, this is lambda_i for r_i,
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

// NOTE to future: A similar design as above can be used for CommonDwellRegions as well...
// class CommonDwellRegionComputation
//{
//    CommonDwellRegionComputation(const Trajectory& T1, const Trajectory& T2, unsigned
//    _window_size, unsigned _num_heaps)
//    {}
//
// private:
//    const Trajectory& trajectory1;
//    const Trajectory& trajectory2;
//    ...
//    const unsigned window_size, num_heaps;
//    std::vector<std::priority_queue<double>> directional_heaps;
//};
}