/// Implements a progam to read a trajectory and print if there exists a dwell region
#include "online_dwell_region_computation.hpp"

#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

void
process_dataset(std::filesystem::path filename,
                unsigned max_window_size,
                unsigned num_heaps,
                double query_radius,
                size_t skip_N_trajectories,
                std::string query_perf_filepath)
{
  using DwellRegions::OnlineDwellRegionComputation;
  using DwellRegions::Point;
  using DwellRegions::Trajectory;

  // Read Points from the file and add_point() to traj here...
  std::ifstream trajFile(filename, std::ios::in);
  std::ofstream query_perf_file;
  bool write_header = true;

  if (query_perf_filepath != "") {
    query_perf_file.open(query_perf_filepath, std::fstream::out);
  }
  std::string line;

  if (!trajFile.is_open()) {
    std::cerr << "ERROR: Could not open file. Exiting." << std::endl;
    exit(1);
  }

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
    // Location reading from a GPS device
    Trajectory traj;
    while (ss.good()) {
      double lon = 0.0, lat = 0.0;
      size_t timestamp = 0;
      ss >> lat >> lon >> timestamp;

      traj.emplace_back(DwellRegions::Point_from_lon_lat(lon, lat), timestamp);
      // std::cout << "Read a point {" << lat << ", " << lon << "} \t";
    }

    if (traj_idx < skip_N_trajectories) {
      continue;
    }

    auto traj_time_window = traj.timestamp_at(traj.size() - 1) - traj.timestamp_at(0);
    if (traj_time_window < max_window_size) {
      continue;
    }

    if (DwellRegions::is_outlier_trajectory(traj.at(0))) {
      // std::cout << "Outlier trajectory id: " << traj_idx << std::endl;
      continue;
    }

    std::cout << "---------------------------------" << std::endl;
    std::cout << "Created trajectory " << traj_idx << std::endl;
    std::cout << "Computing Dwell Regions now." << std::endl;
    OnlineDwellRegionComputation Online_DRC(
      gsl::not_null<Trajectory*>{ &traj }, num_heaps, max_window_size, query_radius);

    bool result = Online_DRC.has_dwell_region();
    auto perf_counters = Online_DRC.get_query_performance_counters();

    if (query_perf_file.is_open()) {
      if (write_header) {
        query_perf_file << perf_counters.get_header_line() << std::endl;
        write_header = false;
      }
      query_perf_file << perf_counters << std::endl;
    }

    std::cout << "Number of times SEC_S was called: " << perf_counters.num_SEC_computed << std::endl;
    std::cout << "Number of queries answered without computing SEC_S: "
              << perf_counters.num_query_answered_without_computing_SEC << std::endl;
    std::cout << "Number of points used to compute actual SEC (filtered points): " << perf_counters.num_filtered_points
              << std::endl;
    std::cout << "Total number of points in trajectory: " << traj.size() << std::endl;
    std::cout << "Number of dwell regions: " << Online_DRC.num_dwell_regions() << std::endl;
    std::cout << "Total time required to update heaps: " << perf_counters.heap_update_time_in_seconds << " seconds"
              << std::endl;
    std::cout << (result ? "Trajectory has a dwell region" : "Trajectory does not have a dwell region") << std::endl;
  }

  query_perf_file.close();
}

int
main(int argc, char** argv)
{
  //     try {
  namespace po = boost::program_options;

  po::options_description desc(
    "This program reads in a trajectory and tells you (true/false) if it has a dwell region");

  unsigned max_window_size = 0, num_heaps = 0;
  double query_radius = 0.0;
  size_t skip_N_trajectories = 0;
  std::string query_perf_filepath = "";

  desc.add_options()("filename,f", po::value<std::string>(), "Reads the trajectory from file")(
    "max_window_size,w", po::value<unsigned>(&max_window_size)->default_value(10), "Max window size")(
    "num_heaps,k", po::value<unsigned>(&num_heaps)->default_value(8), "Number of heaps")(
    "query_radius,q", po::value<double>(&query_radius)->default_value(0.5), "Query Radius (Rq)")(
    "skip_N_trajectories,s",
    po::value<size_t>(&skip_N_trajectories)->default_value(0),
    "Skip N trajectories from the file")("query_perf_filepath,o",
                                         po::value<std::string>(&query_perf_filepath)->default_value(""),
                                         "Output the query performance counters as a separate CSV file")("help,h",
                                                                                                         "Print usage");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

  std::string filename;
  if (vm.count("filename"))
    filename = vm["filename"].as<std::string>();

  auto filename_path = std::filesystem::path(filename);
  std::cout << "Opening " << filename << std::endl;

  process_dataset(filename_path, max_window_size, num_heaps, query_radius, skip_N_trajectories, query_perf_filepath);

  return 0;
}
