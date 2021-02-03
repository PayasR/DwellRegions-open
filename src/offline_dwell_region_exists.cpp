/// Implements a progam to read a trajectory and print if there exists a dwell region
#include "offline_dwell_region_computation.hpp"
#include "gsl-lite.hpp"

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <string>

void
process_offline_dwell_region_queries_using_rho_index(
                                  std::filesystem::path filename,
                                  unsigned max_window_size,
                                  unsigned num_heaps,
                                  double query_radius,
                                  unsigned num_bins,
                                  std::string data_path,
                                  bool create_index,
                                  std::string perf_filepath)
{
  using DwellRegions::OfflineDwellRegionRhoIndex;
  using DwellRegions::OfflineDwellRegionQueryRhoIndex;

  std::ofstream perf_file;
  if (perf_filepath != "") {
    perf_file.open(perf_filepath, std::fstream::app);
  }

  if (create_index) {
    OfflineDwellRegionRhoIndex Offline_DRC_index(num_heaps, num_bins);
    Offline_DRC_index.preprocess_trajectory_dataset(filename, data_path);

    if (perf_file.is_open()) {
      perf_file << num_heaps << "\t" << num_bins << "\t" << Offline_DRC_index.get_preprocess_time() << std::endl;
      perf_file.close();
    }
  }

  // If indexes exist, execute offline dwell region queries
  bool execute_query = create_index ? false : true;
  if (execute_query) {
    OfflineDwellRegionQueryRhoIndex Offline_DRC_query(
              num_heaps, max_window_size, query_radius, num_bins, data_path);
    Offline_DRC_query.execute_dwell_region_query();

    if (perf_file.is_open()) {
      perf_file << num_heaps << "\t" << num_bins << "\t" << query_radius << "\t"
                << Offline_DRC_query.get_num_dwell_regions() << "\t"
                << Offline_DRC_query.get_execution_time() << std::endl;
      perf_file.close();
    }
  }
}

void
process_offline_dwell_region_queries_using_tau_index(
                                  std::filesystem::path filename,
                                  unsigned max_window_size,
                                  unsigned num_heaps,
                                  double query_radius,
                                  unsigned num_bins,
                                  std::string data_path,
                                  bool create_index,
                                  std::string perf_filepath)
{
  using DwellRegions::OfflineDwellRegionTauIndex;
  using DwellRegions::OfflineDwellRegionQueryTauIndex;

  std::ofstream perf_file;
  if (perf_filepath != "") {
    perf_file.open(perf_filepath, std::fstream::app);
  }

  if (create_index) {
    OfflineDwellRegionTauIndex Offline_DRC_index(num_heaps, num_bins);
    Offline_DRC_index.preprocess_trajectory_dataset(filename, data_path);

    if (perf_file.is_open()) {
      perf_file << num_heaps << "\t" << num_bins << "\t" << Offline_DRC_index.get_preprocess_time() << std::endl;
      perf_file.close();
    }
  }

  // If indexes exist, execute offline dwell region queries
  bool execute_query = create_index ? false : true;
  if (execute_query) {
    OfflineDwellRegionQueryTauIndex Offline_DRC_query(
              num_heaps, max_window_size, query_radius, num_bins, data_path);
    Offline_DRC_query.execute_dwell_region_query();

    if (perf_file.is_open()) {
      perf_file << num_heaps << "\t" << num_bins << "\t" << query_radius << "\t"
                << Offline_DRC_query.get_num_dwell_regions() << "\t"
                << Offline_DRC_query.get_execution_time() << std::endl;
      perf_file.close();
    }
  }
}

void
test_offline_dwell_region_queries(unsigned max_window_size,
                                  unsigned num_heaps,
                                  double query_radius,
                                  unsigned num_bins,
                                  std::string data_path,
                                  std::string perf_filepath)
{
  using DwellRegions::OfflineDwellRegionQueryNoIndex;

  std::ofstream perf_file;
  if (perf_filepath != "") {
    perf_file.open(perf_filepath, std::fstream::app);
  }

  size_t num_trajectories = 9;
  OfflineDwellRegionQueryNoIndex Offline_DRC_query(
            num_heaps, max_window_size, query_radius, data_path);
  
  auto start = clock();
  for (size_t i = 0; i < num_trajectories; i++) {
    Offline_DRC_query.test_offline_query_per_traj(i);
  }
  auto end = clock();
  auto test_query_time = (end - start) / (double) CLOCKS_PER_SEC;
  //std::cout << "Execution time wilthout index: " << (end - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
  
  if (perf_file.is_open()) {
    perf_file << num_heaps << "\t" << num_bins << "\t" << query_radius << "\t"
              << test_query_time << std::endl;
    perf_file.close();
  }
}

int
main(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description desc(
    "This program reads trajectories and processes offline dwell region queries");

  unsigned max_window_size = 0, num_heaps = 0, num_bins = 0;
  double query_radius = 0.0;
  bool create_index = false;
  size_t index_type = 2; // 0 = No index, 1 = Rho-Index, 2 = Tau-Index
  std::string perf_filepath = "";
  // std::string filename = "";

  desc.add_options()("filename,f", po::value<std::string>(), "Reads the trajectories from file")(
    "max_window_size,w", po::value<unsigned>(&max_window_size)->default_value(10), "Max window size")(
    "num_heaps,k", po::value<unsigned>(&num_heaps)->default_value(8), "Number of heaps")(
    "query_radius,q", po::value<double>(&query_radius)->default_value(0.5), "Query Radius (Rq)")(
    "num_bins,b", po::value<unsigned>(&num_bins)->default_value(3), "Number of bins")(
    "create_index,i", po::value<bool>(&create_index)->default_value(false), "Want to create index? yes->true")(
    "index_type,t", po::value<size_t>(&index_type)->default_value(2), "Index type: 0 = No index, 1 = Rho-Index, 2 = Tau-Index")(
    "performance_file,o", po::value<std::string>(&perf_filepath)->default_value(""), "Output the offline query performance")(
    "help,h", "Print usage");

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

  // unsigned max_window_size = vm["max_window_size"].as<unsigned>();
  // unsigned num_heaps = vm["num_heaps"].as<unsigned>();
  // double query_radius = vm["query_radius"].as<double>();

  // Data path where the trajectories by id and preprocessed indexes are stored
  std::string data_path = filename.substr(0, filename.find("Ge"));

  auto filename_path = std::filesystem::path(filename);
  std::cout << "Opening " << filename << std::endl;

  if (index_type == 0) { // NO index, i.e., test queries
    test_offline_dwell_region_queries(max_window_size,
                                      num_heaps, query_radius,
                                      num_bins, data_path,
                                      perf_filepath);
  }
  else if (index_type == 1) { // Rho-Index
    process_offline_dwell_region_queries_using_rho_index(
                                      filename_path,
                                      max_window_size,
                                      num_heaps, query_radius,
                                      num_bins, data_path,
                                      create_index, perf_filepath);
  }
  else { // Tau-Index, this is the default one
    process_offline_dwell_region_queries_using_tau_index(
                                      filename_path,
                                      max_window_size,
                                      num_heaps, query_radius,
                                      num_bins, data_path,
                                      create_index, perf_filepath);
  }

  return 0;
}
