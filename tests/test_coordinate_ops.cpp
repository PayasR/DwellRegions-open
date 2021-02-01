#include "../include/xy_coordinate.hpp"

#include <iostream>
#include <random>
#include <vector>

using DwellRegions::Coordinate;

auto
generate_random_coordinates(size_t num_points)
{
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> distribution(0.0, 180.0);
  std::vector<Coordinate> coordinates;
  coordinates.reserve(num_points);

  for (size_t i = 0; i < num_points; i++) {
    // Use dis to transform the random unsigned int generated by gen into a
    // double in [1, 2). Each call to dis(gen) generates a new random double
    coordinates.emplace_back(DwellRegions::Coordinate::from_double(distribution(gen)));
  }

  return coordinates;
}

void
print_coordinates(const std::vector<Coordinate>& coordinates)
{
  for (const auto& c : coordinates) {
    std::cout << c.val() << " " << c.as_float() << std::endl;
  }
}

// Pairwise subtract the coordinates
void
subtraction_test(const std::vector<Coordinate>& coordinates)
{
  std::cout << "---------- Subtraction test--------------" << std::endl;
  for (size_t i = 1; i < coordinates.size(); i++) {
    std::cout << coordinates[i - 1].as_float() << " - " << coordinates[i].as_float() << "  =  ";
    const auto res = (coordinates[i - 1] - coordinates[i]).as_float();
    const auto float_res = coordinates[i-1].as_float() - coordinates[i].as_float();
    if (res != float_res) {
      std::cout << "Error: mismatch! res = " << res << " \t Float res = " << float_res  << std::endl;
    } else {
      std::cout << "Ok." << std::endl;
    }
  }
}

// Pairwise multiply the coordinates
void
multiplication_test(const std::vector<Coordinate>& coordinates)
{
  std::cout << "---------- Multiplication test--------------" << std::endl;
  for (size_t i = 1; i < coordinates.size(); i++) {
    std::cout << coordinates[i - 1].as_float() << " * " << coordinates[i].as_float() << "  =  ";
    const auto res = (coordinates[i - 1] * coordinates[i]).as_float();
    const auto float_res = coordinates[i-1].as_float() * coordinates[i].as_float();
    if (res != float_res) {
      std::cout << "Error: mismatch! res = " << res << " \t Float res = " << float_res  << std::endl;
    } else {
      std::cout << "Ok." << std::endl;
    }
  }
}

// Pairwise divide the coordinates
void
division_test(const std::vector<Coordinate>& coordinates)
{
  std::cout << "---------- Division test--------------" << std::endl;
  for (size_t i = 1; i < coordinates.size(); i++) {
    std::cout << coordinates[i - 1].as_float() << " / " << coordinates[i].as_float() << "  =  ";
    const auto res = (coordinates[i - 1] / coordinates[i]).as_float();
    const auto float_res = coordinates[i-1].as_float() / coordinates[i].as_float();
    if (res != float_res) {
      std::cout << "Error: mismatch! res = " << res << " \t Float res = " << float_res  << std::endl;
    } else {
      std::cout << "Ok." << std::endl;
    }
  }
}

int
main()
{
  const size_t num_coordinates = 10;
  const auto coordinates = generate_random_coordinates(num_coordinates);

  std::cout << "---------- Original Coordinates --------------" << std::endl;
  print_coordinates(coordinates);

  subtraction_test(coordinates);
  multiplication_test(coordinates);
  division_test(coordinates);

  return 0;
}