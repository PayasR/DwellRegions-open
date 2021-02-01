/// This file contains the values of all constants used in the code.
#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

namespace DwellRegions::Constants {
const double tolerance = 0.0001;
const double metres_to_miles = 0.000621;
const double maxdist = 0.0001; // I don't know why this should be any different from tolerance, but
                               // keeping this here for compatibility with Reaz's code

const double earth_radius_in_metres = 6371000;
const double earth_radius_in_miles = 3956.391;

const double central_meridian_in_degrees = 116.3975;   // center longitude of Beijing
const double standard_parallel_in_degrees = 39.906667; // center latitude of Beijing

const double central_meridian_in_radians = 2.031519616;
const double standard_parallel_in_radians = 0.696502732;

// Boundary of Beijing
const double boundary_min_lat = 39.447752;
const double boundary_max_lat = 41.060061;
const double boundary_min_lon = 115.430039;
const double boundary_max_lon = 117.506914;

// See https://www.sciencedirect.com/topics/engineering/statute-mile
// Path calculations are based on an average earth circumference of 24,857
// statute miles. Thus, each degree of arc along that surface represents
// 24,857/360 = 69.047 statute miles.
// In the north–south direction, the distance in miles is just the difference in
// the latitude values times 69.047. 1 degree of latitude = 69.047 miles
const double miles_to_latitude_degrees = 69.047;

// The C++ standard still doesn't have a standard definition of PI, so M_PI it is
// https://stackoverflow.com/questions/49778240/does-c11-14-17-or-20-introduce-a-standard-constant-for-pi
const double PI = 3.141592654;

} // namespace DwellRegions
