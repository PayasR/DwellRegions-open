/// This file contains the utility functions needed through the code.
/// Currently implemented:
/// 1. SphericalLatLonDistance
/// 2. Geodesic Distance
/// To implement later:
/// ChiSquareTest

#pragma once

#include "point.hpp"
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>

namespace Util {

namespace Geo {
using GeographicLib::Geodesic, GeographicLib::GeodesicLine, DwellRegions::Point;

/// Simple wrapper over GeographicLib, ripped off from old code
double
distance_in_metres_from_lon_lat(const Point& p1, const Point& p2);

double
distance_in_miles_from_lon_lat(const Point& p1, const Point& p2);

double
distance_in_miles_from_X_Y(const Point& p1, const Point& p2);
} // namespace WGS84

double
degrees_to_radians(const double& degrees);

double
radians_to_degrees(const double& radians);

// random number generator from Stroustrup:
// http://www.stroustrup.com/C++11FAQ.html#std-random
int
rand_int(const int low, const int high);

size_t
rand_size_t(const size_t low, const size_t high);

double
rand_double(const double low, const double high);
} // namespace Util
