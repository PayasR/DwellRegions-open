#include "../include/util.hpp"
#include <random>

namespace Util {
namespace Geo {

using GeographicLib::Geodesic, DwellRegions::Point;

double
distance_in_metres_from_lon_lat(const Point& p1, const Point& p2)
{
  double dist = 0.0;
  const Geodesic& geod = Geodesic::WGS84();

  geod.Inverse(p1.Lat, p1.Lon, p2.Lat, p2.Lon, dist);

  return dist;
}

double
distance_in_miles_from_lon_lat(const Point& p1, const Point& p2)
{
  double dist = 0.0;
  const Geodesic& geod = Geodesic::WGS84();

  geod.Inverse(p1.Lat, p1.Lon, p2.Lat, p2.Lon, dist);

  return dist * DwellRegions::Constants::metres_to_miles; // Constant for metres
                                                          // to miles conversion
}

double
distance_in_miles_from_X_Y(const Point& p1, const Point& p2)
{
  double dist = std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
  return dist;
}
}

double
degrees_to_radians(const double& degrees)
{
  return degrees * DwellRegions::Constants::PI / 180;
}

double
radians_to_degrees(const double& radians)
{
  return radians * 180 / DwellRegions::Constants::PI;
}

// Returns a random integer in the range [low, high].
// random number generator from Stroustrup:
// http://www.stroustrup.com/C++11FAQ.html#std-random
int
rand_int(const int low, const int high)
{
  static std::default_random_engine re{};
  using dist = std::uniform_int_distribution<int>;
  static dist uid{};
  return uid(re, dist::param_type{ low, high });
}

size_t
rand_size_t(const size_t low, const size_t high)
{
  static std::default_random_engine re{};
  using dist = std::uniform_int_distribution<size_t>;
  static dist uid{};
  return uid(re, dist::param_type{ low, high });
}

double
rand_double(const double low, const double high)
{
  static std::default_random_engine re{};
  using dist = std::uniform_real_distribution<double>;
  static dist uid{};
  return uid(re, dist::param_type{ low, high });
}

}
