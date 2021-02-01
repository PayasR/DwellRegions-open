#include "../include/point.hpp"

namespace DwellRegions {

Point
Point_from_lon_lat(double lon, double lat)
{
  double Lon = lon; // in degrees
  double Lat = lat; // in degrees

  lon = lon * Constants::PI / 180.0; // in radians
  lat = lat * Constants::PI / 180.0; // in radians

  /*double x = Constants::earth_radius_in_miles * std::cos(lon) * std::cos(lat);
  double y = Constants::earth_radius_in_miles * std::sin(lon) * std::cos(lat);
  double z = Constants::earth_radius_in_miles * std::sin(lat);
  */

  // using Equirectangular projection
  // source: https://en.wikipedia.org/wiki/Equirectangular_projection
  double x = Constants::earth_radius_in_miles * (lon - Constants::central_meridian_in_radians) *
             std::cos(Constants::standard_parallel_in_radians);
  double y = Constants::earth_radius_in_miles * (lat - Constants::standard_parallel_in_radians);

  return Point{ Lon, Lat, x, y };
}

Point
Point_from_X_Y(double X, double Y)
{
  /*double Z = 3900;
  double Lon = std::atan2(Y, X) * 180.0 / Constants::PI;  //in degrees
  double Lat = std::asin(Z / Constants::earth_radius_in_miles)
        * 180.0 / Constants::PI;    //in degrees
  */

  // using Equirectangular projection (in radians)
  // source: https://en.wikipedia.org/wiki/Equirectangular_projection
  double lon = (X / (Constants::earth_radius_in_miles * std::cos(Constants::standard_parallel_in_radians))) +
               Constants::central_meridian_in_radians;
  double lat = (Y / Constants::earth_radius_in_miles) + Constants::standard_parallel_in_radians;
  // Convert lon, lat in degrees
  lon = lon * 180.0 / Constants::PI;
  lat = lat * 180.0 / Constants::PI;

  return Point{ lon, lat, X, Y };
}

bool
collinear(const Point& p1, const Point& p2, const Point& p3)
{
  if (p1 == p2 || p2 == p3 || p3 == p1) {
    return true;
  }

  double delta_X_1_2 = p1.x - p2.x;
  double delta_X_2_3 = p2.x - p3.x;

  if (std::fabs(delta_X_1_2) < Constants::tolerance || std::fabs(delta_X_2_3) < Constants::tolerance) {
    return std::fabs(delta_X_1_2 - delta_X_2_3) < Constants::tolerance;
  }

  double delta_Y_1_2 = p1.y - p2.y;
  double delta_Y_2_3 = p2.y - p3.y;
  double m_1_2 = delta_Y_1_2 / delta_X_1_2; // Slope of line joining first two points
  double m_2_3 = delta_Y_2_3 / delta_X_2_3; // Slope of line joining points 2 and 3

  return std::fabs(m_1_2 - m_2_3) < Constants::tolerance;
}
}
