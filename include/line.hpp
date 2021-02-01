/// This file contains the definition for an immutable Line class
#pragma once

#include "constants.hpp"
#include "point.hpp"

#include <optional>
// #include <iostream>

namespace DwellRegions {
class Line
{
public:
  // Defining a custom constructor to delete other default-generated
  // constructors
  Line(double a, double b, double c)
    : a(a)
    , b(b)
    , c(c)
  {}

  // Computes the point where two lines intersect.
  // If lines don't intersect, returns std::nullopt
  // If lines intersect, returns an std::optional<Point>,
  // from which the Point can be extracted with get()
  [[nodiscard]] std::optional<Point> intersects_with(const Line& line) const
  {
    double determinant = (a * line.b - line.a * b);
    //    std::cout << "Determinant: " << determinant << std::endl;

    if (std::fabs(determinant) < Constants::tolerance) {
      return std::nullopt;
    } else {
      /* return std::optional{ Point{
         (b * line.c - line.b * c) / determinant, // XLongitude
         (line.a * c - a * line.c) / determinant  // YLatitude
         */
      return Point_from_X_Y((b * line.c - line.b * c) / determinant, // XLongitude
                            (line.a * c - a * line.c) / determinant  // YLatitude
      );
    }
  }

private:
  // These are constants because a line should never be changed once it is
  // created
  const double a, b, c; // ax + by + c = 0
};
} // namespace DwellRegions
