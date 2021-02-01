#include <iostream>
#include "../include/point.hpp"
#include "../include/disk.hpp"

using DwellRegions::Point;
using DwellRegions::Disk;

std::vector<Point>
oriented(std::vector<Point> points)
{
  double center_x = 0.0, center_y = 0.0;
  for(size_t i = 0; i < points.size(); i++)
  {
    center_x += points[i].x;
    center_y += points[i].y;
  }
  center_x = center_x / points.size();
  center_y = center_y / points.size();

  std::vector<Point> oriented_points;
  std::vector<std::pair<double, Point>> points_with_angle;
  for(size_t i = 0; i < points.size(); i++)
  {
    double angle = std::atan2((points[i].y - center_y), (points[i].x - center_x));
    points_with_angle.push_back(std::make_pair(angle, points[i]));
  }
  std::sort(points_with_angle.begin(), points_with_angle.end(), DwellRegions::PointAnglePairComparator());
  for(size_t i = 0; i < points_with_angle.size(); i++)
    oriented_points.push_back(points_with_angle[i].second);
  
  return oriented_points;
}

std::optional<std::vector<Point>>
form_dwell_region(std::vector<Point> points, double radius)
{
  std::vector<Disk> circles;
  for (auto point: points)
  {
      Disk circle(point, radius);
      circles.push_back(circle);
  }

  if(circles.size() == 2)
  {
    return circles[0].disk_disk_intersection(circles[1]);
  }

  std::vector<Point> intersection_points;
  for (size_t i = 0; i < circles.size(); i++)
  {
      for(size_t j = i + 1; j < circles.size(); j++)
      {
          auto result  = circles[i].disk_disk_intersection(circles[j]);
          if (!result.has_value()) {
            std::cerr
              << "ERROR: Unable to compute disk-disk intersection. Computing dwell region aborted."
              << std::endl;  
          }

          auto intersect_points = result.value();
          for (size_t k = 0; k < intersect_points.size(); k++)
          {
              bool inside_point = true;
              for (size_t l = 0; l < circles.size(); l++)
              {
                  if(l == i || l == j)
                    continue;
                  if(!circles[l].in_disk(intersect_points[k]))
                  {
                      inside_point = false;
                      break;
                  }
              }
              if(inside_point)
                intersection_points.push_back(intersect_points[k]);
          }
      }
  }
  return oriented(intersection_points);
}

int
main()
{
  const double radius = 2.5;

  Point p1 = DwellRegions::Point_from_X_Y(-2.34879, 1.99522);
  Point p2 = DwellRegions::Point_from_X_Y(0.35988, -0.8169);
  Point p3 = DwellRegions::Point_from_X_Y(-0.02336, 3.70937);
  Point p4 = DwellRegions::Point_from_X_Y(2.10195, 1.93356);
  Point p5 = DwellRegions::Point_from_X_Y(-0.94115, -0.72205);
  std::vector<Point> points;
  points.push_back(p1);
  points.push_back(p2);
  
  std::cout << "Two points circle: " << std::endl;
  auto intersection_points = form_dwell_region(points, radius).value();
  for (auto point: intersection_points)
    std::cout << point.x << " " << point.y << std::endl;
  
  points.push_back(p3);
  std::cout << "Three points circle: " << std::endl;
  intersection_points = form_dwell_region(points, radius).value();
  for (auto point: intersection_points)
    std::cout << point.x << " " << point.y << std::endl;

  points.push_back(p4);
  points.push_back(p5);
  std::cout << "Five points on the circle: " << std::endl;
  intersection_points = form_dwell_region(points, radius).value();
  for (auto point: intersection_points)
    std::cout << point.x << " " << point.y << std::endl;

  return 0;
}