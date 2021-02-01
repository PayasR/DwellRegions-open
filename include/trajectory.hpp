#pragma once

#include "gsl-lite.hpp"
#include "point.hpp"

#include <cassert>
#include <vector>

namespace DwellRegions {

// This could probably have been using Trajectory = vector<Point>,
// The point of this wrapper class is to restrict the operations allowed on a Trajectory object.
// e.g. it shouldn't be possible to changes a point once it has been added to a trajectory
class Trajectory
{
  using points_t = std::vector<Point>;

public:
  using const_iterator = points_t::const_iterator;

  Trajectory() = default;
  Trajectory(const Trajectory&) = delete;
  Trajectory(const Trajectory&&) = delete;

  void emplace_back(Point&& point, size_t& timestamp)
  {
    points.emplace_back(point);
    timestamps.emplace_back(timestamp);
  }

  [[nodiscard]] size_t size() const { return points.size(); }
  [[nodiscard]] const_iterator begin() const { return points.begin(); }
  [[nodiscard]] const_iterator end() const { return points.end(); }

  [[nodiscard]] Point at(size_t idx) const
  {
    assert(idx < points.size());
    return points[idx];
  }

  [[nodiscard]] size_t timestamp_at(size_t idx) const
  {
    assert(idx < timestamps.size());
    return timestamps[idx];
  }

private:
  points_t points;
  std::vector<size_t> timestamps;
  friend struct SubTrajectory;
};

// Stores a continuous non-owning range of points in a Trajectory object
// Be CAREFUL when using a Subtrajectory object, since it's lifetime is dependent on the trajectory
// object underneath it So if you continue holding a SubTrajectory even after the Trajectory object
// is deleted, it IS a memory leak
struct SubTrajectory
{
  // Start = index of the first element of the trajectory to be included in the Subtrajectory
  // Size = number of elements to be included after the Start
  SubTrajectory(const Trajectory& trajectory, size_t start, size_t size)
    : points(&(*(trajectory.points.begin() + start)), size)
    , start_idx(start)
    , num_points(size)
  {
    assert(start < trajectory.points.size());
    assert((start + size) <= trajectory.points.size());
  }

  const gsl::span<const Point> points;
  const size_t start_idx;
  const size_t num_points;
};

}