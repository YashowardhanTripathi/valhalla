#ifndef VALHALLA_MIDGARD_POINT2_H_
#define VALHALLA_MIDGARD_POINT2_H_

#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include <valhalla/midgard/util_core.h>
namespace {
constexpr double LL_EPSILON = .00002f;
}

namespace valhalla {
namespace midgard {

// Forward references
template <typename PrecisionT> class VectorXY;

/**
 * 2D Point (cartesian). double x,y components.
 * @author David W. Nesbitt
 */
template <typename PrecisionT> class PointXY {

protected:
  PrecisionT m_x;
  PrecisionT m_y;

public:
  using x_type = PrecisionT;
  using y_type = PrecisionT;

  PointXY()
      : m_x{std::numeric_limits<PrecisionT>::max()}, m_y{std::numeric_limits<PrecisionT>::max()} {};
  PointXY(const PrecisionT x_, const PrecisionT y_) : m_x{x_}, m_y{y_} {
  }

  /**
   * Get the x component.
   * @return  Returns the x component of the point.
   */
  PrecisionT x() const {
    return m_x;
  }

  /**
   * Get the y component.
   * @return  Returns the y component of the point.
   */
  PrecisionT y() const {
    return m_y;
  }

  /**
   * Set the x component.
   * @param  x  x coordinate value.
   */
  void set_x(const PrecisionT x_) {
    m_x = x_;
  }

  /**
   * Set the y component.
   * @param  y  y coordinate value.
   */
  void set_y(const PrecisionT y_) {
    m_y = y_;
  }

  /**
   * Set the coordinate components to the specified values.
   * @param   x   x coordinate position.
   * @param   y   y coordinate position.
   */
  void Set(const PrecisionT x_, const PrecisionT y_) {
    m_x = x_;
    m_y = y_;
  }

  /**
   * Equality approximation.
   * @param   p  Point to compare to the current point.
   * @param   e  An epsilon which determines how close they must be to be considered equal
   * @return  Returns true if two points are approximately equal, false otherwise.
   */
  bool ApproximatelyEqual(const PointXY<PrecisionT>& p, double e = LL_EPSILON) const {
    return equal<PrecisionT>(m_x, p.x(), e) && equal<PrecisionT>(m_y, p.y(), e);
  }

  /**
   * Get the distance squared from this point to point p.
   * @param   p  Other point.
   * @return  Returns the distance squared between this point and p.
   */
  inline PrecisionT DistanceSquared(const PointXY<PrecisionT>& p) const {
    return sqr(m_x - p.x()) + sqr(m_y - p.y());
  }

  /**
   * Get the distance from this point to point p.
   * @param   p  Other point.
   * @return  Returns the distance between this point and p.
   */
  PrecisionT Distance(const PointXY<PrecisionT>& p) const {
    return sqrtf(DistanceSquared(p));
  }

  /**
   * Returns the point a specified percentage along a segment from this point
   * to an end point.
   * @param  end  End point.
   * @param  pct  Percentage along the segment.
   * @return Returns the point along the segment.
   */
  PointXY<PrecisionT> along_segment(const PointXY<PrecisionT>& end, const PrecisionT pct) const {
    return {x() + (end.x() - x()) * pct, y() + (end.y() - y()) * pct};
  }

  /**
   * Affine combination of this point with another point. 2 scalars are
   * provided (a0 and a1) and the must add to 1.
   * @param  a0  Scalar for this point
   * @param  a1  Scalar for p1
   * @param  p1  Point 1
   */
  PointXY<PrecisionT>
  AffineCombination(const PrecisionT a0, const PrecisionT a1, const PointXY<PrecisionT>& p1) const {
    return PointXY<PrecisionT>(a0 * x() + a1 * p1.x(), a0 * y() + a1 * p1.y());
  }
  /**
   * Gets the midpoint on a line segment between this point and point p1.
   * @param   p1  Point
   * @return  Returns the midpoint between this point and p1.
   */
  PointXY<PrecisionT> MidPoint(const PointXY<PrecisionT>& p1) const {
    return PointXY<PrecisionT>(0.5f * (x() + p1.x()), 0.5f * (y() + p1.y()));
  }

  /**
   * Add a vector to the current point.
   * @param   v  Vector to add to the current point.
   * @return  Returns a new point: the result of the current point
   *          plus the specified vector.
   */
  PointXY<PrecisionT> operator+(const VectorXY<PrecisionT>& v) const {
    return PointXY<PrecisionT>(x() + v.x(), y() + v.y());
  }

  /**
   * Subtract a vector from the current point.
   * @param   v  Vector to subtract from the current point.
   * @return  Returns a new point: the result of the current point
   *          minus the specified vector.
   */
  PointXY<PrecisionT> operator-(const VectorXY<PrecisionT>& v) const {
    return {x() - v.x(), y() - v.y()};
  }

  /**
   * Subtraction of a point from the current point.
   * @param p Point to subtract from the current point.
   * @return  Returns a vector.
   */
  VectorXY<PrecisionT> operator-(const PointXY<PrecisionT>& p) const {
    return VectorXY<PrecisionT>(x() - p.x(), y() - p.y());
  }

  /**
   * Equality
   **/
  bool operator==(const PointXY<PrecisionT>& other) const {
    return x() == other.x() && y() == other.y();
  }

  bool operator!=(const PointXY<PrecisionT>& other) const {
    return !(*this == other);
  }

  bool operator<(const PointXY<PrecisionT>& other) const {
    return x() < other.x() ? true : x() > other.x() ? false : y() < other.y();
  }

  /**
   * Finds the closest point to the supplied polyline as well as the distance
   * to that point and the index of the segment where the closest point lies.
   * @param  pts     List of points on the polyline.
   * @return  tuple of <Closest point along the polyline,
   *                    Returns the distance of the closest point,
   *                    Index of the segment of the polyline which contains the closest point
   *                   >
   */
  std::tuple<PointXY<PrecisionT>, PrecisionT, int>
  ClosestPoint(const std::vector<PointXY<PrecisionT>>& pts) const {
    PointXY<PrecisionT> closest;
    float mindist = std::numeric_limits<float>::max();

    // If there are no points we are done
    if (pts.size() == 0) {
      return std::make_tuple(std::move(closest), std::move(mindist), 0);
    }
    // If there is one point we are done
    if (pts.size() == 1) {
      return std::make_tuple(pts.front(), std::sqrt(DistanceSquared(pts.front())), 0);
    }

    // Iterate through the pts
    bool beyond_end = true;     // Need to test past the end point?
    int idx;                    // Index of closest segment so far
    VectorXY<PrecisionT> v1;    // Segment vector (v1)
    VectorXY<PrecisionT> v2;    // Vector from origin to target (v2)
    PointXY<PrecisionT> projpt; // Projected point along v1
    float dot;                  // Dot product of v1 and v2
    float comp;                 // Component of v2 along v1
    float dist;                 // Squared distance from target to closest point on line

    for (size_t index = 0; index < pts.size() - 1; ++index) {
      // Get the current segment
      const PointXY<PrecisionT>& p0 = pts[index];
      const PointXY<PrecisionT>& p1 = pts[index + 1];

      // Construct vector v1 - represents the segment.  Skip 0 length segments
      v1.Set(p0, p1);
      if (v1.x() == 0.0f && v1.y() == 0.0f) {
        continue;
      }

      // Vector v2 from the segment origin to the target point
      v2.Set(p0, *this);

      // Find the dot product of v1 and v2.  If less than 0 the segment
      // origin is the closest point.  Find the distance and continue
      // to the next segment.
      dot = v1.Dot(v2);
      if (dot <= 0.0f) {
        beyond_end = false;
        dist = DistanceSquared(p0);
        if (dist < mindist) {
          mindist = dist;
          closest = p0;
          idx = index;
        }
        continue;
      }

      // Closest point is either beyond the end of the segment or at a point
      // along the segment. Find the component of v2 along v1
      comp = dot / v1.Dot(v1);

      // If component >= 1.0 the segment end is the closest point. A future
      // polyline segment will be closer.  If last segment we need to check
      // distance to the endpoint.  Set flag so this happens.
      if (comp >= 1.0f) {
        beyond_end = true;
      } else {
        // Closest point is along the segment.  The closest point is found
        // by adding the projection of v2 onto v1 to the origin point.
        // The squared distance from this point to the target is then found.
        beyond_end = false;
        projpt = p0 + v1 * comp;
        dist = DistanceSquared(projpt);
        if (dist < mindist) {
          mindist = dist;
          closest = projpt;
          idx = index;
        }
      }
    }

    // Test the end point if flag is set - it may be the closest point
    if (beyond_end) {
      dist = DistanceSquared(pts.back());
      if (dist < mindist) {
        mindist = dist;
        closest = pts.back();
        idx = static_cast<int>(pts.size() - 2);
      }
    }
    return std::make_tuple(std::move(closest), std::move(std::sqrt(mindist)), std::move(idx));
  }

  /**
   * Test whether this point is to the left of a segment from p1 to p2.
   * @param  p1  First point of the segment.
   * @param  p2  End point of the segment.
   * @return  Returns positive if the point is left of segment p1,p2, or negative if it's to the
   * right.
   */
  virtual double IsLeft(const PointXY<PrecisionT>& p1, const PointXY<PrecisionT>& p2) const {
    return (p2.x() - p1.x()) * (y() - p1.y()) - (x() - p1.x()) * (p2.y() - p1.y());
  }

  /**
   * Tests whether this point is within a polygon.
   * @param  poly  List of vertices that form a polygon. Assumes
   *               the following:
   *                  Only the first and last vertices may be duplicated.
   * @return  Returns true if the point is within the polygon, false if not.
   */
  template <class container_t> bool WithinPolygon(const container_t& poly) const {
    auto p1 = poly.front() == poly.back() ? poly.begin() : std::prev(poly.end());
    auto p2 = poly.front() == poly.back() ? std::next(p1) : poly.begin();
    // for each edge
    size_t winding_number = 0;
    for (; p2 != poly.end(); p1 = p2, ++p2) {
      // going upward
      if (p1->y() <= y()) {
        // crosses if its in between on the y and to the left
        winding_number += p2->y() > y() && IsLeft(*p1, *p2) > 0;
      } // going downward maybe
      else {
        // crosses if its in between or on and to the right
        winding_number -= p2->y() <= y() && IsLeft(*p1, *p2) < 0;
      }
    }

    // If it was a full ring we are done otherwise check the last segment
    return winding_number != 0;
  }

  /**
   * Handy for templated functions that use both Point2 or PointLL to know whether or not
   * the classes coordinate system is spherical or in the plane
   *
   * @return true if the system is spherical false if not
   */
  static bool IsSpherical() {
    return false;
  }

  uint64_t hash() const {
    uint64_t h;
    char* b = static_cast<char*>(static_cast<void*>(&h));
    std::memcpy(b, &m_x, 4);
    std::memcpy(b + 4, &m_y, 4);
    return std::hash<uint64_t>()(h);
  }
};

using Point2 = PointXY<float>;

} // namespace midgard
} // namespace valhalla

namespace std {
template <> struct hash<valhalla::midgard::Point2> {
  size_t operator()(const valhalla::midgard::Point2& p) const {

    return p.hash();
  }
};
} // namespace std

#endif // VALHALLA_MIDGARD_POINT2_H_
