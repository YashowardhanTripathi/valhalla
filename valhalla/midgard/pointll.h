#ifndef VALHALLA_MIDGARD_POINTLL_H_
#define VALHALLA_MIDGARD_POINTLL_H_

#include <cmath>
#include <limits>
#include <tuple>
#include <valhalla/midgard/constants.h>
#include <valhalla/midgard/distanceapproximator.h>
#include <valhalla/midgard/point2.h>

namespace valhalla {
namespace midgard {

namespace {
constexpr double RAD_PER_DEG = valhalla::midgard::kPiDouble / 180.0;
constexpr double DEG_PER_RAD = 180.0 / valhalla::midgard::kPiDouble;
} // namespace

/**
 * Longitude, Latitude  point. Derives from Point2 and allows access methods
 * using lng,lat naming. Extends functionality to add heading, curvature,
 * and distance based on spherical geometry. Note that the order in the pair
 * is LONGITUDE first, LATITUDE second.
 */
template <typename PrecisionT> class GeoPoint : public PointXY<PrecisionT> {
public:
  /**
   * Use the constructors provided by pair
   */
  using PointXY<PrecisionT>::PointXY;
  using PointXY<PrecisionT>::x;
  using PointXY<PrecisionT>::y;

  /**
   * Default constructor.  Sets longitude and latitude to INVALID.
   */
  GeoPoint() : PointXY<PrecisionT>(INVALID, INVALID) {
  }

  /**
   * Parent constructor. Forwards to parent.
   */
  GeoPoint(const PointXY<PrecisionT>& p) : PointXY<PrecisionT>(p) {
  }

  /**
   * Get the longitude in degrees.
   * @return  Returns longitude in degrees.
   */
  float lng() const {
    return x();
  }

  /**
   * Get the latitude in degrees.
   * @return  Returns latitude (degrees).
   */
  float lat() const {
    return y();
  }

  /**
   * Checks for validity of the coordinates.
   * @return  Returns the false if lat or lon coordinates are set to INVALID.
   */
  bool IsValid() const {
    // TODO: is a range check appropriate?
    return x() != INVALID && y() != INVALID;
  }

  /**
   * Sets the coordinates to an invalid state
   */
  void Invalidate() {
    set_x(INVALID);
    set_y(INVALID);
  }

  /**
   * Gets the midpoint on a line segment between this point and point p1.
   * @param   p1  Point
   * @return  Returns the midpoint between this point and p1.
   */
  GeoPoint<PrecisionT> MidPoint(const GeoPoint<PrecisionT>& p) const {
    // radians
    const auto lon1 = lng() * -RAD_PER_DEG;
    const auto lat1 = lat() * RAD_PER_DEG;
    const auto lon2 = p.lng() * -RAD_PER_DEG;
    const auto lat2 = p.lat() * RAD_PER_DEG;
    // useful throughout
    const auto sl1 = sin(lat1);
    const auto sl2 = sin(lat2);
    const auto cl1 = cos(lat1);
    const auto cl2 = cos(lat2);
    // fairly accurate distance between points
    const auto d = acos(sl1 * sl2 + cl1 * cl2 * cos(lon1 - lon2));
    // interpolation parameters
    const auto ab = sin(d * .5) / sin(d);
    const auto acs1 = ab * cl1;
    const auto bcs2 = ab * cl2;
    // find the interpolated point along the arc
    const auto x = acs1 * cos(lon1) + bcs2 * cos(lon2);
    const auto y = acs1 * sin(lon1) + bcs2 * sin(lon2);
    const auto z = ab * (sl1 + sl2);
    return GeoPoint<PrecisionT>(atan2(y, x) * -DEG_PER_RAD,
                                atan2(z, sqrt(x * x + y * y)) * DEG_PER_RAD);
  }

  /**
   * Returns the point a specified percentage along a segment from this point
   * to an end point.
   * @param  end  End point.
   * @param  pct  Percentage along the segment.
   * @return Returns the point along the segment.
   */
  GeoPoint along_segment(const GeoPoint& end, const float pct) const {
    return {x() + (end.x() - x()) * pct, y() + (end.y() - y()) * pct};
  }

  /**
   * Calculates the distance between two lng,lat's in meters. Uses spherical
   * geometry (law of cosines).
   * @param   ll2   Second lng,lat position to calculate distance to.
   * @return  Returns the distance in meters.
   */
  PrecisionT Distance(const GeoPoint<PrecisionT>& ll2) const {
    // If points are the same, return 0
    if (*this == ll2) {
      return 0.0;
    }

    // Delta longitude. Don't need to worry about crossing 180
    // since cos(x) = cos(-x)
    double deltalng = (ll2.lng() - lng()) * RAD_PER_DEG;
    double a = lat() * RAD_PER_DEG;
    double c = ll2.lat() * RAD_PER_DEG;

    // Find the angle subtended in radians (law of cosines)
    double cosb = (sin(a) * sin(c)) + (cos(a) * cos(c) * cos(deltalng));

    // Angle subtended * radius of earth (portion of the circumference).
    // Protect against cosb being outside -1 to 1 range.
    if (cosb >= 1.0) {
      return 0.00001;
    } else if (cosb <= -1.0) {
      return kPi * kRadEarthMeters;
    } else {
      return acos(cosb) * kRadEarthMeters;
    }
  }

  /**
   * Approximates the distance squared between two lng,lat points - uses
   * the DistanceApproximator.
   * @param   ll2   Second lng,lat position to calculate distance to.
   * @return  Returns the distance squared in meters.
   */
  float DistanceSquared(const GeoPoint<PrecisionT>& ll2) const {
    DistanceApproximator approx(*this);
    return approx.DistanceSquared(ll2);
  }

  /**
   * Calculates the curvature using this position and 2 others. Found by
   * computing the radius of the circle that circumscribes the 3 positions.
   * @param   ll1   Second lng,lat position
   * @param   ll2   Third lng,lat position
   * @return  Returns the curvature in meters. Returns max float if the points
   *          are collinear.
   */
  PrecisionT Curvature(const GeoPoint<PrecisionT>& ll1, const GeoPoint<PrecisionT>& ll2) const {
    // Find the 3 distances between positions
    float a = Distance(ll1);
    float b = ll1.Distance(ll2);
    float c = Distance(ll2);
    float s = (a + b + c) * 0.5f;
    float k = sqrtf(s * (s - a) * (s - b) * (s - c));
    return (std::isnan(k) || k == 0.0f) ? std::numeric_limits<PrecisionT>::max()
                                        : ((a * b * c) / (4.0f * k));
  }

  /**
   * Calculates the heading or azimuth from the current lng,lat to the
   * specified lng,lat. This uses Haversine method (spherical geometry).
   * @param    ll2   lng,lat position to calculate the heading to.
   * @return   Returns the heading in degrees with range [0,360] where 0 is
   *           due north, 90 is east, 180 is south, and 270 is west.
   */
  float Heading(const GeoPoint<PrecisionT>& ll2) const {
    // If points are the same, return 0
    if (*this == ll2) {
      return 0.0;
    }

    // Convert to radians and compute heading
    double lat1 = lat() * kRadPerDegD;
    double lat2 = ll2.lat() * kRadPerDegD;
    double dlng = (ll2.lng() - lng()) * kRadPerDegD;
    double y = sin(dlng) * cos(lat2);
    double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlng);
    double bearing = atan2(y, x) * kDegPerRadD;
    return (bearing < 0.0) ? bearing + 360.0 : bearing;
  }

  /**
   * Finds the closest point to the supplied polyline as well as the distance
   * to that point and the (floor) index of the segment where the closest
   * point lies. In the case of a tie where the closest point is a point in the
   * linestring, the most extreme index (closest to the end of the linestring
   * in the direction (forward/reverse) of the search) will win
   * @param  pts                  List of points on the polyline.
   * @param  pivot_index          Index where the processing of closest point should start.
   *                              Default value is 0.
   * @param  forward_dist_cutoff  Minimum linear distance along pts that should be considered
   *                              before giving up.
   * @param  reverse_dist_cutoff  Minimum linear distance along pts that should be considered
   *                              before giving up.
   *
   * @return tuple of <Closest point along the polyline,
   *                   Distance in meters of the closest point,
   *                   Index of the segment of the polyline which contains the closest point >
   */
  std::tuple<GeoPoint<PrecisionT>, float, int>
  ClosestPoint(const std::vector<GeoPoint<PrecisionT>>& pts,
               int pivot_index = 0,
               float forward_dist_cutoff = std::numeric_limits<float>::infinity(),
               float reverse_dist_cutoff = 0) const {

    // setup
    if (pts.empty() || pivot_index < 0 || pivot_index > pts.size() - 1)
      return std::make_tuple(GeoPoint<PrecisionT>(), std::numeric_limits<float>::max(), -1);

    int closest_segment = pivot_index;
    GeoPoint<PrecisionT> closest = pts[pivot_index];
    DistanceApproximator approx(*this);
    float mindistsqr = approx.DistanceSquared(closest);

    // start going backwards, then go forwards
    for (bool reverse : {true, false}) {
      // get the range and distance for this direction
      auto dist_cutoff = reverse ? reverse_dist_cutoff : forward_dist_cutoff;
      int increment = reverse ? -1 : 1;
      int indices = reverse ? pivot_index : (pts.size() - 1) - pivot_index;

      GeoPoint<PrecisionT> point;
      for (int index = pivot_index - reverse; indices > 0 && dist_cutoff > 0.f;
           index += increment, --indices) {
        // Get the current segment
        const GeoPoint<PrecisionT>& u = pts[index];
        const GeoPoint<PrecisionT>& v = pts[index + 1];

        // Project a onto b where b is the origin vector representing this segment
        // and a is the origin vector to the point we are projecting, (a.b/b.b)*b
        auto bx = v.lng() - u.lng();
        auto by = v.lat() - u.lat();

        // Scale longitude when finding the projection. Avoid divided-by-zero
        // which gives a NaN scale, otherwise comparisons below will fail
        auto bx2 = bx * approx.GetLngScale();
        auto sq = bx2 * bx2 + by * by;
        auto scale =
            sq > 0 ? (((lng() - u.lng()) * approx.GetLngScale() * bx2 + (lat() - u.lat()) * by) / sq)
                   : 0.f;

        // Projects along the ray before u
        bool right_most = false;
        if (scale <= 0.f) {
          point = {u.lng(), u.lat()};
        } // Projects along the ray after v
        else if (scale >= 1.f) {
          point = {v.lng(), v.lat()};
          right_most = true;
        } // Projects along the ray between u and v
        else {
          point = {u.lng() + bx * scale, u.lat() + by * scale};
        }

        // Check if this point is better
        const auto sq_distance = approx.DistanceSquared(point);
        if (sq_distance < mindistsqr) {
          closest_segment = index + right_most;
          mindistsqr = sq_distance;
          closest = std::move(point);
        }

        // Check if we should bail early because of looking at too much shape
        if (dist_cutoff != std::numeric_limits<float>::infinity())
          dist_cutoff -= u.Distance(v);
      }
    }

    // give back what we found
    return std::make_tuple(std::move(closest), std::sqrt(mindistsqr), closest_segment);
  }

  /**
   * Calculate the heading from the start index within a polyline of lng,lat
   * points to a point at the specified distance from the start.
   * @param  pts   Polyline - list of lng,lat points.
   * @param  dist  Distance in meters from start to find heading to.
   * @param  idx0  Start index within the polyline.
   * @param  idx1  End index within the polyline
   */
  static float HeadingAlongPolyline(const std::vector<GeoPoint<PrecisionT>>& pts,
                                    const float dist,
                                    const uint32_t idx0,
                                    const uint32_t idx1) {
    // Check that at least 2 points exist
    int n = static_cast<int>(idx1) - static_cast<int>(idx0);
    if (n < 1) {
      // LOG_ERROR("PointLL::HeadingAlongPolyline has < 2 vertices");
      return 0.0f;
    }

    // If more than 2 points, walk edges of the polyline until the length
    // is exceeded.
    if (n > 1) {
      double d = 0.0;
      double seglength = 0.0;
      auto pt0 = pts.begin() + idx0;
      auto pt1 = pt0 + 1;
      while (d < dist && pt1 <= pts.begin() + idx1) {
        seglength = pt0->Distance(*pt1);
        if (d + seglength > dist) {
          // Set the extrapolated point along the line.
          float pct = static_cast<float>((dist - d) / seglength);
          GeoPoint<PrecisionT> ll(pt0->lng() + ((pt1->lng() - pt0->lng()) * pct),
                                  pt0->lat() + ((pt1->lat() - pt0->lat()) * pct));
          return pts[idx0].Heading(ll);
        } else {
          d += seglength;
          pt0++;
          pt1++;
        }
      }
    }

    // Only 2 points or the length of polyline is less than the specified
    // distance. Return heading from first to last point.
    return pts[idx0].Heading(pts[idx1]);
  }

  /**
   * Calculate the heading from the start of a polyline of lng,lat points to a
   * point at the specified distance from the start.
   * @param  pts   Polyline - list of lng,lat points.
   * @param  dist  Distance in meters from start to find heading to.
   */
  static float HeadingAlongPolyline(const std::vector<GeoPoint<PrecisionT>>& pts, const float dist) {
    return HeadingAlongPolyline(pts, dist, 0, pts.size() - 1);
  }

  /**
   * Calculate the heading from a point at a specified distance from the end
   * of a polyline of lng,lat points to the end point of the polyline.
   * @param  pts   Polyline - list of lng,lat points.
   * @param  dist  Distance in meters from end. A point that distance is
   *               used to find the heading to the end point.
   * @param  idx0  Start index within the polyline.
   * @param  idx1  End index within the polyline
   */
  static float HeadingAtEndOfPolyline(const std::vector<GeoPoint<PrecisionT>>& pts,
                                      const float dist,
                                      const uint32_t idx0,
                                      const uint32_t idx1) {
    // Check that at least 2 points exist
    int n = static_cast<int>(idx1) - static_cast<int>(idx0);
    if (n < 1) {
      // LOG_ERROR("PointLL::HeadingAtEndOfPolyline has < 2 vertices");
      return 0.0f;
    }

    // If more than 2 points, walk edges of the polyline until the length
    // is exceeded.
    if (n > 1) {
      double d = 0.0;
      double seglength;
      auto pt1 = pts.begin() + idx1;
      auto pt0 = pt1 - 1;
      while (d < dist && pt0 >= pts.begin() + idx0) {
        seglength = pt0->Distance(*pt1);
        if (d + seglength > dist) {
          // Set the extrapolated point along the line.
          float pct = static_cast<float>((dist - d) / seglength);
          GeoPoint<PrecisionT> ll(pt1->lng() + ((pt0->lng() - pt1->lng()) * pct),
                                  pt1->lat() + ((pt0->lat() - pt1->lat()) * pct));
          return ll.Heading(pts[idx1]);
        } else {
          d += seglength;
          pt1--;
          pt0--;
        }
      }
    }

    // Only 2 points or the length of polyline is less than the specified
    // distance. Return heading from first to last point.
    return pts[idx0].Heading(pts[idx1]);
  }

  /**
   * Calculate the heading from a point at a specified distance from the end
   * of a polyline of lng,lat points to the end point of the polyline.
   * @param  pts   Polyline - list of lng,lat points.
   * @param  dist  Distance in meters from end. A point that distance is
   *               used to find the heading to the end point.
   */
  static float HeadingAtEndOfPolyline(const std::vector<GeoPoint<PrecisionT>>& pts,
                                      const float dist) {
    return HeadingAtEndOfPolyline(pts, dist, 0, pts.size() - 1);
  }

  /**
   * Test whether this point is to the left of a segment from p1 to p2.
   * @param  p1  First point of the segment.
   * @param  p2  End point of the segment.
   * @return  Returns true if this point is left of the segment.
   */
  virtual float IsLeft(const GeoPoint<PrecisionT>& p1, const GeoPoint<PrecisionT>& p2) const {
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
      if (p1->lat() <= lat()) {
        // crosses if its in between on the y and to the left
        winding_number += p2->lat() > lat() && IsLeft(*p1, *p2) > 0;
      } // going downward maybe
      else {
        // crosses if its in between or on and to the right
        winding_number -= p2->lat() <= lat() && IsLeft(*p1, *p2) < 0;
      }
    }

    // If it was a full ring we are done otherwise check the last segment
    return winding_number != 0;
  };

  /**
   * Handy for templated functions that use both Point2 or PointLL to know whether or not
   * the classes coordinate system is spherical or in the plane
   *
   * @return true if the system is spherical false if not
   */
  static bool IsSpherical() {
    return true;
  };

  /**
   * Project this point onto the line from u to v
   * @param u          first point of segment
   * @param v          second point of segment
   * @return p  the projected point of this onto the segment uv
   */
  GeoPoint<PrecisionT> Project(const GeoPoint<PrecisionT>& u, const GeoPoint<PrecisionT>& v) const {
    auto lon_scale = cosf(y() * kRadPerDeg);
    return Project(u, v, lon_scale);
  }

  /**
   * Project this point onto the line from u to v
   * @param u          first point of segment
   * @param v          second point of segment
   * @param lon_scale  needed for spherical projections. dont pass this parameter unless
   *                   you cached it and want to avoid trig functions in a tight loop
   * @return p  the projected point of this onto the segment uv
   */
  GeoPoint<PrecisionT>
  Project(const GeoPoint<PrecisionT>& u, const GeoPoint<PrecisionT>& v, float lon_scale) const {
    // we're done if this is a zero length segment
    if (u == v) {
      return u;
    }

    // project a onto b where b is the origin vector representing this segment
    // and a is the origin vector to the point we are projecting, (a.b/b.b)*b
    auto bx = v.lng() - u.lng();
    auto by = v.lat() - u.lat();

    // Scale longitude when finding the projection
    auto bx2 = bx * lon_scale;
    auto sq = bx2 * bx2 + by * by;
    auto scale = (lng() - u.lng()) * lon_scale * bx2 +
                 (lat() - u.lat()) * by; // only need the numerator at first

    // projects along the ray before u
    if (scale <= 0.f) {
      return u;
      // projects along the ray after v
    } else if (scale >= sq) {
      return v;
    }
    // projects along the ray between u and v
    scale /= sq;
    return {u.lng() + bx * scale, u.lat() + by * scale};
  };

  /**
   * Project this point to the supplied polyline as well as the distance
   * to that point and the (floor) index of the segment where the projected location lies.
   * @param  pts                  List of points on the polyline.
   * @return tuple of <Closest point along the polyline,
   *                   Distance in meters of the closest point,
   *                   Index of the segment of the polyline which contains the closest point >
   */
  std::tuple<GeoPoint<PrecisionT>, float, int>
  Project(const std::vector<GeoPoint<PrecisionT>>& pts) const {
    auto u = pts.begin();
    auto v = pts.begin();
    std::advance(v, 1);

    auto min_distance = std::numeric_limits<float>::max();
    auto best = GeoPoint<PrecisionT>{};
    int best_index = 0;
    while (v != pts.end()) {
      auto candidate = Project(*u, *v);
      auto distance = Distance(candidate);
      if (distance < min_distance) {
        min_distance = distance;
        best = candidate;
        best_index = std::distance(pts.begin(), u);
      }
      u++;
      v++;
    }
    return std::make_tuple(best, min_distance, best_index);
  }

  /*
    uint64_t hash() {
      uint64_t h;
      char* b = static_cast<char*>(static_cast<void*>(&h));
      std::memcpy(b, &p.x_, 4);
      std::memcpy(b + 4, &p.y_, 4);
      return std::hash<uint64_t>()(h);
    }
    */

private:
  static constexpr PrecisionT INVALID = 0xBADBADBAD;
};

using PointLL = GeoPoint<float>;

} // namespace midgard
} // namespace valhalla

namespace std {
template <> struct hash<valhalla::midgard::PointLL> {
  size_t operator()(const valhalla::midgard::PointLL& p) const {
    return p.hash();
  }
};
} // namespace std

#endif // VALHALLA_MIDGARD_POINTLL_H_
