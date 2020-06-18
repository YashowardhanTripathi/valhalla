#include "valhalla/midgard/point2.h"

#include <limits>
#include <list>

#include "midgard/util.h"
#include "midgard/vector2.h"

namespace valhalla {
namespace midgard {

// Explicit instantiations
template bool Point2::WithinPolygon(const std::vector<Point2>&) const;
template bool Point2::WithinPolygon(const std::list<Point2>&) const;

} // namespace midgard
} // namespace valhalla
