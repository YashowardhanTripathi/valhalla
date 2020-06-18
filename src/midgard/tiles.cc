#include <array>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>

#include "midgard/distanceapproximator.h"
#include "midgard/ellipse.h"
#include "midgard/polyline2.h"
#include "midgard/tiles.h"
#include "midgard/util.h"

namespace {

// this is modified to include all pixels that are intersected by the floating point line
// at each step it decides to either move in the x or y direction based on which pixels midpoint
// forms a smaller triangle with the line. to avoid edge cases we allow set_pixel to make the
// the loop bail if we leave the valid drawing region
void bresenham_line(double x0,
                    double y0,
                    double x1,
                    double y1,
                    const std::function<bool(int32_t, int32_t)>& set_pixel) {
  // this one for sure
  bool outside = set_pixel(std::floor(x0), std::floor(y0));
  // steps in the proper direction and constants for shoelace formula
  double sx = x0 < x1 ? 1 : -1, dx = x1 - x0, x = std::floor(x0) + .5f;
  double sy = y0 < y1 ? 1 : -1, dy = y1 - y0, y = std::floor(y0) + .5f;
  // keep going until we make it to the ending pixel
  while (std::floor(x) != std::floor(x1) || std::floor(y) != std::floor(y1)) {
    double tx = std::abs(dx * (y - y0) - dy * ((x + sx) - x0));
    double ty = std::abs(dx * ((y + sy) - y0) - dy * (x - x0));
    // less error moving in the x
    if (tx < ty) {
      x += sx;
    }
    // less error moving in the y
    else {
      y += sy;
    }
    // mark this pixel
    bool o = set_pixel(std::floor(x), std::floor(y));
    if (outside == false && o == true) {
      return;
    }
    outside = o;
  }
}

// a functor to generate closest first subdivisions of a set of tiles

} // namespace

namespace valhalla {
namespace midgard {

// Get the list of tiles that lie within the specified bounding box. The method finds the tile
// at the ellipse center. It successively finds neighbors and checks if they are inside or intersect
// with the ellipse.
template <class coord_t>
std::vector<int32_t> Tiles<coord_t>::TileList(const Ellipse<coord_t>& e) const {
  // Create a queue of tiles to check, initialize with the tile at the center of the ellipse
  int32_t tileid = TileId(e.center());
  std::queue<int32_t> check_queue;
  check_queue.push(tileid);

  // Record any tiles added to the check_queue
  std::unordered_set<int32_t> checked_tiles;
  checked_tiles.insert(tileid);

  // Successively check a tile from the queue - if tile bounds are not outside the ellipse then
  // add to the tile list and find its neighbors. Add them to the check list if not in the checked
  // tiles list.
  std::vector<int32_t> tile_list;
  while (!check_queue.empty()) {
    tileid = check_queue.front();
    check_queue.pop();

    // Test if the tile bounds is not outside the ellipse (if not the ellipse is inside the tile,
    // the tile is inside the ellipse, or the tile bounds intersects the ellipse).
    if (e.DoesIntersect(TileBounds(tileid)) != IntersectCase::kOutside) {
      tile_list.push_back(tileid);

      // Add neighboring tiles that have not already been checked
      std::array<int32_t, 4> neighbors = {BottomNeighbor(tileid), TopNeighbor(tileid),
                                          LeftNeighbor(tileid), RightNeighbor(tileid)};
      for (const auto nb : neighbors) {
        if (checked_tiles.find(nb) == checked_tiles.end()) {
          check_queue.push(nb);
          checked_tiles.insert(nb);
        }
      }
    }
  }
  return tile_list;
}

// Color a "connectivity map" starting with a sparse map of uncolored tiles.
// Any 2 tiles that have a connected path between them will have the same
// value in the connectivity map.
template <class coord_t>
void Tiles<coord_t>::ColorMap(std::unordered_map<uint32_t, size_t>& connectivity_map,
                              const std::unordered_map<uint32_t, uint32_t>& not_neighbors) const {
  // Connectivity map - all connected regions will have a unique Id. If any 2
  // tile Ids have a different Id they are judged to be not-connected.

  auto are_feuding = [&not_neighbors](uint32_t a, uint32_t b) {
    auto found = not_neighbors.find(a);
    if (found != not_neighbors.cend() && found->second == b)
      return true;
    found = not_neighbors.find(b);
    return found != not_neighbors.cend() && found->second == a;
  };

  // Iterate through tiles
  size_t color = 1;
  for (auto& tile : connectivity_map) {
    // Continue if already visited
    if (tile.second > 0) {
      continue;
    }

    // Mark this tile Id with the current color and find all its
    // accessible neighbors
    tile.second = color;
    std::unordered_set<uint32_t> checklist{tile.first};
    while (!checklist.empty()) {
      uint32_t next_tile = *checklist.begin();
      checklist.erase(checklist.begin());

      // Check neighbors.
      uint32_t neighbor = LeftNeighbor(next_tile);
      auto neighbor_itr = connectivity_map.find(neighbor);
      if (neighbor_itr != connectivity_map.cend() && neighbor_itr->second == 0 &&
          !are_feuding(next_tile, neighbor)) {
        checklist.emplace(neighbor);
        neighbor_itr->second = color;
      }
      neighbor = RightNeighbor(next_tile);
      neighbor_itr = connectivity_map.find(neighbor);
      if (neighbor_itr != connectivity_map.cend() && neighbor_itr->second == 0 &&
          !are_feuding(next_tile, neighbor)) {
        checklist.emplace(neighbor);
        neighbor_itr->second = color;
      }
      neighbor = TopNeighbor(next_tile);
      neighbor_itr = connectivity_map.find(neighbor);
      if (neighbor_itr != connectivity_map.cend() && neighbor_itr->second == 0 &&
          !are_feuding(next_tile, neighbor)) {
        checklist.emplace(neighbor);
        neighbor_itr->second = color;
      }
      neighbor = BottomNeighbor(next_tile);
      neighbor_itr = connectivity_map.find(neighbor);
      if (neighbor_itr != connectivity_map.cend() && neighbor_itr->second == 0 &&
          !are_feuding(next_tile, neighbor)) {
        checklist.emplace(neighbor);
        neighbor_itr->second = color;
      }
    }

    // Increment color
    color++;
  }
}

template <class coord_t>
template <class container_t>
std::unordered_map<int32_t, std::unordered_set<unsigned short>>
Tiles<coord_t>::Intersect(const container_t& linestring) const {
  std::unordered_map<int32_t, std::unordered_set<unsigned short>> intersection;
  // protect against empty linestring
  if (linestring.empty()) {
    return {};
  }

  // what to do when we want to mark a subdivision as containing a segment of this linestring
  const auto set_pixel = [this, &intersection](int32_t x, int32_t y) {
    // cant mark ones that are outside the valid range of tiles
    // TODO: wrap coordinates around x and y?
    if (x < 0 || y < 0 || x >= nsubdivisions_ * ncolumns_ || y >= nsubdivisions_ * nrows_) {
      return true;
    }
    // find the tile
    int32_t tile_column = x / nsubdivisions_;
    int32_t tile_row = y / nsubdivisions_;
    int32_t tile = tile_row * ncolumns_ + tile_column;
    // find the subdivision
    unsigned short subdivision = (y % nsubdivisions_) * nsubdivisions_ + (x % nsubdivisions_);
    intersection[tile].insert(subdivision);
    return false;
  };

  // if coord_t is spherical and the segment uv is sufficiently long then the geodesic along it
  // cannot be approximated with linear constructs so instead we resample it at a sufficiently
  // small interval so as to approximate the arc with piecewise linear segments
  container_t resampled;
  auto max_meters =
      std::max(1.f, subdivision_size_ * .25f *
                        DistanceApproximator::MetersPerLngDegree(linestring.front().second));
  if (coord_t::IsSpherical() && Polyline2<coord_t>::Length(linestring) > max_meters) {
    resampled = resample_spherical_polyline(linestring, max_meters, true);
  }

  // for each segment
  const auto& line = resampled.size() ? resampled : linestring;
  auto ui = line.cbegin(), vi = line.cbegin();
  while (vi != line.cend()) {
    // figure out what the segment is
    auto u = *ui;
    auto v = u;
    std::advance(vi, 1);
    if (vi != line.cend()) {
      v = *vi;
    } else if (line.size() > 1) {
      return intersection;
    }
    ui = vi;

    // figure out global subdivision start and end points
    auto x0 = (u.first - tilebounds_.minx()) / tilebounds_.Width() * ncolumns_ * nsubdivisions_;
    auto y0 = (u.second - tilebounds_.miny()) / tilebounds_.Height() * nrows_ * nsubdivisions_;
    auto x1 = (v.first - tilebounds_.minx()) / tilebounds_.Width() * ncolumns_ * nsubdivisions_;
    auto y1 = (v.second - tilebounds_.miny()) / tilebounds_.Height() * nrows_ * nsubdivisions_;

    int ix0 = std::floor(x0), ix1 = std::floor(x1);
    int iy0 = std::floor(y0), iy1 = std::floor(y1);
    int dx = ix0 - ix1, dy = iy0 - iy1;
    int ds = dx * dx + dy * dy;
    // its likely for our use case that its all in one cell
    if (ds == 0) {
      set_pixel(ix0, iy0);
    }
    // if not the next most likley thing is adjacent cells
    else if (ds == 1) {
      set_pixel(ix0, iy0);
      set_pixel(ix1, iy1);
    }
    // pretend the subdivisions are pixels and we are doing line rasterization
    else {
      bresenham_line(x0, y0, x1, y1, set_pixel);
    }
  }

  // give them back
  return intersection;
}

template <class coord_t>
std::unordered_map<int32_t, std::unordered_set<uint16_t>>
Tiles<coord_t>::Intersect(const AABB2<coord_t>& box) const {
  std::unordered_map<int32_t, std::unordered_set<uint16_t>> intersection;

  // to calculate the bounds within each tile, we first calculate all the
  // subdivisions (bins) which the bounding box covers in global space, and
  // then iterate over them to fill in each "pixel" or bin.
  const int32_t x_pixels = ncolumns_ * nsubdivisions_;
  const int32_t y_pixels = nrows_ * nsubdivisions_;
  // NOTE: multiply by pixels before dividing to keep as much precision as
  // possible.
  int32_t x0 = std::floor(((box.minx() - tilebounds_.minx()) * x_pixels) / tilebounds_.Width());
  int32_t y0 = std::floor(((box.miny() - tilebounds_.miny()) * y_pixels) / tilebounds_.Height());
  int32_t x1 = std::floor(((box.maxx() - tilebounds_.minx()) * x_pixels) / tilebounds_.Width());
  int32_t y1 = std::floor(((box.maxy() - tilebounds_.miny()) * y_pixels) / tilebounds_.Height());

  // clamp ranges to within the bounds of the tile.
  if (x0 < 0) {
    x0 = 0;
  }
  if (y0 < 0) {
    y0 = 0;
  }
  if (x1 >= x_pixels) {
    x1 = x_pixels - 1;
  }
  if (y1 >= y_pixels) {
    y1 = y_pixels - 1;
  }

  for (auto y = y0; y <= y1; ++y) {
    for (auto x = x0; x <= x1; ++x) {
      auto tile_id = (y / nsubdivisions_) * ncolumns_ + (x / nsubdivisions_);
      auto bin = (y % nsubdivisions_) * nsubdivisions_ + (x % nsubdivisions_);
      intersection[tile_id].insert(bin);
    }
  }

  return intersection;
}

// Explicit instantiation
/*
template class Tiles<Point2>;
template class Tiles<PointLL>;

template class std::unordered_map<int32_t, std::unordered_set<unsigned short>>
Tiles<Point2>::Intersect(const std::list<Point2>&) const;
template class std::unordered_map<int32_t, std::unordered_set<unsigned short>>
Tiles<PointLL>::Intersect(const std::list<PointLL>&) const;
template class std::unordered_map<int32_t, std::unordered_set<unsigned short>>
Tiles<Point2>::Intersect(const std::vector<Point2>&) const;
template class std::unordered_map<int32_t, std::unordered_set<unsigned short>>
Tiles<PointLL>::Intersect(const std::vector<PointLL>&) const;
*/

} // namespace midgard
} // namespace valhalla
