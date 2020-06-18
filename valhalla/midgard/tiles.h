
#ifndef VALHALLA_MIDGARD_TILES_H_
#define VALHALLA_MIDGARD_TILES_H_

#include <cstdint>
#include <functional>
#include <list>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <valhalla/midgard/aabb2.h>
#include <valhalla/midgard/constants.h>
#include <valhalla/midgard/ellipse.h>

namespace valhalla {
namespace midgard {

template <typename coord_t> class Tiles;

namespace {
template <class coord_t> struct closest_first_generator_t {
  coord_t seed;
  valhalla::midgard::Tiles<coord_t> tiles;
  int32_t subcols, subrows;
  std::unordered_set<int32_t> queued;
  using best_t = std::pair<double, int32_t>;
  std::set<best_t, std::function<bool(const best_t&, const best_t&)>> queue;

  closest_first_generator_t(const valhalla::midgard::Tiles<coord_t>& tiles, const coord_t& seed)
      : tiles(tiles), seed(seed), queued(100), queue([](const best_t& a, const best_t& b) {
          return a.first == b.first ? a.second < b.second : a.first < b.first;
        }) {
    // what global subdivision are we starting in
    // TODO: worry about wrapping around valid range
    subcols = tiles.ncolumns() * tiles.nsubdivisions();
    subrows = tiles.nrows() * tiles.nsubdivisions();
    auto x = (seed.x() - tiles.TileBounds().minx()) / tiles.TileBounds().Width() * subcols;
    auto y = (seed.y() - tiles.TileBounds().miny()) / tiles.TileBounds().Height() * subrows;
    auto subdivision = static_cast<int32_t>(y) * subcols + static_cast<int32_t>(x);
    queued.emplace(subdivision);
    queue.emplace(std::make_pair(0, subdivision));
    neighbors(subdivision);
  }

  // something to measure the closest possible point of a subdivision from the given seed point
  double dist(int32_t sub) {
    auto x = sub % subcols;
    auto x0 = static_cast<typename coord_t::x_type>(tiles.TileBounds().minx() +
                                                    x * tiles.SubdivisionSize());
    auto x1 = static_cast<typename coord_t::x_type>(tiles.TileBounds().minx() +
                                                    (x + 1) * tiles.SubdivisionSize());
    auto y = sub / subcols;
    auto y0 = static_cast<typename coord_t::y_type>(tiles.TileBounds().miny() +
                                                    y * tiles.SubdivisionSize());
    auto y1 = static_cast<typename coord_t::y_type>(tiles.TileBounds().miny() +
                                                    (y + 1) * tiles.SubdivisionSize());
    auto distance = std::numeric_limits<double>::max();
    std::list<coord_t> corners{{x0, y0}, {x1, y0}, {x0, y1}, {x1, y1}};
    if (x0 < seed.x() && x1 > seed.x()) {
      corners.emplace_back(seed.x(), y0);
      corners.emplace_back(seed.x(), y1);
    }
    if (y0 < seed.y() && y1 > seed.y()) {
      corners.emplace_back(x0, seed.y());
      corners.emplace_back(x1, seed.y());
    }
    for (const auto& c : corners) {
      auto d = seed.Distance(c);
      if (d < distance) {
        distance = d;
      }
    }
    return distance;
  }

  // something to add the neighbors of a given subdivision
  const std::list<std::pair<int, int>> neighbor_offsets{{0, -1}, {-1, 0}, {1, 0}, {0, 1}};
  void neighbors(int32_t s) {
    // walk over all adjacent subdivisions in row major order
    auto x = s % subcols;
    auto y = s / subcols;
    for (const auto& off : neighbor_offsets) {
      // skip y out of bounds
      auto ny = y + off.second;
      if (ny == -1 || ny == subrows) {
        continue;
      }
      // fix x
      auto nx = x + off.first;
      if (nx == -1 || nx == subcols) {
        if (!coord_t::IsSpherical()) {
          continue;
        }
        nx = (nx + subcols) % subcols;
      }
      // actually add the thing
      auto neighbor = ny * subcols + nx;
      if (queued.find(neighbor) == queued.cend()) {
        queued.emplace(neighbor);
        queue.emplace(std::make_pair(dist(neighbor), neighbor));
      }
    }
  }

  // get the next closest subdivision
  std::tuple<int32_t, unsigned short, float> next() {
    // get the next closest one or bail
    if (!queue.size()) {
      throw std::runtime_error("Subdivisions were exhausted");
    }
    auto best = *queue.cbegin();
    queue.erase(queue.cbegin());
    // add its neighbors
    neighbors(best.second);
    // return it
    auto sx = best.second % subcols;
    auto sy = best.second / subcols;
    auto tile_column = sx / tiles.nsubdivisions();
    auto tile_row = sy / tiles.nsubdivisions();
    auto tile = tile_row * tiles.ncolumns() + tile_column;
    unsigned short subdivision = (sy - tile_row * tiles.nsubdivisions()) * tiles.nsubdivisions() +
                                 (sx - tile_column * tiles.nsubdivisions());
    return std::make_tuple(tile, subdivision, best.first);
  }
};

} // namespace

/**
 * A class that provides a uniform (square) tiling system for a specified
 * bounding box and tile size. This is a template class that works with
 * Point2 (Euclidean x,y) or PointLL (latitude,longitude).
 * A unique tile Id is assigned for each tile based on the following rules:
 *    Tile numbers start at 0 at the min y, x (lower left)
 *    Tile numbers increase by column (x,longitude) then by row (y,latitude)
 *    Tile numbers increase along each row by increasing x,longitude.
 * Contains methods for converting x,y or lat,lng into tile Id and
 * vice-versa.  Methods for relative tiles (using row and column offsets).
 * are also provided. Also includes a method to get a list of tiles covering
 * a bounding box.
 */
template <class coord_t> class Tiles {
public:
  /**
   * Constructor. A bounding box and tile size is specified.
   * Sets class data members and computes the tile size
   * based on the bounding box and rows and columns.
   * @param   bounds       Bounding box
   * @param   tilesize     Size of the tile in both dimensions
   * @param   subdivisions Number of subtiles in both x and y axis of a single tile
   * @param   wrapx        Should neighbor operations wrap around the x axis extents
   */
  Tiles(const AABB2<coord_t>& bounds,
        const double tilesize,
        const unsigned short subdivisions = 1,
        bool wrapx = true)
      : tilebounds_(bounds), tilesize_(tilesize), nsubdivisions_(subdivisions),
        subdivision_size_(tilesize_ / nsubdivisions_), wrapx_(wrapx) {
    auto columns = bounds.Width() / tilesize_;
    auto rows = bounds.Height() / tilesize_;
    // TODO: delete this constructor and force use of the lower one
    // this is not safe because tilesize may not evenly divide into the bounds dimensions
    ncolumns_ = static_cast<int32_t>(std::round(columns));
    nrows_ = static_cast<int32_t>(std::round(rows));
  }

  /**
   * Constructor. A bottom left coord, with tile_size and number of rows and columns
   * @param   min_pt       Bottom left coord of the tileset
   * @param   tile_size    The size of a tile in both dimensions
   * @param   columns      Number of tiles in the x axis
   * @param   rows         Number of tiles in the y axis
   * @param   subdivisions Number of subtiles in both x and y axis of a single tile
   * @param   wrapx        Should neighbor operations wrap around the x axis extents
   */
  Tiles(const coord_t& min_pt,
        const double tile_size,
        const int32_t columns,
        const int32_t rows,
        const unsigned short subdivisions = 1,
        bool wrapx = true)
      : tilebounds_(min_pt,
                    coord_t{min_pt.first + columns * tile_size, min_pt.second + rows * tile_size}),
        tilesize_(tile_size), ncolumns_(columns), nrows_(rows),
        subdivision_size_(tilesize_ / nsubdivisions_), nsubdivisions_(subdivisions), wrapx_(wrapx) {
  }

  /**
   * Get the tile size.
   * @return Tile size.
   */
  double TileSize() const {
    return tilesize_;
  }

  /**
   * Get the tile subdivision size.
   * @return tile subdivision size.
   */
  double SubdivisionSize() const {
    return subdivision_size_;
  }

  /**
   * Get the number of rows in the tiling system.
   * @return  Returns the number of rows.
   */
  int32_t nrows() const {
    return nrows_;
  }

  /**
   * Get the number of columns in the tiling system.
   * @return  Returns the number of columns.
   */
  int32_t ncolumns() const {
    return ncolumns_;
  }

  /**
   * Get the number of subdivisions in a tile in the tiling system.
   * @return  Returns the number of subdivisions.
   */
  unsigned short nsubdivisions() const {
    return nsubdivisions_;
  }

  /**
   * Get the bounding box of the tiling system.
   * @return Bounding box.
   */
  AABB2<coord_t> TileBounds() const {
    return tilebounds_;
  }

  /**
   * Shift the tilebounds - a special method used to nudge the tile bounds
   * so a specific point stays centered in the grid.
   * @param  shift  Amount to shift the bounding box.
   */
  void ShiftTileBounds(const coord_t& shift) {
    tilebounds_ = AABB2<coord_t>(tilebounds_.minx() - shift.x(), tilebounds_.miny() - shift.y(),
                                 tilebounds_.maxx() - shift.x(), tilebounds_.maxy() - shift.y());
  }

  /**
   * Get the "row" based on y.
   * @param   y   y coordinate
   * @return  Returns the tile row. Returns -1 if outside the tile system bounds.
   */
  int32_t Row(const double y) const {
    // Return -1 if outside the tile system bounds
    if (y < tilebounds_.miny() || y > tilebounds_.maxy()) {
      return -1;
    }

    // If equal to the max y return the largest row
    return (y == tilebounds_.maxy()) ? nrows_ - 1
                                     : static_cast<int32_t>((y - tilebounds_.miny()) / tilesize_);
  }

  /**
   * Get the "column" based on x.
   * @param   x   x coordinate
   * @return  Returns the tile column. Returns -1 if outside the tile system bounds.
   */
  int32_t Col(const double x) const {
    // Return -1 if outside the tile system bounds
    if (x < tilebounds_.minx() || x > tilebounds_.maxx()) {
      return -1;
    }

    // If equal to the max x return the largest column
    if (x == tilebounds_.maxx()) {
      return ncolumns_ - 1;
    } else {
      double col = (x - tilebounds_.minx()) / tilesize_;
      return (col >= 0.0) ? static_cast<int32_t>(col) : static_cast<int32_t>(col - 1);
    }
  }

  /**
   * Convert a coordinate into a tile Id. The point is within the tile.
   * @param   c   Coordinate / point.
   * @return  Returns the tile Id. If the coordinate is outside tiling system
   *          extent, an error (-1) is returned.
   */
  int32_t TileId(const coord_t& c) const {
    return TileId(c.y(), c.x());
  }

  /**
   * Convert x,y to a tile Id.
   * @param   y   y (or lat)
   * @param   x   x (or lng)
   * @return  Returns the tile Id. -1 (error is returned if the x,y is
   *          outside the bounding box of the tiling sytem).
   */
  int32_t TileId(const double y, const double x) const {
    // Return -1 if totally outside the extent.
    if (y < tilebounds_.miny() || x < tilebounds_.minx() || y > tilebounds_.maxy() ||
        x > tilebounds_.maxx()) {
      return -1;
    }

    // Find the tileid by finding the latitude row and longitude column
    return (Row(y) * ncolumns_) + Col(x);
  }

  /**
   * Get the tile Id given the row Id and column Id.
   * @param  col  Tile column.
   * @param  row  Tile row.
   * @return  Returns the tile Id.
   */
  int32_t TileId(const int32_t col, const int32_t row) const {
    return (row * ncolumns_) + col;
  }

  /**
   * Get the tile row, col based on tile Id.
   * @param  tileid  Tile Id.
   * @return  Returns a pair indicating {row, col}
   */
  std::pair<int32_t, int32_t> GetRowColumn(const int32_t tileid) const {
    return {tileid / ncolumns_, tileid % ncolumns_};
  }

  /**
   * Get a maximum tileid given a bounds and a tile size.
   * @param bounds      the region for which to compute the maximum tile id
   * @param tile_size   the size of a tile within the region
   * @return the highest tile number within the region
   */
  static uint32_t MaxTileId(const AABB2<coord_t>& bounds, const double tile_size) {
    uint32_t cols = static_cast<uint32_t>(std::ceil(bounds.Width() / tile_size));
    uint32_t rows = static_cast<uint32_t>(std::ceil(bounds.Height() / tile_size));
    return (cols * rows) - 1;
  }

  /**
   * Get the base x,y of a specified tile.
   * @param   tileid   Tile Id.
   * @return  The base x,y of the specified tile.
   */
  coord_t Base(const int32_t tileid) const {
    int32_t row = tileid / ncolumns_;
    int32_t col = tileid - (row * ncolumns_);
    return coord_t(tilebounds_.minx() + (col * tilesize_), tilebounds_.miny() + (row * tilesize_));
  }

  /**
   * Get the bounding box of the specified tile.
   * @param   tileid   Tile Id.
   * @return  The latitude, longitude extent of the specified tile.
   */
  AABB2<coord_t> TileBounds(const int32_t tileid) const {
    Point2 base = Base(tileid);
    return AABB2<coord_t>(base.x(), base.y(), base.x() + tilesize_, base.y() + tilesize_);
  }

  /**
   * Get the bounding box of the tile with specified row, column.
   * @param   col   Tile column.
   * @param   row   Tile row.
   * @return  The latitude, longitude extent of the specified tile.
   */
  AABB2<coord_t> TileBounds(const int32_t col, const int32_t row) const {
    double basex = tilebounds_.minx() + ((double)col * tilesize_);
    double basey = tilebounds_.miny() + ((double)row * tilesize_);
    return AABB2<coord_t>(basex, basey, basex + tilesize_, basey + tilesize_);
  }

  /**
   * Get the center of the specified tile.
   * @param   tileid   Tile Id.
   * @return  The center x,y of the specified tile.
   */
  coord_t Center(const int32_t tileid) const {
    auto base = Base(tileid);
    return coord_t(base.x() + tilesize_ * 0.5, base.y() + tilesize_ * 0.5);
  }

  /**
   * Get the tile offsets (row,column) between the previous tile Id and
   * a new tileid.  The offsets are returned through arguments (references).
   * Offsets can be positive or negative or 0.
   * @param   initial_tileid     Original tile.
   * @param   newtileid      Tile to which relative offset is desired.
   * @param   delta_rows    Return: Relative number of rows.
   * @param   delta_cols    Return: Relative number of columns.
   */
  void TileOffsets(const int32_t initial_tileid,
                   const int32_t newtileid,
                   int32_t& delta_rows,
                   int32_t& delta_cols) const {
    int32_t deltaTile = newtileid - initial_tileid;
    delta_rows = (newtileid / ncolumns_) - (initial_tileid / ncolumns_);
    delta_cols = deltaTile - (delta_rows * ncolumns_);
  }

  /**
   * Get the number of tiles in the tiling system.
   * @return  Number of tiles.
   */
  uint32_t TileCount() const {
    double nrows = (tilebounds_.maxy() - tilebounds_.miny()) / tilesize_;
    return ncolumns_ * static_cast<int32_t>(std::ceil(nrows));
  }

  /**
   * Get the neighboring tileid to the right/east.
   * @param  tileid   Tile Id.
   * @return  Returns the tile Id of the tile to the right/east.
   */
  int32_t RightNeighbor(const int32_t tileid) const {
    return (tileid - ((tileid / ncolumns_) * ncolumns_) < ncolumns_ - 1)
               ? tileid + 1
               : wrapx_ ? tileid - ncolumns_ + 1 : tileid;
  }

  /**
   * Get the neighboring tileid to the left/west.
   * @param  tileid   Tile Id.
   * @return  Returns the tile Id of the tile to the left/west.
   */
  int32_t LeftNeighbor(const int32_t tileid) const {
    return tileid - ((tileid / ncolumns_) * ncolumns_) > 0 ? tileid - 1
                                                           : wrapx_ ? tileid + ncolumns_ - 1 : tileid;
  }

  /**
   * Get the neighboring tileid above or north.
   * @param  tileid   Tile Id.
   * @return  Returns the tile Id of the tile to the north. Return tileid
   *          if tile Id is on the top row (no neighbor to the north).
   */
  int32_t TopNeighbor(const int32_t tileid) const {
    return (tileid < static_cast<int32_t>((TileCount() - ncolumns_))) ? tileid + ncolumns_ : tileid;
  }

  /**
   * Get the neighboring tileid below or south.
   * @param  tileid   Tile Id.
   * @return  Returns the tile Id of the tile to the south. Return tileid
   *          if tile Id is on the bottom row (no neighbor to the south).
   */
  int32_t BottomNeighbor(const int32_t tileid) const {
    return (tileid < ncolumns_) ? tileid : tileid - ncolumns_;
  }

  /**
   * Checks if 2 tiles are neighbors (N,E,S,W). Does not support wrap around 180 longitude.
   * @param  id1  Tile Id 1
   * @param  id2  Tile Id 2
   * @return  Returns true if tile id1 and id2 are neighbors, false if not.
   */
  bool AreNeighbors(const uint32_t id1, const uint32_t id2) const {
    return id2 == id1 - 1 || id2 == id1 + 1 || id2 == id1 + ncolumns_ || id2 == id1 - ncolumns_;
  };

  /**
   * Get the list of tiles that lie within the specified bounding box. Since tiles as well as the
   * bounding box are both aligned to the axes we can simply find tiles by iterating over rows
   * and columns of tiles from the minimum to maximum.
   * @param  boundingbox  Bounding box
   * @return Returns a list of tiles that are within or intersect the bounding box.
   */
  std::vector<int32_t> TileList(const AABB2<coord_t>& bbox) const {
    // Check if x range needs to be split
    std::vector<AABB2<coord_t>> bboxes;
    if (wrapx_) {
      if (bbox.minx() < tilebounds_.minx() && bbox.maxx() > tilebounds_.minx()) {
        // Create 2 bounding boxes
        bboxes.emplace_back(tilebounds_.minx(), bbox.miny(), bbox.maxx(), bbox.maxy());
        bboxes.emplace_back(bbox.minx() + tilebounds_.Width(), bbox.miny(), tilebounds_.maxx(),
                            bbox.maxy());
      } else if (bbox.minx() < tilebounds_.maxx() && bbox.maxx() > tilebounds_.maxx()) {
        // Create 2 bounding boxes
        bboxes.emplace_back(bbox.minx(), bbox.miny(), tilebounds_.maxx(), bbox.maxy());
        bboxes.emplace_back(tilebounds_.minx(), bbox.miny(), bbox.maxx() - tilebounds_.Width(),
                            bbox.maxy());
      } else {
        bboxes.push_back(bbox.Intersection(tilebounds_));
      }
    } else {
      bboxes.push_back(bbox.Intersection(tilebounds_));
    }

    std::vector<int32_t> tilelist;
    for (const auto& bb : bboxes) {
      int32_t minrow = std::max(Row(bb.miny()), 0);
      int32_t maxrow = std::max(Row(bb.maxy()), 0);
      int32_t mincol = std::max(Col(bb.minx()), 0);
      int32_t maxcol = std::max(Col(bb.maxx()), 0);
      for (int32_t row = minrow; row <= maxrow; ++row) {
        int32_t tileid = TileId(mincol, row);
        for (int32_t col = mincol; col <= maxcol; ++col, ++tileid) {
          tilelist.push_back(tileid);
        }
      }
    }
    return tilelist;
  }

  /**
   * Get the list of tiles that lie within the specified ellipse. The method finds the tile
   * at the ellipse center. It successively finds neighbors and checks if they are inside or
   * intersect with the ellipse.
   * @param  ellipse  Ellipse
   * @return Returns a list of tiles that are within or intersect the ellipse.
   */
  std::vector<int32_t> TileList(const Ellipse<coord_t>& ellipse) const;

  /**
   * Color a "connectivity map" starting with a sparse map of uncolored tiles.
   * Any 2 tiles that have a connected path between them will have the same
   * value in the connectivity map.
   * @param  connectivity_map  map of tileid to color value
   */
  void ColorMap(std::unordered_map<uint32_t, size_t>& connectivity_map,
                const std::unordered_map<uint32_t, uint32_t>& not_neighbors = {}) const {
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

  /**
   * Intersect the linestring with the tiles to see which tiles and sub cells it intersects
   * @param line_string  the linestring to be tested against the cells
   * @return             the map of each tile intersected to a list of its intersected sub cell
   * indices
   */
  template <class container_t>
  std::unordered_map<int32_t, std::unordered_set<unsigned short>>
  Intersect(const container_t& line_string) const;

  /**
   * Intersect the bounding box with the tiles to see which tiles and sub-cells
   * (a.k.a bins) it intersects with. This can be used to reduce the number of
   * tiles and bins to search for matching items.
   *
   * @param box the bounding box to be tested.
   * @return    a map of tile IDs to a set of bin IDs with that tile.
   */
  std::unordered_map<int32_t, std::unordered_set<uint16_t>>
  Intersect(const AABB2<coord_t>& box) const;

  /**
   * Returns a functor which returns subdivisions close to the original point on each invocation in
   * a best first fashion If the functor can't expand any further (no more subdivisions) it will
   * throw
   * @param seed   the point at for which we measure 'closeness'
   * @return       the functor to be called
   */
  std::function<std::tuple<int32_t, unsigned short, double>()>
  ClosestFirst(const coord_t& seed) const {
    return std::bind(&closest_first_generator_t<coord_t>::next,
                     closest_first_generator_t<coord_t>(*this, seed));
  }

protected:
  // Does the tile bounds wrap in the x direction (e.g. at longitude = 180)
  bool wrapx_;

  // Bounding box of the tiling system.
  AABB2<coord_t> tilebounds_;

  // Tile size.  Tiles are square (equal y and x size).
  double tilesize_;

  // Number of rows ( y or latitude)
  int32_t nrows_;

  // Number of longitude (x or longitude).
  int32_t ncolumns_;

  // Number of subdivisions within a single tile
  unsigned short nsubdivisions_;

  double subdivision_size_;
};

} // namespace midgard
} // namespace valhalla

#endif // VALHALLA_MIDGARD_TILES_H_
