#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "tesselator.h"

namespace {

constexpr int kCaseCount = 200;
constexpr double kRelativeAreaTolerance = 1e-4;
constexpr uint64_t kSeed = 0x3872c0deULL;
constexpr double kPi = 3.141592653589793238462643383279502884;

struct Point {
  float x;
  float y;
};

struct Contour {
  std::vector<Point> points;
  bool reverse = false;
};

struct PolygonCase {
  int id = 0;
  std::string family;
  std::vector<Contour> contours;
};

struct DPoint {
  double x;
  double y;
};

double Cross(const DPoint& a, const DPoint& b) {
  return a.x * b.y - a.y * b.x;
}

DPoint Sub(const DPoint& a, const DPoint& b) {
  return {a.x - b.x, a.y - b.y};
}

DPoint ToDouble(const Point& p) {
  return {static_cast<double>(p.x), static_cast<double>(p.y)};
}

double SignedArea(const std::vector<Point>& points) {
  double area = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    area += Cross(ToDouble(points[i]),
                  ToDouble(points[(i + 1) % points.size()]));
  }
  return 0.5 * area;
}

std::string FormatDouble(double value, int precision = 10) {
  std::ostringstream out;
  out << std::setprecision(precision) << value;
  return out.str();
}

std::string FormatFloat(float value) {
  std::ostringstream out;
  out << std::setprecision(9) << value;
  return out.str();
}

bool IsFiniteCase(const PolygonCase& polygon_case) {
  for (const Contour& contour : polygon_case.contours) {
    if (contour.points.size() < 3) {
      return false;
    }
    for (const Point& point : contour.points) {
      if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
        return false;
      }
    }
  }
  return !polygon_case.contours.empty();
}

// Computes the area covered with odd-even fill. The input is split into
// vertical slabs at every vertex and segment-intersection x coordinate. Inside
// a slab no two non-vertical edges cross, so their y order is fixed. Two-point
// Gauss quadrature is exact for the linear y(x) difference of each pair.
double OddEvenArea(const std::vector<Contour>& contours) {
  struct Segment {
    DPoint a;
    DPoint b;
  };

  std::vector<Segment> segments;
  std::vector<double> x_events;
  for (const Contour& contour : contours) {
    for (size_t i = 0; i < contour.points.size(); ++i) {
      const DPoint a = ToDouble(contour.points[i]);
      const DPoint b =
          ToDouble(contour.points[(i + 1) % contour.points.size()]);
      segments.push_back({a, b});
      x_events.push_back(a.x);
      x_events.push_back(b.x);
    }
  }

  constexpr double kIntersectionEpsilon = 1e-12;
  for (size_t i = 0; i < segments.size(); ++i) {
    const DPoint r = Sub(segments[i].b, segments[i].a);
    if (std::abs(r.x) <= std::numeric_limits<double>::min()) {
      continue;
    }
    for (size_t j = i + 1; j < segments.size(); ++j) {
      const DPoint s = Sub(segments[j].b, segments[j].a);
      if (std::abs(s.x) <= std::numeric_limits<double>::min()) {
        continue;
      }
      const double denominator = Cross(r, s);
      if (std::abs(denominator) <=
          1e-15 * std::max(1.0, std::abs(r.x * s.y) + std::abs(r.y * s.x))) {
        continue;
      }
      const DPoint qp = Sub(segments[j].a, segments[i].a);
      const double t = Cross(qp, s) / denominator;
      const double u = Cross(qp, r) / denominator;
      if (t >= -kIntersectionEpsilon && t <= 1.0 + kIntersectionEpsilon &&
          u >= -kIntersectionEpsilon && u <= 1.0 + kIntersectionEpsilon) {
        x_events.push_back(segments[i].a.x + t * r.x);
      }
    }
  }

  std::sort(x_events.begin(), x_events.end());
  x_events.erase(std::unique(x_events.begin(), x_events.end()),
                 x_events.end());
  if (x_events.size() < 2) {
    return 0.0;
  }

  constexpr double kGauss1 = 0.21132486540518713;
  constexpr double kGauss2 = 0.78867513459481287;
  double area = 0.0;
  for (size_t slab = 0; slab + 1 < x_events.size(); ++slab) {
    const double left = x_events[slab];
    const double right = x_events[slab + 1];
    const double width = right - left;
    if (!(width > 0.0)) {
      continue;
    }

    double pair_height_sum = 0.0;
    for (const double gauss : {kGauss1, kGauss2}) {
      const double x = left + width * gauss;
      std::vector<double> ys;
      for (const Segment& segment : segments) {
        const double dx = segment.b.x - segment.a.x;
        if (std::abs(dx) <= std::numeric_limits<double>::min()) {
          continue;
        }
        const double lo = std::min(segment.a.x, segment.b.x);
        const double hi = std::max(segment.a.x, segment.b.x);
        if (x < lo || x > hi) {
          continue;
        }
        ys.push_back(segment.a.y +
                     (segment.b.y - segment.a.y) * (x - segment.a.x) / dx);
      }
      std::sort(ys.begin(), ys.end());
      if ((ys.size() & 1U) != 0) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      for (size_t i = 0; i + 1 < ys.size(); i += 2) {
        pair_height_sum += ys[i + 1] - ys[i];
      }
    }
    area += 0.5 * width * pair_height_sum;
  }
  return area;
}

std::vector<DPoint> ConvexHull(std::vector<DPoint> points) {
  std::sort(points.begin(), points.end(), [](const DPoint& a, const DPoint& b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  });
  points.erase(std::unique(points.begin(), points.end(),
                           [](const DPoint& a, const DPoint& b) {
                             return a.x == b.x && a.y == b.y;
                           }),
               points.end());
  if (points.size() <= 1) {
    return points;
  }

  std::vector<DPoint> lower;
  for (const DPoint& point : points) {
    while (lower.size() >= 2 &&
           Cross(Sub(lower.back(), lower[lower.size() - 2]),
                 Sub(point, lower.back())) <= 0.0) {
      lower.pop_back();
    }
    lower.push_back(point);
  }

  std::vector<DPoint> upper;
  for (auto it = points.rbegin(); it != points.rend(); ++it) {
    while (upper.size() >= 2 &&
           Cross(Sub(upper.back(), upper[upper.size() - 2]),
                 Sub(*it, upper.back())) <= 0.0) {
      upper.pop_back();
    }
    upper.push_back(*it);
  }

  lower.pop_back();
  upper.pop_back();
  lower.insert(lower.end(), upper.begin(), upper.end());
  return lower;
}

bool ProperlyIntersect(const DPoint& a, const DPoint& b, const DPoint& c,
                       const DPoint& d) {
  const DPoint r = Sub(b, a);
  const DPoint s = Sub(d, c);
  const double denominator = Cross(r, s);
  if (std::abs(denominator) < 1e-14) {
    return false;
  }
  const DPoint qp = Sub(c, a);
  const double t = Cross(qp, s) / denominator;
  const double u = Cross(qp, r) / denominator;
  return t > 1e-9 && t < 1.0 - 1e-9 && u > 1e-9 && u < 1.0 - 1e-9;
}

int CountProperCrossings(const std::vector<Point>& points) {
  int crossings = 0;
  const int count = static_cast<int>(points.size());
  for (int i = 0; i < count; ++i) {
    const DPoint a = ToDouble(points[i]);
    const DPoint b = ToDouble(points[(i + 1) % count]);
    for (int j = i + 1; j < count; ++j) {
      if (j == i || j == i + 1 || (i == 0 && j == count - 1)) {
        continue;
      }
      const DPoint c = ToDouble(points[j]);
      const DPoint d = ToDouble(points[(j + 1) % count]);
      if (ProperlyIntersect(a, b, c, d)) {
        ++crossings;
      }
    }
  }
  return crossings;
}

double MinimumInteriorAngle(const std::vector<Point>& points) {
  if (points.size() < 3) {
    return std::numeric_limits<double>::infinity();
  }
  double minimum = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < points.size(); ++i) {
    const DPoint previous = ToDouble(points[(i + points.size() - 1) %
                                            points.size()]);
    const DPoint current = ToDouble(points[i]);
    const DPoint next = ToDouble(points[(i + 1) % points.size()]);
    const DPoint incoming = Sub(previous, current);
    const DPoint outgoing = Sub(next, current);
    const double denominator =
        std::hypot(incoming.x, incoming.y) * std::hypot(outgoing.x, outgoing.y);
    if (denominator == 0.0) {
      return 0.0;
    }
    const double cosine =
        std::clamp((incoming.x * outgoing.x + incoming.y * outgoing.y) /
                       denominator,
                   -1.0, 1.0);
    minimum = std::min(minimum, std::acos(cosine));
  }
  return minimum;
}

bool HasCollinearTriple(const std::vector<Point>& points) {
  if (points.size() < 3) {
    return false;
  }
  for (size_t i = 0; i < points.size(); ++i) {
    const DPoint a = ToDouble(points[i]);
    const DPoint b = ToDouble(points[(i + 1) % points.size()]);
    const DPoint c = ToDouble(points[(i + 2) % points.size()]);
    if (Cross(Sub(b, a), Sub(c, b)) == 0.0) {
      return true;
    }
  }
  return false;
}

class Generator {
 public:
  explicit Generator(uint64_t seed) : rng_(seed) {}

  double Uniform(double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(rng_);
  }

  int UniformInt(int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng_);
  }

  Point RandomPolarPoint(const DPoint& center, double angle, double radius) {
    return {static_cast<float>(center.x + radius * std::cos(angle)),
            static_cast<float>(center.y + radius * std::sin(angle))};
  }

  PolygonCase HolesCase(int id) {
    const DPoint center{Uniform(-50.0, 50.0), Uniform(-50.0, 50.0)};
    const double base_radius = Uniform(55.0, 110.0);
    const int outer_count = UniformInt(9, 17);
    const double rotation = Uniform(0.0, 2.0 * kPi);
    std::vector<Point> outer;
    for (int i = 0; i < outer_count; ++i) {
      const double angle =
          rotation + 2.0 * kPi * i / outer_count + Uniform(-0.15, 0.15);
      const double radius = base_radius * Uniform(0.78, 1.0);
      outer.push_back(RandomPolarPoint(center, angle, radius));
    }

    std::vector<Contour> contours{{outer, false}};
    const int hole_count = UniformInt(1, 3);
    for (int hole = 0; hole < hole_count; ++hole) {
      double hole_radius = 0.0;
      DPoint hole_center{};
      bool accepted = false;
      for (int attempt = 0; attempt < 100 && !accepted; ++attempt) {
        hole_radius = Uniform(4.0, 13.0);
        const double distance = Uniform(0.0, 0.34 * base_radius);
        const double angle = Uniform(0.0, 2.0 * kPi);
        hole_center = {center.x + distance * std::cos(angle),
                       center.y + distance * std::sin(angle)};
        accepted = true;
        for (size_t previous = 1; previous < contours.size(); ++previous) {
          double previous_cx = 0.0;
          double previous_cy = 0.0;
          for (const Point& point : contours[previous].points) {
            previous_cx += point.x;
            previous_cy += point.y;
          }
          previous_cx /= contours[previous].points.size();
          previous_cy /= contours[previous].points.size();
          const double dx = hole_center.x - previous_cx;
          const double dy = hole_center.y - previous_cy;
          if (std::hypot(dx, dy) < hole_radius + 16.0) {
            accepted = false;
          }
        }
      }

      const int hole_vertices = UniformInt(5, 9);
      const double hole_rotation = Uniform(0.0, 2.0 * kPi);
      std::vector<Point> hole_points;
      for (int i = 0; i < hole_vertices; ++i) {
        const double angle =
            hole_rotation + 2.0 * kPi * i / hole_vertices +
            Uniform(-0.18, 0.18);
        hole_points.push_back(
            RandomPolarPoint(hole_center, angle, hole_radius * Uniform(0.72, 1.0)));
      }
      contours.push_back({hole_points, false});
    }
    return {id, "holes", contours};
  }

  PolygonCase SelfIntersectCase(int id) {
    for (int attempt = 0; attempt < 500; ++attempt) {
      const DPoint center{Uniform(-30.0, 30.0), Uniform(-30.0, 30.0)};
      const int count = UniformInt(8, 18);
      const double rotation = Uniform(0.0, 2.0 * kPi);
      std::vector<Point> points;
      for (int i = 0; i < count; ++i) {
        const double angle =
            rotation + 2.0 * kPi * i / count + Uniform(-0.12, 0.12);
        points.push_back(
            RandomPolarPoint(center, angle, Uniform(35.0, 120.0)));
      }
      std::shuffle(points.begin(), points.end(), rng_);
      if (CountProperCrossings(points) > 0) {
        return {id, "self-intersect", {{points, false}}};
      }
    }
    return SelfIntersectCase(id);
  }

  PolygonCase CollinearCase(int id) {
    for (int attempt = 0; attempt < 200; ++attempt) {
      std::vector<DPoint> cloud;
      for (int i = 0; i < 48; ++i) {
        cloud.push_back(
            {static_cast<double>(UniformInt(-240, 240)),
             static_cast<double>(UniformInt(-240, 240))});
      }
      const std::vector<DPoint> hull = ConvexHull(cloud);
      if (hull.size() < 5) {
        continue;
      }

      std::vector<Point> points;
      for (size_t i = 0; i < hull.size(); ++i) {
        const DPoint& a = hull[i];
        const DPoint& b = hull[(i + 1) % hull.size()];
        points.push_back({static_cast<float>(a.x), static_cast<float>(a.y)});
        if (Uniform(0.0, 1.0) < 0.45) {
          std::vector<double> fractions{0.25, 0.5, 0.75};
          std::shuffle(fractions.begin(), fractions.end(), rng_);
          const int inserted = UniformInt(1, 3);
          fractions.resize(inserted);
          std::sort(fractions.begin(), fractions.end());
          for (double fraction : fractions) {
            points.push_back(
                {static_cast<float>(a.x + (b.x - a.x) * fraction),
                 static_cast<float>(a.y + (b.y - a.y) * fraction)});
          }
        }
      }
      return {id, "collinear", {{points, false}}};
    }
    return CollinearCase(id);
  }

  PolygonCase TinyAngleCase(int id) {
    for (int attempt = 0; attempt < 200; ++attempt) {
      std::vector<DPoint> cloud;
      for (int i = 0; i < 32; ++i) {
        cloud.push_back({Uniform(-100.0, 100.0), Uniform(-100.0, 100.0)});
      }
      const std::vector<DPoint> hull = ConvexHull(cloud);
      if (hull.size() < 5) {
        continue;
      }

      // Insert one tall, very narrow tooth into a hull edge. The tooth lies in
      // the exterior half-plane of that edge, so it cannot cross the hull.
      const size_t spike_edge = static_cast<size_t>(
          UniformInt(0, static_cast<int>(hull.size()) - 1));
      std::vector<DPoint> spiked;
      for (size_t i = 0; i < hull.size(); ++i) {
        const DPoint& a = hull[i];
        const DPoint& b = hull[(i + 1) % hull.size()];
        spiked.push_back(a);
        if (i != spike_edge) {
          continue;
        }
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double length = std::hypot(dx, dy);
        if (length < 1e-9) {
          continue;
        }
        const double center_fraction = Uniform(0.3, 0.7);
        const double base_length = length * Uniform(1e-3, 4e-3);
        const double height = length * Uniform(5.0, 18.0);
        const DPoint base_center{a.x + dx * center_fraction,
                                 a.y + dy * center_fraction};
        const DPoint direction{dx / length, dy / length};
        // The hull is CCW, so (dy, -dx) points to its exterior.
        const DPoint outward{dy / length, -dx / length};
        const double half_base = 0.5 * base_length;
        spiked.push_back({base_center.x - direction.x * half_base,
                          base_center.y - direction.y * half_base});
        spiked.push_back({base_center.x + outward.x * height,
                          base_center.y + outward.y * height});
        spiked.push_back({base_center.x + direction.x * half_base,
                          base_center.y + direction.y * half_base});
      }

      const double rotation = Uniform(0.0, 2.0 * kPi);
      const double scale = Uniform(0.45, 1.2);
      const double cosine = std::cos(rotation) * scale;
      const double sine = std::sin(rotation) * scale;
      const DPoint translation{Uniform(-30.0, 30.0), Uniform(-30.0, 30.0)};
      std::vector<Point> points;
      for (const DPoint& point : spiked) {
        points.push_back(
            {static_cast<float>(translation.x + cosine * point.x -
                                sine * point.y),
             static_cast<float>(translation.y + sine * point.x +
                                cosine * point.y)});
      }
      if (CountProperCrossings(points) != 0 ||
          MinimumInteriorAngle(points) >= 1e-3) {
        continue;
      }
      return {id, "tiny-angle", {{points, false}}};
    }
    return TinyAngleCase(id);
  }

  PolygonCase Generate(int id) {
    switch (id % 4) {
      case 0:
        return HolesCase(id);
      case 1:
        return SelfIntersectCase(id);
      case 2:
        return CollinearCase(id);
      default:
        return TinyAngleCase(id);
    }
  }

 private:
  std::mt19937_64 rng_;
};

struct TrackingAllocator {
  std::map<void*, size_t> live;
};

void* TrackedAlloc(void* user_data, unsigned int size) {
  auto* allocator = static_cast<TrackingAllocator*>(user_data);
  void* pointer = std::malloc(size == 0 ? 1 : size);
  if (pointer != nullptr) {
    allocator->live[pointer] = size;
  }
  return pointer;
}

void* TrackedRealloc(void* user_data, void* pointer, unsigned int size) {
  auto* allocator = static_cast<TrackingAllocator*>(user_data);
  if (pointer == nullptr) {
    return TrackedAlloc(user_data, size);
  }
  const auto existing = allocator->live.find(pointer);
  void* replacement = std::realloc(pointer, size == 0 ? 1 : size);
  if (replacement != nullptr) {
    if (existing != allocator->live.end()) {
      allocator->live.erase(existing);
    }
    allocator->live[replacement] = size;
  }
  return replacement;
}

void TrackedFree(void* user_data, void* pointer) {
  auto* allocator = static_cast<TrackingAllocator*>(user_data);
  if (pointer == nullptr) {
    return;
  }
  allocator->live.erase(pointer);
  std::free(pointer);
}

struct Output {
  bool succeeded = false;
  int status = TESS_STATUS_INVALID_INPUT;
  int vertex_count = 0;
  int element_count = 0;
  int triangle_count = 0;
  double area = 0.0;
  bool invalid_indices = false;
  bool leaked = false;
  uint64_t hash = 1469598103934665603ULL;
  std::vector<float> vertices;
  std::vector<int> elements;
};

void MixHash(uint64_t* hash, const void* data, size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (size_t i = 0; i < size; ++i) {
    *hash ^= bytes[i];
    *hash *= 1099511628211ULL;
  }
}

void MixHashU64(uint64_t* hash, uint64_t value) {
  MixHash(hash, &value, sizeof(value));
}

Output RunTessellator(const PolygonCase& polygon_case,
                      TrackingAllocator& allocator) {
  Output output;
  const size_t live_before = allocator.live.size();
  TESSalloc tess_alloc{TrackedAlloc,
                       TrackedRealloc,
                       TrackedFree,
                       &allocator,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0};

  TESStesselator* tess = tessNewTess(&tess_alloc);
  if (tess == nullptr) {
    output.leaked = allocator.live.size() != live_before;
    return output;
  }
  for (const Contour& contour : polygon_case.contours) {
    tessSetOption(tess, TESS_REVERSE_CONTOURS, contour.reverse ? 1 : 0);
    tessAddContour(tess, 2, contour.points.data(), 2 * sizeof(TESSreal),
                   static_cast<int>(contour.points.size()));
  }

  const int rc = tessTesselate(tess, TESS_WINDING_ODD, TESS_POLYGONS, 3, 2,
                               nullptr);
  output.succeeded = rc != 0;
  output.status = static_cast<int>(tessGetStatus(tess));
  output.vertex_count = tessGetVertexCount(tess);
  output.element_count = tessGetElementCount(tess);
  output.triangle_count = output.succeeded ? output.element_count : 0;

  const TESSreal* vertices = tessGetVertices(tess);
  const TESSindex* elements = tessGetElements(tess);
  if (output.vertex_count > 0 && vertices != nullptr) {
    output.vertices.assign(vertices, vertices + 2 * output.vertex_count);
  } else if (output.vertex_count > 0) {
    output.invalid_indices = true;
  }
  if (output.element_count > 0 && elements != nullptr) {
    output.elements.assign(elements, elements + 3 * output.element_count);
  } else if (output.element_count > 0) {
    output.invalid_indices = true;
  }

  for (int triangle = 0; triangle < output.element_count; ++triangle) {
    const int i0 = output.elements[3 * triangle + 0];
    const int i1 = output.elements[3 * triangle + 1];
    const int i2 = output.elements[3 * triangle + 2];
    if (i0 < 0 || i1 < 0 || i2 < 0 || i0 >= output.vertex_count ||
        i1 >= output.vertex_count || i2 >= output.vertex_count) {
      output.invalid_indices = true;
      continue;
    }
    const double x0 = output.vertices[2 * i0 + 0];
    const double y0 = output.vertices[2 * i0 + 1];
    const double x1 = output.vertices[2 * i1 + 0];
    const double y1 = output.vertices[2 * i1 + 1];
    const double x2 = output.vertices[2 * i2 + 0];
    const double y2 = output.vertices[2 * i2 + 1];
    output.area +=
        0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
  }

  MixHashU64(&output.hash, static_cast<uint64_t>(output.succeeded));
  MixHashU64(&output.hash, static_cast<uint64_t>(output.status));
  MixHashU64(&output.hash, static_cast<uint64_t>(output.vertex_count));
  MixHashU64(&output.hash, static_cast<uint64_t>(output.element_count));
  MixHash(&output.hash, output.vertices.data(),
          output.vertices.size() * sizeof(float));
  MixHash(&output.hash, output.elements.data(),
          output.elements.size() * sizeof(int));

  tessDeleteTess(tess);
  output.leaked = allocator.live.size() != live_before;
  return output;
}

struct Evaluation {
  bool passed = false;
  std::string kind;
  std::string detail;
  double expected_area = 0.0;
  double actual_area = 0.0;
  double relative_error = 0.0;
  int triangle_count = 0;
  uint64_t first_hash = 0;
  uint64_t second_hash = 0;
  bool leaked = false;
};

Evaluation Evaluate(const PolygonCase& polygon_case) {
  Evaluation result;
  result.expected_area = OddEvenArea(polygon_case.contours);
  if (!std::isfinite(result.expected_area) || result.expected_area <= 0.0) {
    result.kind = "invalid-oracle";
    result.detail = "odd-even area is not finite and positive";
    return result;
  }

  TrackingAllocator allocator;
  const Output first = RunTessellator(polygon_case, allocator);
  const Output second = RunTessellator(polygon_case, allocator);
  result.actual_area = first.area;
  result.triangle_count = first.triangle_count;
  result.relative_error =
      std::abs(first.area - result.expected_area) / result.expected_area;
  result.first_hash = first.hash;
  result.second_hash = second.hash;
  result.leaked = first.leaked || second.leaked;

  std::ostringstream detail;
  detail << "expected=" << FormatDouble(result.expected_area)
         << " actual=" << FormatDouble(first.area)
         << " relerr=" << FormatDouble(result.relative_error)
         << " triangles=" << first.triangle_count;
  if (!first.succeeded) {
    detail << " rc=0 status=" << first.status;
  }
  result.detail = detail.str();

  if (result.leaked) {
    result.kind = "leak";
    result.detail += " allocator leak detected";
  } else if (first.hash != second.hash) {
    result.kind = "nondeterministic";
    result.detail += " repeated tessellation differs";
  } else if (!first.succeeded || first.status != TESS_STATUS_OK) {
    result.kind = "tessellation-failed";
  } else if (first.invalid_indices) {
    result.kind = "invalid-output";
  } else if (result.relative_error > kRelativeAreaTolerance) {
    result.kind = "area-mismatch";
  } else {
    result.passed = true;
  }
  return result;
}

bool SameFailure(const PolygonCase& polygon_case,
                 const std::string& original_kind) {
  const Evaluation evaluation = Evaluate(polygon_case);
  if (evaluation.passed || evaluation.kind != original_kind) {
    return false;
  }
  if (evaluation.kind == "area-mismatch") {
    return evaluation.expected_area > 1e-8;
  }
  return true;
}

PolygonCase MinimizeCase(const PolygonCase& original,
                         const std::string& failure_kind) {
  PolygonCase current = original;
  bool changed = true;
  while (changed) {
    changed = false;

    if (current.contours.size() > 1) {
      for (size_t contour = 0; contour < current.contours.size(); ++contour) {
        PolygonCase candidate = current;
        candidate.contours.erase(candidate.contours.begin() + contour);
        if (SameFailure(candidate, failure_kind)) {
          current = std::move(candidate);
          changed = true;
          break;
        }
      }
      if (changed) {
        continue;
      }
    }

    for (size_t contour = 0; contour < current.contours.size() && !changed;
         ++contour) {
      for (size_t point = 0;
           point < current.contours[contour].points.size() && !changed;
           ++point) {
        if (current.contours[contour].points.size() <= 3) {
          break;
        }
        PolygonCase candidate = current;
        candidate.contours[contour].points.erase(
            candidate.contours[contour].points.begin() + point);
        if (SameFailure(candidate, failure_kind)) {
          current = std::move(candidate);
          changed = true;
        }
      }
    }
  }
  return current;
}

struct Counterexample {
  int source_case = 0;
  std::string family;
  std::string kind;
  int original_point_count = 0;
  PolygonCase minimal;
  Evaluation evaluation;
};

std::string JsonEscape(const std::string& input) {
  std::string output;
  for (char character : input) {
    switch (character) {
      case '\\':
        output += "\\\\";
        break;
      case '"':
        output += "\\\"";
        break;
      case '\n':
        output += "\\n";
        break;
      case '\r':
        output += "\\r";
        break;
      case '\t':
        output += "\\t";
        break;
      default:
        output += character;
    }
  }
  return output;
}

void WriteCounterexampleText(const std::string& path,
                             const std::vector<Counterexample>& examples) {
  std::ofstream out(path);
  if (examples.empty()) {
    out << "No counterexamples: all 200 cases passed.\n";
    return;
  }
  for (size_t index = 0; index < examples.size(); ++index) {
    const Counterexample& example = examples[index];
    int minimal_points = 0;
    for (const Contour& contour : example.minimal.contours) {
      minimal_points += static_cast<int>(contour.points.size());
    }
    out << "failure=" << example.kind << " source_case=" << example.source_case
        << " family=" << example.family << " detail="
        << example.evaluation.detail << "\n";
    out << "minimal: contours=" << example.minimal.contours.size()
        << " points=" << minimal_points
        << " (original points=" << example.original_point_count << ")\n";
    for (size_t contour = 0; contour < example.minimal.contours.size();
         ++contour) {
      const Contour& points = example.minimal.contours[contour];
      out << "  contour[" << contour << "] n=" << points.points.size()
          << " reverse=" << (points.reverse ? 1 : 0)
          << " signed_area=" << FormatDouble(SignedArea(points.points))
          << " pts=\n";
      for (const Point& point : points.points) {
        out << "    {" << FormatFloat(point.x) << "f, " << FormatFloat(point.y)
            << "f},\n";
      }
    }
    out << "  serialized: " << example.minimal.contours.size();
    for (const Contour& contour : example.minimal.contours) {
      out << " | " << (contour.reverse ? 1 : 0) << " "
          << contour.points.size();
      for (const Point& point : contour.points) {
        out << " " << FormatFloat(point.x) << " " << FormatFloat(point.y);
      }
    }
    out << "\n\n";
  }
}

void WriteReproducer(const std::string& path,
                     const std::vector<Counterexample>& examples) {
  std::ofstream out(path);
  out << "// Minimal reproducers generated by Tests/property_test.cc\n"
      << "// Build from the repository root:\n"
      << "//   gcc -std=c11 -IInclude -c Source/*.c\n"
      << "//   g++ -std=c++17 -IInclude this.cc *.o -lm -o repro\n"
      << "#include <cmath>\n"
      << "#include <cstdio>\n"
      << "#include \"tesselator.h\"\n\n"
      << "struct InputContour { const TESSreal* points; int count; int "
         "reverse; };\n\n"
      << "[[maybe_unused]] static void RunCase(const char* name, const "
         "InputContour* contours, "
         "int contour_count, double expected) {\n"
      << "  TESStesselator* tess = tessNewTess(nullptr);\n"
      << "  for (int i = 0; i < contour_count; ++i) {\n"
      << "    tessSetOption(tess, TESS_REVERSE_CONTOURS, "
         "contours[i].reverse);\n"
      << "    tessAddContour(tess, 2, contours[i].points, "
         "2 * sizeof(TESSreal), contours[i].count);\n"
      << "  }\n"
      << "  int rc = tessTesselate(tess, TESS_WINDING_ODD, TESS_POLYGONS, 3, "
         "2, nullptr);\n"
      << "  const TESSreal* vertices = tessGetVertices(tess);\n"
      << "  const TESSindex* elements = tessGetElements(tess);\n"
      << "  int triangles = tessGetElementCount(tess);\n"
      << "  double area = 0.0;\n"
      << "  for (int i = 0; rc && i < triangles; ++i) {\n"
      << "    int a = elements[3*i], b = elements[3*i+1], c = "
         "elements[3*i+2];\n"
      << "    double ax = vertices[2*a], ay = vertices[2*a+1];\n"
      << "    double bx = vertices[2*b], by = vertices[2*b+1];\n"
      << "    double cx = vertices[2*c], cy = vertices[2*c+1];\n"
      << "    area += 0.5 * std::abs((bx-ax)*(cy-ay) - "
         "(cx-ax)*(by-ay));\n"
      << "  }\n"
      << "  double relerr = expected > 0.0 ? std::abs(area - expected) / "
         "expected : 0.0;\n"
      << "  std::printf(\"%s rc=%d status=%d triangles=%d area=%.10g "
         "expected=%.10g relerr=%.6g\\n\", name, rc, (int)tessGetStatus(tess), "
         "triangles, area, expected, relerr);\n"
      << "  tessDeleteTess(tess);\n"
      << "}\n\n"
      << "int main() {\n";

  for (size_t index = 0; index < examples.size(); ++index) {
    const Counterexample& example = examples[index];
    out << "  {\n";
    for (size_t contour = 0; contour < example.minimal.contours.size();
         ++contour) {
      const Contour& points = example.minimal.contours[contour];
      out << "    static const TESSreal c" << index << "_" << contour
          << "[] = {\n";
      for (const Point& point : points.points) {
        out << "      " << FormatFloat(point.x) << "f, " << FormatFloat(point.y)
            << "f,\n";
      }
      out << "    };\n";
    }
    out << "    static const InputContour contours[] = {\n";
    for (size_t contour = 0; contour < example.minimal.contours.size();
         ++contour) {
      out << "      {c" << index << "_" << contour << ", "
          << example.minimal.contours[contour].points.size() << ", "
          << (example.minimal.contours[contour].reverse ? 1 : 0) << "},\n";
    }
    out << "    };\n";
    out << "    RunCase(\"case_" << example.source_case << "_"
        << example.kind << "\", contours, "
        << example.minimal.contours.size() << ", "
        << FormatDouble(example.evaluation.expected_area, 17) << ");\n";
    out << "  }\n";
  }
  out << "  return 0;\n}\n";
}

struct FamilyStats {
  int cases = 0;
  int failures = 0;
  long long triangles = 0;
  double max_relative_error = 0.0;
};

void RunOracleSelfCheck() {
  const std::vector<Contour> square{
      {{{0.0f, 0.0f}, {4.0f, 0.0f}, {4.0f, 4.0f}, {0.0f, 4.0f}}, false}};
  const std::vector<Contour> square_with_hole{
      {{{0.0f, 0.0f}, {4.0f, 0.0f}, {4.0f, 4.0f}, {0.0f, 4.0f}}, false},
      {{{1.0f, 1.0f}, {2.0f, 1.0f}, {2.0f, 2.0f}, {1.0f, 2.0f}}, false}};
  const std::vector<Contour> bow_tie{
      {{{0.0f, 0.0f}, {2.0f, 2.0f}, {2.0f, 0.0f}, {0.0f, 2.0f}}, false}};
  const std::vector<Contour> two_squares{
      {{{0.0f, 0.0f}, {2.0f, 0.0f}, {2.0f, 2.0f}, {0.0f, 2.0f}}, false},
      {{{3.0f, 0.0f}, {5.0f, 0.0f}, {5.0f, 2.0f}, {3.0f, 2.0f}}, false}};
  const std::vector<Contour> identical_twice{
      {{{0.0f, 0.0f}, {2.0f, 0.0f}, {2.0f, 2.0f}, {0.0f, 2.0f}}, false},
      {{{0.0f, 0.0f}, {2.0f, 0.0f}, {2.0f, 2.0f}, {0.0f, 2.0f}}, false}};

  const std::array<std::pair<const std::vector<Contour>*, double>, 5> checks{{
      {&square, 16.0},
      {&square_with_hole, 15.0},
      {&bow_tie, 2.0},
      {&two_squares, 8.0},
      {&identical_twice, 0.0},
  }};
  for (const auto& [contours, expected] : checks) {
    const double actual = OddEvenArea(*contours);
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-9) {
      std::cerr << "OddEvenArea self-check failed: expected " << expected
                << ", got " << actual << "\n";
      std::exit(2);
    }
  }
}

void PrintMinimalCase(const Counterexample& example) {
  std::cout << "FAIL " << example.kind << " case=" << example.source_case
            << " family=" << example.family << " " << example.evaluation.detail
            << "\nminimal contours=" << example.minimal.contours.size();
  for (size_t contour = 0; contour < example.minimal.contours.size();
       ++contour) {
    std::cout << " [" << contour << ":";
    for (const Point& point : example.minimal.contours[contour].points) {
      std::cout << " (" << FormatFloat(point.x) << ","
                << FormatFloat(point.y) << ")";
    }
    std::cout << "]";
  }
  std::cout << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::string output_directory = "artifacts";
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    constexpr const char* kPrefix = "--out-dir=";
    if (argument.rfind(kPrefix, 0) == 0) {
      output_directory = argument.substr(std::strlen(kPrefix));
    } else {
      std::cerr << "Unknown argument: " << argument << "\n";
      return 2;
    }
  }
  std::filesystem::create_directories(output_directory);
  RunOracleSelfCheck();

  Generator generator(kSeed);
  std::array<std::string, 4> family_names{
      "holes", "self-intersect", "collinear", "tiny-angle"};
  std::map<std::string, FamilyStats> family_stats;
  for (const std::string& family : family_names) {
    family_stats[family] = {};
  }

  int failures = 0;
  int leaks = 0;
  int nondeterministic = 0;
  int crashes = 0;
  int hole_feature_cases = 0;
  int self_intersect_feature_cases = 0;
  int collinear_feature_cases = 0;
  int tiny_angle_feature_cases = 0;
  double minimum_tiny_angle = std::numeric_limits<double>::infinity();
  long long total_triangles = 0;
  double max_relative_error = 0.0;
  uint64_t first_run_hash = 1469598103934665603ULL;
  uint64_t second_run_hash = 1469598103934665603ULL;
  std::vector<Counterexample> counterexamples;

  for (int id = 0; id < kCaseCount; ++id) {
    const PolygonCase polygon_case = generator.Generate(id);
    if (!IsFiniteCase(polygon_case)) {
      std::cerr << "Generator produced invalid case " << id << "\n";
      return 2;
    }
    if (polygon_case.family == "holes" &&
        polygon_case.contours.size() > 1) {
      ++hole_feature_cases;
    } else if (polygon_case.family == "self-intersect" &&
               CountProperCrossings(polygon_case.contours.front().points) > 0) {
      ++self_intersect_feature_cases;
    } else if (polygon_case.family == "collinear" &&
               HasCollinearTriple(polygon_case.contours.front().points)) {
      ++collinear_feature_cases;
    } else if (polygon_case.family == "tiny-angle") {
      const double angle =
          MinimumInteriorAngle(polygon_case.contours.front().points);
      minimum_tiny_angle = std::min(minimum_tiny_angle, angle);
      if (angle < 1e-3) {
        ++tiny_angle_feature_cases;
      }
    }
    const Evaluation evaluation = Evaluate(polygon_case);
    FamilyStats& stats = family_stats[polygon_case.family];
    ++stats.cases;
    stats.triangles += evaluation.triangle_count;
    stats.max_relative_error =
        std::max(stats.max_relative_error, evaluation.relative_error);
    total_triangles += evaluation.triangle_count;
    max_relative_error =
        std::max(max_relative_error, evaluation.relative_error);
    MixHashU64(&first_run_hash, evaluation.first_hash);
    MixHashU64(&second_run_hash, evaluation.second_hash);
    leaks += evaluation.leaked ? 1 : 0;
    nondeterministic +=
        evaluation.kind == "nondeterministic" ? 1 : 0;

    if (evaluation.passed) {
      continue;
    }
    ++failures;
    ++stats.failures;
    int original_points = 0;
    for (const Contour& contour : polygon_case.contours) {
      original_points += static_cast<int>(contour.points.size());
    }
    const PolygonCase minimal =
        MinimizeCase(polygon_case, evaluation.kind);
    const Evaluation minimal_evaluation = Evaluate(minimal);
    Counterexample example{polygon_case.id,
                           polygon_case.family,
                           evaluation.kind,
                           original_points,
                           minimal,
                           minimal_evaluation};
    PrintMinimalCase(example);
    counterexamples.push_back(std::move(example));
  }

  const bool hashes_equal = first_run_hash == second_run_hash;
  const bool features_ok =
      hole_feature_cases == 50 && self_intersect_feature_cases == 50 &&
      collinear_feature_cases == 50 && tiny_angle_feature_cases == 50 &&
      minimum_tiny_angle < 1e-3;
  const std::string summary_path =
      output_directory + "/property_test_summary.json";
  const std::string counterexample_path =
      output_directory + "/property_counterexample.txt";
  const std::string reproducer_path =
      output_directory + "/property_counterexample_repro.cc";

  std::ofstream summary(summary_path);
  summary << "{\n"
          << "  \"seed\": " << kSeed << ",\n"
          << "  \"cases\": " << kCaseCount << ",\n"
          << "  \"family_cycle\": [\"holes\", \"self-intersect\", "
             "\"collinear\", \"tiny-angle\"],\n"
          << "  \"winding_rule\": \"TESS_WINDING_ODD\",\n"
          << "  \"tolerance_relative\": " << kRelativeAreaTolerance << ",\n"
          << "  \"crashes\": " << crashes << ",\n"
          << "  \"leaks\": " << leaks << ",\n"
          << "  \"determinism_mismatches\": " << nondeterministic << ",\n"
          << "  \"two_tessellations_equal\": "
          << (hashes_equal ? "true" : "false") << ",\n"
          << "  \"feature_checks\": {\n"
          << "    \"holes_with_inner_contours\": " << hole_feature_cases
          << ",\n"
          << "    \"self_intersecting\": " << self_intersect_feature_cases
          << ",\n"
          << "    \"with_collinear_triples\": " << collinear_feature_cases
          << ",\n"
          << "    \"with_angle_below_1e-3\": " << tiny_angle_feature_cases
          << ",\n"
          << "    \"minimum_tiny_angle_radians\": "
          << FormatDouble(minimum_tiny_angle, 17) << "\n"
          << "  },\n"
          << "  \"failures\": " << failures << ",\n"
          << "  \"passed\": " << (kCaseCount - failures) << ",\n"
          << "  \"total_triangles\": " << total_triangles << ",\n"
          << "  \"max_relative_area_error\": "
          << FormatDouble(max_relative_error, 17) << ",\n"
          << "  \"aggregate_output_hash\": \"" << std::hex
          << first_run_hash << std::dec << "\",\n"
          << "  \"aggregate_output_hash_second_run\": \"" << std::hex
          << second_run_hash << std::dec << "\",\n"
          << "  \"families\": [\n";
  for (size_t i = 0; i < family_names.size(); ++i) {
    const FamilyStats& stats = family_stats[family_names[i]];
    summary << "    {\"name\": \"" << JsonEscape(family_names[i])
            << "\", \"cases\": " << stats.cases
            << ", \"failures\": " << stats.failures
            << ", \"triangles\": " << stats.triangles
            << ", \"max_relerr\": "
            << FormatDouble(stats.max_relative_error, 17) << "}";
    summary << (i + 1 == family_names.size() ? "\n" : ",\n");
  }
  summary << "  ]\n}\n";

  WriteCounterexampleText(counterexample_path, counterexamples);
  WriteReproducer(reproducer_path, counterexamples);

  std::cout << "cases=" << kCaseCount << " failures=" << failures
            << " leaks=" << leaks
            << " determinism_mismatches=" << nondeterministic
            << " max_relerr=" << FormatDouble(max_relative_error)
            << " hash=" << std::hex << first_run_hash << std::dec << "\n";
  std::cout << "features: holes=" << hole_feature_cases
            << " self-intersect=" << self_intersect_feature_cases
            << " collinear=" << collinear_feature_cases
            << " tiny-angle=" << tiny_angle_feature_cases
            << " min_angle=" << FormatDouble(minimum_tiny_angle) << "\n";
  std::cout << "summary=" << summary_path << "\n"
            << "counterexamples=" << counterexample_path << "\n"
            << "reproducer=" << reproducer_path << "\n";
  return failures == 0 && hashes_equal && leaks == 0 && features_ok ? 0 : 1;
}
