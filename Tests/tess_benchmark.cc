#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "tesselator.h"

namespace {

struct Point {
  float x;
  float y;
};

struct Case {
  const char* name;
  int holes_per_side;
  int iterations;
};

struct Triangle {
  std::array<Point, 3> p;

  bool operator<(const Triangle& other) const {
    for (int i = 0; i < 3; ++i) {
      if (p[i].x != other.p[i].x) return p[i].x < other.p[i].x;
      if (p[i].y != other.p[i].y) return p[i].y < other.p[i].y;
    }
    return false;
  }

  bool operator==(const Triangle& other) const {
    for (int i = 0; i < 3; ++i) {
      if (p[i].x != other.p[i].x || p[i].y != other.p[i].y) return false;
    }
    return true;
  }
};

uint64_t Fnv1a(const void* data, size_t size, uint64_t hash = 1469598103934665603ULL) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (size_t i = 0; i < size; ++i) {
    hash ^= bytes[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

std::vector<Point> MakeOuter(int holes_per_side) {
  const float extent = static_cast<float>(holes_per_side);
  return {{-1.0f, -1.0f}, {extent + 1.0f, -1.0f},
          {extent + 1.0f, extent + 1.0f}, {-1.0f, extent + 1.0f}};
}

std::vector<Point> MakeHole(int x, int y) {
  const float x0 = static_cast<float>(x) + 0.25f;
  const float y0 = static_cast<float>(y) + 0.25f;
  const float x1 = static_cast<float>(x) + 0.75f;
  const float y1 = static_cast<float>(y) + 0.75f;
  return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
}

void AddContour(TESStesselator* tess, const std::vector<Point>& points,
                bool reverse) {
  tessSetOption(tess, TESS_REVERSE_CONTOURS, reverse ? 1 : 0);
  tessAddContour(tess, 2, points.data(), sizeof(Point),
                 static_cast<int>(points.size()));
}

TESStesselator* BuildTess(int holes_per_side) {
  TESStesselator* tess = tessNewTess(nullptr);
  if (tess == nullptr) return nullptr;

  AddContour(tess, MakeOuter(holes_per_side), false);
  for (int y = 0; y < holes_per_side; ++y) {
    for (int x = 0; x < holes_per_side; ++x) {
      AddContour(tess, MakeHole(x, y), true);
    }
  }
  return tess;
}

double TriangleArea(const Triangle& tri) {
  const double ax = tri.p[0].x;
  const double ay = tri.p[0].y;
  const double bx = tri.p[1].x;
  const double by = tri.p[1].y;
  const double cx = tri.p[2].x;
  const double cy = tri.p[2].y;
  return std::fabs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax)) * 0.5;
}

std::vector<Triangle> ExtractTriangles(TESStesselator* tess,
                                       int* vertex_count,
                                       int* element_count) {
  *vertex_count = tessGetVertexCount(tess);
  *element_count = tessGetElementCount(tess);
  const TESSreal* vertices = tessGetVertices(tess);
  const TESSindex* elements = tessGetElements(tess);

  std::vector<Triangle> triangles;
  triangles.reserve(*element_count);
  for (int i = 0; i < *element_count; ++i) {
    Triangle tri;
    for (int j = 0; j < 3; ++j) {
      const TESSindex index = elements[i * 3 + j];
      tri.p[j] = {vertices[index * 2], vertices[index * 2 + 1]};
    }
    std::sort(tri.p.begin(), tri.p.end(), [](Point a, Point b) {
      return a.x != b.x ? a.x < b.x : a.y < b.y;
    });
    triangles.push_back(tri);
  }
  std::sort(triangles.begin(), triangles.end());
  return triangles;
}

uint64_t TriangleHash(const std::vector<Triangle>& triangles) {
  uint64_t hash = 1469598103934665603ULL;
  for (const Triangle& tri : triangles) {
    for (const Point& p : tri.p) {
      hash = Fnv1a(&p.x, sizeof(p.x), hash);
      hash = Fnv1a(&p.y, sizeof(p.y), hash);
    }
  }
  return hash;
}

double TotalArea(const std::vector<Triangle>& triangles) {
  double area = 0.0;
  for (const Triangle& tri : triangles) area += TriangleArea(tri);
  return area;
}

bool RunCase(const Case& test_case) {
  TESStesselator* tess = BuildTess(test_case.holes_per_side);
  if (tess == nullptr) {
    std::fprintf(stderr, "failed to create tessellator\n");
    return false;
  }

  if (!tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2,
                     nullptr)) {
    std::fprintf(stderr, "%s: tessellation failed: %d\n", test_case.name,
                 static_cast<int>(tessGetStatus(tess)));
    tessDeleteTess(tess);
    return false;
  }
  int vertex_count = 0;
  int element_count = 0;
  const std::vector<Triangle> reference =
      ExtractTriangles(tess, &vertex_count, &element_count);
  const uint64_t reference_hash = TriangleHash(reference);
  const double reference_area = TotalArea(reference);
  tessDeleteTess(tess);

  std::vector<double> samples_ms;
  samples_ms.reserve(test_case.iterations);
  for (int iteration = 0; iteration < test_case.iterations; ++iteration) {
    tess = BuildTess(test_case.holes_per_side);
    if (tess == nullptr) return false;

    const auto start = std::chrono::steady_clock::now();
    const int ok = tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3,
                                 2, nullptr);
    const auto end = std::chrono::steady_clock::now();
    if (!ok) {
      std::fprintf(stderr, "%s: timed tessellation failed\n", test_case.name);
      tessDeleteTess(tess);
      return false;
    }
    int check_vertex_count = 0;
    int check_element_count = 0;
    const std::vector<Triangle> check =
        ExtractTriangles(tess, &check_vertex_count, &check_element_count);
    if (check != reference || TriangleHash(check) != reference_hash ||
        check_vertex_count != vertex_count ||
        check_element_count != element_count ||
        TotalArea(check) != reference_area) {
      std::fprintf(stderr, "%s: output changed between runs\n",
                   test_case.name);
      tessDeleteTess(tess);
      return false;
    }

    samples_ms.push_back(
        std::chrono::duration<double, std::milli>(end - start).count());
    tessDeleteTess(tess);
  }

  std::sort(samples_ms.begin(), samples_ms.end());
  const double min_ms = samples_ms.front();
  const double median_ms = samples_ms[samples_ms.size() / 2];
  const int input_vertices = test_case.holes_per_side *
                                 test_case.holes_per_side * 4 +
                             4;
  std::printf(
      "%s,%d,%d,%.6f,%.6f,%d,%d,%.6f,%016llx\n", test_case.name,
      input_vertices, test_case.iterations, min_ms, median_ms, element_count,
      vertex_count, reference_area,
      static_cast<unsigned long long>(reference_hash));
  return true;
}

}  // namespace

int main() {
  const Case cases[] = {
      {"small", 96, 8},
      {"medium", 112, 6},
      {"large", 128, 4},
  };
  int scale = 1;
  if (const char* value = std::getenv("BENCH_SCALE")) {
    scale = std::max(1, std::atoi(value));
  }

  std::printf(
      "case,input_vertices,iterations,min_ms,median_ms,triangles,output_"
      "vertices,area,triangle_hash\n");
  for (const Case& test_case : cases) {
    Case scaled = test_case;
    scaled.iterations *= scale;
    if (!RunCase(scaled)) return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
