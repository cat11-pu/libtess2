#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

#include "gtest/gtest.h"
#include "tesselator.h"

namespace {

constexpr int kComponentCount = 2;
constexpr int kNumTriangleVertices = 3;

struct Vector2f {
  float x;
  float y;
};

void* HeapAlloc(void* user_data, unsigned int size) { return malloc(size); }
void* HeapRealloc(void* user_data, void* ptr, unsigned int size) {
  return realloc(ptr, size);
}
void HeapFree(void* user_data, void* ptr) { free(ptr); }

// Adds a polygon with a hole to the tessellator.
//
// a = outer loop
// b = inner loop
// x = interior
// o = exterior
//
// +aaaaaaaaaaaaaa+
// a xx | xx | xx a
// a xx | xx | xx a
// a----+bbbb+----a
// a xx b oo b xx a
// a xx b oo b xx a
// a----+bbbb+----a
// a xx | xx | xx a
// a xx | xx | xx a
// +aaaaaaaaaaaaaa+
//
// This should tessellate to 8 triangles.
void AddPolygonWithHole(TESStesselator* tess) {
  std::vector<Vector2f> outer_loop = {
      {0, 0},
      {3, 0},
      {3, 3},
      {0, 3},
  };
  std::vector<Vector2f> inner_hole = {
      {1, 1},
      {2, 1},
      {2, 2},
      {1, 2},
  };

  tessSetOption(tess, TESS_REVERSE_CONTOURS, 0);
  tessAddContour(tess, kComponentCount, outer_loop.data(),
                 kComponentCount * sizeof(TESSreal), outer_loop.size());
  tessSetOption(tess, TESS_REVERSE_CONTOURS, 1);
  tessAddContour(tess, kComponentCount, inner_hole.data(),
                 kComponentCount * sizeof(TESSreal), inner_hole.size());
}

void AddPolyline(TESStesselator* tess, const std::vector<Vector2f>& polyline) {
  tessAddContour(tess, kComponentCount, polyline.data(),
                 kComponentCount * sizeof(TESSreal), polyline.size());
}

class Libtess2Test : public testing::Test {
 protected:
  Libtess2Test(TESSalloc* alloc) { tess = tessNewTess(alloc); }

  Libtess2Test() : Libtess2Test(nullptr) {}

  ~Libtess2Test() {
    if (tess != nullptr) {
      tessDeleteTess(tess);
    }
  }

  void SetUp() override { ASSERT_NE(tess, nullptr); }

  TESStesselator* tess;
};

static TESSalloc custom_alloc{
    HeapAlloc,
    HeapRealloc,
    HeapFree,
    /*user_data=*/nullptr,
    /*meshEdgeBucketSize=*/512,
    /*meshVertexBucketSize=*/512,
    /*meshFaceBucketSize=*/256,
    /*dictNodeBucketSize=*/512,
    /*regionBucketSize=*/256,
    /*extraVertices=*/0,
};

class Libtess2CustomAllocTest : public Libtess2Test {
 protected:
  Libtess2CustomAllocTest() : Libtess2Test(&custom_alloc) {}
};

// Tests that tessellation succeeds when the default allocator is used.
TEST_F(Libtess2Test, DefaultAllocSuccess) {
  // Add the polygon and tessellate it.
  AddPolygonWithHole(tess);
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);

  // It should take 8 triangles to represent the polygon with the hole.
  EXPECT_EQ(tessGetElementCount(tess), 8);
}

TEST_F(Libtess2CustomAllocTest, CustomAllocSuccess) {
  // Add the polygon and tessellate it.
  AddPolygonWithHole(tess);
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);

  // It should take 8 triangles to represent the polygon with the hole.
  EXPECT_EQ(tessGetElementCount(tess), 8);
}

TEST_F(Libtess2Test, EmptyPolyline) {
  AddPolyline(tess, {});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetElementCount(tess), 0);
}

TEST_F(Libtess2Test, SingleLine) {
  AddPolyline(tess, {{0, 0}, {0, 1}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetElementCount(tess), 0);
}

TEST_F(Libtess2Test, SingleTriangle) {
  AddPolyline(tess, {{0, 0}, {0, 1}, {1, 0}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetElementCount(tess), 1);
}

TEST_F(Libtess2Test, UnitQuad) {
  AddPolyline(tess, {{0, 0}, {0, 1}, {1, 1}, {1, 0}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetElementCount(tess), 2);
}

TEST_F(Libtess2Test, GetStatusInvalidInput) {
  AddPolyline(tess, {{-2e+37f, 0.f}, {0, 5}, {1e37f, -5}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
}

TEST_F(Libtess2Test, GetStatusOk) {
  AddPolyline(tess, {{0, 0}, {0, 1}, {1, 1}, {1, 0}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_OK);
}

TEST_F(Libtess2Test, FloatOverflowQuad) {
  constexpr float kFloatMin = std::numeric_limits<float>::min();
  constexpr float kFloatMax = std::numeric_limits<float>::max();

  // Add the polygon and tessellate it.
  AddPolyline(tess, {{kFloatMin, kFloatMin},
                     {kFloatMin, kFloatMax},
                     {kFloatMax, kFloatMax},
                     {kFloatMax, kFloatMin}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

TEST_F(Libtess2Test, SingularityQuad) {
  AddPolyline(tess, {{0, 0}, {0, 0}, {0, 0}, {0, 0}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetElementCount(tess), 0);
}

TEST_F(Libtess2Test, DegenerateQuad) {
  // A quad that's extremely close to a giant triangle, with an extra sliver.
  // Caused a segfault previously.
  AddPolyline(tess, {{0.f, 3.40282347e+38f},
                     {0.64113313f, -1.f},
                     {-0.f, -0.f},
                     {-3.40282347e+38f, 1.f}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

TEST_F(Libtess2Test, WidthOverflowsTri) {
  AddPolyline(tess, {{-2e+38f, 0}, {0, 0}, {2e+38f, -1}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

TEST_F(Libtess2Test, HeightOverflowsTri) {
  AddPolyline(tess, {{0, 0}, {0, 2e+38f}, {-1, -2e+38f}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

TEST_F(Libtess2Test, AreaOverflowsTri) {
  AddPolyline(tess, {{-2e+37f, 0.f}, {0, 5}, {1e37f, -5}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

TEST_F(Libtess2Test, NanQuad) {
  AddPolyline(tess, {{std::nanf(""), std::nanf("")},
                     {std::nanf(""), std::nanf("")},
                     {std::nanf(""), std::nanf("")},
                     {std::nanf(""), std::nanf("")}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetElementCount(tess), 0);
}

TEST_F(Libtess2Test, AvoidsCrashWhileFindingIntersection) {
  // Previously, this failed an assert while finding an intersection because
  // that fell back to taking a midpoint between two coordinates in a way that
  // could get the wrong answer because of the sum overflowing max float.
  AddPolyline(tess, {{-1.f, 0.f},
                     {0.868218958f, 0.f},
                     {0.902460039f, 0.0649746507f},
                     {-0.f, 0.854620099f},
                     {-1.f, 0.784999669f},
                     {0.f, 0.f},
                     {-1.f, 1.f},
                     {1.f, 1.f},
                     {0.f, -1.f},
                     {3.40282347e+38f, 3.40282347e+38f},
                     {-1.f, -1.f},
                     {-0.f, 0.442898333f},
                     {0.33078745f, -0.f},
                     {-0.f, 1.f},
                     {-1.f, 0.f},
                     {1.f, -0.f},
                     {0.f, 0.186138511f},
                     {0.212649569f, 0.886535764f},
                     {1.f, 0.34795785f},
                     {0.f, 0.788870096f},
                     {0.853441715f, -1.f},
                     {-1.f, 1.f},
                     {1.f, -0.994903505f},
                     {1.f, 0.105880626f},
                     {3.40282347e+38f, 3.40282347e+38f},
                     {-1.f, 3.40282347e+38f},
                     {-0.f, 0.34419331f},
                     {1.f, 1.f}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

TEST_F(Libtess2Test, AvoidsCrashInAddRightEdges) {
  AddPolyline(tess, {{{-0.5f, 1.f},
                      {3.40282347e+38f, 0.f},
                      {0.349171013f, 1.f},
                      {1.f, 0.f},
                      {1.f, -0.f},
                      {0.594775498f, -0.f},
                      {0.f, -0.f},
                      {-0.f, 1.f},
                      {0.f, 1.f},
                      {2.20929384f, 1.f},
                      {1.f, 1.f},
                      {-0.f, -0.f},
                      {3.40282347e+38f, -0.f},
                      {-1.f, 0.f},
                      {1.70141173e+38f, 0.391036272f},
                      {3.40282347e+38f, 0.371295959f},
                      {3.40282347e+38f, -0.f},
                      {0.f, 0.234747186f},
                      {-1.f, 1.f},
                      {-1.f, -0.f},
                      {3.40282347e+38f, 1.f},
                      {-0.f, -0.f},
                      {3.40282347e+38f, 1.f},
                      {0.434241712f, 0.f},
                      {1.f, 0.211511821f},
                      {3.40282347e+38f, 1.f}}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
}

// ---- Narrow per-triangle attribution interface ----

void AddContour2f(TESStesselator* tess, const std::vector<Vector2f>& loop,
                  bool reverse = false) {
  tessSetOption(tess, TESS_REVERSE_CONTOURS, reverse ? 1 : 0);
  tessAddContour(tess, kComponentCount, loop.data(),
                 kComponentCount * sizeof(TESSreal), loop.size());
}

int TriangulateForInfo(TESStesselator* tess) {
  return tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                       kNumTriangleVertices, kComponentCount, nullptr);
}

// Signed area of a triangle using the output vertex array.
double TriangleArea2d(const TESSreal* verts, const TESSindex* tri) {
  auto p = [&](int corner) {
    return reinterpret_cast<const Vector2f*>(
        &verts[tri[corner] * kComponentCount]);
  };
  const Vector2f* a = p(0);
  const Vector2f* b = p(1);
  const Vector2f* c = p(2);
  return 0.5 * ((b->x - a->x) * (c->y - a->y) -
                (c->x - a->x) * (b->y - a->y));
}

// Contour indices line up 1:1 with tessAddContour() calls, and a hole input
// is classified as a hole even though no triangle is produced inside it.
TEST_F(Libtess2Test, TriangleInfoContourOrderAndRoles) {
  std::vector<Vector2f> outer = {{0, 0}, {3, 0}, {3, 3}, {0, 3}};
  std::vector<Vector2f> hole = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
  AddContour2f(tess, outer);
  AddContour2f(tess, hole, /*reverse=*/true);

  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetContourCount(tess), 2);
  const TESSContourInfo* contours = tessGetContours(tess);
  ASSERT_NE(contours, nullptr);
  EXPECT_EQ(contours[0].role, TESS_CONTOUR_OUTER);
  EXPECT_EQ(contours[1].role, TESS_CONTOUR_HOLE);
  EXPECT_EQ(contours[0].vertexCount, 4);
  EXPECT_EQ(contours[1].vertexCount, 4);
  EXPECT_GT(contours[0].signedArea, 0);
}

// Every triangle is attributed to the single outer contour; hole-adjacent
// triangles carry the hole id, the rest report TESS_UNDEF.
TEST_F(Libtess2Test, TriangleInfoSquareWithHole) {
  std::vector<Vector2f> outer = {{0, 0}, {3, 0}, {3, 3}, {0, 3}};
  std::vector<Vector2f> hole = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
  AddContour2f(tess, outer);
  AddContour2f(tess, hole, /*reverse=*/true);

  ASSERT_NE(TriangulateForInfo(tess), 0);
  const int tri_count = tessGetTriangleInfoCount(tess);
  ASSERT_EQ(tri_count, tessGetElementCount(tess));
  ASSERT_EQ(tri_count, 8);
  const TESStriangleInfo* info = tessGetTriangleInfo(tess);
  ASSERT_NE(info, nullptr);

  int touching_hole = 0;
  for (int i = 0; i < tri_count; ++i) {
    EXPECT_EQ(info[i].outerContour, 0);
    EXPECT_TRUE(info[i].holeContour == TESS_UNDEF ||
                info[i].holeContour == 1);
    if (info[i].holeContour == 1) {
      ++touching_hole;
      bool tagged = false;
      for (int edge = 0; edge < 3; ++edge)
        tagged |= info[i].edgeContours[edge] == 1;
      EXPECT_TRUE(tagged);
    }
    for (int edge = 0; edge < 3; ++edge) {
      EXPECT_TRUE(info[i].edgeContours[edge] == TESS_UNDEF ||
                  info[i].edgeContours[edge] == 0 ||
                  info[i].edgeContours[edge] == 1);
    }
  }
  EXPECT_GT(touching_hole, 0);
}

// No triangle centroid may fall inside the hole: the hole interior stays
// empty, which is the collision-filtering guarantee.
TEST_F(Libtess2Test, TriangleInfoNoTriangleCoversHole) {
  std::vector<Vector2f> outer = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
  std::vector<Vector2f> hole = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
  AddContour2f(tess, outer);
  AddContour2f(tess, hole, /*reverse=*/true);

  ASSERT_NE(TriangulateForInfo(tess), 0);
  const TESSreal* verts = tessGetVertices(tess);
  const TESSindex* elems = tessGetElements(tess);
  const int tri_count = tessGetTriangleInfoCount(tess);
  ASSERT_GT(tri_count, 0);

  for (int i = 0; i < tri_count; ++i) {
    const TESSindex* tri = &elems[i * kNumTriangleVertices];
    double cx = 0, cy = 0;
    for (int corner = 0; corner < 3; ++corner) {
      const Vector2f* p = reinterpret_cast<const Vector2f*>(
          &verts[tri[corner] * kComponentCount]);
      cx += p->x;
      cy += p->y;
    }
    cx /= 3.0;
    cy /= 3.0;
    // Strictly interior points of the hole (with margin) are never covered.
    EXPECT_FALSE(cx > 1.25 && cx < 2.75 && cy > 1.25 && cy < 2.75);
  }
}

// All emitted triangles share one winding, matching their geometric
// signed-area sign, so downstream texture mapping sees a consistent front.
TEST_F(Libtess2Test, TriangleInfoUniformWinding) {
  std::vector<Vector2f> outer = {{0, 0}, {3, 0}, {3, 3}, {0, 3}};
  std::vector<Vector2f> hole = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
  AddContour2f(tess, outer);
  AddContour2f(tess, hole, /*reverse=*/true);

  ASSERT_NE(TriangulateForInfo(tess), 0);
  const TESSreal* verts = tessGetVertices(tess);
  const TESSindex* elems = tessGetElements(tess);
  const TESStriangleInfo* info = tessGetTriangleInfo(tess);
  const int tri_count = tessGetTriangleInfoCount(tess);
  ASSERT_GT(tri_count, 0);

  const int expected = info[0].winding;
  EXPECT_NE(expected, TESS_WINDING_FLAT);
  for (int i = 0; i < tri_count; ++i) {
    EXPECT_EQ(info[i].winding, expected);
    const double area = TriangleArea2d(verts, &elems[i * 3]);
    if (expected == TESS_WINDING_CCW)
      EXPECT_GT(area, 0.0);
    else
      EXPECT_LT(area, 0.0);
  }
}

// Two disconnected outer contours are attributed independently.
TEST_F(Libtess2Test, TriangleInfoDisjointOuters) {
  AddContour2f(tess, {{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  AddContour2f(tess, {{5, 5}, {6, 5}, {6, 6}, {5, 6}});

  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetContourCount(tess), 2);
  const TESStriangleInfo* info = tessGetTriangleInfo(tess);
  const int tri_count = tessGetTriangleInfoCount(tess);
  ASSERT_EQ(tri_count, 4);
  bool saw_outer_0 = false, saw_outer_1 = false;
  for (int i = 0; i < tri_count; ++i) {
    EXPECT_TRUE(info[i].outerContour == 0 || info[i].outerContour == 1);
    EXPECT_EQ(info[i].holeContour, TESS_UNDEF);
    saw_outer_0 |= info[i].outerContour == 0;
    saw_outer_1 |= info[i].outerContour == 1;
  }
  EXPECT_TRUE(saw_outer_0 && saw_outer_1);
}

// Two holes in one outer contour: both hole ids appear, never mixed into one.
TEST_F(Libtess2Test, TriangleInfoTwoHoles) {
  AddContour2f(tess, {{0, 0}, {10, 0}, {10, 10}, {0, 10}});
  AddContour2f(tess, {{1, 1}, {3, 1}, {3, 3}, {1, 3}}, /*reverse=*/true);
  AddContour2f(tess, {{7, 7}, {9, 7}, {9, 9}, {7, 9}}, /*reverse=*/true);

  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetContourCount(tess), 3);
  EXPECT_EQ(tessGetContours(tess)[1].role, TESS_CONTOUR_HOLE);
  EXPECT_EQ(tessGetContours(tess)[2].role, TESS_CONTOUR_HOLE);

  const TESStriangleInfo* info = tessGetTriangleInfo(tess);
  const int tri_count = tessGetTriangleInfoCount(tess);
  bool saw_hole_1 = false, saw_hole_2 = false;
  for (int i = 0; i < tri_count; ++i) {
    EXPECT_EQ(info[i].outerContour, 0);
    EXPECT_TRUE(info[i].holeContour == TESS_UNDEF ||
                info[i].holeContour == 1 || info[i].holeContour == 2);
    saw_hole_1 |= info[i].holeContour == 1;
    saw_hole_2 |= info[i].holeContour == 2;
  }
  EXPECT_TRUE(saw_hole_1 && saw_hole_2);
}

// Boundary: no contours / degenerate line input yields zero attribution.
TEST_F(Libtess2Test, TriangleInfoEmptyAndDegenerate) {
  // With no contours there is no mesh; tessellation reports failure and the
  // attribution getters stay empty rather than returning dangling data.
  EXPECT_EQ(TriangulateForInfo(tess), 0);
  EXPECT_EQ(tessGetContourCount(tess), 0);
  EXPECT_EQ(tessGetTriangleInfoCount(tess), 0);
  EXPECT_EQ(tessGetTriangleInfo(tess), nullptr);
  EXPECT_EQ(tessGetContours(tess), nullptr);

  AddContour2f(tess, {{0, 0}, {0, 1}});
  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetContourCount(tess), 1);
  EXPECT_EQ(tessGetContours(tess)[0].role, TESS_CONTOUR_DEGENERATE);
  EXPECT_EQ(tessGetTriangleInfoCount(tess), 0);
}

// Boundary: a single triangle is fully attributed to its only contour and
// all three edges are boundary edges of that contour.
TEST_F(Libtess2Test, TriangleInfoSingleTriangle) {
  AddContour2f(tess, {{0, 0}, {0, 1}, {1, 0}});
  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetTriangleInfoCount(tess), 1);
  const TESStriangleInfo& info = tessGetTriangleInfo(tess)[0];
  EXPECT_EQ(info.outerContour, 0);
  EXPECT_EQ(info.holeContour, TESS_UNDEF);
  for (int edge = 0; edge < 3; ++edge)
    EXPECT_EQ(info.edgeContours[edge], 0);
  EXPECT_EQ(info.winding, TESS_WINDING_CCW);
}

// Boundary: attribution is only produced for TESS_POLYGONS with polySize 3.
TEST_F(Libtess2Test, TriangleInfoOnlyForTriangles) {
  AddContour2f(tess, {{0, 0}, {0, 1}, {1, 1}, {1, 0}});
  ASSERT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          /*polySize=*/4, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetTriangleInfoCount(tess), 0);
  EXPECT_EQ(tessGetTriangleInfo(tess), nullptr);
}

// Reusing a tesselator after a run resets the attribution for the new input.
TEST_F(Libtess2Test, TriangleInfoResetOnReuse) {
  AddContour2f(tess, {{0, 0}, {0, 1}, {1, 0}});
  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetContourCount(tess), 1);

  AddContour2f(tess, {{2, 2}, {3, 2}, {3, 3}});
  ASSERT_NE(TriangulateForInfo(tess), 0);
  EXPECT_EQ(tessGetContourCount(tess), 1);
  EXPECT_EQ(tessGetContours(tess)[0].vertexCount, 3);
  ASSERT_EQ(tessGetTriangleInfoCount(tess), 1);
  EXPECT_EQ(tessGetTriangleInfo(tess)[0].outerContour, 0);
}

// Invalid coordinates fail tessellation and leave no attribution output.
TEST_F(Libtess2Test, TriangleInfoInvalidInput) {
  AddContour2f(tess, {{-2e+37f, 0.f}, {0, 5}, {1e37f, -5}});
  EXPECT_EQ(TriangulateForInfo(tess), 0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  EXPECT_EQ(tessGetTriangleInfoCount(tess), 0);
}

// A rejected contour is not registered, so valid contours keep indices that
// match their tessAddContour() call order one-for-one.
TEST_F(Libtess2Test, TriangleInfoRejectedContourKeepsIndices) {
  AddContour2f(tess, {{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  tessSetOption(tess, TESS_REVERSE_CONTOURS, 0);
  tessAddContour(tess, kComponentCount,
                 std::vector<Vector2f>{{-2e+37f, 0.f}, {0, 5}, {1e37f, -5}}
                     .data(),
                 kComponentCount * sizeof(TESSreal), 3);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  EXPECT_EQ(tessGetContourCount(tess), 1);
}

// Deterministic consistency data for the canonical square-with-one-hole
// fixture. This is the golden snapshot texture mapping and collision code can
// rely on: 8 triangles, all owned by outer contour 0, four bordering hole 1,
// uniform CCW winding, and every boundary edge tagged with its source contour.
TEST_F(Libtess2Test, TriangleInfoGoldenSquareWithHole) {
  std::vector<Vector2f> outer = {{0, 0}, {3, 0}, {3, 3}, {0, 3}};
  std::vector<Vector2f> hole = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
  AddContour2f(tess, outer);
  AddContour2f(tess, hole, /*reverse=*/true);

  ASSERT_NE(TriangulateForInfo(tess), 0);
  ASSERT_EQ(tessGetContourCount(tess), 2);
  ASSERT_EQ(tessGetTriangleInfoCount(tess), 8);

  const TESSContourInfo* contours = tessGetContours(tess);
  EXPECT_EQ(contours[0].role, TESS_CONTOUR_OUTER);
  EXPECT_EQ(contours[1].role, TESS_CONTOUR_HOLE);

  const TESStriangleInfo* info = tessGetTriangleInfo(tess);
  int hole_border = 0;
  int boundary_edges = 0;
  for (int i = 0; i < 8; ++i) {
    // Every triangle belongs to the single outer contour.
    EXPECT_EQ(info[i].outerContour, 0);
    // The convenience hole id agrees with the per-edge tagging.
    if (info[i].holeContour == 1) {
      ++hole_border;
      int hole_edges = 0;
      for (int edge = 0; edge < 3; ++edge)
        if (info[i].edgeContours[edge] == 1) ++hole_edges;
      EXPECT_GE(hole_edges, 1);
    } else {
      EXPECT_EQ(info[i].holeContour, TESS_UNDEF);
    }
    for (int edge = 0; edge < 3; ++edge) {
      const TESSindex c = info[i].edgeContours[edge];
      EXPECT_TRUE(c == TESS_UNDEF || c == 0 || c == 1);
      if (c != TESS_UNDEF) ++boundary_edges;
    }
    // Winding is uniform across the whole triangulation.
    EXPECT_EQ(info[i].winding, TESS_WINDING_CCW);
  }
  EXPECT_EQ(hole_border, 4);
  // Outer boundary (4) plus hole boundary (4) edges are all attributed.
  EXPECT_EQ(boundary_edges, 8);
}

// Zero-area sliver triangles are reported FLAT so collision filtering can
// drop them while the remaining triangles keep a single, uniform winding.
TEST_F(Libtess2Test, TriangleInfoFlagsFlatSliver) {
  // Two overlapping squares produce a collinear sliver under NONZERO winding.
  AddContour2f(tess, {{0, 0}, {4, 0}, {4, 4}, {0, 4}});
  AddContour2f(tess, {{2, 2}, {6, 2}, {6, 6}, {2, 6}});
  ASSERT_NE(tessTesselate(tess, TESS_WINDING_NONZERO, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);

  const TESStriangleInfo* info = tessGetTriangleInfo(tess);
  const int tri_count = tessGetTriangleInfoCount(tess);
  int flat = 0, oriented = 0;
  for (int i = 0; i < tri_count; ++i) {
    if (info[i].winding == TESS_WINDING_FLAT)
      ++flat;
    else
      ++oriented;
  }
  EXPECT_GE(flat, 1);
  EXPECT_GT(oriented, 0);
  // Every non-flat triangle still shares one direction.
  int direction = TESS_WINDING_FLAT;
  for (int i = 0; i < tri_count; ++i) {
    if (info[i].winding == TESS_WINDING_FLAT) continue;
    if (direction == TESS_WINDING_FLAT)
      direction = info[i].winding;
    else
      EXPECT_EQ(info[i].winding, direction);
  }
}

}  // namespace
