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
  // A two-vertex contour has zero area and is now rejected with a location.
  AddPolyline(tess, {{0, 0}, {0, 1}});
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  ASSERT_EQ(tessGetErrorCount(tess), 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error,
            TESS_ERROR_ZERO_AREA_CONTOUR);
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
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  EXPECT_GT(tessGetErrorCount(tess), 0);
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
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  ASSERT_EQ(tessGetErrorCount(tess), 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error,
            TESS_ERROR_INVALID_COORDINATE);
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

// --- Degenerate contour handling with precise error locations ---

static int TesselateAndExpectRejection(TESStesselator* tess,
                                       const std::vector<Vector2f>& polyline) {
  AddPolyline(tess, polyline);
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_NONZERO, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  return tessGetErrorCount(tess);
}

TEST_F(Libtess2Test, RejectsSelfIntersectingBowtie) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {1, 1}, {1, 0}, {0, 1}});
  EXPECT_EQ(count, 2);
  const TESSContourErrorRecord* errors = tessGetErrors(tess);
  EXPECT_EQ(errors[0].error, TESS_ERROR_SELF_INTERSECTION);
  EXPECT_EQ(errors[0].contour, 0);
  EXPECT_EQ(errors[0].vertex, 0);
  EXPECT_EQ(errors[0].vertex2, 1);
}

TEST_F(Libtess2Test, RejectsConsecutiveDuplicateWithLocation) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {1, 0}, {1, 0}, {1, 1}, {0, 1}});
  ASSERT_EQ(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error, TESS_ERROR_DUPLICATE_VERTEX);
  EXPECT_EQ(tessGetErrors(tess)[0].vertex, 2);
}

TEST_F(Libtess2Test, RejectsClosingDuplicateAtFirstVertex) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}});
  ASSERT_EQ(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error, TESS_ERROR_DUPLICATE_VERTEX);
  EXPECT_EQ(tessGetErrors(tess)[0].vertex, 0);
}

TEST_F(Libtess2Test, RejectsCollinearSpike) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}});
  ASSERT_EQ(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error, TESS_ERROR_COLLINEAR_SPIKE);
  EXPECT_EQ(tessGetErrors(tess)[0].vertex, 1);
}

TEST_F(Libtess2Test, AcceptsBenignCollinearMidpoint) {
  AddPolyline(tess, {{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_NONZERO, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_OK);
  EXPECT_EQ(tessGetErrorCount(tess), 0);
  EXPECT_EQ(tessGetElementCount(tess), 3);
}

TEST_F(Libtess2Test, RejectsZeroAreaHoleWithContourIndex) {
  std::vector<Vector2f> outer = {
      {0, 0}, {4, 0}, {4, 4}, {0, 4}};
  std::vector<Vector2f> hole = {
      {1, 1}, {2, 1}, {1, 1}};
  tessAddContour(tess, kComponentCount, outer.data(),
                 kComponentCount * sizeof(TESSreal), outer.size());
  tessSetOption(tess, TESS_REVERSE_CONTOURS, 1);
  tessAddContour(tess, kComponentCount, hole.data(),
                 kComponentCount * sizeof(TESSreal), hole.size());
  EXPECT_EQ(tessTesselate(tess, TESS_WINDING_NONZERO, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_INVALID_INPUT);
  bool foundHoleError = false;
  for (int i = 0; i < tessGetErrorCount(tess); ++i) {
    const TESSContourErrorRecord& record = tessGetErrors(tess)[i];
    if (record.contour == 1) foundHoleError = true;
  }
  EXPECT_TRUE(foundHoleError);
}

TEST_F(Libtess2Test, RejectsTJunction) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {2, 0}, {2, 2}, {1, 0}, {0, 2}});
  EXPECT_GE(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error, TESS_ERROR_SELF_INTERSECTION);
}

TEST_F(Libtess2Test, RejectsCollinearEdgeOverlap) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {4, 0}, {4, 4}, {1, 0}, {0, 4}});
  EXPECT_GE(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error, TESS_ERROR_SELF_INTERSECTION);
}

TEST_F(Libtess2Test, RejectsPentagram) {
  int count = TesselateAndExpectRejection(
      tess, {{1, 2}, {1.2245f, 1.309f}, {1.951f, 1.309f},
             {1.363f, 0.882f}, {1.588f, 0.191f}});
  EXPECT_GE(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error, TESS_ERROR_SELF_INTERSECTION);
}

TEST_F(Libtess2Test, RejectsOutOfRangeCoordinate) {
  int count = TesselateAndExpectRejection(
      tess, {{0, 0}, {9e9f, 0}, {9e9f, 9e9f}});
  ASSERT_EQ(count, 1);
  EXPECT_EQ(tessGetErrors(tess)[0].error,
            TESS_ERROR_INVALID_COORDINATE);
}

TEST_F(Libtess2Test, TesselatorReusableAfterSuccessfulRun) {
  AddPolyline(tess, {{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  ASSERT_NE(tessTesselate(tess, TESS_WINDING_NONZERO, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetErrorCount(tess), 0);

  AddPolyline(tess, {{0, 0}, {0, 1}, {1, 0}});
  EXPECT_NE(tessTesselate(tess, TESS_WINDING_NONZERO, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  EXPECT_EQ(tessGetStatus(tess), TESS_STATUS_OK);
  EXPECT_EQ(tessGetElementCount(tess), 1);
}

TEST_F(Libtess2Test, NoFlippedTrianglesForValidPolygonWithHole) {
  AddPolygonWithHole(tess);
  ASSERT_NE(tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS,
                          kNumTriangleVertices, kComponentCount, nullptr),
            0);
  const TESSindex* elements = tessGetElements(tess);
  const TESSreal* vertices = tessGetVertices(tess);
  for (int i = 0; i < tessGetElementCount(tess); ++i) {
    const TESSindex* triangle = &elements[i * kNumTriangleVertices];
    const TESSreal* a = &vertices[triangle[0] * kComponentCount];
    const TESSreal* b = &vertices[triangle[1] * kComponentCount];
    const TESSreal* c = &vertices[triangle[2] * kComponentCount];
    float twiceArea = (b[0] - a[0]) * (c[1] - a[1]) -
                      (b[1] - a[1]) * (c[0] - a[0]);
    EXPECT_GT(twiceArea, 0.f);
  }
}

}  // namespace
