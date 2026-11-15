#include <stdio.h>
#include <stdlib.h>
#include "tesselator.h"

static void add(TESStesselator* t, const float* p, int n, int reverse) {
  tessSetOption(t, TESS_REVERSE_CONTOURS, reverse);
  tessAddContour(t, 2, p, sizeof(float) * 2, n);
}

static void* pa(void* u, unsigned int s) { (void)u; return malloc(s); }
static void* pr(void* u, void* p, unsigned int s) { (void)u; return realloc(p, s); }
static void pf(void* u, void* p) { (void)u; free(p); }

static void run(const char* name, const float* outer, int no,
                const float* hole, int nh) {
  TESStesselator* t = tessNewTess(NULL);
  add(t, outer, no, 0);
  if (hole) add(t, hole, nh, 1);
  int ok = tessTesselate(t, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2, NULL);
  printf("== %s ok=%d tris=%d contours=%d triInfo=%d\n", name, ok,
         tessGetElementCount(t), tessGetContourCount(t),
         tessGetTriangleInfoCount(t));
  const TESSContourInfo* ci = tessGetContours(t);
  for (int i = 0; i < tessGetContourCount(t); ++i)
    printf("  contour[%d] role=%d area=%g verts=%d\n", i, ci[i].role,
           ci[i].signedArea, ci[i].vertexCount);
  const TESStriangleInfo* ti = tessGetTriangleInfo(t);
  const TESSindex* el = tessGetElements(t);
  for (int i = 0; i < tessGetTriangleInfoCount(t); ++i) {
    printf("  tri[%d] v=(%d,%d,%d) outer=%d hole=%d edge=(%d,%d,%d) wind=%d\n",
           i, el[i*3], el[i*3+1], el[i*3+2], ti[i].outerContour,
           ti[i].holeContour, ti[i].edgeContours[0], ti[i].edgeContours[1],
           ti[i].edgeContours[2], ti[i].winding);
  }
  tessDeleteTess(t);
}

int main(void) {
  const float square[8] = {0,0, 3,0, 3,3, 0,3};
  const float hole[8] = {1,1, 2,1, 2,2, 1,2};
  const float tri[6] = {0,0, 0,1, 1,0};
  const float quad[8] = {0,0, 0,1, 1,1, 1,0};
  run("square+hole", square, 4, hole, 4);
  run("triangle", tri, 3, NULL, 0);
  run("quad", quad, 4, NULL, 0);

  /* two disjoint squares */
  TESStesselator* t = tessNewTess(NULL);
  const float s1[8] = {0,0, 1,0, 1,1, 0,1};
  const float s2[8] = {5,5, 6,5, 6,6, 5,6};
  add(t, s1, 4, 0);
  add(t, s2, 4, 0);
  tessTesselate(t, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2, NULL);
  printf("== two squares tris=%d contours=%d\n", tessGetElementCount(t),
         tessGetContourCount(t));
  const TESStriangleInfo* ti = tessGetTriangleInfo(t);
  for (int i = 0; i < tessGetTriangleInfoCount(t); ++i)
    printf("  tri[%d] outer=%d wind=%d\n", i, ti[i].outerContour,
           ti[i].winding);
  tessDeleteTess(t);

  /* CW outer contour: orientation normalisation still yields a single
   * uniform winding for every triangle. */
  t = tessNewTess(NULL);
  const float cw[8] = {0,0, 0,3, 3,3, 3,0};
  add(t, cw, 4, 0);
  tessTesselate(t, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2, NULL);
  printf("== CW square tris=%d wind0=%d\n", tessGetElementCount(t),
         tessGetTriangleInfo(t)[0].winding);
  tessDeleteTess(t);

  /* overlapping outers under NONZERO: overlap triangles may be unattributed */
  t = tessNewTess(NULL);
  const float a[8] = {0,0, 4,0, 4,4, 0,4};
  const float b[8] = {2,2, 6,2, 6,6, 2,6};
  add(t, a, 4, 0);
  add(t, b, 4, 0);
  tessTesselate(t, TESS_WINDING_NONZERO, TESS_POLYGONS, 3, 2, NULL);
  printf("== overlap NONZERO tris=%d contours=%d\n",
         tessGetElementCount(t), tessGetContourCount(t));
  for (int i = 0; i < tessGetTriangleInfoCount(t); ++i)
    printf("  tri[%d] outer=%d wind=%d\n", i,
           tessGetTriangleInfo(t)[i].outerContour,
           tessGetTriangleInfo(t)[i].winding);
  {
    const TESSreal* v = tessGetVertices(t);
    const TESSindex* el = tessGetElements(t);
    for (int i = 0; i < tessGetElementCount(t); ++i) {
      const float* p0 = &v[el[i*3]*2];
      const float* p1 = &v[el[i*3+1]*2];
      const float* p2 = &v[el[i*3+2]*2];
      double ar = 0.5*((p1[0]-p0[0])*(p2[1]-p0[1]) -
                       (p2[0]-p0[0])*(p1[1]-p0[1]));
      if (ar == 0.0)
        printf("  flat tri[%d] area=%g verts=(%g,%g)(%g,%g)(%g,%g)\n", i, ar,
               p0[0],p0[1],p1[0],p1[1],p2[0],p2[1]);
    }
  }
  tessDeleteTess(t);

  /* custom allocator path */
  TESSalloc ca = {pa, pr, pf, 0,0,0,0,0,0,0};
  t = tessNewTess(&ca);
  add(t, square, 4, 0);
  add(t, hole, 4, 1);
  int ok3 = tessTesselate(t, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2, NULL);
  printf("== custom alloc ok=%d triInfo=%d\n", ok3,
         tessGetTriangleInfoCount(t));
  tessDeleteTess(t);
  return 0;
}
