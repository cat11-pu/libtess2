#include <stdio.h>
#include <stdlib.h>
#include "tesselator.h"

static void* pa(void* u, unsigned int s) { (void)u; return malloc(s); }
static void* pn(void* u, void* p, unsigned int s) {
  (void)u; (void)p; (void)s; return NULL; /* no dynamic growth allowed */
}
static void pf(void* u, void* p) { (void)u; free(p); }

int main(void) {
  TESSalloc a = {pa, pn, pf, 0,0,0,0,0, /*extraVertices=*/64};
  TESStesselator* t = tessNewTess(&a);
  const float outer[8] = {0,0, 3,0, 3,3, 0,3};
  const float hole[8] = {1,1, 2,1, 2,2, 1,2};
  tessSetOption(t, TESS_REVERSE_CONTOURS, 0);
  tessAddContour(t, 2, outer, sizeof(float)*2, 4);
  tessSetOption(t, TESS_REVERSE_CONTOURS, 1);
  tessAddContour(t, 2, hole, sizeof(float)*2, 4);
  int ok = tessTesselate(t, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2, NULL);
  printf("no-realloc ok=%d status=%d tris=%d info=%d contours=%d\n", ok,
         tessGetStatus(t), tessGetElementCount(t),
         tessGetTriangleInfoCount(t), tessGetContourCount(t));
  if (ok && tessGetTriangleInfoCount(t) > 0)
    printf("first outer=%d hole=%d wind=%d\n",
           tessGetTriangleInfo(t)[0].outerContour,
           tessGetTriangleInfo(t)[0].holeContour,
           tessGetTriangleInfo(t)[0].winding);
  tessDeleteTess(t);
  return 0;
}
