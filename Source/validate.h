/*
** SGI FREE SOFTWARE LICENSE B (Version 2.0, Sept. 18, 2008)
** Copyright (C) [dates of first publication] Silicon Graphics, Inc.
** All Rights Reserved.
*/

#ifndef TESS_VALIDATE_H
#define TESS_VALIDATE_H

#include "tesselator.h"

#ifdef __cplusplus
extern "C" {
#endif

// A packed, validated input contour projected into 2D. The points are stored
// tightly so the validator never dereferences caller memory again.
typedef struct TESSContour2D {
	int pointOffset;  // index into the tesselator point pool (x/y pairs).
	int count;        // number of vertices.
	int reverse;      // mirrors TESS_REVERSE_CONTOURS at insertion time.
} TESSContour2D;

// Returns 1 if the finite value 'coord' lies inside the supported input range.
int tessIsValidCoord( TESSreal coord );

// Validates a single contour independently (range, duplicates, collinear
// spikes, self-intersection, zero area). Records every problem found through
// record_fn. 'contour_index' is the tessAddContour() call index.
//
// Returns 1 when the contour is internally valid, 0 when at least one error
// was recorded. record_fn returning 0 indicates an allocation failure and
// aborts validation.
typedef int (*TESSErrorRecorder)( void *user, int error, int contour,
								  int vertex, int vertex2 );

int tessValidateContour( const TESSreal *pointPool,
						 const TESSContour2D *contour, int contour_index,
						 TESSErrorRecorder record_fn, void *user );

// Checks every segment of contour 'b' against every segment of contour 'a'
// (a precedes b in insertion order). Records one SELF_INTERSECTION entry per
// contour for each crossing, T-contact or positive-length collinear overlap.
// record_fn returning 0 aborts the scan. Returns 1 when no interaction was
// found, 0 otherwise.
int tessValidateContourPair( const TESSreal *pointPool,
							 const TESSContour2D *a, int index_a,
							 const TESSContour2D *b, int index_b,
							 TESSErrorRecorder record_fn, void *user );

#ifdef __cplusplus
};
#endif

#endif
