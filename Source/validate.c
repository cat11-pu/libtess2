/*
** SGI FREE SOFTWARE LICENSE B (Version 2.0, Sept. 18, 2008)
** Copyright (C) [dates of first publication] Silicon Graphics, Inc.
** All Rights Reserved.
*/

/*
** Up-front geometric validation of input contours. Every predicate below uses
** exact floating point comparisons on the caller's coordinate domain. Inputs
** are clamped to +-2^23 by tessIsValidCoord(), so differences, cross products
** and dot products cannot overflow a finite float.
*/

#include <math.h>

#include "validate.h"

#define TESS_UNDEF_INDEX (~(TESSindex)0)

typedef struct TESSVec2 {
	TESSreal x, y;
} TESSVec2;

int tessIsValidCoord( TESSreal coord )
{
	return isfinite( coord ) &&
		coord <= TESS_MAX_VALID_INPUT_VALUE &&
		coord >= TESS_MIN_VALID_INPUT_VALUE;
}

static TESSVec2 PointAt( const TESSreal *pool, const TESSContour2D *c, int i )
{
	TESSVec2 p;
	p.x = pool[( c->pointOffset + i ) * 2];
	p.y = pool[( c->pointOffset + i ) * 2 + 1];
	return p;
}

static int PointsEqual( TESSVec2 a, TESSVec2 b )
{
	return a.x == b.x && a.y == b.y;
}

static TESSreal Cross( TESSVec2 o, TESSVec2 a, TESSVec2 b )
{
	return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

/* Dot product of the two edge directions incident on 'at'. A spike goes
* out to 'at' and back, so both edges point the same way (positive dot). */
static TESSreal DotDir( TESSVec2 from, TESSVec2 at, TESSVec2 to )
{
	return (at.x - from.x) * (at.x - to.x) + (at.y - from.y) * (at.y - to.y);
}

static int BBoxReject( TESSVec2 a, TESSVec2 b, TESSVec2 c, TESSVec2 d )
{
	TESSreal abMinX = a.x < b.x ? a.x : b.x;
	TESSreal abMaxX = a.x > b.x ? a.x : b.x;
	TESSreal abMinY = a.y < b.y ? a.y : b.y;
	TESSreal abMaxY = a.y > b.y ? a.y : b.y;
	TESSreal cdMinX = c.x < d.x ? c.x : d.x;
	TESSreal cdMaxX = c.x > d.x ? c.x : d.x;
	TESSreal cdMinY = c.y < d.y ? c.y : d.y;
	TESSreal cdMaxY = c.y > d.y ? c.y : d.y;

	return abMaxX < cdMinX || cdMaxX < abMinX ||
		abMaxY < cdMinY || cdMaxY < abMinY;
}

/* Strictly between the segment endpoints (collinearity already known). */
static int OnSegmentInterior( TESSVec2 p, TESSVec2 a, TESSVec2 b )
{
	int insideX = (p.x > a.x && p.x < b.x) || (p.x > b.x && p.x < a.x);
	int insideY = (p.y > a.y && p.y < b.y) || (p.y > b.y && p.y < a.y);
	return insideX || insideY;
}

/*
* Classifies the interaction of segment (a,b) with segment (c,d).
* Returns 1 for a proper crossing, an interior T-contact or a positive
* length collinear overlap. Shared endpoints alone are not interactions.
*/
static int SegmentsInteract( TESSVec2 a, TESSVec2 b, TESSVec2 c, TESSVec2 d )
{
	TESSreal abc, abd, cda, cdb;

	if( BBoxReject( a, b, c, d ) )
		return 0;

	abc = Cross( a, b, c );
	abd = Cross( a, b, d );
	cda = Cross( c, d, a );
	cdb = Cross( c, d, b );

	if( (abc > 0 && abd < 0) || (abc < 0 && abd > 0) ) {
		if( (cda > 0 && cdb < 0) || (cda < 0 && cdb > 0) )
			return 1; /* proper crossing */
	}

	/* Endpoint of one segment lands in the interior of the other. */
	if( abc == 0 && OnSegmentInterior( c, a, b ) )
		return 1;
	if( abd == 0 && OnSegmentInterior( d, a, b ) )
		return 1;
	if( cda == 0 && OnSegmentInterior( a, c, d ) )
		return 1;
	if( cdb == 0 && OnSegmentInterior( b, c, d ) )
		return 1;

	/* Collinear segments with a positive-length shared run. */
	if( abc == 0 && abd == 0 ) {
		TESSreal overlapX, overlapY;
		TESSreal abMinX = a.x < b.x ? a.x : b.x;
		TESSreal abMaxX = a.x > b.x ? a.x : b.x;
		TESSreal abMinY = a.y < b.y ? a.y : b.y;
		TESSreal abMaxY = a.y > b.y ? a.y : b.y;
		TESSreal cdMinX = c.x < d.x ? c.x : d.x;
		TESSreal cdMaxX = c.x > d.x ? c.x : d.x;
		TESSreal cdMinY = c.y < d.y ? c.y : d.y;
		TESSreal cdMaxY = c.y > d.y ? c.y : d.y;

		overlapX = (abMaxX < cdMaxX ? abMaxX : cdMaxX) -
			(abMinX > cdMinX ? abMinX : cdMinX);
		overlapY = (abMaxY < cdMaxY ? abMaxY : cdMaxY) -
			(abMinY > cdMinY ? abMinY : cdMinY);
		if( overlapX > 0 || overlapY > 0 )
			return 1;
	}

	return 0;
}

static int Record( TESSErrorRecorder fn, void *user, int error,
				   int contour, int vertex, int vertex2 )
{
	return fn( user, error, contour, vertex, vertex2 );
}

int tessValidateContour( const TESSreal *pointPool,
						 const TESSContour2D *contour, int contour_index,
						 TESSErrorRecorder record_fn, void *user )
{
	const int n = contour->count;
	int i, j;
	int errors = 0;
	TESSreal area2 = 0;

	/* 1. Coordinate range and finiteness. */
	for( i = 0; i < n; ++i ) {
		if( !tessIsValidCoord( pointPool[( contour->pointOffset + i ) * 2] ) ||
			!tessIsValidCoord( pointPool[( contour->pointOffset + i ) * 2 + 1] ) ) {
			if( !Record( record_fn, user, TESS_ERROR_INVALID_COORDINATE,
						 contour_index, i, TESS_UNDEF_INDEX ) )
				return 0;
			return 0;
		}
	}

	if( n < 2 ) {
		if( !Record( record_fn, user, TESS_ERROR_ZERO_AREA_CONTOUR,
					 contour_index, 0, TESS_UNDEF_INDEX ) )
			return 0;
		return 0;
	}

	/* 2. Adjacent duplicate vertices, including the closing point. */
	for( i = 0; i < n; ++i ) {
		int next = (i + 1) % n;
		if( PointsEqual( PointAt( pointPool, contour, i ), PointAt( pointPool, contour, next ) ) ) {
			/* A repeated closing point is attributed to the first vertex. */
			int reported = ( next == 0 ) ? 0 : next;
			if( !Record( record_fn, user, TESS_ERROR_DUPLICATE_VERTEX,
						 contour_index, reported, TESS_UNDEF_INDEX ) )
				return 0;
			++errors;
		}
	}
	if( errors > 0 )
		return 0;

	/* Exact shoelace area, computed once for the spike/area disambiguation
	* and as the final zero-area test. */
	for( i = 0; i < n; ++i ) {
		TESSVec2 p = PointAt( pointPool, contour, i );
		TESSVec2 q = PointAt( pointPool, contour, (i + 1) % n );
		area2 += p.x * q.y - q.x * p.y;
	}

	/* 3. Out-and-back collinear spikes. A fully collinear contour has zero
	* total area and is reported as ZERO_AREA instead; a local spur in an
	* otherwise valid polygon keeps its non-zero area and is a real spike. */
	for( i = 0; i < n; ++i ) {
		int prev = (i + n - 1) % n;
		int next = (i + 1) % n;
		TESSVec2 at = PointAt( pointPool, contour, i );
		TESSVec2 before = PointAt( pointPool, contour, prev );
		TESSVec2 after = PointAt( pointPool, contour, next );
		if( Cross( before, at, after ) == 0 &&
			DotDir( before, at, after ) > 0 ) {
			if( area2 == 0 ) {
				if( !Record( record_fn, user, TESS_ERROR_ZERO_AREA_CONTOUR,
							 contour_index, 0, TESS_UNDEF_INDEX ) )
					return 0;
				return 0;
			}
			if( !Record( record_fn, user, TESS_ERROR_COLLINEAR_SPIKE,
						 contour_index, i, TESS_UNDEF_INDEX ) )
				return 0;
			return 0;
		}
	}

	/* 4. Non-adjacent segment interactions within this contour. */
	for( i = 0; i < n; ++i ) {
		TESSVec2 a = PointAt( pointPool, contour, i );
		TESSVec2 b = PointAt( pointPool, contour, (i + 1) % n );
		for( j = i + 1; j < n; ++j ) {
			if( j == i + 1 || (i == 0 && j == n - 1) )
				continue; /* adjacent segments share a legal endpoint */
			if( SegmentsInteract( a, b,
								  PointAt( pointPool, contour, j ),
								  PointAt( pointPool, contour, (j + 1) % n ) ) ) {
				if( !Record( record_fn, user, TESS_ERROR_SELF_INTERSECTION,
							 contour_index, i, (i + 1) % n ) ||
					!Record( record_fn, user, TESS_ERROR_SELF_INTERSECTION,
							 contour_index, j, (j + 1) % n ) )
					return 0;
				return 0;
			}
		}
	}

	/* 5. Zero area with no other defect: point cloud, segment or collinear
	* chain. Self-intersection takes priority, so a bowtie is caught above. */
	if( area2 == 0 ) {
		if( !Record( record_fn, user, TESS_ERROR_ZERO_AREA_CONTOUR,
					 contour_index, 0, TESS_UNDEF_INDEX ) )
			return 0;
		return 0;
	}

	return 1;
}

int tessValidateContourPair( const TESSreal *pointPool,
							 const TESSContour2D *a, int index_a,
							 const TESSContour2D *b, int index_b,
							 TESSErrorRecorder record_fn, void *user )
{
	int i, j;
	int found = 0;

	for( i = 0; i < a->count; ++i ) {
		TESSVec2 pa = PointAt( pointPool, a, i );
		TESSVec2 pb = PointAt( pointPool, a, (i + 1) % a->count );
		for( j = 0; j < b->count; ++j ) {
			if( SegmentsInteract( pa, pb,
								  PointAt( pointPool, b, j ),
								  PointAt( pointPool, b, (j + 1) % b->count ) ) ) {
				if( !Record( record_fn, user, TESS_ERROR_SELF_INTERSECTION,
							 index_a, i, (i + 1) % a->count ) ||
					!Record( record_fn, user, TESS_ERROR_SELF_INTERSECTION,
							 index_b, j, (j + 1) % b->count ) )
					return 0;
				found = 1;
			}
		}
	}

	return !found;
}
