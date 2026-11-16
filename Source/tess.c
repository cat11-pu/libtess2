/*
** SGI FREE SOFTWARE LICENSE B (Version 2.0, Sept. 18, 2008)
** Copyright (C) [dates of first publication] Silicon Graphics, Inc.
** All Rights Reserved.
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
** of the Software, and to permit persons to whom the Software is furnished to do so,
** subject to the following conditions:
**
** The above copyright notice including the dates of first publication and either this
** permission notice or a reference to http://oss.sgi.com/projects/FreeB/ shall be
** included in all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
** INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
** PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL SILICON GRAPHICS, INC.
** BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
** TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
** OR OTHER DEALINGS IN THE SOFTWARE.
**
** Except as contained in this notice, the name of Silicon Graphics, Inc. shall not
** be used in advertising or otherwise to promote the sale, use or other dealings in
** this Software without prior written authorization from Silicon Graphics, Inc.
*/
/*
** Author: Eric Veach, July 1994.
*/

#include <stddef.h>
#include <assert.h>
#include <setjmp.h>
#include "bucketalloc.h"
#include "tess.h"
#include "tesselator.h"
#include "mesh.h"
#include "sweep.h"
#include "geom.h"
#include "validate.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRUE 1
#define FALSE 0

#define Dot(u,v)	(u[0]*v[0] + u[1]*v[1] + u[2]*v[2])

#if defined(FOR_TRITE_TEST_PROGRAM) || defined(TRUE_PROJECT)
static void Normalize( TESSreal v[3] )
{
	TESSreal len = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

	assert( len > 0 );
	len = sqrtf( len );
	v[0] /= len;
	v[1] /= len;
	v[2] /= len;
}
#endif

#define ABS(x)	((x) < 0 ? -(x) : (x))

/* Bounded growth for bookkeeping arrays. The first allocation always goes
* through memalloc so custom allocators without memrealloc are supported.
* Capacity doubles and is checked against an allocation cap so a hostile
* input cannot trigger an unbounded or overflowing allocation. */
#define TESS_MAX_ARRAY_CAPACITY (1 << 28)

static int tessGrowArray( TESSalloc *alloc, void **ptr, int *capacity,
						  unsigned int itemSize, int needed )
{
	void *grown;
	int newCapacity;

	if( *capacity >= needed )
		return 1;
	if( needed > TESS_MAX_ARRAY_CAPACITY )
		return 0;

	if( *ptr == NULL || *capacity == 0 ) {
		newCapacity = needed < 8 ? 8 : needed;
		grown = alloc->memalloc( alloc->userData,
								 (unsigned int)newCapacity * itemSize );
	} else {
		newCapacity = *capacity;
		while( newCapacity < needed ) {
			if( newCapacity > TESS_MAX_ARRAY_CAPACITY / 2 ) {
				newCapacity = TESS_MAX_ARRAY_CAPACITY;
				break;
			}
			newCapacity *= 2;
		}
		if( !alloc->memrealloc )
			return 0;
		grown = alloc->memrealloc( alloc->userData, *ptr,
								   (unsigned int)newCapacity * itemSize );
	}
	if( grown == NULL )
		return 0;
	*ptr = grown;
	*capacity = newCapacity;
	return 1;
}

static int RecordContourError( void *user, int error, int contour,
							   int vertex, int vertex2 )
{
	TESStesselator *tess = (TESStesselator *)user;
	TESSContourErrorRecord record;

	if( !tessGrowArray( &tess->alloc, (void **)&tess->errors,
						&tess->errorCapacity,
						sizeof( TESSContourErrorRecord ),
						tess->errorCount + 1 ) ) {
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return 0;
	}
	record.error = error;
	record.contour = contour;
	record.vertex = vertex;
	record.vertex2 = vertex2;
	tess->errors[tess->errorCount++] = record;
	return 1;
}

static int LongAxis( TESSreal v[3] )
{
	int i = 0;

	if( ABS(v[1]) > ABS(v[0]) ) { i = 1; }
	if( ABS(v[2]) > ABS(v[i]) ) { i = 2; }
	return i;
}

static int ShortAxis( TESSreal v[3] )
{
	int i = 0;

	if( ABS(v[1]) < ABS(v[0]) ) { i = 1; }
	if( ABS(v[2]) < ABS(v[i]) ) { i = 2; }
	return i;
}

static void ComputeNormal( TESStesselator *tess, TESSreal norm[3] )
{
	TESSvertex *v, *v1, *v2;
	TESSreal c, tLen2, maxLen2;
	TESSreal maxVal[3], minVal[3], d1[3], d2[3], tNorm[3];
	TESSvertex *maxVert[3], *minVert[3];
	TESSvertex *vHead = &tess->mesh->vHead;
	int i;

	v = vHead->next;
	if (v == vHead) {
		/* No vertex is initialized -- normal doesn't matter */
		norm[0] = 0; norm[1] = 0; norm[2] = 1;
		return;
	}
	for( i = 0; i < 3; ++i ) {
		c = v->coords[i];
		minVal[i] = c;
		minVert[i] = v;
		maxVal[i] = c;
		maxVert[i] = v;
	}

	for( v = vHead->next; v != vHead; v = v->next ) {
		for( i = 0; i < 3; ++i ) {
			c = v->coords[i];
			if( c < minVal[i] ) { minVal[i] = c; minVert[i] = v; }
			if( c > maxVal[i] ) { maxVal[i] = c; maxVert[i] = v; }
		}
	}

	/* Find two vertices separated by at least 1/sqrt(3) of the maximum
	* distance between any two vertices
	*/
	i = 0;
	if( maxVal[1] - minVal[1] > maxVal[0] - minVal[0] ) { i = 1; }
	if( maxVal[2] - minVal[2] > maxVal[i] - minVal[i] ) { i = 2; }
	if( minVal[i] >= maxVal[i] ) {
		/* All vertices are the same -- normal doesn't matter */
		norm[0] = 0; norm[1] = 0; norm[2] = 1;
		return;
	}

	/* Look for a third vertex which forms the triangle with maximum area
	* (Length of normal == twice the triangle area)
	*/
	maxLen2 = 0;
	v1 = minVert[i];
	v2 = maxVert[i];
	d1[0] = v1->coords[0] - v2->coords[0];
	d1[1] = v1->coords[1] - v2->coords[1];
	d1[2] = v1->coords[2] - v2->coords[2];
	for( v = vHead->next; v != vHead; v = v->next ) {
		d2[0] = v->coords[0] - v2->coords[0];
		d2[1] = v->coords[1] - v2->coords[1];
		d2[2] = v->coords[2] - v2->coords[2];
		tNorm[0] = d1[1]*d2[2] - d1[2]*d2[1];
		tNorm[1] = d1[2]*d2[0] - d1[0]*d2[2];
		tNorm[2] = d1[0]*d2[1] - d1[1]*d2[0];
		tLen2 = tNorm[0]*tNorm[0] + tNorm[1]*tNorm[1] + tNorm[2]*tNorm[2];
		if( tLen2 > maxLen2 ) {
			maxLen2 = tLen2;
			norm[0] = tNorm[0];
			norm[1] = tNorm[1];
			norm[2] = tNorm[2];
		}
	}

	if( maxLen2 <= 0 ) {
		/* All points lie on a single line -- any decent normal will do */
		norm[0] = norm[1] = norm[2] = 0;
		norm[ShortAxis(d1)] = 1;
	}
}


static void CheckOrientation( TESStesselator *tess )
{
	TESSreal area;
	TESSface *f, *fHead = &tess->mesh->fHead;
	TESSvertex *v, *vHead = &tess->mesh->vHead;
	TESShalfEdge *e;

	/* When we compute the normal automatically, we choose the orientation
	* so that the the sum of the signed areas of all contours is non-negative.
	*/
	area = 0;
	for( f = fHead->next; f != fHead; f = f->next ) {
		e = f->anEdge;
		if( e->winding <= 0 ) continue;
		do {
			area += (e->Org->s - e->Dst->s) * (e->Org->t + e->Dst->t);
			e = e->Lnext;
		} while( e != f->anEdge );
	}
	if( area < 0 ) {
		/* Reverse the orientation by flipping all the t-coordinates */
		for( v = vHead->next; v != vHead; v = v->next ) {
			v->t = - v->t;
		}
		tess->tUnit[0] = - tess->tUnit[0];
		tess->tUnit[1] = - tess->tUnit[1];
		tess->tUnit[2] = - tess->tUnit[2];
	}
}

#ifdef FOR_TRITE_TEST_PROGRAM
#include <stdlib.h>
extern int RandomSweep;
#define S_UNIT_X	(RandomSweep ? (2*drand48()-1) : 1.0)
#define S_UNIT_Y	(RandomSweep ? (2*drand48()-1) : 0.0)
#else
#if defined(SLANTED_SWEEP)
/* The "feature merging" is not intended to be complete.  There are
* special cases where edges are nearly parallel to the sweep line
* which are not implemented.  The algorithm should still behave
* robustly (ie. produce a reasonable tesselation) in the presence
* of such edges, however it may miss features which could have been
* merged.  We could minimize this effect by choosing the sweep line
* direction to be something unusual (ie. not parallel to one of the
* coordinate axes).
*/
#define S_UNIT_X	(TESSreal)0.50941539564955385	/* Pre-normalized */
#define S_UNIT_Y	(TESSreal)0.86052074622010633
#else
#define S_UNIT_X	(TESSreal)1.0
#define S_UNIT_Y	(TESSreal)0.0
#endif
#endif

/* Determine the polygon normal and project vertices onto the plane
* of the polygon.
*/
void tessProjectPolygon( TESStesselator *tess )
{
	TESSvertex *v, *vHead = &tess->mesh->vHead;
	TESSreal norm[3];
	TESSreal *sUnit, *tUnit;
	int i, first, computedNormal = FALSE;

	norm[0] = tess->normal[0];
	norm[1] = tess->normal[1];
	norm[2] = tess->normal[2];
	if( norm[0] == 0 && norm[1] == 0 && norm[2] == 0 ) {
		ComputeNormal( tess, norm );
		computedNormal = TRUE;
	}
	sUnit = tess->sUnit;
	tUnit = tess->tUnit;
	i = LongAxis( norm );

#if defined(FOR_TRITE_TEST_PROGRAM) || defined(TRUE_PROJECT)
	/* Choose the initial sUnit vector to be approximately perpendicular
	* to the normal.
	*/
	Normalize( norm );

	sUnit[i] = 0;
	sUnit[(i+1)%3] = S_UNIT_X;
	sUnit[(i+2)%3] = S_UNIT_Y;

	/* Now make it exactly perpendicular */
	w = Dot( sUnit, norm );
	sUnit[0] -= w * norm[0];
	sUnit[1] -= w * norm[1];
	sUnit[2] -= w * norm[2];
	Normalize( sUnit );

	/* Choose tUnit so that (sUnit,tUnit,norm) form a right-handed frame */
	tUnit[0] = norm[1]*sUnit[2] - norm[2]*sUnit[1];
	tUnit[1] = norm[2]*sUnit[0] - norm[0]*sUnit[2];
	tUnit[2] = norm[0]*sUnit[1] - norm[1]*sUnit[0];
	Normalize( tUnit );
#else
	/* Project perpendicular to a coordinate axis -- better numerically */
	sUnit[i] = 0;
	sUnit[(i+1)%3] = S_UNIT_X;
	sUnit[(i+2)%3] = S_UNIT_Y;

	tUnit[i] = 0;
	tUnit[(i+1)%3] = (norm[i] > 0) ? -S_UNIT_Y : S_UNIT_Y;
	tUnit[(i+2)%3] = (norm[i] > 0) ? S_UNIT_X : -S_UNIT_X;
#endif

	/* Project the vertices onto the sweep plane */
	for( v = vHead->next; v != vHead; v = v->next )
	{
		v->s = Dot( v->coords, sUnit );
		v->t = Dot( v->coords, tUnit );
	}
	if( computedNormal ) {
		CheckOrientation( tess );
	}

	/* Compute ST bounds. */
	first = 1;
	for( v = vHead->next; v != vHead; v = v->next )
	{
		if (first)
		{
			tess->bmin[0] = tess->bmax[0] = v->s;
			tess->bmin[1] = tess->bmax[1] = v->t;
			first = 0;
		}
		else
		{
			if (v->s < tess->bmin[0]) tess->bmin[0] = v->s;
			if (v->s > tess->bmax[0]) tess->bmax[0] = v->s;
			if (v->t < tess->bmin[1]) tess->bmin[1] = v->t;
			if (v->t > tess->bmax[1]) tess->bmax[1] = v->t;
		}
	}
}

#define AddWinding(eDst,eSrc)	(eDst->winding += eSrc->winding, \
	eDst->Sym->winding += eSrc->Sym->winding)

/* tessMeshTessellateMonoRegion( face ) tessellates a monotone region
* (what else would it do??)  The region must consist of a single
* loop of half-edges (see mesh.h) oriented CCW.  "Monotone" in this
* case means that any vertical line intersects the interior of the
* region in a single interval.
*
* Tessellation consists of adding interior edges (actually pairs of
* half-edges), to split the region into non-overlapping triangles.
*
* The basic idea is explained in Preparata and Shamos (which I don''t
* have handy right now), although their implementation is more
* complicated than this one.  The are two edge chains, an upper chain
* and a lower chain.  We process all vertices from both chains in order,
* from right to left.
*
* The algorithm ensures that the following invariant holds after each
* vertex is processed: the untessellated region consists of two
* chains, where one chain (say the upper) is a single edge, and
* the other chain is concave.  The left vertex of the single edge
* is always to the left of all vertices in the concave chain.
*
* Each step consists of adding the rightmost unprocessed vertex to one
* of the two chains, and forming a fan of triangles from the rightmost
* of two chain endpoints.  Determining whether we can add each triangle
* to the fan is a simple orientation test.  By making the fan as large
* as possible, we restore the invariant (check it yourself).
*/
int tessMeshTessellateMonoRegion( TESSmesh *mesh, TESSface *face )
{
	TESShalfEdge *up, *lo;

	/* All edges are oriented CCW around the boundary of the region.
	* First, find the half-edge whose origin vertex is rightmost.
	* Since the sweep goes from left to right, face->anEdge should
	* be close to the edge we want.
	*/
	up = face->anEdge;
	assert( up->Lnext != up && up->Lnext->Lnext != up );

	for( ; VertLeq( up->Dst, up->Org ); up = up->Lprev )
		;
	for( ; VertLeq( up->Org, up->Dst ); up = up->Lnext )
		;
	lo = up->Lprev;

	while( up->Lnext != lo ) {
		if( VertLeq( up->Dst, lo->Org )) {
			/* up->Dst is on the left.  It is safe to form triangles from lo->Org.
			* The EdgeGoesLeft test guarantees progress even when some triangles
			* are CW, given that the upper and lower chains are truly monotone.
			*/
			while( lo->Lnext != up && (EdgeGoesLeft( lo->Lnext )
				|| EdgeSign( lo->Org, lo->Dst, lo->Lnext->Dst ) <= 0 )) {
					TESShalfEdge *tempHalfEdge= tessMeshConnect( mesh, lo->Lnext, lo );
					if (tempHalfEdge == NULL) return 0;
					lo = tempHalfEdge->Sym;
			}
			lo = lo->Lprev;
		} else {
			/* lo->Org is on the left.  We can make CCW triangles from up->Dst. */
			while( lo->Lnext != up && (EdgeGoesRight( up->Lprev )
				|| EdgeSign( up->Dst, up->Org, up->Lprev->Org ) >= 0 )) {
					TESShalfEdge *tempHalfEdge= tessMeshConnect( mesh, up, up->Lprev );
					if (tempHalfEdge == NULL) return 0;
					up = tempHalfEdge->Sym;
			}
			up = up->Lnext;
		}
	}

	/* Now lo->Org == up->Dst == the leftmost vertex.  The remaining region
	* can be tessellated in a fan from this leftmost vertex.
	*/
	assert( lo->Lnext != up );
	while( lo->Lnext->Lnext != up ) {
		TESShalfEdge *tempHalfEdge= tessMeshConnect( mesh, lo->Lnext, lo );
		if (tempHalfEdge == NULL) return 0;
		lo = tempHalfEdge->Sym;
	}

	return 1;
}

/* tessMeshTessellateInterior( mesh ) tessellates each region of
* the mesh which is marked "inside" the polygon.  Each such region
* must be monotone.
*/
int tessMeshTessellateInterior( TESSmesh *mesh )
{
	TESSface *f, *next;

	/*LINTED*/
	for( f = mesh->fHead.next; f != &mesh->fHead; f = next ) {
		/* Make sure we don''t try to tessellate the new triangles. */
		next = f->next;
		if( f->inside ) {
			if ( !tessMeshTessellateMonoRegion( mesh, f ) ) return 0;
		}
	}
	return 1;
}


typedef struct EdgeStackNode EdgeStackNode;
typedef struct EdgeStack EdgeStack;

struct EdgeStackNode {
	TESShalfEdge *edge;
	EdgeStackNode *next;
};

struct EdgeStack {
	EdgeStackNode *top;
	struct BucketAlloc *nodeBucket;
};

int stackInit( EdgeStack *stack, TESSalloc *alloc )
{
	stack->top = NULL;
	stack->nodeBucket = createBucketAlloc( alloc, "CDT nodes", sizeof(EdgeStackNode), 512 );
	return stack->nodeBucket != NULL;
}

void stackDelete( EdgeStack *stack )
{
    deleteBucketAlloc( stack->nodeBucket );
}

int stackEmpty( EdgeStack *stack )
{
	return stack->top == NULL;
}

void stackPush( EdgeStack *stack, TESShalfEdge *e )
{
	EdgeStackNode *node = (EdgeStackNode *)bucketAlloc( stack->nodeBucket );
	if ( ! node ) return;
	node->edge = e;
	node->next = stack->top;
	stack->top = node;
}

TESShalfEdge *stackPop( EdgeStack *stack )
{
	TESShalfEdge *e = NULL;
	EdgeStackNode *node = stack->top;
	if (node) {
		stack->top = node->next;
		e = node->edge;
		bucketFree( stack->nodeBucket, node );
	}
	return e;
}


//	Starting with a valid triangulation, uses the Edge Flip algorithm to
//	refine the triangulation into a Constrained Delaunay Triangulation.
void tessMeshRefineDelaunay( TESSmesh *mesh, TESSalloc *alloc )
{
	// At this point, we have a valid, but not optimal, triangulation.
	// We refine the triangulation using the Edge Flip algorithm
	//
	//  1) Find all internal edges
	//	2) Mark all dual edges
	//	3) insert all dual edges into a queue

	TESSface *f;
	EdgeStack stack;
	TESShalfEdge *e;
	int maxFaces = 0, maxIter = 0, iter = 0;

	stackInit(&stack, alloc);

	for( f = mesh->fHead.next; f != &mesh->fHead; f = f->next ) {
		if ( f->inside) {
			e = f->anEdge;
			do {
				e->mark = EdgeIsInternal(e); // Mark internal edges
				if (e->mark && !e->Sym->mark) stackPush(&stack, e); // Insert into queue
				e = e->Lnext;
			} while (e != f->anEdge);
			maxFaces++;
		}
	}

	// The algorithm should converge on O(n^2), since the predicate is not robust,
	// we'll save guard against infinite loop.
	maxIter = maxFaces * maxFaces;

	// Pop stack until we find a reversed edge
	// Flip the reversed edge, and insert any of the four opposite edges
	// which are internal and not already in the stack (!marked)
	while (!stackEmpty(&stack) && iter < maxIter) {
		e = stackPop(&stack);
		e->mark = e->Sym->mark = 0;
		if (!tesedgeIsLocallyDelaunay(e)) {
			TESShalfEdge *edges[4];
			int i;
			tessMeshFlipEdge(mesh, e);
			// for each opposite edge
			edges[0] = e->Lnext;
			edges[1] = e->Lprev;
			edges[2] = e->Sym->Lnext;
			edges[3] = e->Sym->Lprev;
			for (i = 0; i < 4; i++) {
				if (!edges[i]->mark && EdgeIsInternal(edges[i])) {
					edges[i]->mark = edges[i]->Sym->mark = 1;
					stackPush(&stack, edges[i]);
				}
			}
		}
		iter++;
	}

	stackDelete(&stack);
}


/* tessMeshDiscardExterior( mesh ) zaps (ie. sets to NULL) all faces
* which are not marked "inside" the polygon.  Since further mesh operations
* on NULL faces are not allowed, the main purpose is to clean up the
* mesh so that exterior loops are not represented in the data structure.
*/
void tessMeshDiscardExterior( TESSmesh *mesh )
{
	TESSface *f, *next;

	/*LINTED*/
	for( f = mesh->fHead.next; f != &mesh->fHead; f = next ) {
		/* Since f will be destroyed, save its next pointer. */
		next = f->next;
		if( ! f->inside ) {
			tessMeshZapFace( mesh, f );
		}
	}
}

/* tessMeshSetWindingNumber( mesh, value, keepOnlyBoundary ) resets the
* winding numbers on all edges so that regions marked "inside" the
* polygon have a winding number of "value", and regions outside
* have a winding number of 0.
*
* If keepOnlyBoundary is TRUE, it also deletes all edges which do not
* separate an interior region from an exterior one.
*/
int tessMeshSetWindingNumber( TESSmesh *mesh, int value,
							 int keepOnlyBoundary )
{
	TESShalfEdge *e, *eNext;

	for( e = mesh->eHead.next; e != &mesh->eHead; e = eNext ) {
		eNext = e->next;
		if( e->Rface->inside != e->Lface->inside ) {

			/* This is a boundary edge (one side is interior, one is exterior). */
			e->winding = (e->Lface->inside) ? value : -value;
		} else {

			/* Both regions are interior, or both are exterior. */
			if( ! keepOnlyBoundary ) {
				e->winding = 0;
			} else {
				if ( !tessMeshDelete( mesh, e ) ) return 0;
			}
		}
	}
	return 1;
}

void* heapAlloc( void* userData, unsigned int size )
{
	TESS_NOTUSED( userData );
	return malloc( size );
}

void* heapRealloc( void *userData, void* ptr, unsigned int size )
{
	TESS_NOTUSED( userData );
	return realloc( ptr, size );
}

void heapFree( void* userData, void* ptr )
{
	TESS_NOTUSED( userData );
	free( ptr );
}

static TESSalloc defaulAlloc =
{
	heapAlloc,
	heapRealloc,
	heapFree,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
};

TESStesselator* tessNewTess( TESSalloc* alloc )
{
	TESStesselator* tess;

	if (alloc == NULL)
		alloc = &defaulAlloc;

	/* Only initialize fields which can be changed by the api.  Other fields
	* are initialized where they are used.
	*/

	tess = (TESStesselator *)alloc->memalloc( alloc->userData, sizeof( TESStesselator ));
	if ( tess == NULL ) {
		return 0;          /* out of memory */
	}
	tess->alloc = *alloc;
	/* Check and set defaults. */
	if (tess->alloc.meshEdgeBucketSize == 0)
		tess->alloc.meshEdgeBucketSize = 512;
	if (tess->alloc.meshVertexBucketSize == 0)
		tess->alloc.meshVertexBucketSize = 512;
	if (tess->alloc.meshFaceBucketSize == 0)
		tess->alloc.meshFaceBucketSize = 256;
	if (tess->alloc.dictNodeBucketSize == 0)
		tess->alloc.dictNodeBucketSize = 512;
	if (tess->alloc.regionBucketSize == 0)
		tess->alloc.regionBucketSize = 256;

	tess->normal[0] = 0;
	tess->normal[1] = 0;
	tess->normal[2] = 0;

	tess->bmin[0] = 0;
	tess->bmin[1] = 0;
	tess->bmax[0] = 0;
	tess->bmax[1] = 0;

	tess->reverseContours = 0;
    
	tess->windingRule = TESS_WINDING_ODD;
	tess->processCDT = 0;

	if (tess->alloc.regionBucketSize < 16)
		tess->alloc.regionBucketSize = 16;
	if (tess->alloc.regionBucketSize > 4096)
		tess->alloc.regionBucketSize = 4096;
	tess->regionPool = createBucketAlloc( &tess->alloc, "Regions",
										 sizeof(ActiveRegion), tess->alloc.regionBucketSize );
	if (tess->regionPool == NULL) {
		alloc->memfree( alloc->userData, tess );
		return NULL;
	}

	// Initialize to begin polygon.
	tess->mesh = NULL;

	tess->status = TESS_STATUS_OK;
	tess->vertexIndexCounter = 0;

	tess->dict = NULL;
	tess->pq = NULL;

	tess->contours = NULL;
	tess->contourCount = 0;
	tess->contourCapacity = 0;
	tess->contourPoints = NULL;
	tess->pointCount = 0;
	tess->pointCapacity = 0;
	tess->contourPoints3D = NULL;
	tess->pointCapacity3D = 0;
	tess->errors = NULL;
	tess->errorCount = 0;
	tess->errorCapacity = 0;
	tess->validationFailed = 0;

	tess->vertices = 0;
	tess->vertexIndices = 0;
	tess->vertexCount = 0;
	tess->elements = 0;
	tess->elementCount = 0;

	return tess;
}

void tessDeleteTess( TESStesselator *tess )
{

	struct TESSalloc alloc = tess->alloc;

	deleteBucketAlloc( tess->regionPool );

	if (tess->dict != NULL)
		dictDeleteDict( &alloc, tess->dict );
	if (tess->pq != NULL)
		pqDeletePriorityQ( &alloc, tess->pq );

	if (tess->contourPoints != NULL) {
		alloc.memfree( alloc.userData, tess->contourPoints );
		tess->contourPoints = NULL;
	}
	if (tess->contourPoints3D != NULL) {
		alloc.memfree( alloc.userData, tess->contourPoints3D );
		tess->contourPoints3D = NULL;
	}
	if (tess->contours != NULL) {
		alloc.memfree( alloc.userData, tess->contours );
		tess->contours = NULL;
	}
	if (tess->errors != NULL) {
		alloc.memfree( alloc.userData, tess->errors );
		tess->errors = NULL;
	}

	if( tess->mesh != NULL ) {
		tessMeshDeleteMesh( &alloc, tess->mesh );
		tess->mesh = NULL;
	}
	if (tess->vertices != NULL) {
		alloc.memfree( alloc.userData, tess->vertices );
		tess->vertices = 0;
	}
	if (tess->vertexIndices != NULL) {
		alloc.memfree( alloc.userData, tess->vertexIndices );
		tess->vertexIndices = 0;
	}
	if (tess->elements != NULL) {
		alloc.memfree( alloc.userData, tess->elements );
		tess->elements = 0;
	}

	alloc.memfree( alloc.userData, tess );
}


static TESSindex GetNeighbourFace(TESShalfEdge* edge)
{
	if (!edge->Rface)
		return TESS_UNDEF;
	if (!edge->Rface->inside)
		return TESS_UNDEF;
	return edge->Rface->n;
}

void OutputPolymesh( TESStesselator *tess, TESSmesh *mesh, int elementType, int polySize, int vertexSize )
{
	TESSvertex* v = 0;
	TESSface* f = 0;
	TESShalfEdge* edge = 0;
	int maxFaceCount = 0;
	int maxVertexCount = 0;
	int faceVerts, i;
	TESSindex *elements = 0;
	TESSreal *vert;

	// Assume that the input data is triangles now.
	// Try to merge as many polygons as possible
	if (polySize > 3)
	{
		if (!tessMeshMergeConvexFaces( mesh, polySize ))
		{
			tess->status = TESS_STATUS_OUT_OF_MEMORY;
			return;
		}
	}

	// Mark unused
	for ( v = mesh->vHead.next; v != &mesh->vHead; v = v->next )
		v->n = TESS_UNDEF;

	// Create unique IDs for all vertices and faces.
	for ( f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		f->n = TESS_UNDEF;
		if( !f->inside ) continue;

		edge = f->anEdge;
		faceVerts = 0;
		do
		{
			v = edge->Org;
			if ( v->n == TESS_UNDEF )
			{
				v->n = maxVertexCount;
				maxVertexCount++;
			}
			faceVerts++;
			edge = edge->Lnext;
		}
		while (edge != f->anEdge);

		assert( faceVerts <= polySize );

		f->n = maxFaceCount;
		++maxFaceCount;
	}

	tess->elementCount = maxFaceCount;
	if (elementType == TESS_CONNECTED_POLYGONS)
		maxFaceCount *= 2;
	tess->elements = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData,
													  sizeof(TESSindex) * maxFaceCount * polySize );
	if (!tess->elements)
	{
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	tess->vertexCount = maxVertexCount;
	tess->vertices = (TESSreal*)tess->alloc.memalloc( tess->alloc.userData,
													 sizeof(TESSreal) * tess->vertexCount * vertexSize );
	if (!tess->vertices)
	{
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	tess->vertexIndices = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData,
														    sizeof(TESSindex) * tess->vertexCount );
	if (!tess->vertexIndices)
	{
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	// Output vertices.
	for ( v = mesh->vHead.next; v != &mesh->vHead; v = v->next )
	{
		if ( v->n != TESS_UNDEF )
		{
			// Store coordinate
			vert = &tess->vertices[v->n*vertexSize];
			vert[0] = v->coords[0];
			vert[1] = v->coords[1];
			if ( vertexSize > 2 )
				vert[2] = v->coords[2];
			// Store vertex index.
			tess->vertexIndices[v->n] = v->idx;
		}
	}

	// Output indices.
	elements = tess->elements;
	for ( f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		if ( !f->inside ) continue;

		// Store polygon
		edge = f->anEdge;
		faceVerts = 0;
		do
		{
			v = edge->Org;
			*elements++ = v->n;
			faceVerts++;
			edge = edge->Lnext;
		}
		while (edge != f->anEdge);
		// Fill unused.
		for (i = faceVerts; i < polySize; ++i)
			*elements++ = TESS_UNDEF;

		// Store polygon connectivity
		if ( elementType == TESS_CONNECTED_POLYGONS )
		{
			edge = f->anEdge;
			do
			{
				*elements++ = GetNeighbourFace( edge );
				edge = edge->Lnext;
			}
			while (edge != f->anEdge);
			// Fill unused.
			for (i = faceVerts; i < polySize; ++i)
				*elements++ = TESS_UNDEF;
		}
	}
}

void OutputContours( TESStesselator *tess, TESSmesh *mesh, int vertexSize )
{
	TESSface *f = 0;
	TESShalfEdge *edge = 0;
	TESShalfEdge *start = 0;
	TESSreal *verts = 0;
	TESSindex *elements = 0;
	TESSindex *vertInds = 0;
	int startVert = 0;
	int vertCount = 0;

	tess->vertexCount = 0;
	tess->elementCount = 0;

	for ( f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		if ( !f->inside ) continue;

		start = edge = f->anEdge;
		do
		{
			++tess->vertexCount;
			edge = edge->Lnext;
		}
		while ( edge != start );

		++tess->elementCount;
	}

	tess->elements = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData,
													  sizeof(TESSindex) * tess->elementCount * 2 );
	if (!tess->elements)
	{
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	tess->vertices = (TESSreal*)tess->alloc.memalloc( tess->alloc.userData,
													  sizeof(TESSreal) * tess->vertexCount * vertexSize );
	if (!tess->vertices)
	{
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	tess->vertexIndices = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData,
														    sizeof(TESSindex) * tess->vertexCount );
	if (!tess->vertexIndices)
	{
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	verts = tess->vertices;
	elements = tess->elements;
	vertInds = tess->vertexIndices;

	startVert = 0;

	for ( f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		if ( !f->inside ) continue;

		vertCount = 0;
		start = edge = f->anEdge;
		do
		{
			*verts++ = edge->Org->coords[0];
			*verts++ = edge->Org->coords[1];
			if ( vertexSize > 2 )
				*verts++ = edge->Org->coords[2];
			*vertInds++ = edge->Org->idx;
			++vertCount;
			edge = edge->Lnext;
		}
		while ( edge != start );

		elements[0] = startVert;
		elements[1] = vertCount;
		elements += 2;

		startVert += vertCount;
	}
}

int IsValidCoord(TESSreal coord) {
  return tessIsValidCoord( coord );
}

/* For 3D contours pick the two axes with the largest bounding-box spans;
* dropping the smallest span keeps 2D data (constant z) on its native x/y
* plane. 2D contours are copied through unchanged. */
static void ProjectContour( TESStesselator *tess, int size,
							const unsigned char *src, int stride,
							int numVertices, TESSreal *dst, TESSreal *dst3D )
{
	int drop = 2;
	int i, axis;
	TESSreal minVal[3], maxVal[3], span[3];
	TESS_NOTUSED( tess );

	if( size < 3 ) {
		for( i = 0; i < numVertices; ++i ) {
			const TESSreal *coords = (const TESSreal *)( src + (size_t)i * stride );
			dst[i * 2] = coords[0];
			dst[i * 2 + 1] = coords[1];
			dst3D[i * 3] = coords[0];
			dst3D[i * 3 + 1] = coords[1];
			dst3D[i * 3 + 2] = 0;
		}
		return;
	}

	for( axis = 0; axis < 3; ++axis ) {
		const TESSreal *first = (const TESSreal *)src;
		minVal[axis] = maxVal[axis] = first[axis];
	}
	for( i = 1; i < numVertices; ++i ) {
		const TESSreal *coords = (const TESSreal *)( src + (size_t)i * stride );
		for( axis = 0; axis < 3; ++axis ) {
			if( coords[axis] < minVal[axis] ) minVal[axis] = coords[axis];
			if( coords[axis] > maxVal[axis] ) maxVal[axis] = coords[axis];
		}
	}
	for( axis = 0; axis < 3; ++axis )
		span[axis] = maxVal[axis] - minVal[axis];
	if( span[0] <= span[drop] ) drop = 0;
	if( span[1] <= span[drop] ) drop = 1;

	for( i = 0; i < numVertices; ++i ) {
		const TESSreal *coords = (const TESSreal *)( src + (size_t)i * stride );
		dst3D[i * 3] = coords[0];
		dst3D[i * 3 + 1] = coords[1];
		dst3D[i * 3 + 2] = coords[2];
		int out = 0;
		for( axis = 0; axis < 3; ++axis ) {
			if( axis == drop ) continue;
			dst[i * 2 + out++] = coords[axis];
		}
	}
}

/* Builds the quad-edge mesh from the validated contour cache. Mirrors the
* original per-vertex construction, but only runs for accepted contours. */
static int BuildMeshFromContours( TESStesselator *tess )
{
	int c, i;

	tess->mesh = tessMeshNewMesh( &tess->alloc );
	if( tess->mesh == NULL )
		return 0;

	for( c = 0; c < tess->contourCount; ++c ) {
		const TESSContour2D *contour = &tess->contours[c];
		TESShalfEdge *e = NULL;

		for( i = 0; i < contour->count; ++i ) {
			const TESSreal *coords =
				&tess->contourPoints3D[( contour->pointOffset + i ) * 3];
			if( e == NULL ) {
				e = tessMeshMakeEdge( tess->mesh );
				if( e == NULL || !tessMeshSplice( tess->mesh, e, e->Sym ) )
					return 0;
			} else {
				if( tessMeshSplitEdge( tess->mesh, e ) == NULL )
					return 0;
				e = e->Lnext;
			}

			e->Org->coords[0] = coords[0];
			e->Org->coords[1] = coords[1];
			e->Org->coords[2] = coords[2];
			e->Org->idx = tess->vertexIndexCounter++;

			e->winding = contour->reverse ? -1 : 1;
			e->Sym->winding = contour->reverse ? 1 : -1;
		}
	}
	return 1;
}

void tessAddContour( TESStesselator *tess, int size, const void* vertices,
					int stride, int numVertices )
{
	const unsigned char *src = (const unsigned char*)vertices;
	TESSContour2D contour;

	if ( tess->status == TESS_STATUS_OUT_OF_MEMORY )
		return;

	if ( size < 2 )
		size = 2;
	if ( size > 3 )
		size = 3;

	/* An empty contour keeps the historical no-op behavior (an empty
	* tesselator still returns success with zero elements). */
	if ( numVertices <= 0 )
		return;

	if ( !tessGrowArray( &tess->alloc, (void **)&tess->contourPoints,
						&tess->pointCapacity, sizeof( TESSreal ),
						( tess->pointCount + numVertices ) * 2 ) ||
		!tessGrowArray( &tess->alloc, (void **)&tess->contourPoints3D,
						&tess->pointCapacity3D, sizeof( TESSreal ),
						( tess->pointCount + numVertices ) * 3 ) ||
		!tessGrowArray( &tess->alloc, (void **)&tess->contours,
						&tess->contourCapacity, sizeof( TESSContour2D ),
						tess->contourCount + 1 ) ) {
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return;
	}

	contour.pointOffset = tess->pointCount;
	contour.count = numVertices;
	contour.reverse = tess->reverseContours;
	ProjectContour( tess, size, src, stride, numVertices,
					&tess->contourPoints[tess->pointCount * 2],
					&tess->contourPoints3D[tess->pointCount * 3] );
	tess->pointCount += numVertices;
	tess->contours[tess->contourCount] = contour;

	if( !tessValidateContour( tess->contourPoints, &contour,
							  tess->contourCount,
							  RecordContourError, tess ) ) {
		if( tess->status != TESS_STATUS_OUT_OF_MEMORY )
			tess->status = TESS_STATUS_INVALID_INPUT;
		tess->validationFailed = 1;
	}
	++tess->contourCount;
}

static void FreeInputCache( TESStesselator *tess )
{
	if( tess->contourPoints != NULL ) {
		tess->alloc.memfree( tess->alloc.userData, tess->contourPoints );
		tess->contourPoints = NULL;
	}
	if( tess->contourPoints3D != NULL ) {
		tess->alloc.memfree( tess->alloc.userData, tess->contourPoints3D );
		tess->contourPoints3D = NULL;
	}
	if( tess->contours != NULL ) {
		tess->alloc.memfree( tess->alloc.userData, tess->contours );
		tess->contours = NULL;
	}
	tess->contourCount = 0;
	tess->contourCapacity = 0;
	tess->pointCount = 0;
	tess->pointCapacity = 0;
	tess->pointCapacity3D = 0;
	tess->validationFailed = 0;
}

static void FreeErrorRecords( TESStesselator *tess )
{
	if( tess->errors != NULL ) {
		tess->alloc.memfree( tess->alloc.userData, tess->errors );
		tess->errors = NULL;
	}
	tess->errorCount = 0;
	tess->errorCapacity = 0;
}

void tessSetOption( TESStesselator *tess, int option, int value )
{
	switch(option)
	{
	case TESS_CONSTRAINED_DELAUNAY_TRIANGULATION:
		tess->processCDT = value > 0 ? 1 : 0;
		break;
	case TESS_REVERSE_CONTOURS:
		tess->reverseContours = value > 0 ? 1 : 0;
		break;
	}
}


int tessTesselate( TESStesselator *tess, int windingRule, int elementType,
				  int volatile polySize, int volatile vertexSize, const TESSreal* normal )
{
	TESSmesh *mesh;
	int rc = 1;
	int i, j;

	if (tess->vertices != NULL) {
		tess->alloc.memfree( tess->alloc.userData, tess->vertices );
		tess->vertices = 0;
	}
	if (tess->elements != NULL) {
		tess->alloc.memfree( tess->alloc.userData, tess->elements );
		tess->elements = 0;
	}
	if (tess->vertexIndices != NULL) {
		tess->alloc.memfree( tess->alloc.userData, tess->vertexIndices );
		tess->vertexIndices = 0;
	}

	tess->vertexIndexCounter = 0;

	if (normal)
	{
		tess->normal[0] = normal[0];
		tess->normal[1] = normal[1];
		tess->normal[2] = normal[2];
	}

	tess->windingRule = windingRule;

	if (vertexSize < 2)
		vertexSize = 2;
	if (vertexSize > 3)
		vertexSize = 3;

	/* Nothing was ever added: keep the historical success with no output. */
	if( tess->contourCount == 0 && tess->errorCount == 0 ) {
		tess->status = TESS_STATUS_OK;
		return 1;
	}

	/* Cross-contour interactions only make sense once every contour has been
	* supplied. Rejected contours still participate so their call index stays
	* aligned with tessAddContour() ordering. */
	if( tess->status != TESS_STATUS_OUT_OF_MEMORY ) {
		for( i = 0; i < tess->contourCount; ++i ) {
			for( j = i + 1; j < tess->contourCount; ++j ) {
				if( !tessValidateContourPair(
						tess->contourPoints,
						&tess->contours[i], i,
						&tess->contours[j], j,
						RecordContourError, tess ) ) {
					if( tess->status == TESS_STATUS_OUT_OF_MEMORY )
						return 0;
					tess->status = TESS_STATUS_INVALID_INPUT;
					tess->validationFailed = 1;
				}
			}
		}
	}

	if( tess->status != TESS_STATUS_OK )
		return 0;

	if( !BuildMeshFromContours( tess ) ) {
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return 0;
	}
	if (setjmp(tess->env) != 0) {
		/* come back here if out of memory: release every half-built sweep
		* structure and the mesh, leaving the tesselator reusable. */
		tessCleanupSweep( tess );
		if( tess->mesh != NULL ) {
			tessMeshDeleteMesh( &tess->alloc, tess->mesh );
			tess->mesh = NULL;
		}
		tess->status = TESS_STATUS_OUT_OF_MEMORY;
		return 0;
	}

	/* Determine the polygon normal and project vertices onto the plane
	* of the polygon.
	*/
	tessProjectPolygon( tess );

	/* tessComputeInterior( tess ) computes the planar arrangement specified
	* by the given contours, and further subdivides this arrangement
	* into regions.  Each region is marked "inside" if it belongs
	* to the polygon, according to the rule given by tess->windingRule.
	* Each interior region is guaranteed be monotone.
	*/
	if ( !tessComputeInterior( tess ) ) {
		longjmp(tess->env,1);  /* could've used a label */
	}

	mesh = tess->mesh;

	/* If the user wants only the boundary contours, we throw away all edges
	* except those which separate the interior from the exterior.
	* Otherwise we tessellate all the regions marked "inside".
	*/
	if (elementType == TESS_BOUNDARY_CONTOURS) {
		rc = tessMeshSetWindingNumber( mesh, 1, TRUE );
	} else {
		rc = tessMeshTessellateInterior( mesh );
		if (rc != 0 && tess->processCDT != 0)
			tessMeshRefineDelaunay( mesh, &tess->alloc );
	}
	if (rc == 0) longjmp(tess->env,1);  /* could've used a label */

	tessMeshCheckMesh( mesh );

	if (elementType == TESS_BOUNDARY_CONTOURS) {
		OutputContours( tess, mesh, vertexSize );     /* output contours */
	}
	else
	{
		OutputPolymesh( tess, mesh, elementType, polySize, vertexSize );     /* output polygons */
	}
	if( tess->status != TESS_STATUS_OK )
		longjmp( tess->env, 1 );

	tessMeshDeleteMesh( &tess->alloc, mesh );
	tess->mesh = NULL;

	if( tess->status == TESS_STATUS_OK ) {
		/* The cached inputs were consumed; reset so the object can accept a
		* fresh polygon on the next run. */
		FreeInputCache( tess );
		FreeErrorRecords( tess );
		tess->vertexIndexCounter = 0;
	}

	return tess->status == TESS_STATUS_OK;
}

int tessGetVertexCount( TESStesselator *tess )
{
	return tess->vertexCount;
}

const TESSreal* tessGetVertices( TESStesselator *tess )
{
	return tess->vertices;
}

const TESSindex* tessGetVertexIndices( TESStesselator *tess )
{
	return tess->vertexIndices;
}

int tessGetElementCount( TESStesselator *tess )
{
	return tess->elementCount;
}

const int* tessGetElements( TESStesselator *tess )
{
	return tess->elements;
}

TESSstatus tessGetStatus( TESStesselator *tess )
{
	return tess->status;
}

int tessGetErrorCount( TESStesselator *tess )
{
	return tess->errorCount;
}

const TESSContourErrorRecord* tessGetErrors( TESStesselator *tess )
{
	return tess->errors;
}
