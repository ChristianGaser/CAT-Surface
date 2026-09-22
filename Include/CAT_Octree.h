/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 * $Id$
 *
*/

#ifndef _CAT_OCTREE_H_
#define _CAT_OCTREE_H_

#include <bicpl.h>

#define LEVEL 5
#define NBOXES 4096 /* pow(8, LEVEL - 1) */
#define YINC 2*2*2*2 /* pow(2, LEVEL - 1) */
#define XINC YINC*YINC

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

struct polynode {
        int num;           /* number of the triangle in the original mesh */
        int pts[3];        /* points in the triangle */

        double bounds[6];  /* triangle bounds */

        struct polynode *next; /* the next one in a linked list */
};

struct octree {
        double bbox[6]; /* the entire bounding box for the brain */
        int nodeflag[NBOXES]; /* flag for keeping track of what's been done */
        int *polyflag; /* flag for keeping track of what's been done */
        int npoly;

        struct polynode **nodelist; /* raw list of triangles */
        struct polynode *nodes[NBOXES]; /* the triangles in each box */
        double bounds[NBOXES][6]; /* bounds for boxes */
};


/**
 * \brief Compute axis-aligned bounds for a triangle.
 *
 * \param polygons (in)  source polygon mesh
 * \param node     (in/out) triangle node with vertex indices set
 */
void get_triangle_bounds(polygons_struct *polygons, struct polynode *node);
/**
 * \brief Test x-interval overlap between two bounds.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if x intervals overlap, 0 otherwise
 */
unsigned char xintersect(double bounds[6], double bounds2[6]);
/**
 * \brief Test y-interval overlap between two bounds.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if y intervals overlap, 0 otherwise
 */
unsigned char yintersect(double bounds[6], double bounds2[6]);
/**
 * \brief Test z-interval overlap between two bounds.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if z intervals overlap, 0 otherwise
 */
unsigned char zintersect(double bounds[6], double bounds2[6]);
/**
 * \brief Test full 3D bounds overlap between two boxes.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if all axes overlap, 0 otherwise
 */
unsigned char intersect(double bounds[6], double bounds2[6]);
/**
 * \brief Check if a point lies within an axis-aligned bounding box.
 *
 * \param pt     (in) point to test
 * \param bounds (in) bounds array
 * \return 1 if inside, 0 otherwise
 */
unsigned char point_in_bounds(Point pt, double bounds[6]);
/**
 * \brief Build an octree for fast triangle lookup.
 *
 * \param polygons (in) input triangular mesh
 * \return Allocated octree or NULL on error
 */
struct octree * build_octree(polygons_struct *polygons);
/**
 * \brief Free an octree and all associated nodes.
 *
 * \param tree (in/out) octree to delete
 */
void delete_octree(struct octree *tree);

#endif
