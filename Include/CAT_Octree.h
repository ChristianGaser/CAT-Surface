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


void get_triangle_bounds(polygons_struct *, struct polynode *);
unsigned char xintersect(double [6], double [6]);
unsigned char yintersect(double [6], double [6]);
unsigned char zintersect(double [6], double [6]);
unsigned char intersect(double [6], double [6]);
unsigned char point_in_bounds(Point, double [6]);
struct octree * build_octree(polygons_struct *);
void delete_octree(struct octree *);

#endif
