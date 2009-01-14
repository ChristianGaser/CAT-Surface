/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 * $Id$
 *
*/

//#define LEVEL 4
//#define NNODES 512
#define LEVEL 5
#define NNODES 4096 /* pow(8, level - 1) */
//#define LEVEL 6
//#define NNODES 8*8*8*8*8
//#define LEVEL 8
//#define NNODES 8*8*8*8*8*8*8
#define YINC 2*2*2*2 /* pow(2, LEVEL - 1) */
#define XINC YINC*YINC

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

struct polynode {
        int num;           /* number of the triangle in the original mesh */
        int pts[3];        /* points in the triangle */

        struct polynode *next; /* the next one in a linked list */
};

struct octree {
        double bbox[6]; /* the entire bounding box for the brain */
        int nodeflag[NNODES]; /* flag for keeping track of what's been done */
        int *polyflag; /* flag for keeping track of what's been done */
        int npoly;

        struct polynode *nodes[NNODES]; /* the triangles in each node */
        double bounds[NNODES][6]; /* bounding boxes for nodes */
};


struct octree * build_octree(polygons_struct *);
void delete_octree(struct octree *);
void hausdorff_distance(Point, polygons_struct *, struct octree *, double *);
