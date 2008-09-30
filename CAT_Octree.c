/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 */

/* Routines for constructing an octree for finding the nearest polygon
 * to a point.
 */

#include <stdio.h>
#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Surf.h"

//#define LEVEL 4
//#define NNODES 512
//#define LEVEL 5
//#define NNODES 4096 /* pow(8, level - 1) */
//#define LEVEL 6
//#define NNODES 8*8*8*8*8
#define LEVEL 8
#define NNODES 8*8*8*8*8*8*8

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


void
get_polynode_bounds(polygons_struct *polygons, struct polynode *pn,
                    double bounds[6])
{
        bounds[0] = bounds[1] = polygons->points[pn->pts[0]].coords[0];
        bounds[2] = bounds[3] = polygons->points[pn->pts[0]].coords[1];
        bounds[4] = bounds[5] = polygons->points[pn->pts[0]].coords[2];

        /* bounds[0] and bounds[1] */
        if (bounds[0] > polygons->points[pn->pts[1]].coords[0])
                bounds[0] = polygons->points[pn->pts[1]].coords[0];
        if (bounds[0] > polygons->points[pn->pts[2]].coords[0])
                bounds[0] = polygons->points[pn->pts[2]].coords[0];
        if (bounds[1] < polygons->points[pn->pts[1]].coords[0])
                bounds[1] = polygons->points[pn->pts[1]].coords[0];
        if (bounds[1] < polygons->points[pn->pts[2]].coords[0])
                bounds[1] = polygons->points[pn->pts[2]].coords[0];

        /* bounds[2] and bounds[3] */
        if (bounds[2] > polygons->points[pn->pts[1]].coords[1])
                bounds[2] = polygons->points[pn->pts[1]].coords[1];
        if (bounds[2] > polygons->points[pn->pts[2]].coords[1])
                bounds[2] = polygons->points[pn->pts[2]].coords[1];
        if (bounds[3] < polygons->points[pn->pts[1]].coords[1])
                bounds[3] = polygons->points[pn->pts[1]].coords[1];
        if (bounds[3] < polygons->points[pn->pts[2]].coords[1])
                bounds[3] = polygons->points[pn->pts[2]].coords[1];

        /* bounds[4] and bounds[5] */
        if (bounds[4] > polygons->points[pn->pts[1]].coords[2])
                bounds[4] = polygons->points[pn->pts[1]].coords[2];
        if (bounds[4] > polygons->points[pn->pts[2]].coords[2])
                bounds[4] = polygons->points[pn->pts[2]].coords[2];
        if (bounds[5] < polygons->points[pn->pts[1]].coords[2])
                bounds[5] = polygons->points[pn->pts[1]].coords[2];
        if (bounds[5] < polygons->points[pn->pts[2]].coords[2])
                bounds[5] = polygons->points[pn->pts[2]].coords[2];
}

/* x-, y- and z-intersect tests */
unsigned char
xintersect(double bounds1[6], double bounds2[6])
{
        if (bounds1[1] < bounds2[0])
                return 0;
        if (bounds1[0] > bounds2[1])
                return 0;
        return 1;
}

unsigned char
yintersect(double bounds1[6], double bounds2[6])
{
        if (bounds1[3] < bounds2[2])
                return 0;
        if (bounds1[2] > bounds2[3])
                return 0;
        return 1;
}

unsigned char
zintersect(double bounds1[6], double bounds2[6])
{
        if (bounds1[5] < bounds2[4])
                return 0;
        if (bounds1[4] > bounds2[5])
                return 0;
        return 1;
}

void
insert_triangle(struct octree *tree, struct polynode *node, double bounds[6])
{
        int n, i;
        unsigned char insertedflag = 0;
        struct polynode *tmp, *cur;
        int xinc, yinc;

        yinc = pow(2, LEVEL - 1);
        xinc = yinc*yinc;

        for (n = 0; n < NNODES; n++) {
                if (xintersect(bounds, tree->bounds[n]) == 0) {
                        n += xinc - 1;
                        continue;
                }
                if (yintersect(bounds, tree->bounds[n]) == 0) {
                        n += yinc - 1;
                        continue;
                }
                if (zintersect(bounds, tree->bounds[n]) == 1) {
                        if (insertedflag == 0) {
                                tmp = node;
                                insertedflag = 1;
                        } else { /* make a duplicate */
                                tmp = (struct polynode *)
                                      malloc(sizeof(struct polynode));
                                tmp->num = node->num;
                                tmp->next = NULL;
                                for (i = 0; i < 3; i++)
                                        tmp->pts[i] = node->pts[i];
                        }
                        tmp->next = tree->nodes[n];
                        tree->nodes[n] = tmp;
                }
        }
}

struct octree *
build_octree(polygons_struct *polygons)
{
        int t, i, size, nnodes;
        struct polynode *cur;
        struct octree *tree;
        double tbounds[6];
        double xinc, yinc, zinc;
        double xcur, ycur, zcur, xprev, yprev, zprev;
        progress_struct progress;

        /* initialize the tree */
        tree = (struct octree *) malloc(sizeof(struct octree));
        tree->polyflag = (int *) malloc(sizeof(int) * polygons->n_items);
        tree->npoly = polygons->n_items;

        /* set the bounds on the nodes */
        get_bounds(polygons, tree->bbox);
        xinc = (tree->bbox[1] - tree->bbox[0]) / pow(2, LEVEL - 1);
        yinc = (tree->bbox[3] - tree->bbox[2]) / pow(2, LEVEL - 1);
        zinc = (tree->bbox[5] - tree->bbox[4]) / pow(2, LEVEL - 1);
        xprev = tree->bbox[0]; yprev = tree->bbox[2]; zprev = tree->bbox[4];
        xcur = xprev + xinc; ycur = yprev + yinc; zcur = zprev + zinc;

        for (i = 0; i < NNODES; i++) {
                tree->bounds[i][0] = xprev;
                tree->bounds[i][1] = xcur;
                tree->bounds[i][2] = yprev;
                tree->bounds[i][3] = ycur;
                tree->bounds[i][4] = zprev;
                tree->bounds[i][5] = zcur;

                if (zcur >= tree->bbox[5]) {
                     zprev = tree->bbox[4]; zcur = zprev + zinc;
                     yprev = ycur; ycur = yprev + yinc;
                     if (yprev >= tree->bbox[3]) {
                             yprev = tree->bbox[2]; ycur = yprev + yinc;
                             xprev = xcur; xcur = xprev + xinc;
                     } 
                } else {
                     zprev = zcur; zcur = zprev + zinc;
                }
        }

        initialize_progress_report(&progress, FALSE, polygons->n_items,
                                   "InitOctree");

        /* insert the triangles into the tree */
        for (t = 0; t < polygons->n_items; t++) {
                cur = (struct polynode *) malloc(sizeof(struct polynode));
                cur->num = t;
                cur->next = NULL;

                size = GET_OBJECT_SIZE(*polygons, t);
                if (size != 3) {
                        printf("Mesh must only contain triangles. Exiting..\n");
                        return(NULL);
                }

                cur->pts[0] = polygons->indices[
                                      POINT_INDEX(polygons->end_indices, t, 0)];
                cur->pts[1] = polygons->indices[
                                      POINT_INDEX(polygons->end_indices, t, 1)];
                cur->pts[2] = polygons->indices[
                                      POINT_INDEX(polygons->end_indices, t, 2)];
                get_polynode_bounds(polygons, cur, tbounds);
                insert_triangle(tree, cur, tbounds);
                update_progress_report(&progress, t);
        }
        terminate_progress_report(&progress);
        return tree;
}

/* helper function for deleting the octree */
void
recursive_node_delete(struct polynode *n)
{
        if(n->next != NULL)
                recursive_node_delete(n->next);
        free(n);
}

/* delete the octree and all nodes */
void
delete_octree(struct octree *tree) {
        int n;

        for (n = 0; n < NNODES; n++) {
                if (tree->nodes[n] != NULL)
                        recursive_node_delete(tree->nodes[n]);
        }

        free(tree->polyflag);
        free(tree);
}

void
hausdorff_distance(Point p, polygons_struct *polygons, struct octree *tree,
                   double *val)
{
        int n;
        double bounds[6], hbounds[6];
        Point tp[3], closest;
        struct polynode *cur;
        double dist, min_dist = PINF;
        int xinc, yinc;

        yinc = pow(2, LEVEL - 1);
        xinc = yinc*yinc;

        for (n = 0; n < tree->npoly; n++)
                tree->polyflag[n] = 0;

        bounds[0] = bounds[1] = p.coords[0];
        bounds[2] = bounds[3] = p.coords[1];
        bounds[4] = bounds[5] = p.coords[2];

        /* make sure that the point bbox is not outside the total bbox */
        if (bounds[0] > tree->bbox[1])
                bounds[0] = tree->bbox[1];
        if (bounds[1] < tree->bbox[0])
                bounds[1] = tree->bbox[0];
        if (bounds[2] > tree->bbox[3])
                bounds[2] = tree->bbox[3];
        if (bounds[3] < tree->bbox[2])
                bounds[3] = tree->bbox[2];
        if (bounds[4] > tree->bbox[5])
                bounds[4] = tree->bbox[5];
        if (bounds[5] < tree->bbox[4])
                bounds[5] = tree->bbox[4];

        for (n = 0; n < NNODES; n++) {
                tree->nodeflag[n] = 0;
                if (xintersect(bounds, tree->bounds[n]) == 0) {
                        n += xinc - 1;
                        continue;
                }
                if (yintersect(bounds, tree->bounds[n]) == 0) {
                        n += yinc - 1;
                        continue;
                }
                if (zintersect(bounds, tree->bounds[n]) == 1) {
                        tree->nodeflag[n] = 1;
                        cur = tree->nodes[n];
                        while (cur != NULL) {
                                tree->polyflag[cur->num] = 1;
                                tp[0] = polygons->points[cur->pts[0]];
                                tp[1] = polygons->points[cur->pts[1]];
                                tp[2] = polygons->points[cur->pts[2]];
                                dist = find_point_polygon_distance_sq(&p, 3, tp,
                                                                      &closest);
                                if (dist < min_dist)
                                        min_dist = dist;
                                cur = cur->next;
                        }
                }
        }

        /* update the bounding box based on hausdorff distance and check if
         * the nearest point is in a neighboring node */
        dist = min_dist;
        bounds[0] -= dist; bounds[1] += dist;
        bounds[2] -= dist; bounds[3] += dist;
        bounds[4] -= dist; bounds[5] += dist;
        for (n = 0; n < NNODES; n++) {
                if (xintersect(bounds, tree->bounds[n]) == 0) {
                        n += xinc - 1;
                        continue;
                }
                if (yintersect(bounds, tree->bounds[n]) == 0) {
                        n += yinc - 1;
                        continue;
                }
                if (tree->nodeflag[n] == 0 &&
                    zintersect(bounds, tree->bounds[n]) == 1) {
                        cur = tree->nodes[n];
                        while (cur != NULL) {
                                if (tree->polyflag[cur->num] == 0) {
                                        tree->polyflag[cur->num] = 1;
                                        tp[0] = polygons->points[cur->pts[0]];
                                        tp[1] = polygons->points[cur->pts[1]];
                                        tp[2] = polygons->points[cur->pts[2]];
                                        dist = find_point_polygon_distance_sq(
                                               &p, 3, tp, &closest);
                                        if (dist < min_dist) {
                                                min_dist = dist;
                                        }
                                }
                                cur = cur->next;
                        }
                }
        }
        *val = sqrt(min_dist);
}

