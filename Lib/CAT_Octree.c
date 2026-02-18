/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 * $Id$
 *
 */

/* Routines for constructing an octree for finding the nearest polygon
 * to a point.
 */

#include <stdio.h>
#include <math.h>
#include <bicpl.h>

#include "CAT_SurfUtils.h"
#include "CAT_Octree.h"

/**
 * \brief Compute axis-aligned bounds for a triangle.
 *
 * Populates the bounds array of the polynode with min/max coordinates
 * of the triangle vertices.
 *
 * \param polygons (in)  source polygon mesh
 * \param node     (in/out) triangle node with vertex indices set
 * \return void
 */
void get_triangle_bounds(polygons_struct *polygons, struct polynode *node)
{
    node->bounds[0] = polygons->points[node->pts[0]].coords[0];
    node->bounds[1] = polygons->points[node->pts[0]].coords[0];
    node->bounds[2] = polygons->points[node->pts[0]].coords[1];
    node->bounds[3] = polygons->points[node->pts[0]].coords[1];
    node->bounds[4] = polygons->points[node->pts[0]].coords[2];
    node->bounds[5] = polygons->points[node->pts[0]].coords[2];

    /* bounds[0] and bounds[1] */
    if (node->bounds[0] > polygons->points[node->pts[1]].coords[0])
        node->bounds[0] = polygons->points[node->pts[1]].coords[0];
    if (node->bounds[0] > polygons->points[node->pts[2]].coords[0])
        node->bounds[0] = polygons->points[node->pts[2]].coords[0];
    if (node->bounds[1] < polygons->points[node->pts[1]].coords[0])
        node->bounds[1] = polygons->points[node->pts[1]].coords[0];
    if (node->bounds[1] < polygons->points[node->pts[2]].coords[0])
        node->bounds[1] = polygons->points[node->pts[2]].coords[0];

    /* bounds[2] and bounds[3] */
    if (node->bounds[2] > polygons->points[node->pts[1]].coords[1])
        node->bounds[2] = polygons->points[node->pts[1]].coords[1];
    if (node->bounds[2] > polygons->points[node->pts[2]].coords[1])
        node->bounds[2] = polygons->points[node->pts[2]].coords[1];
    if (node->bounds[3] < polygons->points[node->pts[1]].coords[1])
        node->bounds[3] = polygons->points[node->pts[1]].coords[1];
    if (node->bounds[3] < polygons->points[node->pts[2]].coords[1])
        node->bounds[3] = polygons->points[node->pts[2]].coords[1];

    /* bounds[4] and bounds[5] */
    if (node->bounds[4] > polygons->points[node->pts[1]].coords[2])
        node->bounds[4] = polygons->points[node->pts[1]].coords[2];
    if (node->bounds[4] > polygons->points[node->pts[2]].coords[2])
        node->bounds[4] = polygons->points[node->pts[2]].coords[2];
    if (node->bounds[5] < polygons->points[node->pts[1]].coords[2])
        node->bounds[5] = polygons->points[node->pts[1]].coords[2];
    if (node->bounds[5] < polygons->points[node->pts[2]].coords[2])
        node->bounds[5] = polygons->points[node->pts[2]].coords[2];
}

/* x-, y- and z-intersect tests */
/**
 * \brief Test x-interval overlap between two bounds.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if x intervals overlap, 0 otherwise
 */
unsigned char
xintersect(double bounds[6], double bounds2[6])
{
    if (bounds[1] < bounds2[0] || bounds[0] > bounds2[1])
        return 0;
    return 1;
}

/**
 * \brief Test y-interval overlap between two bounds.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if y intervals overlap, 0 otherwise
 */
unsigned char
yintersect(double bounds[6], double bounds2[6])
{
    if (bounds[3] < bounds2[2] || bounds[2] > bounds2[3])
        return 0;
    return 1;
}

/**
 * \brief Test z-interval overlap between two bounds.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if z intervals overlap, 0 otherwise
 */
unsigned char
zintersect(double bounds[6], double bounds2[6])
{
    if (bounds[5] < bounds2[4] || bounds[4] > bounds2[5])
        return 0;
    return 1;
}

/**
 * \brief Test full 3D bounds overlap between two boxes.
 *
 * \param bounds  (in) first bounds array
 * \param bounds2 (in) second bounds array
 * \return 1 if all axes overlap, 0 otherwise
 */
unsigned char
intersect(double bounds[6], double bounds2[6])
{
    if (bounds[1] < bounds2[0] || bounds[0] > bounds2[1] ||
        bounds[3] < bounds2[2] || bounds[2] > bounds2[3] ||
        bounds[5] < bounds2[4] || bounds[4] > bounds2[5])
        return 0;
    return 1;
}

/**
 * \brief Check if a point lies within an axis-aligned bounding box.
 *
 * \param pt     (in) point to test
 * \param bounds (in) bounds array
 * \return 1 if inside, 0 otherwise
 */
unsigned char
point_in_bounds(Point pt, double bounds[6])
{
    if (Point_x(pt) >= bounds[0] && Point_x(pt) <= bounds[1] &&
        Point_y(pt) >= bounds[2] && Point_y(pt) <= bounds[3] &&
        Point_z(pt) >= bounds[4] && Point_z(pt) <= bounds[5])
        return 1;
    return 1;
}

/**
 * \brief Insert a triangle node into all intersecting octree boxes.
 *
 * Duplicates nodes as needed when a triangle overlaps multiple boxes.
 *
 * \param tree (in/out) octree
 * \param node (in) triangle node to insert
 * \return void
 */
void insert_triangle(struct octree *tree, struct polynode *node)
{
    int n, i;
    unsigned char insertedflag = 0;
    struct polynode *tmp, *cur;

    for (n = 0; n < NBOXES; n++)
    {
        if (xintersect(node->bounds, tree->bounds[n]) == 0)
        {
            n += XINC - 1;
            continue;
        }
        if (yintersect(node->bounds, tree->bounds[n]) == 0)
        {
            n += YINC - 1;
            continue;
        }
        if (zintersect(node->bounds, tree->bounds[n]) == 1)
        {
            if (insertedflag == 0)
            {
                tmp = node;
                insertedflag = 1;
            }
            else
            { /* make a duplicate */
                tmp = (struct polynode *)
                    malloc(sizeof(struct polynode));
                tmp->num = node->num;
                for (i = 0; i < 3; i++)
                    tmp->pts[i] = node->pts[i];
                for (i = 0; i < 6; i++)
                    tmp->bounds[i] = node->bounds[i];
                tmp->next = NULL;
            }
            tmp->next = tree->nodes[n];
            tree->nodes[n] = tmp;
        }
    }
}

/**
 * \brief Build an octree for fast triangle lookup.
 *
 * Computes the global bounds, subdivides into fixed boxes, and inserts
 * each triangle into the appropriate nodes based on bounds overlap.
 *
 * \param polygons (in) input triangular mesh
 * \return Allocated octree or NULL on error
 */
struct octree *
build_octree(polygons_struct *polygons)
{
    int t, i, size, nnodes;
    struct polynode *cur;
    struct octree *tree;
    double xinc, yinc, zinc;
    double xcur, ycur, zcur, xprev, yprev, zprev;
    progress_struct progress;

    /* initialize the tree */
    tree = (struct octree *)malloc(sizeof(struct octree));
    tree->polyflag = (int *)malloc(sizeof(int) * polygons->n_items);
    tree->npoly = polygons->n_items;
    tree->nodelist = (struct polynode **)
        malloc(sizeof(struct polynode *) * polygons->n_items);

    /* set the bounds on the boxes */
    get_bounds(polygons, tree->bbox);
    xinc = (tree->bbox[1] - tree->bbox[0]) / pow(2, LEVEL - 1);
    yinc = (tree->bbox[3] - tree->bbox[2]) / pow(2, LEVEL - 1);
    zinc = (tree->bbox[5] - tree->bbox[4]) / pow(2, LEVEL - 1);
    xprev = tree->bbox[0];
    yprev = tree->bbox[2];
    zprev = tree->bbox[4];
    xcur = xprev + xinc;
    ycur = yprev + yinc;
    zcur = zprev + zinc;

    for (i = 0; i < NBOXES; i++)
    {
        tree->bounds[i][0] = xprev;
        tree->bounds[i][1] = xcur;
        tree->bounds[i][2] = yprev;
        tree->bounds[i][3] = ycur;
        tree->bounds[i][4] = zprev;
        tree->bounds[i][5] = zcur;

        if (zcur >= tree->bbox[5])
        {
            zprev = tree->bbox[4];
            zcur = zprev + zinc;
            yprev = ycur;
            ycur = yprev + yinc;
            if (yprev >= tree->bbox[3])
            {
                yprev = tree->bbox[2];
                ycur = yprev + yinc;
                xprev = xcur;
                xcur = xprev + xinc;
            }
        }
        else
        {
            zprev = zcur;
            zcur = zprev + zinc;
        }
    }

    initialize_progress_report(&progress, FALSE, polygons->n_items,
                               "InitOctree");

    /* insert the triangles into the tree */
    for (t = 0; t < polygons->n_items; t++)
    {
        cur = (struct polynode *)malloc(sizeof(struct polynode));
        cur->num = t;
        cur->next = NULL;

        size = GET_OBJECT_SIZE(*polygons, t);
        if (size != 3)
        {
            printf("Mesh must only contain triangles. Exiting..\n");
            return (NULL);
        }

        cur->pts[0] = polygons->indices[POINT_INDEX(polygons->end_indices, t, 0)];
        cur->pts[1] = polygons->indices[POINT_INDEX(polygons->end_indices, t, 1)];
        cur->pts[2] = polygons->indices[POINT_INDEX(polygons->end_indices, t, 2)];
        get_triangle_bounds(polygons, cur);
        tree->nodelist[t] = cur;
        insert_triangle(tree, cur);
        update_progress_report(&progress, t);
    }
    terminate_progress_report(&progress);
    return tree;
}

/* helper function for deleting the octree */
/**
 * \brief Recursively free a linked list of polynode entries.
 *
 * \param n (in) head of list to free
 * \return void
 */
void recursive_node_delete(struct polynode *n)
{
    if (n->next != NULL)
        recursive_node_delete(n->next);
    free(n);
}

/* delete the octree and all nodes */
/**
 * \brief Free an octree and all associated nodes.
 *
 * \param tree (in/out) octree to delete
 * \return void
 */
void delete_octree(struct octree *tree)
{
    int n;

    for (n = 0; n < NBOXES; n++)
    {
        if (tree->nodes[n] != NULL)
            recursive_node_delete(tree->nodes[n]);
    }

    free(tree->polyflag);
    free(tree->nodelist);
    free(tree);
}
