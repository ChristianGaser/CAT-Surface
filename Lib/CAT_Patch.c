/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Patch.h"

/* internal helper routine */
/**
 * \brief Allocate a patch list node for a polygon.
 *
 * Copies the polygon vertex indices into a patchinfo node.
 *
 * \param polygons (in)  source mesh
 * \param poly     (in)  polygon index
 * \return Allocated patchinfo node
 */
struct patchinfo *
newnode(polygons_struct *polygons, int poly)
{
    struct patchinfo *node;
    int size, i, p;

    node = (struct patchinfo *)malloc(sizeof(struct patchinfo));
    if (!node)
    {
        fprintf(stderr, "Memory allocation error in newnode().\n");
        exit(EXIT_FAILURE);
    }
    node->num = poly;
    node->next = NULL;

    size = GET_OBJECT_SIZE(*polygons, poly);
    for (i = 0; i < size; i++)
    {
        p = polygons->indices[POINT_INDEX(polygons->end_indices,
                                          poly, i)];
        node->pts[i] = p;
    }
    return node;
}

/* internal routine for getting the polygons around a point */
/**
 * \brief Collect polygons incident to a vertex.
 *
 * Scans polygons and records those containing the specified point.
 *
 * \param polygons  (in)  source mesh
 * \param ptidx     (in)  point index
 * \param polyidx   (out) polygon index list
 * \param max_polys (in)  capacity of polyidx array
 * \return Number of polygons found
 */
int get_polys_around_point(polygons_struct *polygons, int ptidx, int polyidx[],
                           int max_polys)
{
    int size, poly, i, p, num_polys = 0;

    for (poly = 0; poly < polygons->n_items; poly++)
    {
        size = GET_OBJECT_SIZE(*polygons, poly);
        for (i = 0; i < size; i++)
        {
            p = polygons->indices[POINT_INDEX(polygons->end_indices,
                                              poly, i)];
            if (p == ptidx)
            {
                polyidx[num_polys] = poly;
                num_polys++;
                if (num_polys == max_polys)
                    return num_polys; /* no more room */
                break;
            }
        }
    }
    return num_polys;
}

/* internal iterative routine for getting the polygon info */
/**
 * \brief Recursively gather neighboring polygons up to a depth.
 *
 * Performs a breadth-first expansion from a seed polygon, adding
 * neighbors that share a vertex until the requested level is reached.
 *
 * \param polygons (in)  source mesh
 * \param poly     (in)  starting polygon index
 * \param polyflag (in/out) visited flags per polygon
 * \param head     (in/out) current list tail
 * \param level    (in)  remaining expansion depth
 * \return Updated list tail
 */
struct patchinfo *
get_neighbor_polys(polygons_struct *polygons, int poly, int *polyflag,
                   struct patchinfo *head, int level)
{
    int size, i, n, p, n_polys;
    int polyidx[200];
    struct patchinfo *cur;

    if (level == 0)
        return head;

    cur = head;

    size = GET_OBJECT_SIZE(*polygons, poly);
    for (i = 0; i < size; i++)
    {
        p = polygons->indices[POINT_INDEX(polygons->end_indices,
                                          poly, i)];

        n_polys = get_polys_around_point(polygons, p, polyidx, 200);

        for (n = 0; n < n_polys; n++)
        {
            if (polyflag[polyidx[n]] == 0)
            {
                polyflag[polyidx[n]] = 1;
                cur->next = newnode(polygons, polyidx[n]);
                cur = cur->next;

                cur = get_neighbor_polys(polygons, polyidx[n],
                                         polyflag, cur,
                                         level - 1);
            }
        }
    }
    return cur;
}

/* internal routine for making a patch from a linked list */
/**
 * \brief Build a mesh patch from a list of polygons.
 *
 * Creates a new polygon object containing only the polygons in the
 * patch list, renumbering points and computing normals.
 *
 * \param polygons (in)  source mesh
 * \param head     (in)  list of patch polygons
 * \return Allocated object list containing the patch mesh
 */
object_struct **
make_patch(polygons_struct *polygons, struct patchinfo *head)
{
    struct patchinfo *cur;
    int i, p;
    int *ptmap;
    object_struct **object;
    polygons_struct *patch;

    ptmap = (int *)malloc(sizeof(int) * polygons->n_points);
    if (!ptmap)
    {
        fprintf(stderr, "Memory allocation error in make_patch(): ptmap.\n");
        exit(EXIT_FAILURE);
    }

    /* create the output object */
    object = (object_struct **)malloc(sizeof(object_struct *));
    if (!object)
    {
        fprintf(stderr, "Memory allocation error in make_patch(): object.\n");
        free(ptmap);
        exit(EXIT_FAILURE);
    }
    *object = create_object(POLYGONS);
    patch = get_polygons_ptr(*object);
    initialize_polygons(patch, WHITE, NULL);

    Surfprop_a(patch->surfprop) = 0.3;
    Surfprop_d(patch->surfprop) = 0.6;
    Surfprop_s(patch->surfprop) = 0.6;
    Surfprop_se(patch->surfprop) = 60;
    Surfprop_t(patch->surfprop) = 1.0;

    patch->colour_flag = 0;
    patch->line_thickness = 1.0f;

    /* find the number of points & polygons */
    for (cur = head; cur != NULL; cur = cur->next)
    {
        patch->n_items++;
        for (i = 0; i < 3; i++)
        {
            if (ptmap[cur->pts[i]] == 0)
                patch->n_points++;
            ptmap[cur->pts[i]] = 1;
        }
    }

    patch->points = (Point *)malloc(sizeof(Point) * patch->n_points);
    if (!patch->points)
    {
        fprintf(stderr, "Memory allocation error in make_patch(): points.\n");
        free(ptmap);
        free(object);
        exit(EXIT_FAILURE);
    }

    p = 1;
    for (i = 0; i < polygons->n_points; i++)
    {
        if (ptmap[i] > 0)
        {
            ptmap[i] = p; /* renumber the points */
            patch->points[p - 1] = polygons->points[i];
            p++;
        }
    }

    /* renumber the points in our patch */
    for (cur = head; cur != NULL; cur = cur->next)
    {
        for (i = 0; i < 3; i++)
            cur->pts[i] = ptmap[cur->pts[i]] - 1;
    }

    /* make the indices and end_indices */
    patch->indices = (int *)malloc(sizeof(int) * 3 * patch->n_items);
    patch->end_indices = (int *)malloc(sizeof(int) * patch->n_items);
    patch->normals = (Vector *)malloc(sizeof(Vector) * patch->n_points);
    if (!patch->indices || !patch->end_indices || !patch->normals)
    {
        fprintf(stderr, "Memory allocation error in make_patch(): indices/end_indices/normals.\n");
        free(ptmap);
        free(object);
        free(patch->points);
        if (patch->indices)
            free(patch->indices);
        if (patch->end_indices)
            free(patch->end_indices);
        if (patch->normals)
            free(patch->normals);
        exit(EXIT_FAILURE);
    }

    cur = head;
    for (i = 0; i < patch->n_items; i++)
    {
        patch->end_indices[i] = 3 * (i + 1);
        patch->indices[(i * 3)] = cur->pts[0];
        patch->indices[(i * 3) + 1] = cur->pts[1];
        patch->indices[(i * 3) + 2] = cur->pts[2];
        cur = cur->next;
    }

    compute_polygon_normals(patch);

    free(ptmap);

    return (object);
}

/* a routine for extracting a patch around a triangle in the mesh */
/**
 * \brief Extract a polygon patch around a seed triangle.
 *
 * Gathers neighboring polygons up to the specified level and returns
 * a mesh patch containing the collected polygons.
 *
 * \param polygons (in)  source mesh
 * \param poly     (in)  seed polygon index
 * \param level    (in)  neighborhood depth
 * \return Allocated object list containing the patch mesh
 */
object_struct **
extract_patch_around_polygon(polygons_struct *polygons, int poly, int level)
{
    struct patchinfo *head, *cur;
    int *polyflag;
    object_struct **object;
    polygons_struct *patch;

    polyflag = (int *)malloc(sizeof(int) * polygons->n_items);

    /* start our list of polygons in the neighborhood */
    head = newnode(polygons, poly);
    polyflag[poly] = 1;

    get_neighbor_polys(polygons, poly, polyflag, head, level);

    object = make_patch(polygons, head);

    cur = head;
    while (head != NULL)
    {
        head = cur->next;
        free(cur);
        cur = head;
    }

    free(polyflag);
    return (object);
}

/* a routine for extracting a patch around a point in the mesh */
/**
 * \brief Extract a polygon patch around a seed vertex.
 *
 * Collects polygons incident to the vertex and expands neighbors up to
 * the specified level, returning a patch mesh.
 *
 * \param polygons (in)  source mesh
 * \param point    (in)  seed vertex index
 * \param level    (in)  neighborhood depth
 * \return Allocated object list containing the patch mesh
 */
object_struct **
extract_patch_around_point(polygons_struct *polygons, int point, int level)
{
    struct patchinfo *head, *cur;
    int i, n_polys, polyidx[200];
    int *polyflag, *ptmap;
    object_struct **object;
    polygons_struct *patch;

    polyflag = (int *)malloc(sizeof(int) * polygons->n_items);

    n_polys = get_polys_around_point(polygons, point, polyidx, 200);

    if (n_polys == 0)
    {
        free(polyflag);
        object = (object_struct **)malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        return (object);
    }

    head = newnode(polygons, polyidx[0]);
    polyflag[point] = 1;
    cur = head;

    for (i = 1; i < n_polys; i++)
        cur = get_neighbor_polys(polygons, polyidx[i], polyflag,
                                 cur, level);

    object = make_patch(polygons, head);

    cur = head;
    while (head != NULL)
    {
        head = cur->next;
        free(cur);
        cur = head;
    }

    free(polyflag);

    return (object);
}

/* a routine for extracting a patch from a list of polygons... polys must
 * have size polygons->n_items.
 */
/**
 * \brief Extract a patch from a polygon selection mask.
 *
 * Builds a patch containing polygons flagged in the polys array. If num
 * is non-zero, only polygons with matching value are included.
 *
 * \param polygons (in)  source mesh
 * \param polys    (in)  polygon selection mask
 * \param num      (in)  selection value (0 to include all non-zero)
 * \return Allocated object list containing the patch mesh
 */
object_struct **
extract_patch_polys(polygons_struct *polygons, int *polys, int num)
{
    int poly;
    struct patchinfo *head = NULL, *cur;
    object_struct **object;

    for (poly = 0; poly < polygons->n_items; poly++)
    {
        if (polys[poly] == 0 || (num != 0 && polys[poly] != num))
            continue; /* skip */

        if (head == NULL)
        { /* initialize list */
            head = newnode(polygons, poly);
            cur = head;
        }
        else
        {
            cur->next = newnode(polygons, poly);
            cur = cur->next;
        }
    }

    if (head == NULL)
    {
        object = (object_struct **)malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        return (object);
    }

    object = make_patch(polygons, head);

    cur = head;
    while (head != NULL)
    {
        head = cur->next;
        free(cur);
        cur = head;
    }

    return (object);
}

/* a routine for extracting a patch from a list of points... points must
 * have size polygons->n_points.
 */
/**
 * \brief Extract a patch from a point selection mask.
 *
 * Selects polygons containing any selected points and returns a patch.
 * If num is non-zero, only points with matching value are considered.
 *
 * \param polygons (in)  source mesh
 * \param points   (in)  point selection mask
 * \param num      (in)  selection value (0 to include all non-zero)
 * \return Allocated object list containing the patch mesh
 */
object_struct **
extract_patch_points(polygons_struct *polygons, int *points, int num)
{
    int i, p, poly, size, n_polys;
    int *polys;
    object_struct **object;

    /* get the associated polygons */
    polys = (int *)malloc(sizeof(int) * polygons->n_items);
    memset(polys, 0, sizeof(int) * polygons->n_items);
    n_polys = 0;

    for (poly = 0; poly < polygons->n_items; poly++)
    {
        size = GET_OBJECT_SIZE(*polygons, poly);
        for (i = 0; i < size; i++)
        {
            p = polygons->indices[POINT_INDEX(polygons->end_indices,
                                              poly, i)];
            if (points[p] == 0 || (num != 0 && points[p] != num))
                continue; /* skip */

            polys[poly] = 1; /* not in set, delete */
            n_polys++;
            break;
        }
    }

    if (n_polys == 0)
    {
        object = (object_struct **)malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
    }
    else
    {
        object = extract_patch_polys(polygons, polys, 0);
    }
    free(polys);
    return (object);
}
