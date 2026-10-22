/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>

#include "CAT_Separate.h"
#include "CAT_Smooth.h"

/**
 * \brief Recursively triangulate a polygon using neighbor connectivity.
 *
 * Splits a polygon into triangles by searching for valid internal
 * diagonals based on neighbor connectivity, then recurses on the
 * resulting sub-polygons.
 *
 * \param size         (in)  number of vertices in polygon
 * \param poly         (in)  vertex indices of polygon
 * \param n_neighbours (in)  neighbor counts per vertex
 * \param neighbours   (in)  neighbor index lists
 * \param indices      (out) output triangle indices
 * \return TRUE on success, FALSE on failure
 */
private BOOLEAN recursive_triangulate_one_polygon(
    int size,
    int poly[],
    int n_neighbours[],
    int *neighbours[],
    int indices[])
{
    int *left, *right, n_left, n_right, p1, p2;
    int start_index, end_index, count, i, n;
    BOOLEAN found;

    if (size < 3)
        fprintf(stderr, "recursive_triangulate_one_polygon");

    if (size == 3)
    {
        indices[0] = poly[0];
        indices[1] = poly[1];
        indices[2] = poly[2];
        return (TRUE);
    }

    found = FALSE;

    ALLOC(left, size);
    ALLOC(right, size);

    for_less(start_index, 0, size - 2)
    {
        p1 = poly[start_index];

        for_less(end_index, start_index + 2, size)
        {
            if (start_index == 0 && end_index == size - 1)
                continue;

            p2 = poly[end_index];

            count = 0;
            for_less(n, 0, n_neighbours[p1])
            {
                if (neighbours[p1][n] == p2)
                    ++count;
            }

            if (count != 1)
                continue;

            n_left = 0;
            for_inclusive(i, 0, start_index)
                left[n_left++] = poly[i];
            for_less(i, end_index, size)
                left[n_left++] = poly[i];

            n_right = 0;
            for_inclusive(i, start_index, end_index)
                right[n_right++] = poly[i];

            if (n_left + n_right != size + 2)
            {
                fprintf(stderr, "n_left\n");
            }

            if (recursive_triangulate_one_polygon(n_left, left,
                                                  n_neighbours, neighbours, indices) &&
                recursive_triangulate_one_polygon(n_right, right,
                                                  n_neighbours, neighbours,
                                                  &indices[3 * (n_left - 2)]))
            {
                found = TRUE;
                break;
            }
        }

        if (found)
            break;
    }

    FREE(left);
    FREE(right);

    return (found);
}

/**
 * \brief Triangulate a polygon using recursive splitting.
 *
 * \param size         (in)  number of vertices in polygon
 * \param poly         (in)  vertex indices of polygon
 * \param n_neighbours (in)  neighbor counts per vertex
 * \param neighbours   (in)  neighbor index lists
 * \param indices      (out) output triangle indices
 * \return TRUE on success, FALSE on failure
 */
private BOOLEAN triangulate_one_polygon(
    int size,
    int poly[],
    int n_neighbours[],
    int *neighbours[],
    int indices[])
{
    return (recursive_triangulate_one_polygon(size, poly,
                                              n_neighbours, neighbours, indices));
}

/**
 * \brief Triangulate all polygons in a mesh.
 *
 * Converts polygons with arbitrary vertex counts into triangles while
 * preserving topology and vertex ordering. Result is written into the
 * output triangles mesh.
 *
 * \param polygons  (in)  input mesh with polygons
 * \param triangles (out) triangulated mesh
 */
void triangulate_polygons(
    polygons_struct *polygons,
    polygons_struct *triangles)
{
    int poly, size, index, ind, n_matches, n;
    int *n_neighbours, **neighbours, *indices, max_size;
    progress_struct progress;
    BOOLEAN done;

    create_polygon_point_neighbours(polygons, TRUE,
                                    &n_neighbours, &neighbours,
                                    NULL, NULL);

    *triangles = *polygons;

    triangles->colour_flag = ONE_COLOUR;
    ALLOC(triangles->colours, 1);
    triangles->colours[0] = polygons->colours[0];

    triangles->points = polygons->points;
    triangles->normals = polygons->normals;

    triangles->n_items = 0;
    for_less(poly, 0, polygons->n_items)
        triangles->n_items += GET_OBJECT_SIZE(*polygons, poly) - 2;

    ALLOC(triangles->indices, 3 * triangles->n_items);

    max_size = 0;
    for_less(poly, 0, polygons->n_items)
        max_size = MAX(max_size, GET_OBJECT_SIZE(*polygons, poly));

    ALLOC(indices, max_size);

    initialize_progress_report(&progress, FALSE, polygons->n_items,
                               "Triangulating");

    ind = 0;
    for_less(poly, 0, polygons->n_items)
    {
        size = GET_OBJECT_SIZE(*polygons, poly);
        for_less(index, 0, size)
        {
            indices[index] = polygons->indices[POINT_INDEX(
                polygons->end_indices, poly, index)];
        }

        for_less(index, 2, size - 1)
        {
            n_matches = 0;
            for_less(n, 0, n_neighbours[indices[0]])
            {
                if (neighbours[indices[0]][n] == indices[index])
                    ++n_matches;
            }

            if (n_matches != 1)
                break;
        }

        if (size > 3 && index < size - 1)
        {
            done = triangulate_one_polygon(size, indices,
                                           n_neighbours, neighbours,
                                           &triangles->indices[ind]);
            if (!done)
                print("Could not find good triangulation: %d\n", poly);
            else
                ind += 3 * (size - 2);
        }
        else
            done = FALSE;

        if (!done)
        {
            for_less(index, 1, size - 1)
            {
                triangles->indices[ind] = indices[0];
                ++ind;
                triangles->indices[ind] = indices[index];
                ++ind;
                triangles->indices[ind] = indices[index + 1];
                ++ind;
            }
        }

        update_progress_report(&progress, poly + 1);
    }

    terminate_progress_report(&progress);

    if (ind != 3 * triangles->n_items)
        fprintf(stderr, "Summation of ind\n");

    FREE(polygons->end_indices);
    FREE(polygons->indices);
    FREE(indices);

    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours,
                                    NULL, NULL);

    ALLOC(triangles->end_indices, triangles->n_items);
    for_less(poly, 0, triangles->n_items)
        triangles->end_indices[poly] = 3 * (poly + 1);
}

/**
 * \brief Label connected polygon components in a mesh.
 *
 * Performs a BFS over polygon adjacency to assign component labels and
 * counts the number of polygons in each component.
 *
 * \param polygons        (in)  input mesh
 * \param polygon_classes (out) component label per polygon
 * \param n_in_class      (out) counts per component
 * \return Number of connected components
 */
private int make_connected_components(
    polygons_struct *polygons,
    int polygon_classes[],
    int n_in_class[])
{
    int poly, current_poly, edge, size;
    int neigh;
    int n_components;
    int not_done;
    QUEUE_STRUCT(int)
    queue;

    n_components = 0;

    not_done = 9999;

    for_less(poly, 0, polygons->n_items)
        polygon_classes[poly] = not_done;

    for_less(poly, 0, polygons->n_items)
    {
        if (polygon_classes[poly] != not_done)
            continue;

        if (n_components == 9999)
        {
            ++n_components;
            break;
        }

        INITIALIZE_QUEUE(queue);
        INSERT_IN_QUEUE(queue, poly);
        polygon_classes[poly] = n_components;
        n_in_class[n_components] = 1;

        while (!IS_QUEUE_EMPTY(queue))
        {
            REMOVE_FROM_QUEUE(queue, current_poly);
            size = GET_OBJECT_SIZE(*polygons, current_poly);

            for_less(edge, 0, size)
            {
                neigh = polygons->neighbours[POINT_INDEX(polygons->end_indices, current_poly, edge)];
                if (neigh >= 0 &&
                    polygon_classes[neigh] == not_done)
                {
                    polygon_classes[neigh] = n_components;
                    ++n_in_class[n_components];
                    INSERT_IN_QUEUE(queue, neigh);
                }
            }
        }

        DELETE_QUEUE(queue);

        ++n_components;
    }

    return (n_components);
}

/**
 * \brief Separate a mesh into connected components.
 *
 * Extracts connected components into separate polygon objects. If
 * desired_index is non-negative, only that component is returned.
 * Components are ordered by size (largest first).
 *
 * \param polygons      (in)  input mesh
 * \param desired_index (in)  component index to extract, or -1 for all
 * \param out           (out) array of output objects (allocated)
 * \return Number of output objects
 */
int separate_polygons(
    polygons_struct *polygons,
    int desired_index,
    object_struct **out[])
{
    int point, ind, p_ind, poly, vertex, size, i, j, tmp;
    int point_index, *new_point_ids, n_objects, comp, c;
    int biggest;
    int *poly_classes;
    int n_components, *n_in_class, *ordered;
    polygons_struct *new_poly;

    ALLOC(poly_classes, polygons->n_items);
    ALLOC(n_in_class, 10009);
    ALLOC(ordered, 10000);

    n_components = make_connected_components(polygons, poly_classes,
                                             n_in_class);

    for_less(i, 0, n_components)
        ordered[i] = i;

    for_less(i, 0, n_components - 1)
    {
        biggest = i;
        for_less(j, i + 1, n_components)
        {
            if (n_in_class[ordered[j]] > n_in_class[ordered[biggest]])
                biggest = j;
        }

        tmp = ordered[i];
        ordered[i] = ordered[biggest];
        ordered[biggest] = tmp;
    }

    ALLOC(new_point_ids, polygons->n_points);

    n_objects = 0;

    for_less(c, 0, n_components)
    {
        if (desired_index >= 0 && c != desired_index)
            continue;

        comp = ordered[c];

        for_less(point, 0, polygons->n_points)
            new_point_ids[point] = -1;

        SET_ARRAY_SIZE(*out, n_objects, n_objects + 1,
                       DEFAULT_CHUNK_SIZE);
        (*out)[n_objects] = create_object(POLYGONS);
        new_poly = get_polygons_ptr((*out)[n_objects]);
        ++n_objects;
        initialize_polygons(new_poly, WHITE, NULL);
        if (desired_index >= 0)
        {
            new_poly->points = polygons->points;
            new_poly->normals = polygons->normals;
            new_poly->indices = polygons->indices;
            new_poly->n_items = 0;

            for_less(poly, 0, polygons->n_items)
            {
                if (poly_classes[poly] != comp)
                    continue;
                size = GET_OBJECT_SIZE(*polygons, poly);
                ++new_poly->n_items;
                for_less(vertex, 0, size)
                {
                    point_index = polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)];
                    if (new_point_ids[point_index] < 0)
                        new_point_ids[point_index] = 0;
                }
            }

            ALLOC(new_poly->end_indices, new_poly->n_items);

            ind = 0;
            for_less(point, 0, polygons->n_points)
            {
                if (new_point_ids[point] >= 0)
                {
                    new_point_ids[point] = ind;
                    new_poly->points[ind] = new_poly->points[point];
                    new_poly->normals[ind] = new_poly->normals[point];
                    ++ind;
                }
            }

            new_poly->n_points = ind;

            p_ind = 0;
            ind = 0;
            for_less(poly, 0, polygons->n_items)
            {
                if (poly_classes[poly] != comp)
                    continue;

                size = GET_OBJECT_SIZE(*polygons, poly);

                for_less(vertex, 0, size)
                {
                    point_index = polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)];
                    new_poly->indices[ind] = new_point_ids[point_index];
                    ++ind;
                }

                new_poly->end_indices[p_ind] = ind;
                ++p_ind;
            }
        }
        else
        {
            ind = 0;
            for_less(poly, 0, polygons->n_items)
            {
                if (poly_classes[poly] != comp)
                    continue;

                size = GET_OBJECT_SIZE(*polygons, poly);
                for_less(vertex, 0, size)
                {
                    point_index = polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)];

                    if (new_point_ids[point_index] < 0)
                    {
                        new_point_ids[point_index] = new_poly->n_points;
                        ADD_ELEMENT_TO_ARRAY(new_poly->points,
                                             new_poly->n_points,
                                             polygons->points[point_index],
                                             DEFAULT_CHUNK_SIZE);
                        --new_poly->n_points;
                        ADD_ELEMENT_TO_ARRAY(new_poly->normals,
                                             new_poly->n_points,
                                             polygons->normals[point_index],
                                             DEFAULT_CHUNK_SIZE);
                    }

                    ADD_ELEMENT_TO_ARRAY(new_poly->indices, ind,
                                         new_point_ids[point_index],
                                         DEFAULT_CHUNK_SIZE);
                }

                ADD_ELEMENT_TO_ARRAY(new_poly->end_indices, new_poly->n_items,
                                     ind, DEFAULT_CHUNK_SIZE);
            }

            REALLOC(new_poly->points, new_poly->n_points);
            REALLOC(new_poly->normals, new_poly->n_points);
            REALLOC(new_poly->end_indices, new_poly->n_items);
            REALLOC(new_poly->indices, ind);
        }
    }

    FREE(poly_classes);
    FREE(new_point_ids);
    FREE(n_in_class);
    FREE(ordered);

    return (n_objects);
}

/* size of a cluster and the order in which it was found */
typedef struct
{
    int size;
    int id;
} cluster_size;

/**
 * \brief qsort() order of clusters: largest first, then in the order found.
 *
 * \param a (in) first cluster_size
 * \param b (in) second cluster_size
 * \return negative if a comes first, positive if b does
 */
private int compare_cluster_size(const void *a, const void *b)
{
    const cluster_size *ca = (const cluster_size *)a;
    const cluster_size *cb = (const cluster_size *)b;

    if (ca->size != cb->size)
        return (ca->size < cb->size) ? 1 : -1;
    return (ca->id > cb->id) - (ca->id < cb->id);
}

/**
 * \brief Label the clusters of positive values on a mesh.
 *
 * Each cluster is grown from its lowest vertex by a breadth-first search
 * along the mesh edges, so a cluster is found whole however far it extends.
 * The clusters are then renumbered by size, 1 being the largest; clusters of
 * equal size keep the order of their lowest vertex index.
 *
 * \param polygons (in)  mesh
 * \param values   (in)  double[n_points]; value per vertex
 * \param labels   (out) int[n_points]; cluster number per vertex, 0 where the
 *                       value is <= 0 or NaN
 * \return number of clusters
 */
int label_value_clusters(
    polygons_struct *polygons,
    const double values[],
    int labels[])
{
    int point, current, n, neigh, n_clusters;
    int *n_neighbours, **neighbours, *rank;
    cluster_size *clusters;
    QUEUE_STRUCT(int) queue;

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    /* -1 marks a vertex of a cluster that has not been grown yet; NaN fails
       the comparison and stays outside like every value <= 0 */
    for_less(point, 0, polygons->n_points)
        labels[point] = (values[point] > 0.0) ? -1 : 0;

    ALLOC(clusters, polygons->n_points);
    n_clusters = 0;

    for_less(point, 0, polygons->n_points)
    {
        if (labels[point] != -1)
            continue;

        clusters[n_clusters].id = n_clusters + 1;
        clusters[n_clusters].size = 0;
        ++n_clusters;

        INITIALIZE_QUEUE(queue);
        INSERT_IN_QUEUE(queue, point);
        labels[point] = n_clusters;

        while (!IS_QUEUE_EMPTY(queue))
        {
            REMOVE_FROM_QUEUE(queue, current);
            ++clusters[n_clusters - 1].size;

            for_less(n, 0, n_neighbours[current])
            {
                neigh = neighbours[current][n];
                if (labels[neigh] == -1)
                {
                    labels[neigh] = n_clusters;
                    INSERT_IN_QUEUE(queue, neigh);
                }
            }
        }

        DELETE_QUEUE(queue);
    }

    delete_all_polygon_point_neighbours(n_neighbours, neighbours);

    /* renumber the clusters by size */
    if (n_clusters > 1)
        qsort(clusters, n_clusters, sizeof(cluster_size), compare_cluster_size);
    ALLOC(rank, n_clusters + 1);
    for_less(n, 0, n_clusters)
        rank[clusters[n].id] = n + 1;
    for_less(point, 0, polygons->n_points)
    {
        if (labels[point] > 0)
            labels[point] = rank[labels[point]];
    }

    FREE(rank);
    FREE(clusters);

    return (n_clusters);
}

/**
 * \brief Whether every vertex of a polygon carries the given label.
 *
 * \param polygons (in) mesh
 * \param poly     (in) polygon index
 * \param labels   (in) int[n_points]; label per vertex
 * \param label    (in) label to test for
 * \return TRUE if the polygon lies entirely in the label
 */
private BOOLEAN polygon_has_label(
    polygons_struct *polygons,
    int poly,
    const int labels[],
    int label)
{
    int vertex, size;

    size = GET_OBJECT_SIZE(*polygons, poly);
    for_less(vertex, 0, size)
    {
        if (labels[polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)]] != label)
            return (FALSE);
    }
    return (TRUE);
}

/**
 * \brief Extract the polygons of one label as a new mesh.
 *
 * A polygon belongs to the label if all its vertices do, so polygons that
 * straddle the border of a cluster are left out, and a cluster without a
 * complete polygon gives an empty mesh.  The vertices keep their order.
 *
 * \param polygons (in)  mesh
 * \param labels   (in)  int[n_points]; label per vertex
 * \param label    (in)  label to extract
 * \param out      (out) the polygons whose vertices all carry the label, with
 *                       the vertices renumbered; free with delete_polygons()
 * \return number of polygons in out
 */
int extract_label_polygons(
    polygons_struct *polygons,
    const int labels[],
    int label,
    polygons_struct *out)
{
    int poly, size, vertex, point, n_indices, ind, p_ind;
    int *new_point_ids;

    initialize_polygons(out, WHITE, NULL);

    /* -2 marks the vertices of the polygons that are kept */
    ALLOC(new_point_ids, polygons->n_points);
    for_less(point, 0, polygons->n_points)
        new_point_ids[point] = -1;

    n_indices = 0;
    for_less(poly, 0, polygons->n_items)
    {
        if (!polygon_has_label(polygons, poly, labels, label))
            continue;

        size = GET_OBJECT_SIZE(*polygons, poly);
        ++out->n_items;
        n_indices += size;
        for_less(vertex, 0, size)
            new_point_ids[polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)]] = -2;
    }

    if (out->n_items > 0)
    {
        for_less(point, 0, polygons->n_points)
        {
            if (new_point_ids[point] == -2)
                new_point_ids[point] = out->n_points++;
        }

        ALLOC(out->points, out->n_points);
        ALLOC(out->normals, out->n_points);
        for_less(point, 0, polygons->n_points)
        {
            if (new_point_ids[point] >= 0)
                out->points[new_point_ids[point]] = polygons->points[point];
        }

        ALLOC(out->end_indices, out->n_items);
        ALLOC(out->indices, n_indices);
        ind = 0;
        p_ind = 0;
        for_less(poly, 0, polygons->n_items)
        {
            if (!polygon_has_label(polygons, poly, labels, label))
                continue;

            size = GET_OBJECT_SIZE(*polygons, poly);
            for_less(vertex, 0, size)
                out->indices[ind++] = new_point_ids[polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)]];
            out->end_indices[p_ind++] = ind;
        }

        compute_polygon_normals(out);
    }

    FREE(new_point_ids);

    return (out->n_items);
}
