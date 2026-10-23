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
 * \brief Copy a subset of the polygons of a mesh into a new mesh.
 *
 * The cost is proportional to the subset, not to the mesh, so that a mesh
 * can be split into many small parts in linear time.
 *
 * \param polygons      (in)     mesh
 * \param poly_list     (in)     int[n_polys]; polygons to copy, in this order
 * \param n_polys       (in)     number of polygons
 * \param point_list    (in)     int[n_list]; ascending vertex indices that
 *                               include every vertex of the polygons
 * \param n_list        (in)     number of those vertices
 * \param new_point_ids (in/out) int[n_points]; all -1, and so again on return
 * \param out           (out)    mesh set up by initialize_polygons()
 */
private void copy_polygon_subset(
    polygons_struct *polygons,
    const int poly_list[],
    int n_polys,
    const int point_list[],
    int n_list,
    int new_point_ids[],
    polygons_struct *out)
{
    int i, poly, size, vertex, point, n_indices, ind;

    if (n_polys == 0)
        return;

    /* -2 marks the vertices the polygons use; they are numbered in order */
    n_indices = 0;
    for_less(i, 0, n_polys)
    {
        poly = poly_list[i];
        size = GET_OBJECT_SIZE(*polygons, poly);
        n_indices += size;
        for_less(vertex, 0, size)
            new_point_ids[polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)]] = -2;
    }
    for_less(i, 0, n_list)
    {
        if (new_point_ids[point_list[i]] == -2)
            new_point_ids[point_list[i]] = out->n_points++;
    }

    ALLOC(out->points, out->n_points);
    ALLOC(out->normals, out->n_points);
    for_less(i, 0, n_list)
    {
        point = point_list[i];
        if (new_point_ids[point] >= 0)
            out->points[new_point_ids[point]] = polygons->points[point];
    }

    out->n_items = n_polys;
    ALLOC(out->end_indices, n_polys);
    ALLOC(out->indices, n_indices);
    ind = 0;
    for_less(i, 0, n_polys)
    {
        poly = poly_list[i];
        size = GET_OBJECT_SIZE(*polygons, poly);
        for_less(vertex, 0, size)
            out->indices[ind++] = new_point_ids[polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)]];
        out->end_indices[i] = ind;
    }

    compute_polygon_normals(out);

    for_less(i, 0, n_list)
        new_point_ids[point_list[i]] = -1;
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
    int poly, point, n_polys, n_list;
    int *poly_list, *point_list, *new_point_ids;

    initialize_polygons(out, WHITE, NULL);

    ALLOC(poly_list, polygons->n_items);
    n_polys = 0;
    for_less(poly, 0, polygons->n_items)
    {
        if (polygon_has_label(polygons, poly, labels, label))
            poly_list[n_polys++] = poly;
    }

    ALLOC(point_list, polygons->n_points);
    ALLOC(new_point_ids, polygons->n_points);
    n_list = 0;
    for_less(point, 0, polygons->n_points)
    {
        new_point_ids[point] = -1;
        if (labels[point] == label)
            point_list[n_list++] = point;
    }

    copy_polygon_subset(polygons, poly_list, n_polys, point_list, n_list,
                        new_point_ids, out);

    FREE(poly_list);
    FREE(point_list);
    FREE(new_point_ids);

    return (out->n_items);
}

/**
 * \brief A vertex label, or 0 if it lies outside 1..n_labels.
 *
 * \param label    (in) label of a vertex
 * \param n_labels (in) highest valid label
 * \return label, or 0 if it is out of range
 */
private int label_in_range(int label, int n_labels)
{
    return (label > 0 && label <= n_labels) ? label : 0;
}

/**
 * \brief Split a labelled mesh into one mesh per label.
 *
 * Does for every label 1..n_labels what extract_label_polygons() does for
 * one, but sorts the polygons and vertices by label in a single pass, so the
 * cost is linear in the size of the mesh however many labels there are.
 *
 * \param polygons (in)  mesh
 * \param labels   (in)  int[n_points]; label per vertex, 0 or 1..n_labels
 * \param n_labels (in)  highest label
 * \param out      (out) object_struct*[n_labels], allocated; (*out)[k - 1]
 *                       holds the polygons of label k and is empty if no
 *                       polygon lies entirely in it; free with
 *                       delete_object_list()
 * \return number of labels with at least one polygon
 */
int split_label_polygons(
    polygons_struct *polygons,
    const int labels[],
    int n_labels,
    object_struct **out[])
{
    int poly, point, k, label, n_nonempty;
    int *poly_label, *poly_start, *poly_list, *point_start, *point_list;
    int *new_point_ids, *fill;

    /* label of every polygon that lies entirely in one label, else 0; labels
       outside 1..n_labels count as 0 so that they cannot index past the
       arrays below */
    ALLOC(poly_label, polygons->n_items);
    for_less(poly, 0, polygons->n_items)
    {
        label = labels[polygons->indices[POINT_INDEX(polygons->end_indices, poly, 0)]];
        poly_label[poly] = (label > 0 && label <= n_labels &&
                            polygon_has_label(polygons, poly, labels, label)) ? label : 0;
    }

    /* counting sort of the polygons and the vertices by label; both keep
       their order within a label */
    ALLOC(poly_start, n_labels + 2);
    ALLOC(point_start, n_labels + 2);
    for_less(k, 0, n_labels + 2)
    {
        poly_start[k] = 0;
        point_start[k] = 0;
    }
    for_less(poly, 0, polygons->n_items)
        poly_start[poly_label[poly] + 1]++;
    for_less(point, 0, polygons->n_points)
        point_start[label_in_range(labels[point], n_labels) + 1]++;
    for_less(k, 1, n_labels + 2)
    {
        poly_start[k] += poly_start[k - 1];
        point_start[k] += point_start[k - 1];
    }

    ALLOC(fill, n_labels + 1);
    ALLOC(poly_list, polygons->n_items);
    for_less(k, 0, n_labels + 1)
        fill[k] = poly_start[k];
    for_less(poly, 0, polygons->n_items)
        poly_list[fill[poly_label[poly]]++] = poly;
    ALLOC(point_list, polygons->n_points);
    for_less(k, 0, n_labels + 1)
        fill[k] = point_start[k];
    for_less(point, 0, polygons->n_points)
        point_list[fill[label_in_range(labels[point], n_labels)]++] = point;

    ALLOC(new_point_ids, polygons->n_points);
    for_less(point, 0, polygons->n_points)
        new_point_ids[point] = -1;

    ALLOC(*out, n_labels);
    n_nonempty = 0;
    for_less(k, 1, n_labels + 1)
    {
        polygons_struct *part;

        (*out)[k - 1] = create_object(POLYGONS);
        part = get_polygons_ptr((*out)[k - 1]);
        initialize_polygons(part, WHITE, NULL);
        copy_polygon_subset(polygons, &poly_list[poly_start[k]],
                            poly_start[k + 1] - poly_start[k],
                            &point_list[point_start[k]],
                            point_start[k + 1] - point_start[k],
                            new_point_ids, part);
        if (part->n_items > 0)
            ++n_nonempty;
    }

    FREE(poly_label);
    FREE(poly_start);
    FREE(point_start);
    FREE(fill);
    FREE(poly_list);
    FREE(point_list);
    FREE(new_point_ids);

    return (n_nonempty);
}
