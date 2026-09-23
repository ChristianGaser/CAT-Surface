/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Interpolate.h"
#include "CAT_Map.h"

/**
 * \brief Interpolate scalar value at a 3D point on a unit sphere mesh.
 *
 * Uses barycentric interpolation within the closest polygon to compute
 * a scalar value at an arbitrary 3D coordinate. The sphere is assumed
 * to have unit radius. Builds a spatial index (bintree) on first call.
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param pt        (in)  3D point at which to interpolate (should be near sphere)
 * \return               Interpolated scalar value using barycentric weighting
 */
double
interp_point_unit_sphere(polygons_struct *sphere, double *values, Point pt)
{
    Point on_sphere_pt, poly_pts[1000];
    double value;
    double weights[1000];
    int i, poly, size, ind;

    if (sphere->bintree == NULL)
        create_polygons_bintree(sphere, ROUND((double)sphere->n_items *
                                              BINTREE_FACTOR));

    poly = find_closest_polygon_point(&pt, sphere, &on_sphere_pt);
    size = get_polygon_points(sphere, poly, poly_pts);
    get_polygon_interpolation_weights(&on_sphere_pt, size, poly_pts,
                                      weights);
    value = 0.0;
    for (i = 0; i < size; i++)
    {
        ind = sphere->indices[POINT_INDEX(sphere->end_indices,
                                          poly, i)];
        value += (double)weights[i] * values[ind];
    }

    return value;
}
