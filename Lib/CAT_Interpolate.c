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

double
interp_point_unit_sphere(polygons_struct *sphere, double *values, Point pt)
{
        Point             on_sphere_pt, poly_pts[1000];
        double            value;
        double            weights[1000];
        int               i, poly, size, ind;

        if (sphere->bintree == NULL)
                create_polygons_bintree(sphere, ROUND((double) sphere->n_items *
                                                      BINTREE_FACTOR));

        poly = find_closest_polygon_point(&pt, sphere, &on_sphere_pt);
        size = get_polygon_points(sphere, poly, poly_pts);
        get_polygon_interpolation_weights(&on_sphere_pt, size, poly_pts,
                                          weights);
        value = 0.0;
        for (i = 0; i < size; i++) {
                ind = sphere->indices[POINT_INDEX(sphere->end_indices,
                                                  poly, i)];
                value += weights[i] * values[ind];
        }
        
        return value;
}

double
interp_point_sphere(polygons_struct *sphere, double *values, Point pt)
{
        polygons_struct   *unit_sphere;
        double            value;
        int               i;

        unit_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));
        copy_polygons(sphere, unit_sphere);

        /* set radius to 1 */
        for (i = 0; i < unit_sphere->n_points; i++) 
                set_vector_length(&unit_sphere->points[i], 1.0);

        value = interp_point_unit_sphere(unit_sphere, values, pt);
        free(unit_sphere);
        return(value);
}


double
interp_uv_unit_sphere(polygons_struct *sphere, double *values,
                      double u, double v)
{
        Point             pt;

        uv_to_point(u, v, &pt);
        return(interp_point_unit_sphere(sphere, values, pt));
}


double
interp_uv_sphere(polygons_struct *sphere, double *values, double u, double v)
{
        Point pt;

        uv_to_point(u, v, &pt);
        return(interp_point_sphere(sphere, values, pt));
}


double
interp_xyz_sphere(polygons_struct *sphere, double *values, double x, double y,
                  double z)
{
        Point pt;

        fill_Point(pt, x, y, z);
        return(interp_point_unit_sphere(sphere, values, pt));
}


double
interp_xyz_unit_sphere(polygons_struct *sphere, double *values,
                       double x, double y, double z)
{
        Point pt;

        fill_Point(pt, x, y, z);
        return(interp_point_unit_sphere(sphere, values, pt));
}

