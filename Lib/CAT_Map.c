/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_SurfUtils.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"

#define BINTREE_FACTOR 0.5
#define NEW_COORDINATE_SYSTEM 1

/* to print out fsavg.index2D_256x128.txt */
// #define DEBUG 1

static double
compute_clockwise_rotation2(double x, double y)
{
    double radians;

    if (x == 0.0)
    {
        if (y < 0.0)
            return (PI / 2.0);
        else if (y > 0.0)
            return (3.0 * PI / 2.0);
        else
            return (0.0);
    }
    else if (y == 0.0)
    {
        if (x > 0.0)
            return (0.0);
        else
            return (PI);
    }
    else
    {
        radians = -(double)atan2(y, x);

        if (radians < 0.0)
            radians += 2.0 * PI;

        return (radians);
    }
}

/**
 * \brief Convert 3D Cartesian point on unit sphere to 2D (u,v) spherical coordinates.
 *
 * Maps a 3D point on the unit sphere to 2D texture coordinates where u,v ∈ [0,1].
 * u corresponds to azimuthal angle (longitude), v to polar angle (latitude).
 * Handles both standard and alternative coordinate systems (controlled by NEW_COORDINATE_SYSTEM macro).
 *
 * \param point (in)  3D point on unit sphere (x² + y² + z² = 1)
 * \param u     (out) horizontal coordinate [0, 1] (azimuth/longitude)
 * \param v     (out) vertical coordinate [0, 1] (elevation/latitude)
 */
void point_to_uv(Point *point, double *u, double *v)
{
    double x, y, z, phi, theta;

    x = (double)Point_x(*point);
    y = (double)Point_y(*point);
    z = (double)Point_z(*point);

    theta = acos(z);
    if (NEW_COORDINATE_SYSTEM)
    {
        phi = compute_clockwise_rotation(y, x) + PI / 2.0;
        *u = phi / (PI * 2.0);
        if (*u > 1.0)
            *u -= 1.0;
    }
    else
    {
        phi = compute_clockwise_rotation2(y, x);
        *u = phi / (PI * 2.0);
    }
    *v = theta / PI;
}

/**
 * \brief Convert 2D (u,v) spherical coordinates to 3D point on unit sphere.
 *
 * Inverse operation of point_to_uv: maps 2D texture coordinates (u,v) ∈ [0,1]²
 * to a 3D point on the unit sphere. u controls azimuthal angle, v controls polar angle.
 * Consistent with point_to_uv regarding coordinate system selection.
 *
 * \param u     (in)  horizontal coordinate [0, 1] (azimuth/longitude)
 * \param v     (in)  vertical coordinate [0, 1] (elevation/latitude)
 * \param point (out) 3D point on unit sphere with x² + y² + z² = 1
 */
void uv_to_point(double u, double v, Point *point)
{
    double x, y, z, phi, theta;

    if (NEW_COORDINATE_SYSTEM)
    {
        phi = u * PI * 2.0;
    }
    else
    {
        /* shift phi by 90 degree to obtain correct position of midline */
        phi = u * PI * 2.0 + PI / 2.0;
    }

    theta = v * PI;

    z = cos(theta);
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);

    fill_Point(*point, x, y, z);
}

/**
 * \brief Map per-vertex values of a surface onto a 2D (u,v) sheet via its sphere.
 *
 * Every sheet pixel is converted to a point on the unit sphere, the closest
 * triangle of the (normalized) sphere is found and the values of its vertices
 * are interpolated with barycentric weights. The result is rescaled linearly
 * to [0, 1]. Used to build the flat maps that DARTEL surface registration works
 * on.
 *
 * \param polygons      (in)     surface the values belong to
 * \param sphere        (in)     its spherical mapping; NULL creates a
 *                               tetrahedral unit sphere with the same number of
 *                               triangles
 * \param sphere_values (in/out) per-vertex values; NULL computes smoothed
 *                               curvatures of type curvtype, otherwise they are
 *                               smoothed in-place when fwhm > 0
 * \param mapped_data   (out)    dm[0]*dm[1] sheet values in [0, 1]
 * \param fwhm          (in)     heat-kernel FWHM in mm applied first (0: none)
 * \param dm            (in)     sheet dimensions {nu, nv}
 * \param curvtype      (in)     curvature type used when sphere_values is NULL
 */
void map_sphere_values_to_sheet(polygons_struct *polygons,
                                polygons_struct *sphere, double *sphere_values,
                                double *mapped_data, double fwhm, int *dm,
                                int curvtype)
{
    polygons_struct unit_sphere;
    Point unit_point, on_sphere_point, centre;
    Point poly_points[1000];
    double value;
    double u, v;
    double mn, mx;
    double weights[1000];
    int i;
    int x, y;
    int poly, size, ind;

    /* if sphere_values is empty calculate curvature */
    if (sphere_values == (double *)0)
    {
        sphere_values = (double *)malloc(sizeof(double) * polygons->n_points);
        if (!sphere_values)
        {
            fprintf(stderr, "Memory allocation error in map_sphere_values_to_sheet().\n");
            exit(EXIT_FAILURE);
        }
        get_smoothed_curvatures(polygons, sphere_values, fwhm,
                                curvtype);
    }
    else if (fwhm > 0)
    {
        smooth_heatkernel(polygons, sphere_values, fwhm);
    }

    if (sphere == NULL)
    {
        /* create unit sphere w/ same # of triangles as skin surface */
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                                  polygons->n_items, &unit_sphere);
    }
    else
    {
        copy_polygons(sphere, &unit_sphere);
        /* set radius to 1 */
        for (i = 0; i < unit_sphere.n_points; i++)
            set_vector_length(&unit_sphere.points[i], 1.0);
    }

    create_polygons_bintree(&unit_sphere,
                            ROUND((double)unit_sphere.n_items *
                                  BINTREE_FACTOR));

#ifdef DEBUG
    for (y = 0; y < dm[1]; y++)
    {
        for (x = 0; x < dm[0]; x++)
        {
#else
    for (x = 0; x < dm[0]; x++)
    {
        for (y = 0; y < dm[1]; y++)
        {
#endif
            u = ((double)x + 0.5) / (double)(dm[0]);
            v = ((double)y + 0.5) / (double)(dm[1]);
            uv_to_point(u, v, &unit_point);

            poly = find_closest_polygon_point(&unit_point,
                                              &unit_sphere,
                                              &on_sphere_point);

            size = get_polygon_points(&unit_sphere, poly,
                                      poly_points);

            get_polygon_interpolation_weights(&on_sphere_point,
                                              size, poly_points,
                                              weights);

#ifdef DEBUG
            fprintf(stdout, "%d\n", unit_sphere.indices[POINT_INDEX(unit_sphere.end_indices, poly, 0)] + 1);
#endif
            value = 0.0;
            for (i = 0; i < size; i++)
            {
                ind = unit_sphere.indices[POINT_INDEX(unit_sphere.end_indices,
                                                      poly, i)];
                value += (double)weights[i] * sphere_values[ind];
            }
            mapped_data[x + (dm[0] * y)] = value;
        }
    }

    /* scale mapped_data to uint8 range */
    mn = FLT_MAX;
    mx = -FLT_MAX;
    for (i = 0; i < dm[0] * dm[1]; i++)
    {
        if (mapped_data[i] > mx)
            mx = mapped_data[i];
        if (mapped_data[i] < mn)
            mn = mapped_data[i];
    }

    for (i = 0; i < dm[0] * dm[1]; i++)
        mapped_data[i] = (mapped_data[i] - mn) / (mx - mn);

    delete_the_bintree(&unit_sphere.bintree);
    delete_polygons(&unit_sphere);
}

/**
 * \brief Sample a 2D (u,v) sheet at the vertices of a sphere.
 *
 * Inverse of map_sphere_values_to_sheet(): the vertices of a tetrahedral unit
 * sphere with as many triangles as polygons are converted to (u,v) and the
 * sheet is read there, bilinearly interpolated or at the nearest pixel.
 *
 * \param sheet2d     (in)  dm[0]*dm[1] sheet values
 * \param values      (out) polygons->n_points values, allocated by the caller
 * \param polygons    (in)  mesh whose triangle count defines the sphere
 * \param interpolate (in)  non-zero for bilinear interpolation, zero for the
 *                          nearest pixel
 * \param dm          (in)  sheet dimensions {nu, nv}
 */
void map_sheet2d_to_sphere(double *sheet2d, double *values,
                           polygons_struct *polygons, int interpolate, int *dm)
{
    double tmp_x, tmp_y;
    double u, v;
    Point centre;
    polygons_struct unit_sphere;
    int i, x, y;
    double xp, yp, xm, ym;
    double H00, H01, H10, H11;

    /* create a unit sphere with same number of triangles as skin surface */
    fill_Point(centre, 0.0, 0.0, 0.0);

    create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                              polygons->n_items, &unit_sphere);

    create_polygons_bintree(&unit_sphere,
                            ROUND((double)unit_sphere.n_items *
                                  BINTREE_FACTOR));

    tmp_x = (double)dm[0];
    tmp_y = (double)dm[1];

    for (i = 0; i < polygons->n_points; i++)
    {
        point_to_uv(&unit_sphere.points[i], &u, &v);

        x = (int)(u * tmp_x - 0.5);
        y = (int)(v * tmp_y - 0.5);

        if (interpolate)
        {
            xp = u * tmp_x - 0.5 - x;
            yp = v * tmp_y - 0.5 - y;
            xm = 1.0 - xp;
            ym = 1.0 - yp;
            H00 = sheet2d[bound(x, y, dm)];
            H01 = sheet2d[bound(x, y + 1, dm)];
            H10 = sheet2d[bound(x + 1, y, dm)];
            H11 = sheet2d[bound(x + 1, y + 1, dm)];

            values[i] = (ym * (xm * H00 + xp * H10) +
                         yp * (xm * H01 + xp * H11));
        }
        else
            values[i] = sheet2d[x + y * dm[0]];

        /* prevent unlikely values at the poles */
        if (values[i] < -1e15 || values[i] > 1e15)
            values[i] = 0.0;
    }
    delete_the_bintree(&unit_sphere.bintree);
    delete_polygons(&unit_sphere);
}
