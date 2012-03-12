/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_Blur2d.h"
#include "CAT_Curvature.h"

#define  BINTREE_FACTOR   0.5
#define  NEW_COORDINATE_SYSTEM   1


double
compute_clockwise_rotation2(double x, double y)
{
        double radians;

        if (x == 0.0) {
                if (y < 0.0)
                        return(PI / 2.0);
                else if (y > 0.0)
                        return(3.0 * PI / 2.0);
                else
                        return(0.0);
        } else if (y == 0.0) {
                if (x > 0.0)
                        return(0.0);
                else
                        return(PI);
        } else {
                radians = -(double) atan2(y, x);

                if (radians < 0.0)
                        radians += 2.0 * PI;

                return(radians);
        }
}


void
point_to_uv(Point *point, double *u, double *v)
{
        double x, y, z, theta, phi;

        x = (double) Point_x(*point);
        y = (double) Point_y(*point);
        z = (double) Point_z(*point);

        phi = acos(z);
        if (NEW_COORDINATE_SYSTEM) {
                theta = compute_clockwise_rotation(y, x) + PI/2.0;
                *u = theta / (PI * 2.0);
                if (*u > 1.0) *u -= 1.0;
        } else {
                theta = compute_clockwise_rotation2(y, x);
                *u = theta / (PI * 2.0);
        }
        *v = phi / PI;
}

void
uv_to_point(double u, double v, Point *point)
{
        double x, y, z, theta, phi, cos_u;
        double sin_u, cos_v, sin_v;
    
        if (NEW_COORDINATE_SYSTEM) {
                theta = u * PI * 2.0;
        } else {
                /* shift theta by 90 degree to obtain correct position of midline */    
                theta = u * PI * 2.0 + PI/2.0;
        }
        phi = v * PI;
 
        cos_u = cos(theta);
        sin_u = sin(theta);
        cos_v = cos(phi);
        sin_v = sin(phi);

        z = cos_v;
        x = sin_v * cos_u;
        y = sin_v * sin_u;

        fill_Point(*point, x, y, z);
}


void
wrap_sheet2d(double *data, int *dm, int wrap, double *wdata)
{
        int x, y, x0, y0, wdm[2];

        wdm[0] = dm[0] + 2*wrap;
        wdm[1] = dm[1] + 2*wrap;

        for (x = 0; x < dm[0]; x++) {
                for (y = 0; y < dm[1]; y++) {
                        x0 = x + wrap;
                        y0 = y + wrap;
                        wdata[x0 + (wdm[0]*y0)] = data[x + (dm[0]*y)];
                }
        }

        for (y = 0; y < wrap; y++) {
                for (x = wrap; x < wdm[0] - wrap; x++) {
                        x0 = (int) ((x - wrap + dm[0]/2) % dm[0]) + wrap;
                        y0 = 2*wrap - y - 1;
                        wdata[x + wdm[0]*y] = wdata[x0 + wdm[0]*y0];
                        y0 = wdm[1] - 2*wrap + y;
                        wdata[x + wdm[0]*(wdm[1] - y - 1)] =
                                                          wdata[x0 + wdm[0]*y0];
                }
        }

        for (x = 0; x < wrap; x++) {
                for (y = 0; y < wdm[1]; y++) {
                        x0 = wdm[0] - 2*wrap + x;
                        wdata[x + wdm[0]*y] = wdata[x0 + wdm[0]*y];
                        x0 = wrap + x;
                        wdata[x + wdm[0] - wrap + wdm[0]*y] =
                                                          wdata[x0 + wdm[0]*y];
                }
        }
}


void
unwrap_sheet2d(double *wdata, int *wdm, int wrap, double *data)
{
        int x, y, x0, y0, dm[2];

        dm[0] = wdm[0] - 2*wrap;
        dm[1] = wdm[1] - 2*wrap;

        for (x = 0; x < dm[0]; x++) {
                for (y = 0; y < dm[1]; y++) {
                        x0 = x + wrap;
                        y0 = y + wrap;
                        data[x + dm[0]*y] = wdata[x0 + wdm[0]*y0];
                }
        }

        for (y = 0; y < wrap; y++) {
                for (x = wrap; x < wdm[0] - wrap; x++) {
                        x0 = (int) ((x - wrap + dm[0]/2) % dm[0]);
                        y0 = wrap - y - 1;
                        data[x0 + dm[0]*y0] += wdata[x + wdm[0]*y];
                        if (x0 > wrap - 1 && x0 < dm[0] - wrap)
                                data[x0 + dm[0]*y0] /= 2;
                        y0 = dm[1] - wrap + y;
                        data[x0 + dm[0]*y0] += wdata[x + wdm[0]*(wdm[1] - y-1)];
                        if (x0 > wrap - 1 && x0 < dm[0] - wrap)
                                data[x0 + dm[0]*y0] /= 2;
                }
        }

        for (x = 0; x < wrap; x++) {
                for (y = 0; y < dm[1]; y++) {
                        data[x + dm[0]*y] += wdata[x + wdm[0] - wrap +
                                                   wdm[0]*(y + wrap)];
                        x0 = dm[0] - x - 1;
                        data[x0 + dm[0]*y] += wdata[wrap - x  - 1 +
                                                    wdm[0]*(y + wrap)];
                        if (y > wrap - 1 && y < dm[1] - wrap) {
                                data[x + dm[0]*y] /= 2;
                                data[x0 + dm[0]*y] /= 2;
                        }
                }
        }


        for (x = 0; x < wrap; x++) {
                for (y = 0; y < wrap; y++) {
                        x0 = dm[0] - x - 1;
                        y0 = dm[1] - y - 1;
                        data[x  + dm[0]*y ] /= 3;
                        data[x0 + dm[0]*y ] /= 3;
                        data[x  + dm[0]*y0] /= 3;
                        data[x0 + dm[0]*y0] /= 3;
                }
        }
}


void
map_smoothed_curvature_to_sphere(polygons_struct *polygons,
                                 polygons_struct *sphere, double *values,
                                 double *data, double fwhm, int *dm,
                                 int curvtype)
{
        polygons_struct   unit_sphere;
        Point             unit_point, on_sphere_point, centre;
        Point             poly_points[1000];
        double            value;
        double            u, v;
        double            weights[1000], mn, mx;
        int               i;
        int               x, y;
        int               poly, size, ind;

        /* if values is empty calculate curvature */
        if (values == (double *)0) {
                values = (double *) malloc(sizeof(double) * polygons->n_points);
                get_smoothed_curvatures(polygons, values, fwhm,
                                        curvtype);
        } else if (fwhm > 0) {
                smooth_heatkernel(polygons, values, fwhm);
        }
        
        if (sphere == NULL) {
                /* create unit sphere w/ same # of triangles as skin surface */
                fill_Point(centre, 0.0, 0.0, 0.0);
                create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                                  polygons->n_items, &unit_sphere);
        } else {
                copy_polygons(sphere, &unit_sphere);
                /* set radius to 1 */
                for (i = 0; i < unit_sphere.n_points; i++) 
                        set_vector_length(&unit_sphere.points[i], 1.0);
        }

        create_polygons_bintree(&unit_sphere,
                                ROUND((double) unit_sphere.n_items *
                                      BINTREE_FACTOR));

        for (x = 0; x < dm[0]; x++) {
                for (y = 0; y < dm[1]; y++) {
                        u = ((double) x + 0.5)/(double) (dm[0]);
                        v = ((double) y + 0.5)/(double) (dm[1]);
                        uv_to_point(u, v, &unit_point);
            
                        poly = find_closest_polygon_point(&unit_point,
                                                          &unit_sphere,
                                                          &on_sphere_point);
            
                        size = get_polygon_points(&unit_sphere, poly,
                                                  poly_points);

                        get_polygon_interpolation_weights(&on_sphere_point,
                                                          size, poly_points,
                                                          weights);

                        value = 0.0;
                        for (i = 0; i < size; i++) {
                                ind = unit_sphere.indices[
                                      POINT_INDEX(unit_sphere.end_indices,
                                                  poly, i)];
                                value += weights[i] * values[ind];
                        }
                        data[x + (dm[0]*y)] = value;
                }
        }

        /* scale data to uint8 range */
        mn = FLT_MAX; mx = -FLT_MAX;
        for (i = 0; i < dm[0]*dm[1]; i++) {
            if (data[i] > mx) mx = data[i];
            if (data[i] < mn) mn = data[i];
        }
    
        for (i = 0; i < dm[0]*dm[1]; i++) 
                data[i] = (data[i] - mn)/(mx - mn);
        
        delete_the_bintree(&unit_sphere.bintree);
        delete_polygons(&unit_sphere);
}

void
map_sheet2d_to_sphere(double *sheet2d, double *values,
                      polygons_struct *polygons, int interpolate, int *dm)
{
        double               tmp_x, tmp_y;
        double               u, v;
        Point                centre;
        polygons_struct      unit_sphere;
        int                  i, x, y;
        double               xp, yp, xm, ym;
        double               H00, H01, H10, H11;

        /* create a unit sphere with same number of triangles as skin surface */
        fill_Point(centre, 0.0, 0.0, 0.0);

        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                                  polygons->n_items, &unit_sphere);

        create_polygons_bintree(&unit_sphere,
                                ROUND((double) unit_sphere.n_items *
                                      BINTREE_FACTOR));

        tmp_x = (double)dm[0];
        tmp_y = (double)dm[1];

        for (i = 0; i < polygons->n_points; i++) {
                point_to_uv(&unit_sphere.points[i], &u, &v);

                x = (int) (u*tmp_x - 0.5);
                y = (int) (v*tmp_y - 0.5);
            
                if (interpolate) {
                        xp = u*tmp_x - 0.5 - x;
                        yp = v*tmp_y - 0.5 - y;
                        xm = 1.0 - xp;
                        ym = 1.0 - yp;
                        H00 = sheet2d[bound(x,  y,  dm)];
                        H01 = sheet2d[bound(x,  y+1,dm)];
                        H10 = sheet2d[bound(x+1,y,  dm)];
                        H11 = sheet2d[bound(x+1,y+1,dm)];
            
                        values[i] = (ym * ( xm * H00 + xp * H10) + 
                                     yp * ( xm * H01 + xp * H11));
                } else values[i] = sheet2d[x + y*dm[0]];
                
                /* prevent unlikely values at the poles */
                if (values[i] < -1e15 || values[i] > 1e15) values[i] = 0.0;
        
        }
        delete_the_bintree(&unit_sphere.bintree);
        delete_polygons(&unit_sphere);
}

void
map_sheet2d_to_unit_sphere(double *sheet2d, double *values,
                           polygons_struct *sphere, int interpolate, int *dm)
{
        double               tmp_x, tmp_y;
        double               u, v;
        int                  i, x, y;
        double               xp, yp, xm, ym;
        double               H00, H01, H10, H11;
        Point                unit_pt;

        create_polygons_bintree(sphere, ROUND((double) sphere->n_items *
                                      BINTREE_FACTOR));

        tmp_x = (double) dm[0];
        tmp_y = (double) dm[1];

        for (i = 0; i < sphere->n_points; i++) {
                unit_pt = sphere->points[i];
                set_vector_length(&unit_pt, 1.0);
                point_to_uv(&unit_pt, &u, &v);

                x = (int) (u*tmp_x - 0.5);
                y = (int) (v*tmp_y - 0.5);

                if (interpolate) {
                        xp = u*tmp_x - 0.5 - x;
                        yp = v*tmp_y - 0.5 - y;
                        xm = 1.0 - xp;
                        ym = 1.0 - yp;
                        H00 = sheet2d[bound(x,  y,  dm)];
                        H01 = sheet2d[bound(x,  y+1,dm)];
                        H10 = sheet2d[bound(x+1,y,  dm)];
                        H11 = sheet2d[bound(x+1,y+1,dm)];
            
                        values[i] = (ym * ( xm * H00 + xp * H10) + 
                                     yp * ( xm * H01 + xp * H11));
                } else values[i] = sheet2d[x + y*dm[0]];
                
                /* prevent unlikely values at the poles */
                if (values[i] < -1e15 || values[i] > 1e15) values[i] = 0.0;
        
        }
        delete_the_bintree(&sphere->bintree);
}

/* smooth the sheet such that diff between neighbours <= tol */
void
smooth_sheet2d(double *data, int dm[2], double tol)
{
        int x, y, i, done;
        double x0, x1, y0, y1, *buf;

        buf = (double *) malloc(sizeof(double) * dm[0] * dm[1]);

        done = 0;
        while (!done) {
                done = 1;
                for (x = 0; x < dm[0]; x++) {
                        for (y = 0; y < dm[1]; y++) {
                                i = x + dm[0]*y;
                                buf[i] = data[i];

                                x0 = data[bound(x-1,y-1,dm)];
                                x1 = data[bound(x+1,y-1,dm)];
                                y0 = data[bound(x-1,y+1,dm)];
                                y1 = data[bound(x+1,y+1,dm)];

                                if (fabs(x0 - data[i]) > tol ||
                                    fabs(x1 - data[i]) > tol ||
                                    fabs(y0 - data[i]) > tol ||
                                    fabs(y1 - data[i]) > tol) {
                                        buf[i] = data[i]/2.0 + (x0+x1)/8.0 +
                                                 (y0+y1)/8.0;
                                        done = 0;
                                }
                        }
                }
                for (x = 0; x < dm[0]*dm[1]; x++) {
                        data[x] = buf[x];
                }
        }
        free(buf);
}
