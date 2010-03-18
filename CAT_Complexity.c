/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Complexity.c 139 2009-10-02 15:39:28Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "Cat_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Complexity.h"


/* get the slope */
double
slope(double *x, double *y, int len)
{
        int i;
        double slope = 0;

        for (i = 1; i < len; i++) {
                slope += (y[i] - y[i-1]) / (x[i] - x[i-1]);
        }

        slope /= len;

        return(slope);
}

/* get global FD values ... only use area values below the threshold */
double
get_globalfd(double *x, double *y, int len, double threshold)
{
        int i, tlen;
        double *logx, *logy, fd;

        for (tlen = len-1; tlen > 0; tlen--) {
                if (y[tlen] < threshold) break;
        }
        tlen++;

        logx = (double *) malloc(sizeof(double) * tlen);
        logy = (double *) malloc(sizeof(double) * tlen);
        for (i = 0; i < tlen; i++) {
                logx[i] = log(x[i]);
                logy[i] = log(y[i]);
        }

        fd = slope(logx, logy, tlen);

        free(logx);
        free(logy);

        return fd;
}

/* get local FD values ... x is bandwidth or dimension.
 * only uses area values below the threshold */
void
get_localfd(double *x, double **areas, int x_len, int n_points, double *fd,
            double threshold)
{
        int xx, p, len;
        double *logx, *logy;

        logx = (double *) malloc(sizeof(double) * x_len);
        for (xx = 0; xx < x_len; xx++)
                logx[xx] = log(x[xx]);

        logy = (double *) malloc(sizeof(double) * x_len);
        for (p = 0; p < n_points; p++) {
                for (len = x_len-1; len > 0; len--) {
                        if (areas[len][p] < threshold) break;
                }
                len++;

                for (xx = 0; xx < len; xx++)
                        logy[xx] = log(areas[xx][p]);

                fd[p] = slope(logx, logy, len);

        }

        free(logx);
        free(logy);
}

object_struct **
create_resampling_sphere(double radius, int n_triangles)
{
        printf("radius = %f, n_triangles = %d\n", radius, n_triangles);
        int iterations, n_points, n_items, n_p, n_i, *indices, base;
        polygons_struct *polygons;
        object_struct **object;
        int n, n2, i, a, b, c, p, p2, count, offset, iter;
        Point *points, v;
        double len, dist;

        /* overestimate the # of points initially */
        n_points = n_triangles;
        n_items = n_triangles;

        points = (Point *) malloc(sizeof(Point) * n_points);
        indices = (int *) malloc(sizeof(int) * n_items * 3);

        /* figure out the base and the iterations first */
        iterations = 0;
        base = n_triangles;
        while (base > 8 && base != 20) {
                base /= 4;
                iterations++;
        }
        printf("base = %d, iters = %d\n", base, iterations);

        /* create the base shape */
        if (base == 4) { /* tetrahedral */
                n_points = 4; n_items = 4;
                fill_Point(points[0], 1.0, 1.0, 1.0);
                fill_Point(points[1], -1.0, -1.0, 1.0);
                fill_Point(points[2], 1.0, -1.0, -1.0);
                fill_Point(points[3], -1.0, 1.0, -1.0);
                indices[0] = 0; indices[1] = 1; indices[2] = 2;
                indices[3] = 0; indices[4] = 2; indices[5] = 3;
                indices[6] = 0; indices[7] = 3; indices[8] = 1;
                indices[9] = 1; indices[10] = 3; indices[11] = 2;
        } else if (base == 6) { /* hexahedral */
                n_points = 5; n_items = 6;
                fill_Point(points[0], 1.0, 0.0, 0.0);
                fill_Point(points[1], -1.0, 0.0, 0.0);
                fill_Point(points[2], 0.0, 1.0, 0.0);
                fill_Point(points[3], 0.0, -0.5, sqrt(3)/2.0);
                fill_Point(points[4], 0.0, -0.5, -sqrt(3)/2.0);
                indices[0] = 0; indices[1] = 2; indices[2] = 3;
                indices[3] = 0; indices[4] = 3; indices[5] = 4;
                indices[6] = 0; indices[7] = 4; indices[8] = 2;
                indices[9] = 1; indices[10] = 3; indices[11] = 2;
                indices[12] = 1; indices[13] = 4; indices[14] = 3;
                indices[15] = 1; indices[16] = 2; indices[17] = 4;
        } else if (base == 8) { /* octahedral */
                n_points = 6; n_items = 8;
                fill_Point(points[0], 1.0, 0.0, 0.0);
                fill_Point(points[1], -1.0, 0.0, 0.0);
                fill_Point(points[2], 0.0, 1.0, 0.0);
                fill_Point(points[3], 0.0, 0.0, 1.0);
                fill_Point(points[4], 0.0, -1.0, 0.0);
                fill_Point(points[5], 0.0, 0.0, -1.0);
                indices[0] = 0; indices[1] = 2; indices[2] = 3;
                indices[3] = 0; indices[4] = 3; indices[5] = 4;
                indices[6] = 0; indices[7] = 4; indices[8] = 5;
                indices[9] = 0; indices[10] = 5; indices[11] = 2;
                indices[12] = 1; indices[13] = 3; indices[14] = 2;
                indices[15] = 1; indices[16] = 4; indices[17] = 3;
                indices[18] = 1; indices[19] = 5; indices[20] = 4;
                indices[21] = 1; indices[22] = 2; indices[23] = 5;
        } else if (base == 20) { /* icosahedral */
                n_points = 12; n_items = 20;
                fill_Point(points[0], (sqrt(5)+1)/2, 1.0, 0.0);
                fill_Point(points[1], (sqrt(5)+1)/2, -1.0, 0.0);
                fill_Point(points[2], -(sqrt(5)+1)/2, 1.0, 0.0);
                fill_Point(points[3], -(sqrt(5)+1)/2, -1.0, 0.0);
                fill_Point(points[4], 0.0, (sqrt(5)+1)/2, 1.0);
                fill_Point(points[5], 0.0, (sqrt(5)+1)/2, -1.0);
                fill_Point(points[6], 0.0, -(sqrt(5)+1)/2, 1.0);
                fill_Point(points[7], 0.0, -(sqrt(5)+1)/2, -1.0);
                fill_Point(points[8], 1.0, 0.0, (sqrt(5)+1)/2);
                fill_Point(points[9], -1.0, 0.0, (sqrt(5)+1)/2);
                fill_Point(points[10], 1.0, 0.0, -(sqrt(5)+1)/2);
                fill_Point(points[11], -1.0, 0.0, -(sqrt(5)+1)/2);
                indices[0] = 2; indices[1] = 11; indices[2] = 3; 
                indices[3] = 2; indices[4] = 5; indices[5] = 11; 
                indices[6] = 11; indices[7] = 5; indices[8] = 10; 
                indices[9] = 5; indices[10] = 0; indices[11] = 10; 
                indices[12] = 0; indices[13] = 5; indices[14] = 4; 
                indices[15] = 4; indices[16] = 5; indices[17] = 2; 
                indices[18] = 4; indices[19] = 2; indices[20] = 9; 
                indices[21] = 9; indices[22] = 2; indices[23] = 3; 
                indices[24] = 9; indices[25] = 3; indices[26] = 6; 
                indices[27] = 9; indices[28] = 6; indices[29] = 8; 
                indices[30] = 4; indices[31] = 9; indices[32] = 8; 
                indices[33] = 4; indices[34] = 8; indices[35] = 0; 
                indices[36] = 0; indices[37] = 8; indices[38] = 1; 
                indices[39] = 1; indices[40] = 8; indices[41] = 6; 
                indices[42] = 1; indices[43] = 6; indices[44] = 7; 
                indices[45] = 7; indices[46] = 6; indices[47] = 3; 
                indices[48] = 11; indices[49] = 7; indices[50] = 3; 
                indices[51] = 11; indices[52] = 10; indices[53] = 7; 
                indices[54] = 10; indices[55] = 1; indices[56] = 7; 
                indices[57] = 10; indices[58] = 0; indices[59] = 1; 
        } else {
                fprintf(stderr, "Error: base of %d is invalid!\n", base);
                exit(EXIT_FAILURE);
        }

        SUB_POINTS(v, points[indices[0]], points[indices[1]]);
        len = sqrt(DOT_POINTS(v, v));

        /* cut each triangle into 4 triangles */
        for (iter = 0; iter < iterations; iter++) {
                len /= 2;

                n_i = n_items;
                for (n = 0; n < n_i; n++) {
                        a = indices[n*3];
                        b = indices[n*3 + 1];
                        c = indices[n*3 + 2];

                        INTERPOLATE_POINTS(points[n_points],
                                           points[a], points[b], 0.5);
                        INTERPOLATE_POINTS(points[n_points+1],
                                           points[b], points[c], 0.5);
                        INTERPOLATE_POINTS(points[n_points+2],
                                           points[a], points[c], 0.5);
                        n_points += 3;

                        /* add the triangle indices */
                        indices[n*3 + 1] = n_points - 3;
                        indices[n*3 + 2] = n_points - 1;

                        indices[n_items*3] = n_points - 3;
                        indices[n_items*3 + 1] = b;
                        indices[n_items*3 + 2] = n_points - 2;
                        n_items++;
                        indices[n_items*3] = n_points - 2;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points - 1;
                        n_items++;
                        indices[n_items*3 + 2] = n_points - 3;
                        indices[n_items*3] = n_points - 2;
                        indices[n_items*3 + 1] = n_points - 1;
                        n_items++;
                }
                /* remove duplicate points */
                for (p = 0; p < n_points; p++) {
                        for (p2 = p+1; p2 < n_points; p2++) {
                                dist = distance_between_points(&points[p],
                                                               &points[p2]);
                                if (dist < 0.01*len) { /* delete it */
                                        n_points--;
                                        for (n = 0; n < n_items * 3; n++) {
                                                if (indices[n] == p2)
                                                        indices[n] = p;
                                                else if (indices[n] == n_points)
                                                        indices[n] = p2;
                                        }
                                        points[p2] = points[n_points];
                                        break;
                                }
                        }
                }
                /* remove duplicate triangles */
                offset = 0;
                for (n = 0; n < n_items; n++) {
                        for (n2 = n+1; n2 < n_items; n2++) {
                                if ( (indices[n*3  ] == indices[n2*3  ] ||
                                      indices[n*3  ] == indices[n2*3+1] ||
                                      indices[n*3  ] == indices[n2*3+2]) &&
                                     (indices[n*3+1] == indices[n2*3  ] ||
                                      indices[n*3+1] == indices[n2*3+1] ||
                                      indices[n*3+1] == indices[n2*3+2]) &&
                                     (indices[n*3+2] == indices[n2*3  ] ||
                                      indices[n*3+2] == indices[n2*3+1] ||
                                      indices[n*3+2] == indices[n2*3+2])) {
                                        n_items--;
                                        offset++;
                                }
                        }
                        indices[n*3] = indices[(n+offset)*3];
                        indices[n*3+1] = indices[(n+offset)*3+1];
                        indices[n*3+2] = indices[(n+offset)*3+2];
                }
                for (i = 0; i < n_points; i++) {
                        set_vector_length(&points[i], radius);
                }
        }

        /* build surface */
        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        polygons = get_polygons_ptr(*object);
        initialize_polygons(polygons, WHITE, NULL);

        Surfprop_a(polygons->surfprop) = 0.3;
        Surfprop_d(polygons->surfprop) = 0.6;
        Surfprop_s(polygons->surfprop) = 0.6;
        Surfprop_se(polygons->surfprop) = 60;
        Surfprop_t(polygons->surfprop) = 1.0;

        polygons->colour_flag = 0;
        polygons->line_thickness = 1.0f;

        polygons->points = (Point *) malloc(sizeof(Point) * n_points);
        polygons->indices = (int *) malloc(sizeof(int) * 3 * n_items);
        polygons->end_indices = (int *) malloc(sizeof(int) * n_items);
        polygons->normals = (Vector *) malloc(sizeof(Vector) * n_points);
        polygons->n_points = n_points;
        polygons->n_items = n_items;

        for (i = 0; i < polygons->n_points; i++) {
                polygons->points[i] = points[i];
                set_vector_length(&polygons->points[i], radius);
        }

        for (i = 0; i < n_items; i++) {
                polygons->end_indices[i] = 3 * (i + 1);
                polygons->indices[i] = indices[i];
        }

        for (i = n_items; i < n_items * 3; i++) 
                polygons->indices[i] = indices[i];

        compute_polygon_normals(polygons);

        free(points);
        free(indices);

        return(object);
}

int
min_triangles_update(int *base4, int *base6, int *base8, int *base20)
{
        int retval;

        if (*base4 < *base6 && *base4 < *base8 && *base4 < *base20) {
                retval = *base4;
                *base4 *= 4;
        } else if (*base6 < *base8 && *base6 < *base20) {
                retval = *base6;
                *base6 *= 4;
        } else if (*base8 < *base20) {
                retval = *base8;
                *base8 *= 4;
        } else {
                retval = *base20;
                *base20 *= 4;
        }

        return(retval);
}

object_struct **
resample_surface(polygons_struct *surface, polygons_struct *sphere,
                 int n_triangles, double *invals, double *outvals)
{
        int i, j, t, poly, n_points;
        Point point, centre, scaled_point;
        Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        object_struct **objects;
        polygons_struct *platonic_solid;
        double radius, r, bounds[6];
        Real weights[MAX_POINTS_PER_POLYGON];

        /* Check tetrahedral topology. Best areal distribution of triangles
         * is achieved for 20 edges
         */ /*
        t = n_triangles;
        while (t != 20 && t > 8 && t % 4 == 0)
                t /= 4;

        if (t != 20) {
                fprintf(stderr, "Warning: Number of triangles %d", n_triangles);
                fprintf(stderr," is not recommended because\ntetrahedral ");
                fprintf(stderr,"topology is not optimal.\n");
                fprintf(stderr,"Please try 20*(4*x) triangles (e.g. 81920).\n");
        } */
	
        /* Determine radius for the output sphere.  The sphere is not always
         * perfectly spherical, thus use average radius
         */
        radius = 0.0;
        for (i = 0; i < sphere->n_points; i++) {
                r = 0.0;
                for (j = 0; j < 3; j++) 
                        r += Point_coord(sphere->points[i], j) *
                             Point_coord(sphere->points[i], j);
                radius += sqrt(r);
        }
        radius /= sphere->n_points;

        /* Calc. sphere center based on bounds of input (correct for shifts) */
        get_bounds(sphere, bounds);
        fill_Point(centre, (bounds[0]+bounds[1]) / 2,
                           (bounds[2]+bounds[3]) / 2,
                           (bounds[4]+bounds[5]) / 2);
    
        objects = create_resampling_sphere(radius, n_triangles);
        platonic_solid = get_polygons_ptr(*objects);

        char str[80];
        sprintf(str, "base_%d.obj", n_triangles);
        if (output_graphics_any_format(str, ASCII_FORMAT, 1, objects) != OK)
                        exit(EXIT_FAILURE);

        create_polygons_bintree(sphere, ROUND((Real) sphere->n_items * 0.5));

        new_points = (Point *) malloc(sizeof(Point) * platonic_solid->n_points);
        if (invals != NULL) {
                outvals = (double *) malloc(sizeof(double) *
                                            platonic_solid->n_points);
        }

        for (i = 0; i < platonic_solid->n_points; i++) {
                poly = find_closest_polygon_point(&platonic_solid->points[i],
                                                  sphere, &point);
		
                n_points = get_polygon_points(sphere, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                if (get_polygon_points(surface, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                fill_Point(new_points[i], 0.0, 0.0, 0.0);
                if (invals != NULL)
                        outvals[i] = 0.0;

                for (j = 0; j < n_points; j++) {
                        SCALE_POINT(scaled_point, poly_points[j], weights[j]);
                        ADD_POINTS(new_points[i], new_points[i], scaled_point);
                        if (invals != NULL) {
                                outvals[i] += weights[j] *
                                               invals[surface->indices[
                                     POINT_INDEX(surface->end_indices,poly,j)]];
                        }
                }
        }

        free(platonic_solid->points);
        platonic_solid->points = new_points;

        compute_polygon_normals(platonic_solid);

        return(objects);
}

/*
 * Compute the direct fractal dimension by re-parametrizing the surface
 * using progressively fewer triangles & computing the area
 */
double
fractal_dimension(polygons_struct *surface, polygons_struct *sphere,
                  int maxiters, char *file)
{
        object_struct **object;
        polygons_struct *polygons;
        double orig_area;
        char str[80];
        int base4 = 64, base6 = 96, base8 = 128, base20 = 80;
        int iter, n, n_triangles;
        double *areas, *dimension, fd;
        
        orig_area = get_polygons_surface_area(surface);
        n_triangles = base4;
        base4 *= 4;

        areas = (double *) malloc(sizeof(double) * maxiters);
        dimension = (double *) malloc(sizeof(double) * maxiters);

        for (iter = 0; iter < maxiters; iter++) {
                dimension[iter] = 1/sqrt(n_triangles);
                object = resample_surface(surface, sphere, n_triangles, NULL,
                                           NULL);

                sprintf(str, "resamp_%d.obj", n_triangles);
                if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                               object) != OK)
                        exit(EXIT_FAILURE);

                polygons = get_polygons_ptr(*object);
                areas[iter] = get_polygons_surface_area(polygons) / orig_area;

                n_triangles = min_triangles_update(&base4, &base6, &base8,
                                                   &base20);
                delete_object_list(1, object);
        }

        for (n = 0; n < maxiters; n++) {
                if (areas[iter] > 0.9)
                        break;
        }

        get_localfd(dimension, &areas, n, 1, &fd, 0.9);

        if (output_values_any_format(file, maxiters, areas, TYPE_DOUBLE) != OK)
                exit(EXIT_FAILURE);
        
        return fd;
}

double
get_smoothed_areas(polygons_struct *polygons, double *orig_areas, double *areas)
{
        int i, j, n_iter;
        double fwhm, sigma, sa, *sm_areas;
        int *n_neighbours, **neighbours;
        Point point;

        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        sa = get_area_of_points(polygons, areas);
        sm_areas = (double *) malloc(sizeof(double) * polygons->n_points);

        /* normalize the areas */
        if (orig_areas != NULL) {
                for (i = 0; i < polygons->n_points; i++)
                        areas[i] = areas[i] / orig_areas[i];
        }

        /* smooth the areas */
        fwhm = 25.0;
        sigma = 2.0;
        n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482/sigma);
    
        /* diffusion smoothing using heat kernel */
        for (j = 0; j < n_iter; j++) {
                for (i = 0; i < polygons->n_points; i++) {
                        heatkernel_blur_points(polygons->n_points,
                                               polygons->points, areas,
                                               n_neighbours[i], neighbours[i],
                                               i, sigma, &point, &sm_areas[i]);
                }
                for (i = 0; i < polygons->n_points; i++)
                        areas[i] = sm_areas[i];
        }

        free(sm_areas);
        return sa;
}


double bws[10] = {8.0, 10.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 48.0, 64.0};

/*
 * Compute the fractal dimension using spherical harmonics: progressively
 * lower bandwidths & area computation.
 */
double
fractal_dimension_sph(polygons_struct *surface, polygons_struct *sphere,
                      char *file, int n_triangles)
{
        polygons_struct *polygons;
        object_struct **object;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz, *spectral_power;
        double *areas, *orig_areas, fd, orig_area, **sph_areas, *local_fd;
        int it, bw2, l, m, i;
        char str[80];

        bw2 = BW * BW;

        rdatax = (double *) malloc(sizeof(double) * 4 * bw2);
        rdatay = (double *) malloc(sizeof(double) * 4 * bw2);
        rdataz = (double *) malloc(sizeof(double) * 4 * bw2);
        rcx    = (double *) malloc(sizeof(double) * bw2);
        rcy    = (double *) malloc(sizeof(double) * bw2);
        rcz    = (double *) malloc(sizeof(double) * bw2);
        icx    = (double *) malloc(sizeof(double) * bw2);
        icy    = (double *) malloc(sizeof(double) * bw2);
        icz    = (double *) malloc(sizeof(double) * bw2);
        lrcx   = (double *) malloc(sizeof(double) * bw2);
        lrcy   = (double *) malloc(sizeof(double) * bw2);
        lrcz   = (double *) malloc(sizeof(double) * bw2);
        licx   = (double *) malloc(sizeof(double) * bw2);
        licy   = (double *) malloc(sizeof(double) * bw2);
        licz   = (double *) malloc(sizeof(double) * bw2);

        fprintf(stderr,"Samp SPH coords..");
        get_equally_sampled_coords_of_polygon(surface, sphere, BW,
                                              rdatax, rdatay, rdataz);

        fprintf(stderr,"Fwd SPH xform..");
        get_sph_coeffs_of_realdata(rdatax, BW, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, BW, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, BW, DATAFORMAT, rcz, icz);

        /* output the spectral power */
        spectral_power = (double *) malloc(sizeof(double) * BW);
        memset(spectral_power, 0.0, sizeof(double) * BW);
        for (l = 0; l < BW; l++) {
                for (m = -l; m < l+1; m++) {
                        i = seanindex(m, l, BW);
                        spectral_power[l] += (rcx[i]*rcx[i]) + (rcy[i]*rcy[i]) +
                                             (rcz[i]*rcz[i]) + (icx[i]*icx[i]) +
                                             (icy[i]*icy[i]) + (icz[i]*icz[i]);
                }
        }
        output_values_any_format("psd.txt", BW, spectral_power, TYPE_DOUBLE);
        free(spectral_power);

        fprintf(stderr,"Inv SPH xform..");
        get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT, rcz, icz);

        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        polygons = get_polygons_ptr(*object);

        fprintf(stderr,"Resamp surf.\n");
        sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons,
                               n_triangles, BW);

        orig_areas = (double *) malloc(sizeof(double) * polygons->n_points);
        get_smoothed_areas(polygons, NULL, orig_areas);
        orig_area = get_polygons_surface_area(polygons);

        sprintf(str, "sarea_%d.txt", BW);
        output_values_any_format(str, polygons->n_points, orig_areas,
                                 TYPE_DOUBLE);
        sprintf(str, "sph_%d.obj", BW);
        if (output_graphics_any_format(str, ASCII_FORMAT, 1, object) != OK)
                exit(EXIT_FAILURE);
        delete_polygons(polygons);

        areas = (double *) malloc(sizeof(double) * SPH_ITERS);
        sph_areas = (double **) malloc(sizeof(double *) * SPH_ITERS);
        for (i = 0; i < SPH_ITERS; i++) {
                sph_areas[i] = (double *) malloc(sizeof(double) *
                                                 polygons->n_points);
        }

        for (it = 0; it < SPH_ITERS; it++) {
                fprintf(stderr, "BW %d: ", (int) bws[it]);

                fprintf(stderr,"Samp SPH coords..");
                limit_bandwidth(BW, (int) bws[it], rcx, lrcx);
                limit_bandwidth(BW, (int) bws[it], rcy, lrcy);
                limit_bandwidth(BW, (int) bws[it], rcz, lrcz);
                limit_bandwidth(BW, (int) bws[it], icx, licx);
                limit_bandwidth(BW, (int) bws[it], icy, licy);
                limit_bandwidth(BW, (int) bws[it], icz, licz);

                fprintf(stderr,"Inv SPH xform..");
                get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT,
                                             lrcx, licx);
                get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT,
                                             lrcy, licy);
                get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT,
                                             lrcz, licz);

                fprintf(stderr,"Resamp surf.\n");
                sample_sphere_from_sph(rdatax, rdatay, rdataz,
                                       polygons, n_triangles, BW);

                areas[it] = get_polygons_surface_area(polygons) / orig_area;
                get_smoothed_areas(polygons, orig_areas, sph_areas[it]);
                printf("bw = %d, area = %f\n", (int) bws[it], areas[it]);

                sprintf(str, "sph_%d.obj", (int) bws[it]);
                if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                               object) != OK)
                        exit(EXIT_FAILURE);

                sprintf(str, "sarea_%d.txt", (int) bws[it]);
                output_values_any_format(str, polygons->n_points,
                                         sph_areas[it], TYPE_DOUBLE);

                delete_polygons(polygons);
        }

        local_fd = (double *) malloc(sizeof(double) * polygons->n_points);
        get_localfd(bws, sph_areas, SPH_ITERS, polygons->n_points,
                    local_fd, 0.9);
        sprintf(str, "fd_local.txt");
        output_values_any_format(str, polygons->n_points, local_fd,
                                 TYPE_DOUBLE);
 
        if (output_values_any_format(file, SPH_ITERS, areas, TYPE_DOUBLE) != OK)
                exit(EXIT_FAILURE);

        fd = get_globalfd(bws, areas, SPH_ITERS, 0.8);

	free(rcx); free(rcy); free(rcz); 
        free(icx); free(icy); free(icz); 
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);
        
        return fd;
}

