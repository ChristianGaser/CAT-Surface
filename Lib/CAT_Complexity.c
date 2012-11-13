/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "Cat_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Complexity.h"
#include "CAT_Resample.h"


/* get the slope */
double
slope(double *x, double *y, int len)
{
        int i;
        double slope = 0;

        for (i = 1; i < len; i++) 
                slope += (y[i] - y[i-1]) / (x[i] - x[i-1]);

        if (len > 1) 
                slope /= (len-1);

        return(slope);
}

/* get global FD values ... only use area values below the threshold */
double
get_globalfd(double *x, double *y, int len)
{
        int i;
        double *logx, *logy, fd;

        logx = (double *) malloc(sizeof(double) * len);
        logy = (double *) malloc(sizeof(double) * len);
        for (i = 0; i < len; i++) {
                logx[i] = log(x[i]);
                logy[i] = log(y[i]);
        }

        fd = 2 + slope(logx, logy, len);

        free(logx);
        free(logy);

        return fd;
}

/* get local FD values ... x is bandwidth or dimension. */
void
get_localfd(polygons_struct *polygons, double *x, double **areas, int x_len,
            double *fd, int smoothflag)
{
        int xx, p, offset;
        double *logx, *logy;
        int *pcount;
        double *polyfd, poly_size, area;
        Point points[MAX_POINTS_PER_POLYGON];
        int ptidx, poly, vertidx, size;

        polyfd = (double *) malloc(sizeof(double) * polygons->n_items);
        logx = (double *) malloc(sizeof(double) * x_len);
        logy = (double *) malloc(sizeof(double) * x_len);

        /* calculate the FD for each polygon */
        for (xx = 0; xx < x_len; xx++)
                logx[xx] = log(x[xx]);

        for (p = 0; p < polygons->n_items; p++) {
                offset = 0;
                for (xx = 0; xx < x_len; xx++) {
                        if (areas[xx][p] == 0)
                                offset = xx+1;
                        logy[xx] = log(areas[xx][p]);
                }

                polyfd[p] = slope(logx + offset, logy + offset, x_len - offset);
        }

        /* calculate a point-wise FD value based on neighboring polygons */
        pcount = (int *) malloc(sizeof(int) * polygons->n_points);
        memset(pcount, 0, sizeof(int) * polygons->n_points);
        memset(fd, 0.0, sizeof(double) * polygons->n_points);

        for (poly = 0; poly < polygons->n_items; poly++) {    
                size = GET_OBJECT_SIZE(*polygons, poly);

                for (p = 0; p < size; p++) {
                        ptidx = polygons->indices[
                                              POINT_INDEX(polygons->end_indices,
                                              poly, p)];
                        pcount[ptidx]++;
                        fd[ptidx] += polyfd[poly];
                }
        }

        for (p = 0; p < polygons->n_points; p++) {
                if (pcount[p] > 0)
                        fd[p] /= pcount[p];
                fd[p] += 2;
        }

        if (smoothflag) {
                /* smooth the FD values with FWHM = 30.0 mm */
                get_smoothed_values(polygons, fd, FWHM);
        }

        free(logx);
        free(logy);
        free(pcount);
}


int
min_triangles_update(int *base6, int *base8, int *base20)
{
        int retval;

        if (*base6 < *base8 && *base6 < *base20) {
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

/*
 * Compute the direct fractal dimension by re-parametrizing the surface
 * using progressively fewer triangles & computing the area
 */
double
fractal_dimension(polygons_struct *surface, polygons_struct *sphere,
                  int maxiters, char *file, int smoothflag, int debugflag)
{
        object_struct **object, **object2;
        polygons_struct *polygons, *resampled;
        double orig_area, *orig_areas;
        char str[80];
        int base6 = 96, base8 = 128, base20 = 80;
        int i, iter, n, n_triangles;
        double **bc_areas, *areas, *dimension, fd, *local_fd;
        Point centre;
        
        translate_to_center_of_mass(sphere);

        n_triangles = base20;
        base20 *= 4;

        orig_areas = (double *) malloc(sizeof(double) * surface->n_items);
        orig_area = get_area_of_polygons(surface, orig_areas);

        dimension = (double *) malloc(sizeof(double) * maxiters);
        areas = (double *) malloc(sizeof(double) * maxiters);

        bc_areas = (double **) malloc(sizeof(double *) * maxiters);
        for (i = 0; i < maxiters; i++) {
                bc_areas[i] = (double *) malloc(sizeof(double) *
                                                surface->n_items);
        }

        for (iter = 0; iter < maxiters; iter++) {
                dimension[iter] = sqrt(n_triangles);
                object = resample_surface(surface, sphere, n_triangles, NULL,
                                           NULL);

                polygons = get_polygons_ptr(*object);
                areas[iter] = get_polygons_surface_area(polygons) /
                              orig_area;

                /* resample back into the original object space */
                object2 = resample_surface_sphere(polygons, sphere);
                resampled = get_polygons_ptr(*object2);

                get_area_of_polygons(resampled, bc_areas[iter]);

                printf("n_tri = %d, areas = %f\n", n_triangles, areas[iter]);

                for (i = 0; i < sphere->n_items; i++) {
                        if (orig_areas[i] != 0)
                                bc_areas[iter][i] /= orig_areas[i];
                }

                if (debugflag) {
                        sprintf(str, "resamp_%d.obj", n_triangles);
                        if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                                       object) != OK)
                                exit(EXIT_FAILURE);
                        sprintf(str, "reresamp_%d.obj", n_triangles);
                        if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                                       object2) != OK)
                                exit(EXIT_FAILURE);
                        sprintf(str, "areas_%d.txt", n_triangles);
                        if (output_values_any_format(str, sphere->n_items,
                                                     bc_areas[iter],
                                                     TYPE_DOUBLE) != OK)
                                exit(EXIT_FAILURE);
                }

                n_triangles = min_triangles_update(&base6, &base8, &base20);

                delete_object_list(1, object);
                delete_object_list(1, object2);

                //if (areas[iter] > 0.9)
                        //break; /* stop here */
        }

        if (iter < maxiters) iter--; /* skip >90% values */

        local_fd = (double *) malloc(sizeof(double) * surface->n_points);
        get_localfd(surface, dimension, bc_areas, iter, local_fd, smoothflag);

        if (output_values_any_format(file, surface->n_points,
                                     local_fd, TYPE_DOUBLE) != OK)
                exit(EXIT_FAILURE);

        fd = get_globalfd(dimension, areas, iter);

        if (1) {
                if (output_values_any_format("fd_global.txt", maxiters,
                                     areas, TYPE_DOUBLE) != OK)
                                exit(EXIT_FAILURE);
        }

        free(orig_areas);
        free(dimension);
        free(areas);
        free(local_fd);

        for (i = 0; i < maxiters; i++)
                free(bc_areas[i]);
        free(bc_areas);

        return fd;
}

void
get_smoothed_values(polygons_struct *polygons, double *values, double fwhm)
{
        int i, j, n_iter;
        double sigma, *sm_values;
        int *n_neighbours, **neighbours;
        Point point;

        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
        smooth_heatkernel(polygons, values, fwhm);
}


double bws[SPH_ITERS] = {11.0, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0, 23.0, 26.0, 29.0};

//double bws[SPH_ITERS] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0}; // for testing
//double bws[SPH_ITERS] = {5.0, 6.0, 7.0, 8.0}; // for von Koch surfaces

/*
 * Compute the fractal dimension using spherical harmonics: progressively
 * lower bandwidths & area computation.
 */
double
fractal_dimension_sph(polygons_struct *surface, polygons_struct *sphere,
                      char *file, int n_triangles, polygons_struct *reparam,
                      int smoothflag, int debugflag)
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

        if (debugflag) fprintf(stderr,"Samp SPH coords..");
        get_equally_sampled_coords_of_polygon(surface, sphere, BW,
                                              rdatax, rdatay, rdataz);

        if (debugflag) fprintf(stderr,"Fwd SPH xform..");
        get_sph_coeffs_of_realdata(rdatax, BW, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, BW, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, BW, DATAFORMAT, rcz, icz);

        if (debugflag) {
                /* output the spectral power */
                spectral_power = (double *) malloc(sizeof(double) * BW);
                memset(spectral_power, 0.0, sizeof(double) * BW);
                for (l = 0; l < BW; l++) {
                        for (m = -l; m < l+1; m++) {
                                i = seanindex(m, l, BW);
                                spectral_power[l] += (rcx[i]*rcx[i]) +
                                                     (rcy[i]*rcy[i]) +
                                                     (rcz[i]*rcz[i]) +
                                                     (icx[i]*icx[i]) +
                                                     (icy[i]*icy[i]) +
                                                     (icz[i]*icz[i]);
                        }
                }
                output_values_any_format("psd.txt", BW, spectral_power,
                                         TYPE_DOUBLE);
                free(spectral_power);
        }

        if (debugflag) fprintf(stderr,"Inv SPH xform..");
        get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT, rcz, icz);

        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        polygons = get_polygons_ptr(*object);

        if (debugflag) fprintf(stderr,"Resamp surf.\n");
        sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons,
                               n_triangles, reparam, BW);

        orig_areas = (double *) malloc(sizeof(double) * polygons->n_items);
        orig_area = get_area_of_polygons(polygons, orig_areas);

        if (debugflag) {
                sprintf(str, "area_%d.txt", BW);
                output_values_any_format(str, polygons->n_points, orig_areas,
                                         TYPE_DOUBLE);
                sprintf(str, "sph_%d.obj", BW);
                if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                               object) != OK)
                        exit(EXIT_FAILURE);
        }
        delete_polygons(polygons);

        areas = (double *) malloc(sizeof(double) * SPH_ITERS);
        sph_areas = (double **) malloc(sizeof(double *) * SPH_ITERS);
        for (i = 0; i < SPH_ITERS; i++) {
                sph_areas[i] = (double *) malloc(sizeof(double) *
                                                 polygons->n_items);
        }

        for (it = 0; it < SPH_ITERS; it++) {
                if (debugflag) fprintf(stderr, "BW %d: ", (int) bws[it]);

                if (debugflag) fprintf(stderr,"Samp SPH coords..");
                limit_bandwidth(BW, (int) bws[it], rcx, lrcx);
                limit_bandwidth(BW, (int) bws[it], rcy, lrcy);
                limit_bandwidth(BW, (int) bws[it], rcz, lrcz);
                limit_bandwidth(BW, (int) bws[it], icx, licx);
                limit_bandwidth(BW, (int) bws[it], icy, licy);
                limit_bandwidth(BW, (int) bws[it], icz, licz);

                if (debugflag) fprintf(stderr,"Inv SPH xform..");
                get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT,
                                             lrcx, licx);
                get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT,
                                             lrcy, licy);
                get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT,
                                             lrcz, licz);

                if (debugflag) fprintf(stderr,"Resamp surf.\n");
                sample_sphere_from_sph(rdatax, rdatay, rdataz,
                                       polygons, n_triangles, reparam, BW);

                areas[it] = get_area_of_polygons(polygons, sph_areas[it]) /
                            orig_area;
                for (i = 0; i < polygons->n_items; i++) {
                        if (orig_areas[i] != 0)
                                sph_areas[it][i] /= orig_areas[i];
                }
                printf("bw = %d, area = %0.7f\n", (int) bws[it], areas[it]);

                if (debugflag) {
                        sprintf(str, "sph_%d.obj", (int) bws[it]);
                        if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                                       object) != OK)
                                exit(EXIT_FAILURE);

                        sprintf(str, "area_%d.txt", (int) bws[it]);
                        output_values_any_format(str, polygons->n_points,
                                                 sph_areas[it], TYPE_DOUBLE);
                }

                delete_polygons(polygons);
        }

        local_fd = (double *) malloc(sizeof(double) * reparam->n_points);
        get_localfd(reparam, bws, sph_areas, SPH_ITERS, local_fd, smoothflag);
        output_values_any_format(file, reparam->n_points, local_fd,
                                 TYPE_DOUBLE);
 
        if (debugflag) {
                if (output_values_any_format("fd_global.txt", SPH_ITERS, areas,
                                             TYPE_DOUBLE) != OK)
                        exit(EXIT_FAILURE);
        }

        fd = get_globalfd(bws, areas, SPH_ITERS);

	    free(rcx); free(rcy); free(rcz); 
        free(icx); free(icy); free(icz); 
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);
        
        return fd;
}

