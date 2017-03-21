/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */
#include <float.h>

#include "CAT_SurfaceIO.h"
#include "CAT_SPH.h"
#include "CAT_Intersect.h"
#include "CAT_Defect.h"
#include "CAT_Map.h"
#include "CAT_Surf.h"

#define DUMP_FILES 0


// not yet working (pointers!!!)
int
read_SPHxyz(char *file, int bandwidth, double **rcx, double **rcy, double **rcz,
            double **icx, double **icy, double **icz)
{
        FILE  *fp;
        char  line[256];
        int   bw2, n_dims, i;
        int   dataformat;

        /* read coefficients */
        if ((fp = fopen(file, "rb"))  == NULL) {
                fprintf(stderr, "Error opening file %s.\n", file);
                return(1);
        }

        fgets(line, 256, fp);
        if (strncmp(line, "SPH", 3)) {
                fprintf(stderr, "Wrong file format.\n");
                return(1);
        }
        fgets(line, 256, fp);
        sscanf(line, "%d %d %d", &bandwidth, &n_dims, &dataformat);

        if (dataformat !=  1) {
                fprintf(stderr, "Wrong dataformat: Data should be real and not complex.\n");
                return(1);
        }

        if (n_dims !=  3) {
                fprintf(stderr, "Data dimension should be 3.\n");
                return(1);
        }

        bw2 = bandwidth * bandwidth;

        *rcx = (double *) malloc(sizeof(double) * bw2);
        *rcy = (double *) malloc(sizeof(double) * bw2);
        *rcz = (double *) malloc(sizeof(double) * bw2);
        *icx = (double *) malloc(sizeof(double) * bw2);
        *icy = (double *) malloc(sizeof(double) * bw2);
        *icz = (double *) malloc(sizeof(double) * bw2);

        if (fread(*rcx, sizeof(double), bw2, fp) !=  bw2) {
                fprintf(stderr, "Error reading data.\n");
                return(1);
        }
        if (fread(*icx, sizeof(double), bw2, fp) !=  bw2) {
                fprintf(stderr, "Error reading data.\n");
                return(1);
        }
        if (fread(*rcy, sizeof(double), bw2, fp) !=  bw2) {
                fprintf(stderr, "Error reading data.\n");
                return(1);
        }
        if (fread(*icy, sizeof(double), bw2, fp) !=  bw2) {
                fprintf(stderr, "Error reading data.\n");
                return(1);
        }
        if (fread(*rcz, sizeof(double), bw2, fp) !=  bw2) {
                fprintf(stderr, "Error reading data.\n");
                return(1);
        }
        if (fread(*icz, sizeof(double), bw2, fp) !=  bw2) {
                fprintf(stderr, "Error reading data.\n");
                return(1);
        }
        
        return(fclose(fp));
}

int
write_SPHxyz(char *file, int bandwidth, double *rcx, double *rcy, double *rcz,
             double *icx, double *icy, double *icz)
{
        FILE  *fp;
        int i, bw2 = bandwidth * bandwidth;

        /* output coefficients */
        fp = fopen(file, "w");
        fprintf(fp, "SPH\n%d %d %d\n", bandwidth, 3, 1);
        for (i = 0; i < bw2; i++) fprintf(fp,"%g %g %g %g %g %g\n",rcx[i],icx[i],rcy[i],icy[i],rcz[i],icz[i]);
        fprintf(fp,"\n");
        
        return(fclose(fp));
}

void
sample_sphere_from_sph(double *rdatax, double *rdatay, double *rdataz,
                       polygons_struct *sphere, int n_triangles,
                       polygons_struct *reparam, int bandwidth)
{
        Point     centre, new_point;    
        double    H00, H01, H10, H11, valuex, valuey, valuez;
        double    u, v;
        int       i, size_map[2], bandwidth2, n_polygons;

        bandwidth2 = bandwidth*2;
    
        if (reparam  == NULL) {
                /* check tetrahedral topology */
                /* best areal distribution is achieved for 20 edges */
                n_polygons = n_triangles;
        
                while (n_polygons % 4  == 0)
                        n_polygons /=  4;

                if (n_polygons !=  5 && n_polygons !=  2 &&
                    n_polygons !=  6 && n_polygons !=  1) { /* 20, 8, 6, 4 */
                        fprintf(stderr,"Warning: Number of triangles %d is",
                                       n_triangles);
                        fprintf(stderr," not recommend because\ntetrahedral");
                        fprintf(stderr," topology is not optimal.\nPlease try");
                        fprintf(stderr," 20*(4*x) triangles (e.g. 81920).\n");
                }
    
                fill_Point(centre, 0.0, 0.0, 0.0);
                create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, n_triangles,
                                          sphere);
        } else {
                copy_polygons(reparam, sphere);
        }

        compute_polygon_normals(sphere);
        create_polygons_bintree(sphere, round((double) sphere->n_items *
                                              BINTREE_FACTOR));

        size_map[0] = bandwidth2;
        size_map[1] = bandwidth2;

        for (i = 0; i < sphere->n_points; i++) {
                point_to_uv(&sphere->points[i], &u, &v);

                /* interpolate points */
                double x0 = (double) (bandwidth2 - 1) * u;
                double y0 = (double) (bandwidth2 - 1) * v;
                int x = (int) x0;
                int y = (int) y0;
                double xp = x0 - x;
                double yp = y0 - y;
                double xm = 1.0 - xp;
                double ym = 1.0 - yp;
        
                H00 = rdatax[bound(x,   y,   size_map)];
                H01 = rdatax[bound(x,   y+1, size_map)];
                H10 = rdatax[bound(x+1, y,   size_map)];
                H11 = rdatax[bound(x+1, y+1, size_map)];

                valuex = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);
 
                H00 = rdatay[bound(x,   y,   size_map)];
                H01 = rdatay[bound(x,   y+1, size_map)];
                H10 = rdatay[bound(x+1, y,   size_map)];
                H11 = rdatay[bound(x+1, y+1, size_map)];

                valuey = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                H00 = rdataz[bound(x,   y,   size_map)];
                H01 = rdataz[bound(x,   y+1, size_map)];
                H10 = rdataz[bound(x+1, y,   size_map)];
                H11 = rdataz[bound(x+1, y+1, size_map)];

                valuez = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                fill_Point(new_point, valuex, valuey, valuez);
                sphere->points[i] = new_point;
        }

        compute_polygon_normals(sphere);
}


void
replaceSPH(int bandwidth, int bandwidth_limited,
           double *coeffs, double *coeffs_filter)
{
        int l, m, i;

        for (l = 0; l < bandwidth; l++) {
                for (m = -l; m < l+1; m++) {
                        i = seanindex(m, l, bandwidth);
                        if (l < bandwidth_limited) 
                                coeffs[i] = coeffs_filter[i];
                }
        }
}


void
shape_description(int bandwidth, double *rcx, double *rcy, double *rcz,
                  double *icx, double *icy, double *icz, double *shape_desc)
{

        int k, l, m;

        for (l = 0; l < bandwidth; l++) {
                shape_desc[l] = 0;
                for (m = -l; m < l+1; m++) {
                        k = seanindex(m, l, bandwidth);
                        shape_desc[l] +=  (rcx[k] * rcx[k]) +
                                         (rcy[k] * rcy[k]) +
                                         (rcz[k] * rcz[k]);
                        shape_desc[l] +=  (icx[k] * icx[k]) +
                                         (icy[k] * icy[k]) +
                                         (icz[k] * icz[k]);
                }
        }
}


void
butterworth_filter(int bandwidth, int bandwidth_limited,
                   double *coeffs_old, double *coeffs_new)
{
        int l, m, i;
        const double order = 128;
        double *coeffs;

        coeffs = (double *) malloc(sizeof(double) * bandwidth);

        for (l = 0; l < bandwidth; l++) {
                coeffs[l] = l;
                coeffs[l] /=  bandwidth_limited;
                coeffs[l] = 1/sqrt(1 + pow(coeffs[l], order));
       }

        for (l = 0; l < bandwidth; l++) {
                for (m = -l; m < l+1; m++) {
                        i = seanindex(m, l, bandwidth);
                        coeffs_new[i] = coeffs_old[i] * coeffs[l];
                }
        }
        free(coeffs);
}

/* bandpass boxcar filter */
void
bandpass_bandwidth(int bandwidth, int bw_lo, int bw_hi,
                   double *coeffs_old, double *coeffs_new)
{
        int l, m, i;

        for (l = 0; l < bandwidth; l++) {
                for (m = -l; m < l+1; m++) {
                        i = seanindex(m, l, bandwidth);
                        if (l >=  bw_lo && l <=  bw_hi && l < bandwidth) {
                                coeffs_new[i] = coeffs_old[i];
                        } else coeffs_new[i] = 0;
                }
        }
}


/* boxcar filter */
void
limit_bandwidth(int bandwidth, int bandwidth_limited,
                double *coeffs_old, double *coeffs_new)
{
        int l, m, i;

        for (l = 0; l < bandwidth; l++) {
                for (m = -l; m < l+1; m++) {
                        i = seanindex(m, l, bandwidth);
                        if (l < bandwidth_limited && l < bandwidth) {
                                coeffs_new[i] = coeffs_old[i];
                        } else coeffs_new[i] = 0;
                }
        }
}


void
get_sph_coeffs_of_realdata(double *rdata, int bandwidth, int dataformat,
                           double *rc, double *ic)
{
        double               *workspace, *weights, *idata;
        int                  rank, howmany_rank, n_objects, bandwidth2;
        fftw_plan            dctPlan;
        fftw_plan            fftPlan;
        fftw_iodim           dims[1], howmany_dimx[1];

        bandwidth2 = bandwidth*2;
    
        idata = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        weights = (double *) malloc(sizeof(double) * 4 * bandwidth);
        workspace = (double *) malloc(sizeof(double) * 
                    ((10 * (bandwidth*bandwidth)) + (24 * bandwidth)));

        /* set imag data to zero */
        memset(idata, 0, sizeof(double) * bandwidth2*bandwidth2);

        /* forward DCT */
        dctPlan = fftw_plan_r2r_1d(bandwidth2, weights, rdata,
                                   FFTW_REDFT10, FFTW_ESTIMATE);

        rank = 1;
        dims[0].n = bandwidth2;
        dims[0].is = 1;
        dims[0].os = bandwidth2;
        howmany_rank = 1;
        howmany_dimx[0].n = bandwidth2;
        howmany_dimx[0].is = bandwidth2;
        howmany_dimx[0].os = 1;

        /* forward fft */
        fftPlan = fftw_plan_guru_split_dft(rank, dims,
                                           howmany_rank, howmany_dimx,
                                           rdata, idata,
                                           workspace,
                                           workspace+(4*bandwidth*bandwidth),
                                           FFTW_ESTIMATE);

        /* now make the weights */
        makeweights(bandwidth, weights);

        FST_semi_fly(rdata, idata,
                     rc, ic,
                     bandwidth,
                     workspace,
                     dataformat,
                     bandwidth, /* Use seminaive for all orders */
                     &dctPlan,
                     &fftPlan,
                     weights);

        fftw_destroy_plan(fftPlan);
        fftw_destroy_plan(dctPlan);
        free(idata);
        free(workspace);
        free(weights);
}


void
get_realdata_from_sph_coeffs(double *rdata, int bandwidth, int dataformat,
                             double *rc, double *ic)
{
        double       *workspace, *weights, *idata;
        int          rank, howmany_rank, n_objects, bandwidth2;
        fftw_plan    idctPlan;
        fftw_plan    ifftPlan;
        fftw_iodim   dims[1], howmany_dimx[1];

        bandwidth2 = bandwidth*2;
    
        idata =     (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
        weights =   (double *) malloc(sizeof(double) * 4 * bandwidth);
        workspace = (double *) malloc(sizeof(double) * 
                    ((10 * (bandwidth*bandwidth)) + (24 * bandwidth)));

        /* inverse DCT */
        idctPlan = fftw_plan_r2r_1d(bandwidth2, weights, rc,
                                    FFTW_REDFT01, FFTW_ESTIMATE);

        rank = 1;
        dims[0].n = bandwidth2;
        dims[0].is = bandwidth2;
        dims[0].os = 1;
        howmany_rank = 1;
        howmany_dimx[0].n = bandwidth2;
        howmany_dimx[0].is = 1;
        howmany_dimx[0].os = bandwidth2;

        /* forward fft */
        ifftPlan = fftw_plan_guru_split_dft(rank, dims,
                                            howmany_rank, howmany_dimx,
                                            rc, ic, workspace,
                                            workspace+(4*bandwidth*bandwidth),
                                            FFTW_ESTIMATE);

        /* now make the weights */
        makeweights(bandwidth, weights);

        InvFST_semi_fly(rc, ic,
                        rdata, idata,
                        bandwidth,
                        workspace,
                        dataformat,
                        bandwidth, /* Use seminaive for all orders */
                        &idctPlan,
                        &ifftPlan);

        fftw_destroy_plan(ifftPlan);
        fftw_destroy_plan(idctPlan);
        free(idata);
        free(workspace);
        free(weights);
}


void
get_equally_sampled_coords_of_polygon(polygons_struct *polygons,
                                      polygons_struct *sphere,
                                      int bandwidth, double xcoord[],
                                      double ycoord[], double zcoord[])
{
        int i, j, x, y;
        double value, u, v, r;
        Point unit_point, on_sphere_point, new_point;
        Point poly_points[1000], poly_points_src[1000], scaled_point;
        int poly, size, ind, bandwidth2;
        double weights[1000], radius;
        object_struct *objects;
        polygons_struct *scaled_sphere;

        objects = create_object(POLYGONS);
        scaled_sphere = get_polygons_ptr(objects);
        copy_polygons(sphere, scaled_sphere);

        /*
         * Set centre and radius.  Make radius slightly smaller to get sure
         * that the inner side of handles will be found as nearest point
         * on the surface
         */
        translate_to_center_of_mass(scaled_sphere);
        for (i = 0; i < scaled_sphere->n_points; i++) 
                set_vector_length(&scaled_sphere->points[i], 1.0);

        create_polygons_bintree(scaled_sphere,
                                round((double) scaled_sphere->n_items *
                                      BINTREE_FACTOR));

        bandwidth2 = bandwidth*2;
    
        for (x = 0; x < bandwidth2; x++) {
                for (y = 0; y < bandwidth2; y++) {
                        u = ((double) x) / (double) (bandwidth2 - 1);
                        v = ((double) y) / (double) (bandwidth2 - 1);
  
                        uv_to_point(u, v, &unit_point);
            
                        poly = find_closest_polygon_point(&unit_point,
                                                          scaled_sphere,
                                                          &on_sphere_point);

                        size = get_polygon_points(scaled_sphere, poly,
                                                  poly_points_src);

                        get_polygon_interpolation_weights(&on_sphere_point,
                                                          size,
                                                          poly_points_src,
                                                          weights);

                        if (get_polygon_points(polygons, poly,
                                               poly_points) !=  size)
                                handle_internal_error("map_point_between_polygons");
                
                        fill_Point(new_point, 0.0, 0.0, 0.0);
                        for (i = 0; i < size; i++) {
                                SCALE_POINT(scaled_point, poly_points[i],
                                            weights[i]);
                                ADD_POINTS(new_point, new_point, scaled_point);
                        }
                        xcoord[x + (bandwidth2*y)] = Point_x(new_point);
                        ycoord[x + (bandwidth2*y)] = Point_y(new_point);
                        zcoord[x + (bandwidth2*y)] = Point_z(new_point);
                }
        }

}


void
get_equally_sampled_coords_holes(polygons_struct *polygons,
                                 polygons_struct *sphere, int *defects,
                                 int n_defects, int *holes, int bandwidth,
                                 double xcoord[], double ycoord[],
                                 double zcoord[], int force)
{
        int i, j, x, y, *bisected;
        double value, u, v, r;
        Point unit_point, on_sphere_point, new_point;
        Point poly_points[1000], poly_points_src[1000], scaled_point;
        int poly, size, ind, bandwidth2;
        double weights[1000], radius;
        object_struct **objects;
        polygons_struct *scaled_sphere;

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        scaled_sphere = get_polygons_ptr(*objects);
        initialize_polygons(scaled_sphere, WHITE, NULL);
        copy_polygons(sphere, scaled_sphere);
        translate_to_center_of_mass(scaled_sphere);

        /* cut holes and handles in half, saving the correct half 
         * for cutting/filling */
        bisected = (int *) malloc(sizeof(int) * polygons->n_points);
        bisect_defects(polygons, sphere, defects, n_defects, holes, bisected, (!force));
        if (DUMP_FILES) 
                output_values_any_format("orig_bisected.txt", polygons->n_points, 
                    bisected, TYPE_INTEGER);

        /* Set centre and radius */
        for (i = 0; i < scaled_sphere->n_points; i++) {
                switch (bisected[i]) {
                        case HANDLE: /* handle, cut */
                                r = 1.05;
                                break;
                        case HOLE: /* hole, fill */
                                r = 0.95;
                                break;
                        default:
                                r = 1.0;
                }

                set_vector_length(&scaled_sphere->points[i], r);
        }

        if (DUMP_FILES) {
                output_graphics_any_format("holesphere.obj", ASCII_FORMAT,
                                           1, objects, NULL);
        }

        create_polygons_bintree(scaled_sphere,
                                round((double) scaled_sphere->n_items *
                                      BINTREE_FACTOR));

        bandwidth2 = bandwidth*2;
    
        for (x = 0; x < bandwidth2; x++) {
                for (y = 0; y < bandwidth2; y++) {
                        u = ((double) x) / (double) (bandwidth2 - 1);
                        v = ((double) y) / (double) (bandwidth2 - 1);
                        uv_to_point(u, v, &unit_point);
            
                        poly = find_closest_polygon_point(&unit_point,
                                                          scaled_sphere,
                                                          &on_sphere_point);

                        size = get_polygon_points(scaled_sphere, poly,
                                                  poly_points_src);

                        get_polygon_interpolation_weights(&on_sphere_point,
                                                          size,
                                                          poly_points_src,
                                                          weights);

                        if (get_polygon_points(polygons, poly,
                                               poly_points) !=  size)
                                handle_internal_error("map_point_between_polygons");
                
                        fill_Point(new_point, 0.0, 0.0, 0.0);
                        for (i = 0; i < size; i++) {
                                SCALE_POINT(scaled_point, poly_points[i],
                                            weights[i]);
                                ADD_POINTS(new_point, new_point, scaled_point);
                        }
                        xcoord[x + (bandwidth2*y)] = Point_x(new_point);
                        ycoord[x + (bandwidth2*y)] = Point_y(new_point);
                        zcoord[x + (bandwidth2*y)] = Point_z(new_point);
                }

        }
        free(bisected);

}

object_struct **
create_equally_sampled_unit_sphere(int n_theta, int n_phi)
{
        object_struct **object;
        polygons_struct *sphere;
        Point unit_point;
        int p, theta, phi;
        double u, v, x, y;

        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        sphere = get_polygons_ptr(*object);
        initialize_polygons(sphere, WHITE, NULL);

        sphere->n_items = 2*n_theta*(n_phi - 1);
        sphere->n_points = (n_theta * (n_phi-1)) + 2;

        sphere->points = (Point *) malloc(sizeof(Point) * sphere->n_points);
        sphere->normals = (Vector *) malloc(sizeof(Vector) * sphere->n_points);
        sphere->indices = (int *) malloc(sizeof(int) * 3 * sphere->n_items);
        sphere->end_indices = (int *) malloc(sizeof(int) * sphere->n_items);

        p = 0;
        for (phi = 0; phi <=  n_phi; phi++) {
                v = (double) phi / (double) n_phi;
                for (theta = 0; theta < n_theta; theta++) {
                        u = (double) theta / (double) n_theta;

                        uv_to_point(u, v, &sphere->points[p]);
                        p++;

                        if (phi  == 0 || phi  == n_phi) break; // pole
                }
        }

        for (p = 0; p < sphere->n_items; p++)
                sphere->end_indices[p] = 3 * (p + 1); // all triangles

        /* north pole */
        p = 0;
        for (theta = 0; theta < n_theta; theta++) {
                sphere->indices[p++] = 0;
                sphere->indices[p++] = theta + 1;
                sphere->indices[p++] = theta + 2;
                if (theta  == n_theta-1) sphere->indices[p-1] -=  n_theta;
        }

        /* longitudinal rows */
        for (phi = 1; phi < n_phi-1; phi++) {
                for (theta = 0; theta < n_theta; theta++) {
                        sphere->indices[p++] = (phi-1)*n_theta + theta + 1;
                        sphere->indices[p++] = phi*n_theta + theta + 1;
                        sphere->indices[p++] = phi*n_theta + theta + 2;
                        if (theta  == n_theta-1) sphere->indices[p-1] -=  n_theta;

                        sphere->indices[p++] = (phi-1)*n_theta + theta + 1;
                        sphere->indices[p++] = phi*n_theta + theta + 2;
                        if (theta  == n_theta-1) sphere->indices[p-1] -=  n_theta;
                        sphere->indices[p++] = (phi-1)*n_theta + theta + 2;
                        if (theta  == n_theta-1) sphere->indices[p-1] -=  n_theta;
                }
        }

        /* south pole */
        for (theta = 0; theta < n_theta; theta++) {
                sphere->indices[p++] = (n_phi-2)*n_theta + theta + 1;
                sphere->indices[p++] = (n_phi-2)*n_theta + theta + 2;
                if (theta  == n_theta-1) sphere->indices[p-1] -=  n_theta;
                sphere->indices[p++] = (n_theta * (n_phi-1)) + 1;
        }
        compute_polygon_normals(sphere);
  
        return object;
}

/* estimate x,y,z position of index i in an array size sx,sxy = sx*sy... */
void 
ind2sub(int i, int *x, int *y, int sxy, int sy) {
           i = i % (sxy);
          *y = (int)floor( i / (double)sy );        
          *x = i % sy;
}

/* 2D laplace filter for continuous boundaries (e.g. longitude and latitude coordinates) */
double *
laplace2d(double *im, unsigned char *msk, int dimx, int dimy, double TH)
{

        int xy = dimx*dimy;
        double *L1, *L2;
        unsigned char *LN;
        const int NI[]  = { -1, 1, -dimx, dimx};  
        const int sN = sizeof(NI)/4;    
        int i, n;
        unsigned char * msk2;

        int u,v,nu,nv,ni,iter = 0, maxiter = 2000;
        double Nn, diff, maxdiffi, maxdiff = 1;
  
        if ( TH>= 0.5 || TH<0.000001 ) {
                fprintf(stderr,"ERROR:laplace2d: threshold must be >0.000001 and smaller than 0.5\n");
                return(NULL);
        }
    
        /* output data */
        L1  = (double *)malloc(xy * sizeof(double));
        L2  = (double *)malloc(xy * sizeof(double));
        LN  = (unsigned char *)malloc(xy * sizeof(unsigned char));
  
        msk2  = (unsigned char *)malloc(xy * sizeof(unsigned char));
  
        if (msk  == NULL)
                for (i = 0; i<xy; i++) msk2[i] = 1;
        else    for (i = 0; i<xy; i++) msk2[i] = msk[i];
        
        /* intitialisiation */
        for (i = 0; i<xy; i++) {
                if ( isnan(im[i]) ) L1[i] = FLT_MAX; else L1[i] = im[i];
                L2[i] = L1[i];
                LN[i] = msk2[i];
        }
                                        
        while (maxdiff > TH && iter < maxiter) {
                maxdiffi = 0; iter++;
                for (i = 0; i<xy; i++) {
                        if ( msk2[i] && LN[i] ) {  
                                ind2sub(i,&u,&v,xy,dimx);

                                /* read neighbor values */
                                L2[i] = 0; Nn = 0;
                                for (n = 0; n<sN; n++) {
                                        ni = i + NI[n]; 
                                        ind2sub(ni,&nu,&nv,xy,dimx);
                                        if (ni < 0) ni += dimx;
                                        if (ni >= xy) ni -= dimx;
                                        if (((ni<0) || (ni>= xy) || (L1[ni] == -FLT_MAX) || (L1[ni] == FLT_MAX)) == 0) {
                                                L2[i] = L2[i] + L1[ni];
                                                Nn++;
                                        }
                                }
                                if (Nn>0) L2[i]/= Nn; else L2[i] = L1[i];
        
                                diff  = fabs( L1[i] - L2[i] );
                                if ( diff>(TH/10) ) { 
                                        for (n = 0; n<sN; n++) {
                                                ni = i + NI[n]; 
                                                ind2sub(ni,&nu,&nv,xy,dimx);
                                                if (ni < 0) ni += dimx;
                                                if (ni >= xy) ni -= dimx;
                                                if (((ni<0) || (ni>= xy) || (L1[ni] == -FLT_MAX) || (L1[ni] == FLT_MAX) ) == 0) 
                                                        LN[ni] = 1; /* if i changes his neigbors has to be recalculated */
                                        }
                                }
      
                                LN[i] = 0;
                                if (maxdiffi<diff) maxdiffi = diff; 
                        }
                }
                maxdiff = maxdiffi;

                /* update of L1 */
                for (i = 0; i<xy; i++) L1[i] = L2[i];
        }
  
        free(L2);
        free(LN);
        free(msk2);

        return(L1);
}

/* 2D gradient magnitude for continuous boundaries (e.g. longitude and latitude coordinates) */
double *
gradient_magnitude(double *im, int dimx, int dimy)
{
        int xy = dimx*dimy;
        int x, y, x2, y2, ind;
        double * mag;
        double gx, gy;
        
        /* output data */
        mag  = (double *)malloc(xy * sizeof(double));
        memset(mag, 0, sizeof(double) * xy);

        for (y = 1; y<dimy; y++) {
                for (x = 1; x<dimx; x++) {
                        
                        /* deal with borders */
                        if (x>0) x2 = x; else x2 = x + dimx;
                        if (y>0) y2 = y; else y2 = y + dimy;
                        if (x<dimx-1) x2 = x; else x2 = x - dimx;
                        if (y<dimy-1) y2 = y; else y2 = y - dimy;
                        
		  		        ind = (y2*dimx) + x2;
		  		        if ((ind < 1) | (ind > xy-1)) continue;
		  		        gx = (im[ind+1] - im[ind-1])/2;
		  		        gy = (im[((y2+1)*dimx)+x2] - im[((y2-1)*dimx)+x2])/2;
		  		        mag[ind] = sqrt(gx*gx + gy*gy);       
                }
        }
        return(mag);
}

/* threshold 2D image */
unsigned char *
threshold_image(double *im, int dimx, int dimy, double threshold)
{

        int xy = dimx*dimy;
        int i;
        unsigned char * im_thresholded;
        
        /* output data */
        im_thresholded  = (unsigned char *)malloc(xy * sizeof(unsigned char));
        
        for (i = 0; i<xy; i++) {
                im_thresholded[i] = (im[i] < threshold) ? 0 : 1;
        }

        return(im_thresholded);
}
