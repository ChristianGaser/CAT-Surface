/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <fftw3.h>
#include  "s2kit10/makeweights.h"
#include  "s2kit10/FST_semi_fly.h"
#define BINTREE_FACTOR   0.5

int read_SPHxyz(
    char *filename,
    int *bandwidth,
    double *rcoeffsx,
    double *rcoeffsy,
    double *rcoeffsz,
    double *icoeffsx,
    double *icoeffsy,
    double *icoeffsz);

int write_SPHxyz(
    char *filename,
    int bandwidth,
    double *rcoeffsx,
    double *rcoeffsy,
    double *rcoeffsz,
    double *icoeffsx,
    double *icoeffsy,
    double *icoeffsz);

void sampe_sphere_from_sph(
    double *rdatax,
    double *rdatay,
    double *rdataz,
    polygons_struct *polygons_sphere,
    int n_triangles,
    int bandwidth);

void  replaceSPH(
    int bandwidth,
    int bandwidth_limited,
    double *coeffs,
    double *coeffs_filter);
    
void shape_description(
    int bandwidth, 
    double *rcoeffsx, 
    double *rcoeffsy,
    double *rcoeffsz,
    double *icoeffsx, 
    double *icoeffsy, 
    double *icoeffsz, 
    double *shape_descriptor);

void  limit_bandwidth(
    int bandwidth,
    int bandwidth_limited,
    double *coeffs_old,
    double *coeffs_new);

void get_sph_coeffs_of_realdata(
    double     *rdata,
    int        bandwidth,
    int        dataformat,
    double     *rcoeffs,
    double     *icoeffs);

void get_realdata_from_sph_coeffs(
    double     *rdata,
    int        bandwidth,
    int        dataformat,
    double     *rcoeffs,
    double     *icoeffs);

void get_equally_sampled_coords_of_polygon(
    polygons_struct *polygons,
    polygons_struct *polygons_sphere,
    int             bandwidth,
    double          xcoord[],
    double          ycoord[],
    double          zcoord[]);
