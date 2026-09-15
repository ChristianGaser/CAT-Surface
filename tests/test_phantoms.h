#ifndef TEST_PHANTOMS_H
#define TEST_PHANTOMS_H

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_NiftiLib.h"

/*
 * Analytic phantoms on a 0.5 mm grid.  Tissue boundaries are linear ramps one
 * voxel wide, so a label of 1.5 marks the CSF/GM boundary exactly and every
 * expected surface position has a closed-form value.
 */

#define VX 0.5

/**
 * \brief Linear partial-volume ramp: 0 outside, 1 inside, 0.5 at d = 0.
 *
 * \param d (in) signed distance to the boundary in mm, positive inside
 * \return ramp value in [0, 1]
 */
static inline double
ramp(double d)
{
    double v = d / VX + 0.5;
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}

/**
 * \brief Create a header-only NIfTI image with a diagonal sform.
 *
 * \param n    (in) dimensions
 * \param sign (in) +1 stores the axes RAS, -1 stores them LAS-like (all flipped)
 * \return new image; free with nifti_image_free()
 */
static inline nifti_image *
make_nim(const int n[3], int sign)
{
    nifti_image *nim = nifti_simple_init_nim();
    int k;

    nim->nx = nim->dim[1] = n[0];
    nim->ny = nim->dim[2] = n[1];
    nim->nz = nim->dim[3] = n[2];
    nim->nvox = (size_t)n[0] * n[1] * n[2];
    nim->dx = nim->pixdim[1] = VX;
    nim->dy = nim->pixdim[2] = VX;
    nim->dz = nim->pixdim[3] = VX;
    memset(&nim->sto_xyz, 0, sizeof(mat44));
    for (k = 0; k < 3; k++)
    {
        /* centre the grid on the origin */
        nim->sto_xyz.m[k][k] = (float)(sign * VX);
        nim->sto_xyz.m[k][3] = (float)(-sign * VX * (n[k] - 1) / 2.0);
    }
    nim->sto_xyz.m[3][3] = 1.0f;
    nim->sform_code = 1;
    return nim;
}

/**
 * \brief Fill a volume by evaluating f at the world position of each voxel.
 *
 * \param nim (in)  header defining the grid
 * \param f   (in)  phantom function of the world position
 * \return newly allocated volume
 */
static inline float *
fill_volume(const nifti_image *nim, double (*f)(double, double, double))
{
    float *vol = (float *)malloc(sizeof(float) * nim->nvox);
    int i, j, k;

    for (k = 0; k < nim->nz; k++)
        for (j = 0; j < nim->ny; j++)
            for (i = 0; i < nim->nx; i++)
            {
                double x = nim->sto_xyz.m[0][0] * i + nim->sto_xyz.m[0][3];
                double y = nim->sto_xyz.m[1][1] * j + nim->sto_xyz.m[1][3];
                double z = nim->sto_xyz.m[2][2] * k + nim->sto_xyz.m[2][3];
                vol[i + nim->nx * (j + nim->ny * k)] = (float)f(x, y, z);
            }
    return vol;
}

/* WM inside r = 6, pial boundary at r = 9 */
static inline double
sphere_phantom(double x, double y, double z)
{
    double r = sqrt(x * x + y * y + z * z);
    return 1.0 + ramp(9.0 - r) + ramp(6.0 - r);
}

/* Glued sulcus: GM continues to r = 11, with a valley down to 1.8 at r = 9 */
static inline double
valley_phantom(double x, double y, double z)
{
    double r = sqrt(x * x + y * y + z * z);
    return 1.0 + ramp(11.0 - r) + ramp(6.0 - r) - 0.2 * exp(-(r - 9.0) * (r - 9.0) / (2.0 * 0.16));
}

static inline void
make_sphere(polygons_struct *p, double radius)
{
    Point centre;
    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, radius, radius, radius, 5120, p);
    compute_polygon_normals(p);
}

static inline void
radius_stats(const polygons_struct *p, double *mean, double *min, double *max)
{
    int i;
    *mean = 0.0;
    *min = 1e30;
    *max = -1e30;
    for (i = 0; i < p->n_points; i++)
    {
        double r = sqrt(Point_x(p->points[i]) * Point_x(p->points[i]) +
                        Point_y(p->points[i]) * Point_y(p->points[i]) +
                        Point_z(p->points[i]) * Point_z(p->points[i]));
        *mean += r;
        *min = r < *min ? r : *min;
        *max = r > *max ? r : *max;
    }
    *mean /= p->n_points;
}

#endif /* TEST_PHANTOMS_H */
