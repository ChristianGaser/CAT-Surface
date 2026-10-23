#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bicpl.h>
#include "CAT_NiftiLib.h"
#include "CAT_BBReg.h"

/*
 * Boundary-based registration on a synthetic white surface.
 *
 * The surface is a folded, off-centre ellipsoid without symmetry.  The
 * moving volume is a BOLD-like image of it -- white matter 6% darker than
 * grey matter -- computed at the positions the true transform maps each
 * voxel back to, so no interpolation of test data enters, plus noise.  The
 * error of an estimate is the mean distance between the positions the
 * estimated and the true transform give the surface's vertices, which does
 * not depend on how the transform is parametrised.
 *
 * WM and GM are sampled 1.5 mm from the surface on either side, and there is
 * no third tissue beyond the grey matter, so that the true transform is the
 * minimum of the cost.  Otherwise it is not for scales: the contrast peaks
 * where the boundary lies midway between the two samples, which a rigid
 * transform cannot reach around a closed surface but a scale can -- with 2
 * and 1 mm the scales came out 1% too large, and tissue brighter than grey
 * matter 3 mm outside the boundary added another 1.5%.
 */

#define NX 64
#define NY 72
#define NZ 60
#define VOX 2.5

static const double shape_centre[3] = {4.0, -10.0, 12.0};

/** \brief Radius of the surface in the unit direction (ux, uy, uz). */
static double
surface_radius(double ux, double uy, double uz)
{
    double bx = ux - 0.6, by = uy - 0.3, bz = uz - 0.74;
    double r = 1.0 / sqrt(ux * ux / (55.0 * 55.0) + uy * uy / (68.0 * 68.0) +
                          uz * uz / (48.0 * 48.0));

    r += 3.0 * cos(7.0 * atan2(uy, ux)) * (1.0 - uz * uz);
    r += 2.0 * cos(5.0 * acos(uz));
    r += 4.0 * exp(-(bx * bx + by * by + bz * bz) / 0.1);
    return r;
}

/** \brief BOLD-like intensity at a fixed-space position: WM 94, GM 100. */
static double
intensity(double x, double y, double z)
{
    double dx = x - shape_centre[0], dy = y - shape_centre[1], dz = z - shape_centre[2];
    double len = sqrt(dx * dx + dy * dy + dz * dz), s;

    if (len < 1e-9)
        return 60.0;
    /* radial distance to the surface, negative inside */
    s = len - surface_radius(dx / len, dy / len, dz / len);
    return 94.0 + 6.0 / (1.0 + exp(-s / 0.6));
}

static double
gauss(unsigned long *state)
{
    double u1, u2;

    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    u1 = ((*state >> 11) + 0.5) / 9007199254740992.0;
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    u2 = ((*state >> 11) + 0.5) / 9007199254740992.0;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

/** \brief The white surface: a tetrahedral sphere pushed out to surface_radius(). */
static void
make_surface(polygons_struct *p)
{
    Point origin;
    int i;

    fill_Point(origin, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&origin, 1.0, 1.0, 1.0, 20480, p);
    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]), y = Point_y(p->points[i]);
        double z = Point_z(p->points[i]), len = sqrt(x * x + y * y + z * z), r;

        x /= len;
        y /= len;
        z /= len;
        r = surface_radius(x, y, z);
        fill_Point(p->points[i], shape_centre[0] + r * x, shape_centre[1] + r * y,
                   shape_centre[2] + r * z);
    }
    compute_polygon_normals(p);
}

/** \brief Header of the moving grid: VOX mm voxels centred on the origin. */
static nifti_image *
make_header(void)
{
    int dims[8] = {3, NX, NY, NZ, 1, 1, 1, 1}, n[3] = {NX, NY, NZ}, i, j;
    nifti_image *nii = nifti_make_new_nim(dims, DT_FLOAT32, 0);

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            nii->sto_xyz.m[i][j] = (i == j) ? (i < 3 ? (float)VOX : 1.0f) : 0.0f;
    for (i = 0; i < 3; i++)
        nii->sto_xyz.m[i][3] = (float)(-0.5 * (n[i] - 1) * VOX);
    nii->sto_ijk = nifti_mat44_inverse(nii->sto_xyz);
    nii->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    nii->dx = nii->dy = nii->dz = VOX;
    return nii;
}

/** \brief The moving image of a fixed -> moving transform m, with noise of sd 1. */
static float *
make_volume(const nifti_image *nii, const double m[16], unsigned long seed)
{
    float *vol = (float *)malloc(sizeof(float) * NX * NY * NZ);
    unsigned long state = seed;
    double inv[16];
    int x, y, z, r;

    CAT_BBReg_invert_affine(m, inv);
    for (z = 0; z < NZ; z++)
        for (y = 0; y < NY; y++)
            for (x = 0; x < NX; x++)
            {
                double w[3], f[3];

                for (r = 0; r < 3; r++)
                    w[r] = nii->sto_xyz.m[r][0] * x + nii->sto_xyz.m[r][1] * y +
                           nii->sto_xyz.m[r][2] * z + nii->sto_xyz.m[r][3];
                for (r = 0; r < 3; r++)
                    f[r] = inv[4 * r] * w[0] + inv[4 * r + 1] * w[1] +
                           inv[4 * r + 2] * w[2] + inv[4 * r + 3];
                vol[x + NX * (y + NY * z)] =
                    (float)(intensity(f[0], f[1], f[2]) + gauss(&state));
            }
    return vol;
}

/** \brief Mean distance (mm) between the vertex positions of two transforms. */
static double
vertex_error(const polygons_struct *p, const double a[16], const double b[16])
{
    double sum = 0.0;
    int i, r;

    for (i = 0; i < p->n_points; i++)
    {
        double v[3] = {Point_x(p->points[i]), Point_y(p->points[i]),
                       Point_z(p->points[i])}, d2 = 0.0;

        for (r = 0; r < 3; r++)
        {
            double e = (a[4 * r] - b[4 * r]) * v[0] + (a[4 * r + 1] - b[4 * r + 1]) * v[1] +
                       (a[4 * r + 2] - b[4 * r + 2]) * v[2] + (a[4 * r + 3] - b[4 * r + 3]);
            d2 += e * e;
        }
        sum += sqrt(d2);
    }
    return sum / p->n_points;
}

/**
 * \brief Register the volume of transform truth with dof degrees of freedom.
 *
 * As CAT_SurfBBReg does: the rigid grid search and Powell refinement, then
 * for 9 and 12 DOF the affine refinement.  Returns the mean vertex error.
 */
static double
register_dof(const polygons_struct *p, const double truth[16], int dof)
{
    nifti_image *nii = make_header();
    int dims[3] = {NX, NY, NZ};
    float *vol = make_volume(nii, truth, 7);
    CAT_SurfData sd;
    CAT_RigidParams p0 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, pr;
    double m[16], err;

    sd.surface = (polygons_struct *)p;
    sd.cortex_mask = NULL;
    sd.thickness = NULL;
    sd.gm_proj_frac = 0.0;

    CAT_BBReg_optimise(&p0, &pr, &sd, 1, vol, nii, dims, 1.5, 1.5, 0.5, 1,
                       4.0, 0.07, 2, 200, 1e-6, 0);
    if (dof > 6)
    {
        CAT_AffineParams a0, a;

        CAT_BBReg_affine_init(&pr, &sd, 1, &a0);
        CAT_BBReg_optimise_affine(&a0, &a, dof, &sd, 1, vol, nii, dims, 1.5, 1.5,
                                  0.5, 1, 200, 1e-6, 0);
        CAT_BBReg_affine_to_matrix(&a, m);
    }
    else
        CAT_BBReg_params_to_matrix(&pr, m);
    err = vertex_error(p, m, truth);

    free(vol);
    nifti_image_free(nii);
    return err;
}

/** \brief Without scales and shears the affine matrix is the rigid one; the inverse inverts. */
static void
test_affine_matrix(void)
{
    CAT_RigidParams r = {1.5, -2.0, 0.7, 0.03, -0.02, 0.05};
    CAT_AffineParams a;
    double mr[16], ma[16], inv[16], prod[16], dev = 0.0;
    int i, j, k;

    a.rigid = r;
    a.sx = a.sy = a.sz = 1.0;
    a.hxy = a.hxz = a.hyz = 0.0;
    a.centre[0] = 12.0;
    a.centre[1] = -30.0;
    a.centre[2] = 25.0;
    CAT_BBReg_params_to_matrix(&r, mr);
    CAT_BBReg_affine_to_matrix(&a, ma);
    MU_ASSERT("affine: unit scales give exactly the rigid matrix",
              memcmp(mr, ma, sizeof(mr)) == 0);

    /* the centre goes where the rigid core takes it */
    a.sx = 1.05;
    a.sy = 0.96;
    a.sz = 1.02;
    a.hxy = 0.03;
    a.hxz = -0.02;
    a.hyz = 0.01;
    CAT_BBReg_affine_to_matrix(&a, ma);
    for (i = 0; i < 3; i++)
    {
        double ea = ma[4 * i + 3], er = mr[4 * i + 3];
        for (j = 0; j < 3; j++)
        {
            ea += ma[4 * i + j] * a.centre[j];
            er += mr[4 * i + j] * a.centre[j];
        }
        dev = fmax(dev, fabs(ea - er));
    }
    MU_ASSERT("affine: scales and shears keep the centre", dev < 1e-9);

    CAT_BBReg_invert_affine(ma, inv);
    dev = 0.0;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
        {
            double s = 0.0;
            for (k = 0; k < 4; k++)
                s += ma[4 * i + k] * inv[4 * k + j];
            prod[4 * i + j] = s;
            dev = fmax(dev, fabs(prod[4 * i + j] - (i == j ? 1.0 : 0.0)));
        }
    MU_ASSERT("affine: the inverse inverts", dev < 1e-12);
}

/** \brief A rigid motion is recovered with 6 DOF. */
static void
test_rigid(void)
{
    polygons_struct p;
    CAT_RigidParams truth = {1.5, -1.0, 0.8, 0.03, -0.02, 0.025};
    double m[16], err;

    make_surface(&p);
    CAT_BBReg_params_to_matrix(&truth, m);
    err = register_dof(&p, m, 6);
    if (err >= 0.2)
        printf("  rigid: %.3f mm\n", err);
    MU_ASSERT("rigid: recovered within 0.2 mm", err < 0.2);
    delete_polygons(&p);
}

/**
 * \brief Scales need 9 DOF, and shears 12.
 *
 * Scales of 3-4% displace the cortex by up to 2.5 mm, shears of 2-3% by up
 * to 2 mm: more than a rigid transform can absorb.
 */
static void
test_affine(void)
{
    polygons_struct p;
    CAT_AffineParams a;
    double m[16], e6, e9, e12;

    make_surface(&p);
    a.rigid.tx = 1.0;
    a.rigid.ty = -0.8;
    a.rigid.tz = 0.6;
    a.rigid.rx = 0.02;
    a.rigid.ry = -0.015;
    a.rigid.rz = 0.01;
    a.sx = 1.04;
    a.sy = 0.97;
    a.sz = 1.03;
    a.hxy = a.hxz = a.hyz = 0.0;
    memcpy(a.centre, shape_centre, sizeof(a.centre));

    CAT_BBReg_affine_to_matrix(&a, m);
    e6 = register_dof(&p, m, 6);
    e9 = register_dof(&p, m, 9);
    if (e9 >= 0.3 || e6 <= 1.0)
        printf("  scales: 6 DOF %.3f mm, 9 DOF %.3f mm\n", e6, e9);
    MU_ASSERT("scales: 6 DOF cannot fit them", e6 > 1.0);
    MU_ASSERT("scales: 9 DOF recover them within 0.3 mm", e9 < 0.3);

    a.hxy = 0.03;
    a.hxz = -0.02;
    a.hyz = 0.025;
    CAT_BBReg_affine_to_matrix(&a, m);
    e9 = register_dof(&p, m, 9);
    e12 = register_dof(&p, m, 12);
    if (e12 >= 0.3 || e9 <= 2.0 * e12)
        printf("  shears: 9 DOF %.3f mm, 12 DOF %.3f mm\n", e9, e12);
    MU_ASSERT("shears: 9 DOF fit them worse than 12", e9 > 2.0 * e12);
    MU_ASSERT("shears: 12 DOF recover them within 0.3 mm", e12 < 0.3);

    delete_polygons(&p);
}

int main(void)
{
    MU_RUN_TEST(test_affine_matrix);
    MU_RUN_TEST(test_rigid);
    MU_RUN_TEST(test_affine);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
