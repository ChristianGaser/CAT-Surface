#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bicpl.h>
#include "CAT_WarpDemons.h"

/*
 * Spherical Demons registration (CAT_WarpDemonsRegister).
 *
 * A bumpy sphere is registered to itself through a copy of its spherical map
 * rotated by a known angle, so the correct output is the unrotated map and the
 * error of every vertex is known.  One cap of the sphere carries no bumps:
 * there the data term is weak, as it is around the medial wall of a brain
 * template, and only the global rotation decides where the cap ends up.
 */

#define RADIUS 100.0
#define N_TRIANGLES 20480
#define N_BUMPS 24
#define CAP_DEG 60.0 /* no bump centre within this angle of the cap centre */

static double bump_dir[N_BUMPS][3], bump_sigma[N_BUMPS], bump_height[N_BUMPS];

/** \brief Deterministic uniform numbers in [0, 1), the same on every platform. */
static double
next_uniform(unsigned int *state)
{
    *state = *state * 1103515245u + 12345u;
    return (double)((*state >> 8) & 0xFFFFFF) / (double)0x1000000;
}

/** \brief Bump centres spread over the sphere except the cap around -z. */
static void
make_bumps(void)
{
    unsigned int state = 2024;
    int n = 0;

    while (n < N_BUMPS)
    {
        double z = 2.0 * next_uniform(&state) - 1.0;
        double phi = 2.0 * PI * next_uniform(&state), s = sqrt(1.0 - z * z);

        if (acos(-z) < CAP_DEG * PI / 180.0)
            continue;
        bump_dir[n][0] = s * cos(phi);
        bump_dir[n][1] = s * sin(phi);
        bump_dir[n][2] = z;
        bump_sigma[n] = (8.0 + 8.0 * next_uniform(&state)) * PI / 180.0;
        bump_height[n] = 4.0 + 6.0 * next_uniform(&state);
        n++;
    }
}

/** \brief Regular sphere mesh of RADIUS around the origin. */
static void
make_sphere(polygons_struct *p)
{
    Point centre;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, N_TRIANGLES, p);
}

/** \brief The sphere with every vertex moved out by the bumps above it. */
static void
make_surface(const polygons_struct *sphere, polygons_struct *surface)
{
    int i, k;

    copy_polygons((polygons_struct *)sphere, surface);
    for (i = 0; i < surface->n_points; i++)
    {
        double u[3], r = RADIUS;

        for (k = 0; k < 3; k++)
            u[k] = Point_coord(sphere->points[i], k) / RADIUS;
        for (k = 0; k < N_BUMPS; k++)
        {
            double c = u[0] * bump_dir[k][0] + u[1] * bump_dir[k][1] + u[2] * bump_dir[k][2];
            double t = acos(fmin(1.0, fmax(-1.0, c)));
            r += bump_height[k] * exp(-0.5 * t * t / (bump_sigma[k] * bump_sigma[k]));
        }
        fill_Point(surface->points[i], r * u[0], r * u[1], r * u[2]);
    }
    compute_polygon_normals(surface);
}

/** \brief Rotate every vertex by deg degrees about a fixed oblique axis. */
static void
rotate_sphere(polygons_struct *p, double deg)
{
    double k[3] = {0.3, 0.8, 0.5}, n = sqrt(0.98), t = deg * PI / 180.0;
    double c = cos(t), s = sin(t);
    int i, j;

    for (j = 0; j < 3; j++)
        k[j] /= n;
    for (i = 0; i < p->n_points; i++)
    {
        double v[3], kv, kxv[3], w[3];

        for (j = 0; j < 3; j++)
            v[j] = Point_coord(p->points[i], j);
        kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
        kxv[0] = k[1] * v[2] - k[2] * v[1];
        kxv[1] = k[2] * v[0] - k[0] * v[2];
        kxv[2] = k[0] * v[1] - k[1] * v[0];
        for (j = 0; j < 3; j++)
            w[j] = v[j] * c + kxv[j] * s + k[j] * kv * (1.0 - c);
        fill_Point(p->points[i], w[0], w[1], w[2]);
    }
}

static int
cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/**
 * \brief Register the bumpy sphere to itself through a map rotated by deg.
 *
 * Returns the median angular error (degrees) of the vertices inside the cap
 * (within 40 deg of -z) and of the rest.
 */
static int
register_rotated(double deg, double *err_cap, double *err_rest)
{
    polygons_struct sphere, surface, src_sphere, trg_surface, trg_sphere, out;
    CAT_WarpDemonsOptions opt;
    double *cap, *rest;
    int i, n_cap = 0, n_rest = 0;

    make_sphere(&sphere);
    make_surface(&sphere, &surface);
    copy_polygons(&sphere, &src_sphere);
    rotate_sphere(&src_sphere, deg);
    copy_polygons(&surface, &trg_surface);
    copy_polygons(&sphere, &trg_sphere);

    CAT_WarpDemonsDefaults(&opt);
    opt.rot_max_degrees = 32.0;  /* covers the rotation; the search dominates the run time */
    if (CAT_WarpDemonsRegister(&surface, &src_sphere, &trg_surface, &trg_sphere,
                               &out, &opt) != OK)
        return 0;

    cap = (double *)malloc(sizeof(double) * sphere.n_points);
    rest = (double *)malloc(sizeof(double) * sphere.n_points);
    for (i = 0; i < sphere.n_points; i++)
    {
        double a[3], b[3], la = 0.0, lb = 0.0, d = 0.0;
        int k;

        for (k = 0; k < 3; k++)
        {
            a[k] = Point_coord(out.points[i], k);
            b[k] = Point_coord(sphere.points[i], k);
            la += a[k] * a[k];
            lb += b[k] * b[k];
            d += a[k] * b[k];
        }
        d = acos(fmin(1.0, fmax(-1.0, d / sqrt(la * lb)))) * 180.0 / PI;
        if (acos(-b[2] / sqrt(lb)) < 40.0 * PI / 180.0)
            cap[n_cap++] = d;
        else
            rest[n_rest++] = d;
    }
    qsort(cap, n_cap, sizeof(double), cmp_double);
    qsort(rest, n_rest, sizeof(double), cmp_double);
    *err_cap = cap[n_cap / 2];
    *err_rest = rest[n_rest / 2];

    free(cap);
    free(rest);
    delete_polygons(&sphere);
    delete_polygons(&surface);
    delete_polygons(&src_sphere);
    delete_polygons(&trg_surface);
    delete_polygons(&trg_sphere);
    delete_polygons(&out);
    return 1;
}

/**
 * \brief A rotated map is registered as accurately as an unrotated one.
 *
 * The rotation search finds the rotation; the warp only has to remove what
 * the search leaves.  If the rotation were part of the warp, the elastic prior
 * would shrink it wherever the data term is weak, and the featureless cap
 * would drift back towards its rotated start: 1.2 deg median error there, and
 * 0.5 deg over the bumps, against 0.12 / 0.08 deg without a rotation.
 */
static void
test_rotation_is_not_regularized(void)
{
    double cap, rest;

    make_bumps();
    MU_ASSERT("registration with a 20 deg rotation runs", register_rotated(20.0, &cap, &rest));
    if (cap >= 0.4 || rest >= 0.3)
        printf("  median error: cap %.3f deg, bumps %.3f deg\n", cap, rest);
    MU_ASSERT("rotated map: the bumps are registered", rest < 0.3);
    MU_ASSERT("rotated map: the featureless cap keeps the rotation", cap < 0.4);
}

int main(void)
{
    MU_RUN_TEST(test_rotation_is_not_regularized);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
