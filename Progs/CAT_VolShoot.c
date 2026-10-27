/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/* Geodesic shooting of GM/WM maps onto a shooting template series, as CAT12/
 * CAT26 runs it (cat_main_registration.m, run_Shooting).  This is a thin CLI
 * wrapper around CAT_VolShoot().
 *
 * Input  : GM and WM maps (probabilities) and the five shooting templates.
 * Output : with -y the deformation in SPM's y_ format, with -iy its inverse,
 *          with -jac the Jacobian determinants.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_ArgvCheck.h"
#include "CAT_GeodesicShoot.h"

static CAT_VolShootOpts opts; /* defaults set in main() before parsing */

static char *y_file = NULL;
static char *iy_file = NULL;
static char *jac_file = NULL;
static char *affine_file = NULL;
static char *lesion_file = NULL;
static char *csf_file = NULL;
static double reg_vx = 0.0;
static int no_hydro = 0;
static int no_hotspot = 0;
static char *t1_file = NULL;
static char *t1_tpl_file = NULL;
static char *t1_mask_file = NULL;
static double t1_weight = 1.0;
static int t1_degree = 3;
static int t1_from = CAT_SHOOT_NTEMPLATES - 1;

static ArgvInfo argTable[] = {
    {"-y", ARGV_STRING, (char *)1, (char *)&y_file,
     "Write the deformation in SPM's y_ format: on the template grid, the input-image\n\t\t\tposition in mm of every template voxel (5-D float)."},
    {"-iy", ARGV_STRING, (char *)1, (char *)&iy_file,
     "Write the inverse deformation: on the input grid, the template position in mm."},
    {"-jac", ARGV_STRING, (char *)1, (char *)&jac_file,
     "Write the Jacobian determinant of the deformation (template grid)."},
    {"-regstr", ARGV_FLOAT, (char *)1, (char *)&opts.regstr,
     "CAT's optimised shooting: 1 (default) runs every iteration, (0,1) switches\n\t\t\tlevels early (CAT: 0.5); 4: SPM's default schedule."},
    {"-nits", ARGV_INT, (char *)1, (char *)&opts.nits,
     "Maximum number of Gauss-Newton iterations."},
    {"-lowres", ARGV_FLOAT, (char *)1, (char *)&opts.lowres,
     "Coarsest registration resolution in mm."},
    {"-reg-scale", ARGV_FLOAT, (char *)1, (char *)&opts.reg_scale,
     "Factor on CAT's regularisation (default 0.5, CAT: 1); smaller values allow\n\t\t\tmore local deformation."},
    {"-affine", ARGV_STRING, (char *)1, (char *)&affine_file,
     "Text file with the 4x4 input-to-template world matrix (CAT's res.Affine);\n\t\t\tdefault: identity, i.e. maps already in template space."},
    {"-csf", ARGV_STRING, (char *)1, (char *)&csf_file,
     "Experimental: CSF map as a third class; the templates must hold CSF as their\n\t\t\tthird volume."},
    {"-lesion", ARGV_STRING, (char *)1, (char *)&lesion_file,
     "Lesion mask on the input grid; the template replaces the tissue there."},
    {"-reg-vx", ARGV_FLOAT, (char *)1, (char *)&reg_vx,
     "Voxel size of the regulariser in mm; 0 (default) uses the template's\n\t\t\t(CAT always takes its TPM's 1.5 mm)."},
    {"-no-hydro", ARGV_CONSTANT, (char *)1, (char *)&no_hydro,
     "Do not model very large ventricles as a class of their own."},
    {"-no-hotspot", ARGV_CONSTANT, (char *)1, (char *)&no_hotspot,
     "Do not filter local Jacobian hot spots of the velocity."},
    {"-t1", ARGV_STRING, (char *)1, (char *)&t1_file,
     "Experimental: intensity image on the input grid, matched by sum of squares."},
    {"-t1-template", ARGV_STRING, (char *)1, (char *)&t1_tpl_file,
     "Template intensity on the template grid."},
    {"-t1-mask", ARGV_STRING, (char *)1, (char *)&t1_mask_file,
     "Weight of the intensity term on the template grid."},
    {"-t1-weight", ARGV_FLOAT, (char *)1, (char *)&t1_weight,
     "Curvature of the intensity term relative to the tissue term."},
    {"-t1-degree", ARGV_INT, (char *)1, (char *)&t1_degree,
     "Degree of the polynomial intensity mapping."},
    {"-t1-from", ARGV_INT, (char *)1, (char *)&t1_from,
     "First level (0..4) with the intensity term; default: the last only.  A level\n\t\t\tk < 4 reads <template>_L<k>.nii when that file exists."},
    {"-threads", ARGV_INT, (char *)1, (char *)&opts.n_threads,
     "Number of worker threads (default: the library's)."},
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&opts.verbose,
     "Print CAT's iteration log."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(const char *executable)
{
    CAT_VolShootOpts d;

    CAT_VolShootOptionsInit(&d);
    fprintf(stderr,
            "\nUsage: %s [options] <GM> <WM> <Template_0> [<Template_1> .. <Template_4>]\n"
            "          -y <y.nii> [-iy <iy.nii>] [-jac <jac.nii>]\n\n"
            "Geodesic shooting of GM and WM maps onto a series of five shooting\n"
            "templates, as CAT12/CAT26 runs it:\n\n"
            "    %s -y y_sub.nii rp1sub.nii rp2sub.nii Template_0_GS.nii\n\n"
            "With a single template its name must contain _0_; the other four are\n"
            "found by replacing that with _1_ .. _4_.  The maps are probabilities in\n"
            "[0, 1]; GM and WM are the first two volumes of every template.\n\n"
            "Options (defaults in brackets):\n"
            "  -regstr <v>     CAT's optimised shooting: 1 all iterations, 0.5 CAT's\n"
            "                  early level switching; 4 SPM's default schedule [%g]\n"
            "  -nits <n>       maximum Gauss-Newton iterations [%d]\n"
            "  -lowres <mm>    coarsest registration resolution [%g]\n"
            "  -reg-scale <f>  factor on CAT's regularisation (CAT: 1) [%g]\n"
            "  -affine <file>  4x4 input-to-template world matrix [identity]\n"
            "  -lesion <file>  lesion mask on the input grid\n"
            "  -reg-vx <mm>    voxel size of the regulariser [template's]\n"
            "  -no-hydro       no separate class for very large ventricles\n"
            "  -no-hotspot     no filtering of Jacobian hot spots\n"
            "  -threads <n>    worker threads\n"
            "  -verbose        print the iteration log\n\n",
            executable, executable, d.regstr, d.nits, d.lowres, d.reg_scale);
}

/* sto_xyz of a header as a row-major double matrix */
static void
header_matrix(const nifti_image *nii, double M[16])
{
    int r, c;
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            M[4 * r + c] = (double)nii->sto_xyz.m[r][c];
}

static int
read_matrix(const char *file, double M[16])
{
    FILE *fp = fopen(file, "r");
    int i;
    if (!fp) return 1;
    for (i = 0; i < 16; i++)
        if (fscanf(fp, "%lf", &M[i]) != 1)
        {
            fclose(fp);
            return 1;
        }
    fclose(fp);
    return 0;
}

/* A float map on the grid of ref with the given number of volumes. */
static float *
read_map(const char *file, int dims[3], int *nvol, nifti_image **hdr)
{
    float *data = NULL;
    nifti_image *nii = read_nifti_float(file, &data, 0);
    if (!nii)
    {
        fprintf(stderr, "Error reading %s.\n", file);
        exit(EXIT_FAILURE);
    }
    dims[0] = nii->nx;
    dims[1] = nii->ny;
    dims[2] = nii->nz;
    /* read_nifti_float() leaves nvox at one volume; count them from the dims */
    if (nvol)
        *nvol = (nii->nt > 1 ? nii->nt : 1) * (nii->nu > 1 ? nii->nu : 1) *
                (nii->nv > 1 ? nii->nv : 1) * (nii->nw > 1 ? nii->nw : 1);
    if (hdr) *hdr = nii;
    else nifti_image_free(nii);
    return data;
}

/* Write a deformation in SPM's 5-D layout (x, y, z, 1, 3) on the grid of
 * grid: the field holds voxel coordinates (0-based) of the other grid, which
 * vmat turns into mm. */
static int
write_deformation(const char *file, const float *vox, const int dims[3],
                  const double vmat[16], const nifti_image *grid)
{
    nifti_image *nim = nifti_copy_nim_info(grid);
    size_t n = (size_t)dims[0] * dims[1] * dims[2], i;
    float *mm;
    int j;

    if (!nim) return 1;
    mm = (float *)malloc(sizeof(float) * 3 * n);
    if (!mm)
    {
        nifti_image_free(nim);
        return 1;
    }
    for (i = 0; i < n; i++)
    {
        double a = vox[i], b = vox[n + i], c = vox[2 * n + i];
        for (j = 0; j < 3; j++)
            mm[j * n + i] = (float)(vmat[4 * j] * a + vmat[4 * j + 1] * b +
                                    vmat[4 * j + 2] * c + vmat[4 * j + 3]);
    }
    nim->ndim = nim->dim[0] = 5;
    nim->nx = nim->dim[1] = dims[0];
    nim->ny = nim->dim[2] = dims[1];
    nim->nz = nim->dim[3] = dims[2];
    nim->nt = nim->dim[4] = 1;
    nim->nu = nim->dim[5] = 3;
    nim->nv = nim->dim[6] = 1;
    nim->nw = nim->dim[7] = 1;
    nim->nvox = 3 * n;
    nim->datatype = DT_FLOAT32;
    nim->nbyper = 4;
    nim->scl_slope = 1.0f;
    nim->scl_inter = 0.0f;
    nim->intent_code = NIFTI_INTENT_VECTOR;
    if (nifti_set_filenames(nim, file, 0, 1))
    {
        free(mm);
        nifti_image_free(nim);
        return 1;
    }
    nim->data = mm;
    nifti_image_write(nim);
    nifti_image_free(nim); /* frees mm as nim->data */
    return 0;
}

/* Template_0 -> Template_k by replacing the last "_0_" */
static char *
template_name(const char *first, int k)
{
    const char *p = NULL, *q = first;
    char *name;
    while ((q = strstr(q, "_0_")) != NULL) p = q++;
    if (!p) return NULL;
    name = (char *)malloc(strlen(first) + 1);
    strcpy(name, first);
    name[p - first + 1] = (char)('0' + k);
    return name;
}

int
main(int argc, char *argv[])
{
    const float *cls[3], *tpl[CAT_SHOOT_NTEMPLATES];
    float *gm, *wm, *lesion = NULL, *y, *iy = NULL, *jac = NULL;
    nifti_image *hdr_in, *hdr_tpl = NULL;
    double mat[16], tmat[16], Aff[16];
    int dims[3], d2[3], tdims[3], k, nvol;
    char *tnames[CAT_SHOOT_NTEMPLATES];
    CAT_VolShootInfo info;
    CAT_VolShootIntensity ich;
    size_t nt;

    CAT_VolShootOptionsInit(&opts);
    if (ParseArgv(&argc, argv, argTable, 0) ||
        cat_check_unknown_options(argc, argv, argv[0]) ||
        (argc != 4 && argc != 8) ||
        (!y_file && !iy_file && !jac_file))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    opts.hydro = !no_hydro;
    opts.hotspot = !no_hotspot;

    if (argc == 4)
    {
        for (k = 0; k < CAT_SHOOT_NTEMPLATES; k++)
            if (!(tnames[k] = template_name(argv[3], k)))
            {
                fprintf(stderr, "Error: %s has no _0_ to derive the other templates from.\n",
                        argv[3]);
                exit(EXIT_FAILURE);
            }
    }
    else
        for (k = 0; k < CAT_SHOOT_NTEMPLATES; k++) tnames[k] = argv[3 + k];

    gm = read_map(argv[1], dims, NULL, &hdr_in);
    wm = read_map(argv[2], d2, NULL, NULL);
    if (d2[0] != dims[0] || d2[1] != dims[1] || d2[2] != dims[2])
    {
        fprintf(stderr, "Error: %s and %s differ in size.\n", argv[1], argv[2]);
        exit(EXIT_FAILURE);
    }
    header_matrix(hdr_in, mat);
    cls[0] = gm;
    cls[1] = wm;
    if (csf_file)
    {
        cls[2] = read_map(csf_file, d2, NULL, NULL);
        if (d2[0] != dims[0] || d2[1] != dims[1] || d2[2] != dims[2])
        {
            fprintf(stderr, "Error: %s and %s differ in size.\n", argv[1], csf_file);
            exit(EXIT_FAILURE);
        }
    }

    for (k = 0; k < CAT_SHOOT_NTEMPLATES; k++)
    {
        nifti_image *h;
        tpl[k] = read_map(tnames[k], d2, &nvol, &h);
        if (nvol < (csf_file ? 3 : 2))
        {
            fprintf(stderr, "Error: %s needs GM, WM%s as its first volumes.\n", tnames[k],
                    csf_file ? " and CSF" : "");
            exit(EXIT_FAILURE);
        }
        if (k == 0)
        {
            memcpy(tdims, d2, sizeof(tdims));
            header_matrix(h, tmat);
            hdr_tpl = h;
        }
        else
        {
            if (d2[0] != tdims[0] || d2[1] != tdims[1] || d2[2] != tdims[2])
            {
                fprintf(stderr, "Error: the templates differ in size.\n");
                exit(EXIT_FAILURE);
            }
            nifti_image_free(h);
        }
    }

    if (affine_file && read_matrix(affine_file, Aff))
    {
        fprintf(stderr, "Error reading the 4x4 matrix in %s.\n", affine_file);
        exit(EXIT_FAILURE);
    }
    if (lesion_file)
    {
        lesion = read_map(lesion_file, d2, NULL, NULL);
        if (d2[0] != dims[0] || d2[1] != dims[1] || d2[2] != dims[2])
        {
            fprintf(stderr, "Error: %s is not on the grid of %s.\n", lesion_file, argv[1]);
            exit(EXIT_FAILURE);
        }
    }

    nt = (size_t)tdims[0] * tdims[1] * tdims[2];
    y = (float *)malloc(sizeof(float) * 3 * nt);
    if (iy_file) iy = (float *)malloc(sizeof(float) * 3 * (size_t)dims[0] * dims[1] * dims[2]);
    if (jac_file) jac = (float *)malloc(sizeof(float) * nt);
    if (!y || (iy_file && !iy) || (jac_file && !jac))
    {
        fprintf(stderr, "Error: out of memory.\n");
        exit(EXIT_FAILURE);
    }

    if (t1_file)
    {
        float *img, *timg, *tmask = NULL;
        if (!t1_tpl_file)
        {
            fprintf(stderr, "Error: -t1 needs -t1-template.\n");
            exit(EXIT_FAILURE);
        }
        img = read_map(t1_file, d2, NULL, NULL);
        if (d2[0] != dims[0] || d2[1] != dims[1] || d2[2] != dims[2])
        {
            fprintf(stderr, "Error: %s is not on the grid of %s.\n", t1_file, argv[1]);
            exit(EXIT_FAILURE);
        }
        timg = read_map(t1_tpl_file, d2, NULL, NULL);
        if (d2[0] != tdims[0] || d2[1] != tdims[1] || d2[2] != tdims[2])
        {
            fprintf(stderr, "Error: %s is not on the template grid.\n", t1_tpl_file);
            exit(EXIT_FAILURE);
        }
        if (t1_mask_file)
        {
            tmask = read_map(t1_mask_file, d2, NULL, NULL);
            if (d2[0] != tdims[0] || d2[1] != tdims[1] || d2[2] != tdims[2])
            {
                fprintf(stderr, "Error: %s is not on the template grid.\n", t1_mask_file);
                exit(EXIT_FAILURE);
            }
        }
        memset(&ich, 0, sizeof(ich));
        ich.img = img;
        for (k = 0; k < CAT_SHOOT_NTEMPLATES; k++)
        {
            char *lname, *dot;
            FILE *fp;
            ich.tpl[k] = (k >= t1_from) ? timg : NULL;
            if (k < t1_from || k == CAT_SHOOT_NTEMPLATES - 1) continue;
            /* <template>_L<k>.nii: the template of a coarser level, if given */
            lname = (char *)malloc(strlen(t1_tpl_file) + 8);
            strcpy(lname, t1_tpl_file);
            dot = strstr(lname, ".nii");
            if (dot)
            {
                char ext[16];
                strncpy(ext, strstr(t1_tpl_file, ".nii"), sizeof(ext) - 1);
                ext[sizeof(ext) - 1] = '\0';
                sprintf(dot, "_L%d%s", k, ext);
                if ((fp = fopen(lname, "r")) != NULL)
                {
                    fclose(fp);
                    ich.tpl[k] = read_map(lname, d2, NULL, NULL);
                    if (d2[0] != tdims[0] || d2[1] != tdims[1] || d2[2] != tdims[2])
                    {
                        fprintf(stderr, "Error: %s is not on the template grid.\n", lname);
                        exit(EXIT_FAILURE);
                    }
                    if (opts.verbose) printf("level %d: intensity template %s\n", k, lname);
                }
            }
            free(lname);
        }
        ich.mask = tmask;
        ich.weight = t1_weight;
        ich.degree = t1_degree;
    }

    k = CAT_VolShootEx(cls, csf_file ? 3 : 2, dims, mat, affine_file ? Aff : NULL, tpl, tdims, tmat,
                       reg_vx, lesion, t1_file ? &ich : NULL, &opts, y, iy, jac, &info);
    if (k)
    {
        fprintf(stderr, "Error: shooting failed (%s).\n",
                (k == 2) ? "out of memory" : "invalid input");
        exit(EXIT_FAILURE);
    }
    if (opts.verbose)
        printf("%d iterations, Jacobian determinants %.3f .. %.3f\n", info.iterations,
               info.det_min, info.det_max);

    if (y_file && write_deformation(y_file, y, tdims, mat, hdr_tpl))
        exit(EXIT_FAILURE);
    if (iy_file && write_deformation(iy_file, iy, dims, tmat, hdr_in))
        exit(EXIT_FAILURE);
    if (jac_file)
    {
        double vox[3];
        vox[0] = hdr_tpl->dx;
        vox[1] = hdr_tpl->dy;
        vox[2] = hdr_tpl->dz;
        if (!write_nifti_float(jac_file, jac, DT_FLOAT32, 0.0, tdims, vox, hdr_tpl))
            exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
