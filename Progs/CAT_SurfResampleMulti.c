/*
 * CAT_SurfResampleMulti.c
 *
 * Flexible multi-surface resampler with optional smoothing:
 *  - Any number of units via repeated -unit key=value lists
 *  - Each unit: surf, src_sphere, trg_sphere, vals, (optional) fwhm
 *  - Global -fwhm applies to units that did not specify fwhm
 *  - Per-unit resampling (src_sphere -> trg_sphere)
 *  - Optional per-unit smoothing (smooth_heatkernel)
 *  - Concatenate all resampled meshes and values into ONE output GIfTI
 *  - No labels/annotations
 */

#include <bicpl.h>
#include <ParseArgv.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Resample.h"
#include "CAT_Smooth.h"
#include "CAT_Math.h"

/* ------------------ Unit description ------------------ */

typedef struct {
    char   *surf;
    char   *src_sphere;
    char   *trg_sphere;
    char   *vals;
    char   *mask;
    double  fwhm;          /* < 0.0  => not specified, take global default */

    object_struct   **src_obj;
    object_struct   **src_sphere_obj;
    object_struct   **trg_sphere_obj;
    object_struct   **resampled_obj;
    polygons_struct *src_poly;
    polygons_struct *src_sphere_poly;
    polygons_struct *trg_sphere_poly;

    double *vals_mask;
    double *vals_in;
    double *vals_out;
} unit_t;

static unit_t  *g_units = NULL;
static int      g_nunits = 0;

static char    *g_out = NULL;
static double   g_fwhm_default = 0.0;   /* NEW: global default FWHM */

/* ------------------ helpers ------------------ */

static void trim(char *s)
{
    char *p = s;
    size_t l = strlen(p);
    while (l && isspace((unsigned char)p[l-1])) p[--l] = 0;
    while (*p && isspace((unsigned char)*p)) ++p, --l;
    memmove(s, p, l + 1);
}

static char *strdup0(const char *s)
{
    if (!s) return NULL;
    char *d;
    ALLOC(d, strlen(s) + 1);
    strcpy(d, s);
    return d;
}

/* Parse "k1=v1,k2=v2,..." for -unit */
static int parse_unit_kv(char *arg, unit_t *u)
{
    u->surf = u->src_sphere = u->trg_sphere = u->vals = u->mask = NULL;
    u->fwhm = -1.0;   /* NEW: sentinel => take global default if still <0 */

    char *buf = strdup(arg);
    if (!buf) return 0;
    char *saveptr = NULL;
    for (char *tok = strtok_r(buf, ",", &saveptr); tok; tok = strtok_r(NULL, ",", &saveptr)) {
        trim(tok);
        char *eq = strchr(tok, '=');
        if (!eq) { fprintf(stderr, "Bad -unit token: '%s'\n", tok); free(buf); return 0; }
        *eq = 0;
        const char *key = tok;
        const char *val = eq + 1;
        if (!strcmp(key, "surf"))            u->surf        = strdup0(val);
        else if (!strcmp(key, "src_sphere")) u->src_sphere  = strdup0(val);
        else if (!strcmp(key, "trg_sphere")) u->trg_sphere  = strdup0(val);
        else if (!strcmp(key, "vals"))       u->vals        = strdup0(val);
        else if (!strcmp(key, "mask"))       u->mask        = strdup0(val);
        else if (!strcmp(key, "fwhm"))       u->fwhm        = atof(val);
        else {
            fprintf(stderr, "Unknown key in -unit: '%s'\n", key);
            free(buf);
            return 0;
        }
    }
    free(buf);

    if (!u->surf || !u->src_sphere || !u->trg_sphere) {
        fprintf(stderr, "unit is missing one of: surf, src_sphere, trg_sphere\n");
        return 0;
    }

    return 1;
}

/* ARGV_FUNC to push one unit */
static int add_unit_cb(char *dst, char *key, char *nextArg)
{
    (void)dst; (void)key;
    if (!nextArg) {
        fprintf(stderr, "-unit requires a key=value list.\n");
        return 1;
    }

    g_units = (unit_t *) realloc(g_units, sizeof(unit_t) * (g_nunits + 1));
    if (!g_units) {
        fprintf(stderr, "Out of memory.\n");
        return 1;
    }
    memset(&g_units[g_nunits], 0, sizeof(unit_t));

    if (!parse_unit_kv(nextArg, &g_units[g_nunits])) {
        return 1;
    }

    g_nunits++;
    return 0;
}

static void usage(const char *p)
{
    fprintf(stderr,
      "\nUsage:\n"
      "  %s \\\n"
      "     -unit surf=...,src_sphere=...,trg_sphere=...,vals=...,mask=...,fwhm=... \\\n"
      "     [-unit ...] \\\n"
      "     [-fwhm <float>] \\\n"
      "     -out combined.gii\n\n"
      "Repeat -unit as many times as you like (LH, RH, Cerebellum, ...).\n"
      "Per-unit fwhm overrides global -fwhm; if a unit omits fwhm, it takes the global value.\n"
      "All units are resampled independently (own target spheres), smoothed (optional),\n"
      "then concatenated and written as one GIfTI with values.\n\n",
      p);
}

/* Concatenate many polygons + values */
static object_struct **concat_many(object_struct **objs, double **vals, int nobjs,
                                   double **out_vals, int *out_n_points)
{
    int i, j;
    int total_pts = 0, total_items = 0, total_idx = 0;

    for (i = 0; i < nobjs; i++) {
        polygons_struct *p = get_polygons_ptr(objs[i]);
        total_pts   += p->n_points;
        total_items += p->n_items;
        total_idx   += (p->n_items > 0) ? p->end_indices[p->n_items - 1] : 0;
    }

    object_struct **out_list = (object_struct **) malloc(sizeof(object_struct *));
    out_list[0] = create_object(POLYGONS);
    polygons_struct *out = get_polygons_ptr(out_list[0]);
    out->normals = NULL;

    ALLOC(out->points, total_pts);
    ALLOC(out->end_indices, total_items);
    ALLOC(out->indices, total_idx);

    int p_off = 0;
    int it_off = 0;
    int idx_off = 0;
    int v_off = 0;

    if (out_vals) {
        ALLOC(*out_vals, total_pts);
    }

    for (i = 0; i < nobjs; i++) {
        polygons_struct *p = get_polygons_ptr(objs[i]);
        int npts   = p->n_points;
        int nitems = p->n_items;
        int nidx   = (nitems > 0) ? p->end_indices[nitems - 1] : 0;

        for (j = 0; j < npts; j++) {
            out->points[p_off + j]  = p->points[j];
        }

        for (j = 0; j < nitems; j++) {
            int in_start = (j == 0) ? 0 : p->end_indices[j - 1];
            int in_end   = p->end_indices[j];
            int nin      = in_end - in_start;

            int out_item = it_off + j;
            int out_start = (out_item == 0) ? 0 : out->end_indices[out_item - 1];
            out->end_indices[out_item] = out_start + nin;

            for (int k = 0; k < nin; k++) {
                out->indices[out_start + k] = p->indices[in_start + k] + p_off;
            }
        }

        if (out_vals && vals && vals[i]) {
            for (j = 0; j < npts; j++) {
                (*out_vals)[v_off + j] = vals[i][j];
            }
        }

        p_off  += npts;
        it_off += nitems;
        idx_off += nidx;
        v_off  += npts;
    }

    out->n_points = total_pts;
    out->n_items  = total_items;
    out->colour_flag = 0;
    out->bintree = NULL;

    if (out_n_points) *out_n_points = total_pts;
    return out_list;
}

/* ------------------ ArgvInfo table ------------------ */

static ArgvInfo argTable[] = {
    { "-unit", ARGV_FUNC,  (char *)add_unit_cb, NULL,
      "Add a unit: surf=...,src_sphere=...,trg_sphere=...,vals=...,mask=...,fwhm=..." },

    { "-fwhm", ARGV_FLOAT, (char *)1, (char *)&g_fwhm_default,
      "Global default FWHM for smoothing (applies to units without fwhm=). Default 0.0 (off)." },

    { "-out",  ARGV_STRING, (char *)1, (char *)&g_out,
      "Combined output GIfTI (.gii)" },

    { NULL, ARGV_END, NULL, NULL, NULL }
};

/* ------------------ main ------------------ */

int main(int argc, char *argv[])
{
    File_formats format = BINARY_FORMAT;
    int i, j, n_objects;
    double **vals_per_unit = NULL;

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (g_nunits < 1 || !g_out) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (filename_extension_matches(g_out, "gii")==0 && filename_extension_matches(g_out, "dat")==0) {
        fprintf(stderr, "Output must be a .gii or .dat file.\n");
        return EXIT_FAILURE;
    }

    /* Apply global default FWHM where per-unit not specified */
    for (i = 0; i < g_nunits; i++) {
        if (g_units[i].fwhm < 0.0)
            g_units[i].fwhm = g_fwhm_default;
    }

    ALLOC(vals_per_unit, g_nunits);

    /* read & resample each unit */
    for (i = 0; i < g_nunits; i++) {
        unit_t *u = &g_units[i];

        /* input */
        if (input_graphics_any_format(u->surf, &format, &n_objects, &u->src_obj) != OK ||
            n_objects != 1 || get_object_type(u->src_obj[0]) != POLYGONS) {
            fprintf(stderr, "[%d] %s must contain 1 polygons object.\n", i, u->surf);
            return EXIT_FAILURE;
        }
        u->src_poly = get_polygons_ptr(u->src_obj[0]);

        if (input_graphics_any_format(u->src_sphere, &format, &n_objects, &u->src_sphere_obj) != OK ||
            n_objects != 1 || get_object_type(u->src_sphere_obj[0]) != POLYGONS) {
            fprintf(stderr, "[%d] %s must contain 1 polygons object.\n", i, u->src_sphere);
            return EXIT_FAILURE;
        }
        u->src_sphere_poly = get_polygons_ptr(u->src_sphere_obj[0]);

        if (input_graphics_any_format(u->trg_sphere, &format, &n_objects, &u->trg_sphere_obj) != OK ||
            n_objects != 1 || get_object_type(u->trg_sphere_obj[0]) != POLYGONS) {
            fprintf(stderr, "[%d] %s must contain 1 polygons object.\n", i, u->trg_sphere);
            return EXIT_FAILURE;
        }
        u->trg_sphere_poly = get_polygons_ptr(u->trg_sphere_obj[0]);

        /* values */
        int n_vals = 0;
        if (u->vals) {
            if (input_values_any_format(u->vals, &n_vals, &u->vals_in) != OK) {
                fprintf(stderr, "[%d] Cannot read values in %s.\n", i, u->vals);
                return EXIT_FAILURE;
            }
            if (n_vals != u->src_sphere_poly->n_points) {
                fprintf(stderr, "[%d] #values (%d) != #points in source sphere (%d).\n",
                        i, n_vals, u->src_sphere_poly->n_points);
                return EXIT_FAILURE;
            }
            ALLOC(u->vals_out, u->trg_sphere_poly->n_points);
        }

        int n_vals_mask = 0;
        if (u->mask) {
            if (input_values_any_format(u->mask, &n_vals_mask, &u->vals_mask) != OK) {
                fprintf(stderr, "[%d] Cannot read values in %s.\n", i, u->mask);
                return EXIT_FAILURE;
            }
            if (n_vals_mask != u->trg_sphere_poly->n_points) {
                fprintf(stderr, "[%d] #values (%d) in mask != #points in target sphere (%d).\n",
                        i, n_vals_mask, u->trg_sphere_poly->n_points);
                return EXIT_FAILURE;
            }
        }

        /* resample */
        u->resampled_obj = resample_surface_to_target_sphere(
            u->src_poly, u->src_sphere_poly, u->trg_sphere_poly,
            u->vals_in, u->vals_out, 0 /* no labels */);

        /* smoothing (optional) */
        if (u->vals && u->fwhm > 0.0) {
            polygons_struct *p = get_polygons_ptr(u->resampled_obj[0]);
            
            /* Additionally set masked values to NaN */
            if (u->mask) {
                for (j = 0; j < n_vals_mask; j++) {
                    if (u->vals_mask[j] == 0) u->vals_out[j] = FNAN;
                }
                FREE(u->mask);
            }

            smooth_heatkernel(p, u->vals_out, u->fwhm);
        }

        vals_per_unit[i] = u->vals_out;
        
    }

    /* concatenate */
    object_struct **objs_only = (object_struct **) malloc(sizeof(object_struct *) * g_nunits);
    for (i = 0; i < g_nunits; i++) objs_only[i] = g_units[i].resampled_obj[0];

    double *vals_concat = NULL;
    int out_n_points = 0;
    object_struct **combined = concat_many(objs_only, vals_per_unit, g_nunits,
                                           &vals_concat, &out_n_points);

    if (output_graphics_any_format(g_out, format, 1, combined, vals_concat) != OK) {
        fprintf(stderr, "Failed to write %s.\n", g_out);
        return EXIT_FAILURE;
    }

    /* cleanup */
    if (vals_concat) FREE(vals_concat);
    FREE(objs_only);
    FREE(vals_per_unit);
    for (i = 0; i < g_nunits; i++) {
        if (g_units[i].vals_in)   FREE(g_units[i].vals_in);
        if (g_units[i].vals_out)  FREE(g_units[i].vals_out);
        if (g_units[i].vals_mask) FREE(g_units[i].vals_mask);
    }
    FREE(g_units);

    return EXIT_SUCCESS;
}
