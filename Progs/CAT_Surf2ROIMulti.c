/*
 * CAT_Surf2ROIMulti.c
 *
 * Resample annotation files onto target spheres (label interpolation) and
 * accumulate per-ROI means of a values file defined on the target sphere.
 *
 * Usage (repeat -unit per hemisphere/region):
 *   CAT_Surf2ROIMulti \
 *     -unit src_sphere=...,trg_sphere=...,annot=...,vals=... \
 *     [-unit ...] \
 *     -out output.json
 *
 * Inputs mirror CAT_SurfResampleMulti but specialize on annotations:
 *   src_sphere  : source sphere matching the annotation vertices
 *   trg_sphere  : target sphere to resample onto
 *   annot       : annotation file (e.g., *.annot)
 *   vals        : values on the target sphere (e.g., lh.thickness)
 *   hemi        : optional (lh|rh|left|right); guessed from paths if omitted
 *
 * Output JSON structure:
 * {
 *   "lh": [ { "id": <int>, "name": "<label>", "mean": <double>, "n": <int> }, ... ],
 *   "rh": [ ... ],
 *   "unknown": [ ... ]
 * }
 */

#include <bicpl.h>
#include <ParseArgv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Resample.h"
#include "CAT_SafeAlloc.h"

typedef enum {
    HEMI_UNKNOWN = 0,
    HEMI_LH      = 1,
    HEMI_RH      = 2
} hemi_t;

typedef struct {
    char *src_sphere;
    char *trg_sphere;
    char *annot;
    char *vals;
    hemi_t hemi;
} unit_t;

typedef struct {
    int id;
    char *name;
    double sum;
    int count;
} roi_stat_t;

static unit_t *g_units = NULL;
static int g_nunits = 0;
static char *g_out = NULL;

static void trim(char *s)
{
    char *p = s;
    size_t l = strlen(p);
    while (l && isspace((unsigned char)p[l-1])) p[--l] = 0;
    while (*p && isspace((unsigned char)*p)) ++p, --l;
    memmove(s, p, l + 1);
}

static char *dup_or_die(const char *s)
{
    if (!s) return NULL;
    char *d = strdup(s);
    if (!d) {
        fprintf(stderr, "Out of memory.\n");
        exit(EXIT_FAILURE);
    }
    return d;
}

/* Parse "k=v,..." for -unit */
static hemi_t parse_hemi_value(const char *val)
{
    if (!val) return HEMI_UNKNOWN;
    if (!strcasecmp(val, "lh") || !strcasecmp(val, "left")) return HEMI_LH;
    if (!strcasecmp(val, "rh") || !strcasecmp(val, "right")) return HEMI_RH;
    return HEMI_UNKNOWN;
}

static hemi_t guess_hemi_from_path(const char *p)
{
    if (!p) return HEMI_UNKNOWN;
    if (strstr(p, "lh.")) return HEMI_LH;
    if (strstr(p, "rh.")) return HEMI_RH;
    if (strstr(p, "_lh")) return HEMI_LH;
    if (strstr(p, "_rh")) return HEMI_RH;
    if (strstr(p, "/lh")) return HEMI_LH;
    if (strstr(p, "/rh")) return HEMI_RH;
    return HEMI_UNKNOWN;
}

static hemi_t guess_hemi_from_unit(const unit_t *u)
{
    hemi_t h;
    if (!u) return HEMI_UNKNOWN;
    if (u->hemi != HEMI_UNKNOWN) return u->hemi;
    h = guess_hemi_from_path(u->annot);
    if (h != HEMI_UNKNOWN) return h;
    h = guess_hemi_from_path(u->vals);
    if (h != HEMI_UNKNOWN) return h;
    h = guess_hemi_from_path(u->src_sphere);
    if (h != HEMI_UNKNOWN) return h;
    h = guess_hemi_from_path(u->trg_sphere);
    return h;
}

static int parse_unit(char *arg, unit_t *u)
{
    u->src_sphere = u->trg_sphere = u->annot = u->vals = NULL;
    u->hemi = HEMI_UNKNOWN;

    char *buf = dup_or_die(arg);
    char *tok = strtok(buf, ",");
    while (tok) {
        trim(tok);
        char *eq = strchr(tok, '=');
        if (!eq) {
            fprintf(stderr, "Bad -unit token: '%s'\n", tok);
            free(buf);
            return 0;
        }
        *eq = 0;
        const char *key = tok;
        const char *val = eq + 1;

        if      (!strcmp(key, "src_sphere")) u->src_sphere = dup_or_die(val);
        else if (!strcmp(key, "trg_sphere")) u->trg_sphere = dup_or_die(val);
        else if (!strcmp(key, "annot"))      u->annot      = dup_or_die(val);
        else if (!strcmp(key, "vals"))       u->vals       = dup_or_die(val);
        else if (!strcmp(key, "hemi"))       u->hemi       = parse_hemi_value(val);
        else {
            fprintf(stderr, "Unknown key in -unit: '%s'\n", key);
            free(buf);
            return 0;
        }

        tok = strtok(NULL, ",");
    }

    free(buf);

    if (!u->src_sphere || !u->trg_sphere || !u->annot || !u->vals) {
        fprintf(stderr, "unit is missing one of: src_sphere, trg_sphere, annot, vals\n");
        return 0;
    }
    return 1;
}

/* ARGV callback to append a unit */
static int add_unit_cb(char *dst, char *key, char *nextArg)
{
    (void)dst; (void)key;
    if (!nextArg) {
        fprintf(stderr, "-unit requires a key=value list.\n");
        return 1;
    }

    unit_t *new_units = (unit_t *) realloc(g_units, sizeof(unit_t) * (g_nunits + 1));
    if (!new_units) {
        fprintf(stderr, "Out of memory.\n");
        return 1;
    }
    g_units = new_units;
    memset(&g_units[g_nunits], 0, sizeof(unit_t));

    if (!parse_unit(nextArg, &g_units[g_nunits])) {
        return 1;
    }

    g_nunits++;
    return 0;
}

static ArgvInfo argTable[] = {
            { "-unit", ARGV_FUNC,  (char *)add_unit_cb, NULL,
                "Add a unit: src_sphere=...,trg_sphere=...,annot=...,vals=...[,hemi=lh|rh]" },
    { "-out",  ARGV_STRING, (char *)1, (char *)&g_out,
      "Output JSON file" },
    { NULL, ARGV_END, NULL, NULL, NULL }
};

static const char *name_for_id(ATABLE *atable, int n_labels, int id)
{
    int i;
    for (i = 0; i < n_labels; i++) {
        if (atable[i].annotation == id) return atable[i].name;
    }
    return "unknown";
}

static void add_sum(roi_stat_t **arr, int *count, int *cap, int id, const char *name, double value)
{
    int i;
    for (i = 0; i < *count; i++) {
        if ((*arr)[i].id == id) {
            (*arr)[i].sum += value;
            (*arr)[i].count += 1;
            return;
        }
    }

    if (*count == *cap) {
        int new_cap = (*cap == 0) ? 8 : (*cap * 2);
        roi_stat_t *tmp = (roi_stat_t *) realloc(*arr, sizeof(roi_stat_t) * new_cap);
        if (!tmp) {
            fprintf(stderr, "Out of memory.\n");
            exit(EXIT_FAILURE);
        }
        *arr = tmp;
        *cap = new_cap;
    }

    (*arr)[*count].id = id;
    (*arr)[*count].name = dup_or_die(name);
    (*arr)[*count].sum = value;
    (*arr)[*count].count = 1;
    (*count)++;
}

static void write_roi_array(FILE *fp, const char *key, roi_stat_t *arr, int count, int indent, int trailing_comma)
{
    int i;
    fprintf(fp, "%*s\"%s\": [\n", indent, "", key);
    for (i = 0; i < count; i++) {
        double mean = (arr[i].count > 0) ? (arr[i].sum / (double) arr[i].count) : NAN;
        fprintf(fp, "%*s{ \"id\": %d, \"name\": \"%s\", \"mean\": %.10g, \"n\": %d }%s\n",
                indent + 2, "", arr[i].id,
                arr[i].name ? arr[i].name : "unknown",
                mean, arr[i].count,
                (i + 1 == count) ? "" : ",");
    }
    fprintf(fp, "%*s]%s\n", indent, "", trailing_comma ? "," : "");
}

int main(int argc, char *argv[])
{
    File_formats format = BINARY_FORMAT;
    int i, j;
    roi_stat_t *rois[3] = { NULL, NULL, NULL };
    int roi_count[3] = { 0, 0, 0 };
    int roi_cap[3]   = { 0, 0, 0 };

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0)) {
        return EXIT_FAILURE;
    }

    if (g_nunits < 1 || !g_out) {
        fprintf(stderr, "Usage: CAT_Surf2ROIMulti -unit src_sphere=...,trg_sphere=...,annot=...,vals=...[,hemi=lh|rh] [-unit ...] -out output.json\n");
        return EXIT_FAILURE;
    }

    for (i = 0; i < g_nunits; i++) {
        unit_t *u = &g_units[i];
        int n_objects = 0;
        object_struct **obj_src_sph = NULL, **obj_trg_sph = NULL, **obj_resampled = NULL;
        polygons_struct *poly = NULL, *src_sph = NULL, *trg_sph = NULL;

        /* Input geometry */
        if (input_graphics_any_format(u->src_sphere, &format, &n_objects, &obj_src_sph) != OK ||
            n_objects != 1 || get_object_type(obj_src_sph[0]) != POLYGONS) {
            fprintf(stderr, "[%d] %s must contain 1 polygons object.\n", i, u->src_sphere);
            return EXIT_FAILURE;
        }
        src_sph = get_polygons_ptr(obj_src_sph[0]);
        poly = src_sph;

        if (input_graphics_any_format(u->trg_sphere, &format, &n_objects, &obj_trg_sph) != OK ||
            n_objects != 1 || get_object_type(obj_trg_sph[0]) != POLYGONS) {
            fprintf(stderr, "[%d] %s must contain 1 polygons object.\n", i, u->trg_sphere);
            return EXIT_FAILURE;
        }
        trg_sph = get_polygons_ptr(obj_trg_sph[0]);

        /* Annotation */
        int *in_annot = NULL;
        int n_array = 0, n_labels = 0;
        ATABLE *atable = NULL;
        if (read_annotation_table(u->annot, &n_array, &in_annot, &n_labels, &atable) != OK) {
            fprintf(stderr, "[%d] Cannot read annotation %s.\n", i, u->annot);
            return EXIT_FAILURE;
        }
        if (n_array != src_sph->n_points) {
            fprintf(stderr, "[%d] #annot (%d) != #points in source sphere (%d).\n", i, n_array, src_sph->n_points);
            return EXIT_FAILURE;
        }

        double *ann_in = (double *) SAFE_MALLOC(double, n_array);
        double *ann_out = (double *) SAFE_MALLOC(double, trg_sph->n_points);
        for (j = 0; j < n_array; j++) ann_in[j] = (double) in_annot[j];

        obj_resampled = resample_surface_to_target_sphere(poly, src_sph, trg_sph,
                                  ann_in, ann_out, 1 /* label_interp */, 0);

        /* Values on target sphere */
        double *vals = NULL;
        int n_vals = 0;
        if (input_values_any_format(u->vals, &n_vals, &vals) != OK) {
            fprintf(stderr, "[%d] Cannot read values in %s.\n", i, u->vals);
            return EXIT_FAILURE;
        }
        if (n_vals != trg_sph->n_points) {
            fprintf(stderr, "[%d] #values (%d) != #points in target sphere (%d).\n", i, n_vals, trg_sph->n_points);
            return EXIT_FAILURE;
        }

        hemi_t hemi = guess_hemi_from_unit(u);
        int hemi_idx = (int) hemi;
        if (hemi_idx < 0 || hemi_idx > 2) hemi_idx = 0;

        for (j = 0; j < trg_sph->n_points; j++) {
            int lbl = (int) lround(ann_out[j]);
            double v = vals[j];
            if (isnan(v)) continue;
            const char *nm = name_for_id(atable, n_labels, lbl);
            add_sum(&rois[hemi_idx], &roi_count[hemi_idx], &roi_cap[hemi_idx], lbl, nm, v);
        }

        if (ann_in) free(ann_in);
        if (ann_out) free(ann_out);
        if (vals) FREE(vals);
        if (in_annot) free(in_annot);
        if (atable) free(atable);
        if (obj_resampled) delete_object_list(1, obj_resampled);
        if (obj_src_sph) delete_object_list(1, obj_src_sph);
        if (obj_trg_sph) delete_object_list(1, obj_trg_sph);
    }

    /* Write JSON */
    FILE *fp = SAFE_FOPEN(g_out, "w");
    fprintf(fp, "{\n");
    int has_lh = roi_count[HEMI_LH] > 0;
    int has_rh = roi_count[HEMI_RH] > 0;
    int has_unknown = roi_count[HEMI_UNKNOWN] > 0;
    int remaining = has_lh + has_rh + has_unknown;

    if (has_lh) {
        write_roi_array(fp, "lh", rois[HEMI_LH], roi_count[HEMI_LH], 2, --remaining > 0);
    }
    if (has_rh) {
        write_roi_array(fp, "rh", rois[HEMI_RH], roi_count[HEMI_RH], 2, --remaining > 0);
    }
    if (has_unknown) {
        write_roi_array(fp, "unknown", rois[HEMI_UNKNOWN], roi_count[HEMI_UNKNOWN], 2, --remaining > 0);
    }
    fprintf(fp, "}\n");
    fclose(fp);

    for (i = 0; i < 3; i++) {
        int n = roi_count[i];
        int k;
        for (k = 0; k < n; k++) {
            if (rois[i][k].name) free(rois[i][k].name);
        }
        if (rois[i]) free(rois[i]);
    }
    if (g_units) {
        for (i = 0; i < g_nunits; i++) {
            free(g_units[i].src_sphere);
            free(g_units[i].trg_sphere);
            free(g_units[i].annot);
            free(g_units[i].vals);
        }
        free(g_units);
    }

    return EXIT_SUCCESS;
}
