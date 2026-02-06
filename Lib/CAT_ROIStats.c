/*
 * CAT_ROIStats.c
 */

#include <bicpl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CAT_ROIStats.h"
#include "CAT_Resample.h"
#include "CAT_SafeAlloc.h"

static const char *cat_name_for_id(const ATABLE *atable, int n_labels, int id)
{
    int i;

    if (!atable || n_labels <= 0) return "unknown";

    for (i = 0; i < n_labels; i++) {
        if (atable[i].annotation == id) return atable[i].name;
    }

    return "unknown";
}

Status CAT_ResampleAnnotationLabels(
    polygons_struct *src_sphere,
    polygons_struct *trg_sphere,
    const int *labels_src,
    int **labels_trg)
{
    int i;
    object_struct **obj_resampled;
    double *ann_in;
    double *ann_out;
    int *out;

    if (!src_sphere || !trg_sphere || !labels_src || !labels_trg) {
        fprintf(stderr, "CAT_ResampleAnnotationLabels: invalid arguments.\n");
        return ERROR;
    }

    if (src_sphere->n_points <= 0 || trg_sphere->n_points <= 0) {
        fprintf(stderr, "CAT_ResampleAnnotationLabels: spheres have no points.\n");
        return ERROR;
    }

    ann_in = (double *) SAFE_MALLOC(double, src_sphere->n_points);
    ann_out = (double *) SAFE_MALLOC(double, trg_sphere->n_points);

    for (i = 0; i < src_sphere->n_points; i++) {
        ann_in[i] = (double) labels_src[i];
    }

    obj_resampled = resample_surface_to_target_sphere(src_sphere, src_sphere, trg_sphere,
                                                      ann_in, ann_out,
                                                      1 /* label_interp */, 0);

    out = (int *) SAFE_MALLOC(int, trg_sphere->n_points);
    for (i = 0; i < trg_sphere->n_points; i++) {
        out[i] = (int) lround(ann_out[i]);
    }

    if (obj_resampled) delete_object_list(1, obj_resampled);
    if (ann_in) free(ann_in);
    if (ann_out) free(ann_out);

    *labels_trg = out;

    return OK;
}

Status CAT_ComputeROIMeansFromLabels(
    const int *labels,
    const double *vals,
    int n_points,
    const ATABLE *atable,
    int n_labels,
    CAT_ROIStat **out_stats,
    int *out_n_stats)
{
    CAT_ROIStat *stats;
    int n_stats;
    int cap;
    int i;

    if (!labels || !vals || n_points <= 0 || !out_stats || !out_n_stats) {
        fprintf(stderr, "CAT_ComputeROIMeansFromLabels: invalid arguments.\n");
        return ERROR;
    }

    stats = NULL;
    n_stats = 0;
    cap = 0;

    for (i = 0; i < n_points; i++) {
        int id;
        double v;
        int k;

        v = vals[i];
        if (isnan(v)) continue;

        id = labels[i];

        for (k = 0; k < n_stats; k++) {
            if (stats[k].id == id) break;
        }

        if (k == n_stats) {
            if (n_stats == cap) {
                int new_cap = (cap == 0) ? 16 : (cap * 2);
                CAT_ROIStat *tmp = (CAT_ROIStat *) realloc(stats, sizeof(CAT_ROIStat) * new_cap);
                if (!tmp) {
                    fprintf(stderr, "CAT_ComputeROIMeansFromLabels: out of memory.\n");
                    if (stats) free(stats);
                    return ERROR;
                }
                stats = tmp;
                cap = new_cap;
            }

            stats[n_stats].id = id;
            stats[n_stats].name = cat_name_for_id(atable, n_labels, id);
            stats[n_stats].sum = 0.0;
            stats[n_stats].n = 0;
            k = n_stats;
            n_stats++;
        }

        stats[k].sum += v;
        stats[k].n += 1;
    }

    *out_stats = stats;
    *out_n_stats = n_stats;

    return OK;
}
