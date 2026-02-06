/*
 * CAT_ROIStats.h
 *
 * Helpers to resample annotation labels onto a target sphere and
 * compute per-label ROI statistics from a values array.
 */

#ifndef CAT_ROI_STATS_H
#define CAT_ROI_STATS_H

#include "CAT_SurfaceIO.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int id;               /* annotation id */
    const char *name;     /* pointer into ATABLE (owned by caller) or "unknown" */
    double sum;           /* sum over all non-NaN values for this id */
    int n;                /* number of samples contributing to sum */
} CAT_ROIStat;

/*
 * Resample integer annotation labels defined on the source sphere vertices
 * onto the target sphere vertices, using label interpolation.
 *
 * labels_trg is allocated by this function and must be freed by the caller.
 */
Status CAT_ResampleAnnotationLabels(
    polygons_struct *src_sphere,
    polygons_struct *trg_sphere,
    const int *labels_src,
    int **labels_trg);

/*
 * Compute per-label sums and sample counts (skipping NaN values) from labels
 * and values arrays defined on the same target sphere.
 *
 * The returned array is allocated by this function and must be freed by the
 * caller. The mean for each label is stat[i].sum / (double) stat[i].n.
 */
Status CAT_ComputeROIMeansFromLabels(
    const int *labels,
    const double *vals,
    int n_points,
    const ATABLE *atable,
    int n_labels,
    CAT_ROIStat **out_stats,
    int *out_n_stats);

#ifdef __cplusplus
}
#endif

#endif
