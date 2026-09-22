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

/**
 * \brief Resample integer annotation labels onto a target sphere.
 *
 * \param src_sphere (in)  source sphere mesh
 * \param trg_sphere (in)  target sphere mesh
 * \param labels_src (in)  source labels per vertex
 * \param labels_trg (out) allocated target labels per vertex
 * \return OK on success, ERROR otherwise
 */
Status CAT_ResampleAnnotationLabels(polygons_struct *src_sphere,
                                    polygons_struct *trg_sphere,
                                    const int *labels_src, int **labels_trg);

/**
 * \brief Compute per-label sums and sample counts for ROI statistics.
 *
 * \param labels     (in)  integer labels per vertex
 * \param vals       (in)  values per vertex
 * \param n_points   (in)  number of vertices
 * \param atable     (in)  annotation table (can be NULL)
 * \param n_labels   (in)  number of annotation entries
 * \param out_stats  (out) allocated stats array
 * \param out_n_stats (out) number of stats entries
 * \return OK on success, ERROR otherwise
 */
Status CAT_ComputeROIMeansFromLabels(const int *labels, const double *vals,
                                     int n_points, const ATABLE *atable,
                                     int n_labels, CAT_ROIStat **out_stats,
                                     int *out_n_stats);

#ifdef __cplusplus
}
#endif

#endif
