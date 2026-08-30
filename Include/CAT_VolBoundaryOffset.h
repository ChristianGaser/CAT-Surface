/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/**
 * \file CAT_VolBoundaryOffset.h
 * \brief Measure how far the GM/WM boundary has been displaced by myelination.
 *
 * In the primary motor and somatosensory strip, and along the line of Gennari
 * in V1, the deep cortical layers carry enough myelin to approach white-matter
 * intensity on T1w.  The classifier follows the intensity, so the GM/WM
 * boundary is placed too far out and the cortex comes back too thin.
 *
 * Correcting this by relabelling voxels cannot work, because the displacement
 * is a fraction of a voxel and a label can only move in whole isovalue
 * crossings.  What the displacement *can* be expressed in is PBT's distance
 * field, which already carries the boundary to sub-voxel precision (see the
 * PVE refinement of `dist_WM` in CAT_VolPbt.c).  This module estimates the
 * displacement in millimetres so it can be added there.
 *
 * **The observable is the width of the intensity transition, not its
 * position.**  The label map is derived from the intensity, so the label
 * boundary and the intensity boundary agree by construction and their
 * difference measures nothing.  What distinguishes myelinated cortex is the
 * *shape* of the profile across the boundary: healthy cortex steps from grey
 * to white over about the partial-volume width, while a myelinated deep layer
 * turns that step into a ramp one to three millimetres wide.  The excess width
 * over what this brain shows in its healthy majority is the signal.
 *
 * The profile is read along the cortical normal, taken from the gradient of
 * the label map, and in a contrast-normalized coordinate
 * `t = (T1w - GM_local) / (WM_local - GM_local)` built from local tissue
 * references.  Because both ends of that ratio are local, a residual bias
 * field cancels out of it and the level crossings are comparable across the
 * brain.
 *
 * The result is deliberately *diagnostic first*: write the map, look at
 * whether it picks out the central sulcus and V1, and only then feed it to
 * PBT.
 */

#ifndef _CAT_VOL_BOUNDARY_OFFSET_H_
#define _CAT_VOL_BOUNDARY_OFFSET_H_

#include "CAT_Vol.h"

/**
 * \brief Options controlling the boundary-offset estimate.
 */
typedef struct
{
    double ref_fwhm;      /**< FWHM (mm) of the local GM and WM intensity references (default: 10.0). Only has to track the residual bias field; noise is averaged out by the number of control points at any of these scales. */
    double erosion_mm;    /**< Erosion (mm) applied to the WM mask before it is used as an intensity control set (default: 3.0), so the myelinated band does not calibrate the reference it is measured against. */
    double t_lo;          /**< Lower level crossing of the normalized profile (default: 0.25). */
    double t_hi;          /**< Upper level crossing (default: 0.75). The distance between the two is the transition width. */
    double search_mm;     /**< Half-length (mm) of the profile searched along the normal (default: 4.0). */
    double step_mm;       /**< Sampling step (mm) along the profile (default: 0.25). */
    double width_pct;     /**< Upper percentile of the measured widths above which a transition is read as myelinated (default: 88.0), i.e. the share of the boundary the correction is allowed to touch. A *central* location is no use here: whatever one is chosen the healthy population straddles it, every voxel above acquires a positive offset, and the whole cortex ends up displaced -- at p25 that was 99.8% of the boundary on a real subject. Nor does a location plus a multiple of the spread transfer between subjects, because the spread follows resolution and noise rather than anatomy: the same setting selected 12.0% of the boundary on a 0.75 mm scan and 1.4% on a 1x1x1.25 mm one. A fixed share transfers by construction and matches the anatomy, heavily myelinated cortex being a roughly fixed tenth of the sheet. Fixing the share does not fix the correction: the displacement is the excess *beyond* the threshold, so a brain whose widths are tightly grouped is barely corrected however large the share. */
    double gain;          /**< Displacement per unit excess width (default: 0.5). For a symmetric ramp the cytoarchitectonic boundary sits about half the excess beyond the intensity midpoint. */
    double max_offset_mm; /**< Clamp on the returned displacement (default: 1.5). */
    double smooth_fwhm;   /**< FWHM (mm) of the smoothing applied to the transition width *within* the boundary sheet (default: 8.0), before any threshold is applied. Myelination is a centimetre-scale, stereotyped phenomenon while the per-voxel profile measurement is noisy, and the band is a thin surface, so a normalized convolution restricted to it averages along the boundary and not across it. Doing this before the threshold rather than after is what makes the threshold meaningful, and it removes any reason to worry about the non-negativity clamp rectifying noise into a systematic inflation. */
    int fill_ribbon;      /**< Carry the displacement off the boundary sheet and through the tissue, each voxel taking the value measured at the nearest point of the sheet (default: 1). The measurement is only defined where the boundary is, but the quantity it corrects -- cortical thickness -- is a property of the whole column and is read at the central surface, half a ribbon away. Seeds include sheet voxels whose profile was rejected, so an unmeasurable spot takes its neighbours' value rather than a spurious zero. */
    int verbose;          /**< Print progress and summary statistics (default: 0) */
} CAT_BoundaryOffsetOpts;

/**
 * \brief Initialize boundary-offset options with defaults.
 *
 * \param opts (out) options structure to fill
 */
void CAT_BoundaryOffsetOptionsInit(CAT_BoundaryOffsetOpts *opts);

/**
 * \brief Local tissue-intensity reference by normalized convolution.
 *
 * Estimates the intensity a given tissue has *at each voxel* by smoothing
 * `values` over the control voxels and dividing by the smoothed control mask.
 * The quotient is the control-point average weighted by distance, so it is
 * defined everywhere -- including away from the control set, where it
 * extrapolates smoothly -- and it follows a residual bias field instead of
 * averaging it away.
 *
 * This is what lets an intensity test be a fixed number.  Comparing against
 * one global mean and standard deviation makes the threshold move with the
 * artifact being looked for: a bias field inflates the standard deviation and
 * loosens it, and the myelinated territory itself drags the mean down.
 *
 * \param values    (in)  volume to average, typically the T1w
 * \param control   (in)  1 where the voxel is a control point for this tissue
 * \param dims      (in)  volume dimensions {nx, ny, nz}
 * \param voxelsize (in)  voxel sizes in mm {dx, dy, dz}
 * \param fwhm_mm   (in)  smoothing FWHM in mm
 * \param fallback  (in)  value used where no control point carries any weight
 * \param ref       (out) local reference, one value per voxel
 * \return 0 on success, non-zero on allocation failure
 */
int CAT_VolLocalTissueReference(const float *values,
                                const unsigned char *control,
                                int dims[3], double voxelsize[3],
                                double fwhm_mm, double fallback, float *ref);

/**
 * \brief Estimate the myelination-induced GM/WM boundary displacement.
 *
 * For every voxel in the GM/WM transition band, reads the T1w profile along
 * the cortical normal, measures the width of the transition between the
 * `t_lo` and `t_hi` level crossings, and converts the excess over this brain's
 * healthy reference width into a displacement in millimetres.  The result is
 * smoothed within the boundary sheet, which is the direction the estimate is
 * coherent in.
 *
 * The returned offset is non-negative: myelination can only cause the boundary
 * to be placed too far out, never too far in.
 *
 * \param label     (in)  PVE label volume (float, values in [0..3])
 * \param t1w       (in)  T1w intensity volume on the same grid
 * \param offset    (out) displacement in mm, one value per voxel.  With
 *                        `fill_ribbon` it is carried through the tissue from
 *                        the boundary sheet, so it can be added directly to a
 *                        thickness map; otherwise it is 0 off the sheet.
 * \param width     (out) optional transition width in mm, smoothed within the
 *                        sheet and before the threshold is applied -- the
 *                        quantity the decision is actually made on, and the
 *                        map to inspect first; NULL to skip
 * \param dims      (in)  volume dimensions {nx, ny, nz}
 * \param voxelsize (in)  voxel sizes in mm {dx, dy, dz}
 * \param opts      (in)  algorithm options (NULL for defaults)
 * \return number of band voxels with a valid profile, or negative on error
 */
long CAT_VolBoundaryOffset(const float *label, const float *t1w,
                           float *offset, float *width,
                           int dims[3], double voxelsize[3],
                           const CAT_BoundaryOffsetOpts *opts);

#endif /* _CAT_VOL_BOUNDARY_OFFSET_H_ */
