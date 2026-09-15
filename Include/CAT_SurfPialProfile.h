/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURF_PIAL_PROFILE_H_
#define _CAT_SURF_PIAL_PROFILE_H_

/**
 * @file CAT_SurfPialProfile.h
 * @brief Profile-based placement of the pial surface.
 *
 * Each vertex searches the intensity profile along its normal for its own
 * target: the isovalue crossing, or - where a glued sulcus never reaches the
 * isovalue - the bottom of the valley the two banks form.  Facing walls of a
 * sulcus meet in the middle instead of crossing, and the regularization acts
 * on the per-iteration update only, so it does not pull convex gyral crowns
 * inwards the way smoothing the accumulated displacement does.
 */

#include <bicpl.h>
#include "CAT_NiftiLib.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * Options for profile-based pial placement.
     */
    typedef struct
    {
        double isovalue;          /**< Boundary value in the volume (default: 1.5, CSF/GM) */
        double search_out;        /**< Outward search distance along the normal in mm (default: 2.0) */
        double search_in;         /**< Inward search distance along the normal in mm (default: 1.0) */
        double sample_step;       /**< Profile sampling step in mm (default: 0.1) */
        double valley_depth;      /**< Rise above the running minimum that ends a valley (default: 0.05) */
        double step_fraction;     /**< Fraction of the target offset applied per iteration (default: 0.5) */
        double max_step;          /**< Maximum normal displacement per iteration in mm (default: 0.25) */
        double max_offset;        /**< Maximum normal offset from the start surface in mm (default: 2.0) */
        double smooth_lambda;     /**< Blend weight of the neighbour average when smoothing the update (default: 0.5) */
        int smooth_passes_start;  /**< Update smoothing passes in the first iteration (default: 5) */
        int smooth_passes_end;    /**< Update smoothing passes in the last iteration (default: 1) */
        double tangential_weight; /**< Tangential relaxation weight (default: 0.3) */
        double concave_weight;    /**< Normal Laplacian weight, applied in concave regions only (default: 0.1) */
        double contact_margin;    /**< Gap kept between facing sulcal walls in mm (default: 0.1) */
        double fold_angle;        /**< Maximum rotation of a face per iteration in degrees (default: 72.5) */
        int iterations;           /**< Number of iterations (default: 60) */
        int verbose;              /**< Verbose output (default: 0) */
    } CAT_PialProfileOptions;

    /**
     * \brief Initialize profile-based pial placement options with defaults.
     *
     * \param opts (out) options structure to initialize
     * \return void
     */
    void CAT_PialProfileOptionsInit(CAT_PialProfileOptions *opts);

    /**
     * \brief Move a pial surface onto the CSF/GM boundary by profile search.
     *
     * Iteratively places each vertex on the isovalue crossing along its normal
     * or, where the profile turns upwards again before reaching the isovalue,
     * on the bottom of that valley.  Vertices that see neither within the
     * search range are held in place.  The update is smoothed over the mesh,
     * clamped so that facing sulcal walls stop half-way, limited to
     * max_offset from the start surface, and reverted locally where a face
     * would fold over.
     *
     * The volume is sampled in world coordinates through nii->sto_xyz, so the
     * result does not depend on the image orientation.
     *
     * \param pial (in/out) pial surface; start positions in, placed surface out
     * \param vol  (in)     volume data (e.g. label map with CSF=1, GM=2, WM=3)
     * \param nii  (in)     NIfTI header of vol
     * \param opts (in)     options, see CAT_PialProfileOptionsInit()
     * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
     */
    int CAT_SurfDeformPialProfile(polygons_struct *pial, const float *vol,
                                  nifti_image *nii,
                                  const CAT_PialProfileOptions *opts);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_SURF_PIAL_PROFILE_H_ */
