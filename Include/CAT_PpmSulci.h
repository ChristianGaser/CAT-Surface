/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_PPMSULCI_H_
#define _CAT_PPMSULCI_H_

#include "CAT_Sheetness.h"

/*
 * Opening buried sulci directly in a percentage position map (PPM).
 *
 * The pre-PBT repairs above need the intensity image. By the time the central
 * surface is extracted the T1 is gone -- marching cubes sees only the PPM --
 * but the PPM carries the geometry itself, so no intensity is required:
 *
 *   crossing a sulcus:  1 (WM) -> 0.5 -> ~0 (pial) -> 0.5 -> 1 (other bank)
 *   crossing a gyrus:   0 (pial) -> 0.5 -> ~1 (WM blade) -> 0.5 -> 0
 *
 * A sulcus is therefore a *valley* in the PPM and a gyral blade a *ridge*, and
 * the polarity guard of the sheetness filter separates them exactly. A buried
 * sulcus is a valley whose floor never drops below the isovalue: the geometry
 * is still there, only the amplitude is missing, which is why an isovalue-based
 * surface fuses the two banks while the valley remains plainly visible to a
 * shape filter.
 */

/** \brief Parameters for opening buried sulci in a PPM. */
typedef struct
{
    double sigma_min; /**< smallest sheetness scale in mm */
    double sigma_max; /**< largest sheetness scale in mm */
    int n_scales;     /**< number of log-spaced scales */
    double sheet_normalize; /**< value the response's p99.9 is scaled to (default
                                 CAT_SHEETNESS_NORMALIZE); <= 0 keeps it raw */
    int sheet_skeleton;    /**< 1 to thin the response to its medial sheet (default 0);
                                see sheet_skeleton above */
    double sheet_strength; /**< overall gain on the response (default 1.0), passed
                                through as CAT_SheetnessOpts::gain.  The automatic
                                noise scale of the sheetness filter is half the
                                largest Hessian norm in the volume, which on real
                                data is set by the cortical ribbon itself -- a thin
                                sulcal valley is far weaker than that, so the raw
                                response usually sits an order of magnitude below
                                the thresholds here and nothing happens at a gain
                                of 1.  With the p99.9 anchor enabled this is
                                normally unnecessary; measure with
                                CAT_VolSheetness -polarity -1 on the PPM and choose
                                a gain that puts p99 near twice thresh. */
    double thresh;    /**< sheetness below this is ignored */
    double band;      /**< only act on values in (isovalue, isovalue + band) */
    double margin;    /**< how far below the isovalue an opened valley is pushed */
    double strength;  /**< 0..1 blend towards that target */
    double offset;    /**< signed sheetness offset applied to the map before the
                           surface is extracted, in map units (default 0 = off).
                           A sulcus is a valley in a PPM and a gyral blade a
                           ridge, so a *signed* sheetness is negative on one and
                           positive on the other, and adding it lowers the map
                           along sulci while raising it along blades in a single
                           pass.  This is what a global isovalue shift cannot do:
                           lowering the isovalue opens glued sulci and severs
                           thin white-matter fingers with the same stroke,
                           because it moves every voxel the same way regardless
                           of what is there.  The offset is scaled by the
                           response, so it is strongest where the structure is
                           clearest and vanishes where nothing was found.
                           sheet_skeleton applies here too, and confines the
                           offset to the medial line instead of a band as wide
                           as the Gaussian that produced the response -- the more
                           defensible place to displace a surface from.  It is
                           also a weaker correction at the same offset (on a real
                           PPM it roughly halves the responding field and cuts
                           isovalue crossings by about two thirds), so raise the
                           offset to match if you enable it. */
    double cutoff;    /**< admission cutoff of the oriented median the caller may run
                           alongside; <= 0 selects CAT_ORIENTED_MEDIAN_CUTOFF */
    int verbose;
} CAT_PpmSulciOpts;

/**
 * \brief Fill a CAT_PpmSulciOpts with defaults for a 0..1 PPM at 0.5 mm.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_PpmSulciOptionsInit(CAT_PpmSulciOpts *opts);

/**
 * \brief Push buried sulcal valleys in a PPM below the isovalue.
 *
 * Finds thin low-value sheets -- valleys -- in the map and lowers the ones
 * that sit just above the isovalue until they cross it, so the isosurface
 * separates the two banks instead of bridging them.
 *
 * Three conditions must hold together, and it is the conjunction that makes
 * the operation safe rather than just another filter:
 *
 *   - the value lies in (isovalue, isovalue + band), so only marginal cases
 *     are touched and solid white matter (PPM near 1) never is;
 *   - the dark-sheet response exceeds opts->thresh, so there really is a thin
 *     planar structure and not a blob or noise;
 *   - the value is a genuine local minimum *along the sheet normal*, which a
 *     gyral crown can never satisfy -- crossing a blade the PPM has a maximum,
 *     not a minimum, so the operation cannot cut a gyrus open.
 *
 * The correction is one-sided: it only ever lowers a value.
 *
 * \param ppm       (in/out) percentage position map in [0,1], corrected in place
 * \param sheetness (in)     precomputed dark-sheet response, or NULL to compute
 * \param normal    (in)     precomputed unit sheet normals, or NULL to compute
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param isovalue  (in)     threshold the surface will be extracted at
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolOpenPpmSulci(float *ppm, const float *sheetness, const float *normal,
                        int dims[3], double voxelsize[3], double isovalue,
                        const CAT_PpmSulciOpts *opts);

#endif /* _CAT_PPMSULCI_H_ */
