/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_THICKNESSQC_H_
#define _CAT_THICKNESSQC_H_

/*
 * Shape triage of implausibly thick cortex.
 *
 * A thickness map from CAT_VolThicknessPbt usually carries a population of
 * voxels above any anatomically defensible value -- human cortex does not
 * exceed roughly 4.5 mm. Two quite different faults produce them, and the
 * repair for one is actively harmful applied to the other:
 *
 *   glued sulcus   two banks with no CSF between them, measured as one band.
 *                  The fix is to separate them, e.g. CAT_VolSulcusRepair.
 *
 *   solid mass     cortex merged with subcortical grey matter, or a genuinely
 *                  thick region such as a temporal or orbitofrontal pole.
 *                  There is no sulcus to recover; carving one invents anatomy.
 *
 * Thickness alone cannot tell them apart, which makes it easy to attribute the
 * whole population to gluing and build a separator for a problem that is not
 * there. Shape can tell them apart, and cheaply. A glued sulcus is two cortical
 * banks back to back, so the band it forms is at most about 5 mm across
 * whatever its extent along the sulcus; the largest sphere that fits inside it
 * therefore has a radius of about 2.5 mm. A solid mass has no such bound. So
 * the maximum inscribed radius of a connected component -- the maximum of the
 * Euclidean distance transform taken inside it -- separates the two cases with
 * one number, and one that does not depend on where the component sits or how
 * large it is.
 *
 * For scale: a normal 2.5 mm cortical ribbon gives an inscribed radius near
 * 1.25 mm, a glued pair of banks stays under 2.5 mm, and a subcortical mass or
 * a pole runs to 3.5 mm and beyond.
 */

/** \brief Parameters of the thickness triage. */
typedef struct
{
    double thresh;       /**< thickness above this is implausible, in mm (default 4.5) */
    double plate_radius; /**< largest inscribed radius still called a plate, in mm (default 2.5) */
    double min_volume;   /**< discard components below this volume, in mm^3 (default 20) */
    int conn;            /**< voxel connectivity: 6, 18 or 26 (default 26) */
    int verbose;         /**< 1 to report progress */
} CAT_ThicknessQCOpts;

/** \brief Component shapes reported by CAT_VolThicknessQC(). */
#define CAT_QC_PLATE 1 /**< thin band: consistent with a glued sulcus */
#define CAT_QC_SOLID 2 /**< solid core: not a sulcus */

/** \brief One connected component of implausibly thick cortex. */
typedef struct
{
    long n_voxels;      /**< size of the component */
    double volume_mm3;  /**< its volume */
    double max_radius;  /**< maximum inscribed radius, in mm; the shape criterion */
    double centroid[3]; /**< centre of mass, in voxel coordinates */
    double gmt_mean;    /**< mean thickness inside the component, in mm */
    double gmt_max;     /**< maximum thickness inside the component, in mm */
    int shape;          /**< CAT_QC_PLATE or CAT_QC_SOLID */
} CAT_ThicknessComponent;

/**
 * \brief Fill a CAT_ThicknessQCOpts with defaults.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_ThicknessQCOptionsInit(CAT_ThicknessQCOpts *opts);

/**
 * \brief Find implausibly thick components and classify them by shape.
 *
 * Thresholds the thickness map, optionally restricted to the cortical band of a
 * PVE label map, groups the surviving voxels into connected components, and
 * measures the maximum inscribed radius of each from the Euclidean distance
 * transform taken inside the component. A component whose radius stays at or
 * below opts->plate_radius is reported as CAT_QC_PLATE and is consistent with a
 * glued sulcus; anything thicker is CAT_QC_SOLID and is not a sulcus at all.
 *
 * The caller owns the returned array and must free() it. Passing NULL for
 * `comps` and `n_comps` computes only the class map, and vice versa.
 *
 * \param gmt       (in)  thickness map in mm
 * \param label     (in)  PVE label map in [0..3] restricting the search to the
 *                        cortical band, or NULL to use the whole thickness map
 * \param dims      (in)  {nx, ny, nz}
 * \param voxelsize (in)  voxel spacing in mm
 * \param opts      (in)  parameters; NULL selects the defaults
 * \param comps     (out) newly allocated array of components, or NULL to skip
 * \param n_comps   (out) number of components returned, or NULL to skip
 * \param classmap  (out) float[nvox] tagged with CAT_QC_PLATE / CAT_QC_SOLID
 *                        and 0 elsewhere, or NULL to skip
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolThicknessQC(const float *gmt, const float *label, int dims[3],
                       double voxelsize[3], const CAT_ThicknessQCOpts *opts,
                       CAT_ThicknessComponent **comps, int *n_comps,
                       float *classmap);

#endif /* _CAT_THICKNESSQC_H_ */
