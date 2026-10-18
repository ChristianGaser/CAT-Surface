/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef CAT_REALIGN_H
#define CAT_REALIGN_H

/**
 * \file CAT_Realign.h
 * \brief Head-motion correction of a 4D series, in mcflirt's conventions.
 *
 * Every volume is registered rigidly to a reference with the robust
 * Gauss-Newton registration of CAT_VolumeReg, global intensity scale
 * included, so a reference brighter than the series (the median of its
 * non-steady-state volumes, a single-band reference) is no problem.  The
 * volumes are processed outward from the middle one, each starting from the
 * estimate of its neighbour, as mcflirt does.
 *
 * The results follow mcflirt (FSL 6.0.7.23) closely enough to replace it in
 * fMRIPrep's head-motion workflow, which calls
 *   mcflirt -in bold -out bold_mcf.nii.gz -reffile ref -mats -plots -rmsabs -rmsrel
 * and reads the matrices, the parameters and the relative RMS:
 *
 * - MAT_%04d: volume -> reference in FSL's scaled-voxel coordinates
 *   (ras_to_fsl_matrix()), the reference volume itself the identity.
 * - .par: rx ry rz (rad) tx ty tz (mm) per volume.  R = Rx Ry Rz with FSL's
 *   elementary rotations, the translation taken after rotating about the
 *   centre of mass of the reference (FSL mm, image minimum subtracted first).
 * - _abs.rms / _rel.rms: Jenkinson's RMS deviation over a sphere of 80 mm
 *   around the centre of the reference grid, (n - 1) / 2 * pixdim -- not the
 *   centre of mass -- against the reference (nt values) and between
 *   consecutive volumes (nt - 1); _abs_mean.rms and _rel_mean.rms hold their
 *   means.
 *
 * Each convention was checked against mcflirt's own output.  The estimates
 * themselves differ from mcflirt's (a different cost and optimiser); on
 * synthetic series with known motion they are closer to the truth.
 */

#include "CAT_NiftiLib.h"
#include "CAT_VolumeReg.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /** \brief Options of the motion correction. */
    typedef struct
    {
        int refvol;            /**< reference volume of the series when no
                                    separate reference is given; -1 (default)
                                    is nt / 2, as in mcflirt */
        CAT_VolumeRegOpts reg; /**< options of each rigid registration; the
                                    defaults are CAT_VolumeRegOptionsInit()'s
                                    with 3 pyramid levels */
        int verbose;           /**< non-zero: one line per volume */
    } CAT_VolRealignOpts;

    /**
     * \brief Motion estimates of a series, in world space and in mcflirt's terms.
     *
     * Filled by CAT_VolRealignEstimate(); release with CAT_VolRealignFree().
     * All matrices are 4x4 row-major and map a volume to the reference.
     */
    typedef struct
    {
        int     nt;           /**< number of volumes */
        int     ref_index;    /**< volume used as reference, -1 for a separate one */
        double *ras;          /**< nt world (RAS) transforms, volume -> reference */
        double *fsl;          /**< nt transforms in FSL's convention (MAT_%04d) */
        double *par;          /**< nt x 6 parameters rx ry rz (rad) tx ty tz (mm) */
        double *rms_abs;      /**< nt RMS deviations from the reference (mm) */
        double *rms_rel;      /**< nt - 1 RMS deviations between volumes (mm) */
        double *scale;        /**< nt intensity scales s, reference = s * volume */
        double  rms_abs_mean; /**< mean of rms_abs */
        double  rms_rel_mean; /**< mean of rms_rel (0 for a single volume) */
    } CAT_VolRealignResult;

    /**
     * \brief Set the default options of the motion correction.
     *
     * The reference is the middle volume, and each registration uses the
     * defaults of CAT_VolumeRegOptionsInit() -- intensity scale included --
     * with 3 pyramid levels: an EPI of 3-4 mm voxels then sees scales of
     * about 4, 8 and 16 mm, close to mcflirt's 4 and 8 mm.
     *
     * \param opts (out) options to initialise
     */
    void CAT_VolRealignOptionsInit(CAT_VolRealignOpts *opts);

    /**
     * \brief Estimate the motion of every volume of a series.
     *
     * The volumes are registered to the reference starting with the middle
     * one (the reference volume itself is the identity), then outward in both
     * directions on two threads, each volume starting from the estimate and
     * intensity scale of its neighbour.  The reference pyramid is built once.
     *
     * \param series     (in)  nt volumes of dims[0]*dims[1]*dims[2] floats, as
     *                         read_nifti_float() returns them
     * \param series_nii (in)  header of the series (sto_xyz, voxel size)
     * \param dims       (in)  [nx, ny, nz] of one volume
     * \param nt         (in)  number of volumes, at least 1
     * \param ref        (in)  separate reference volume, or NULL to use volume
     *                         opts->refvol of the series
     * \param ref_nii    (in)  header of the separate reference (ignored for NULL)
     * \param ref_dims   (in)  [nx, ny, nz] of the separate reference (ignored for NULL)
     * \param opts       (in)  options; NULL for CAT_VolRealignOptionsInit()
     * \param res        (out) estimates; release with CAT_VolRealignFree()
     * \return 0 on success, -1 on error (a message is written to stderr)
     */
    int CAT_VolRealignEstimate(const float *series, const nifti_image *series_nii,
                               const int dims[3], int nt,
                               const float *ref, const nifti_image *ref_nii,
                               const int ref_dims[3],
                               const CAT_VolRealignOpts *opts,
                               CAT_VolRealignResult *res);

    /**
     * \brief Resample every volume into the reference grid.
     *
     * Trilinear interpolation, as mcflirt's default; positions outside a
     * volume are 0.  Runs on up to MAX_NTHREADS threads.
     *
     * \param series     (in)  the series given to CAT_VolRealignEstimate()
     * \param series_nii (in)  its header
     * \param dims       (in)  [nx, ny, nz] of one volume
     * \param nt         (in)  number of volumes
     * \param ref_nii    (in)  header of the reference grid (the series header
     *                         when a volume of the series was the reference)
     * \param ref_dims   (in)  [nx, ny, nz] of the reference grid
     * \param res        (in)  estimates from CAT_VolRealignEstimate()
     * \param out        (out) nt volumes of ref_dims[0]*ref_dims[1]*ref_dims[2]
     *                         floats, allocated by the caller
     * \return 0 on success, -1 if out of memory
     */
    int CAT_VolRealignResample(const float *series, const nifti_image *series_nii,
                               const int dims[3], int nt,
                               const nifti_image *ref_nii, const int ref_dims[3],
                               const CAT_VolRealignResult *res, float *out);

    /**
     * \brief Write mcflirt's text outputs next to an output base name.
     *
     * Suffixes are appended to the whole base, as mcflirt does: with the base
     * bold_mcf.nii.gz the parameters go to bold_mcf.nii.gz.par.  Files and
     * formats:
     *   mats:   <base>.mat/MAT_%04d, FSL matrices ("%f " values)
     *   plots:  <base>.par, one line of six "%g" values per volume
     *   rmsabs: <base>_abs.rms (nt lines) and <base>_abs_mean.rms
     *   rmsrel: <base>_rel.rms (nt - 1 lines) and <base>_rel_mean.rms
     * An existing .mat directory is reused.
     *
     * \param base   (in) output base name (mcflirt's -out)
     * \param res    (in) estimates from CAT_VolRealignEstimate()
     * \param mats   (in) non-zero: write the matrices
     * \param plots  (in) non-zero: write the parameters
     * \param rmsabs (in) non-zero: write the RMS against the reference
     * \param rmsrel (in) non-zero: write the RMS between consecutive volumes
     * \return 0 on success, -1 if a file could not be written
     */
    int CAT_VolRealignWrite(const char *base, const CAT_VolRealignResult *res,
                            int mats, int plots, int rmsabs, int rmsrel);

    /**
     * \brief File name of the realigned series for an output base name.
     *
     * The base itself when it ends in .nii or .nii.gz, otherwise the base
     * with .nii.gz appended (.nii when FSLOUTPUTTYPE is NIFTI), as mcflirt
     * names it.
     *
     * \param base  (in)  output base name
     * \param image (out) file name
     * \param size  (in)  size of image in bytes
     */
    void CAT_VolRealignImageName(const char *base, char *image, size_t size);

    /**
     * \brief Release the arrays of a result and zero it.
     *
     * \param res (in/out) result of CAT_VolRealignEstimate(); may be zeroed already
     */
    void CAT_VolRealignFree(CAT_VolRealignResult *res);

    /**
     * \brief Centre of mass of a volume in FSL coordinates, as FSL computes it.
     *
     * Intensities minus the minimum of the volume are the weights; the result
     * is mapped with vox2fsl_matrix().  The centre of the grid if the volume is
     * constant.
     *
     * \param vol  (in)  volume (float, x fastest)
     * \param nii  (in)  its header
     * \param dims (in)  [nx, ny, nz]
     * \param com  (out) centre of mass, FSL mm
     */
    void CAT_VolRealignCentreOfMass(const float *vol, const nifti_image *nii,
                                    const int dims[3], double com[3]);

    /**
     * \brief mcflirt's six parameters of an FSL transform.
     *
     * Euler angles of R = Rx Ry Rz (FSL's rotmat2euler), and the translation
     * that remains after rotating about centre: R * centre + t - centre.
     *
     * \param fsl    (in)  4x4 row-major rigid transform in FSL coordinates
     * \param centre (in)  rotation centre, FSL mm (the reference's centre of mass)
     * \param par    (out) rx ry rz (rad) tx ty tz (mm)
     */
    void CAT_VolRealignFslParams(const double fsl[16], const double centre[3],
                                 double par[6]);

    /**
     * \brief Jenkinson's RMS deviation between two transforms.
     *
     * The root mean square displacement of the points of a sphere of 80 mm
     * radius around centre: with D = a * b^-1 - I = [A t],
     * sqrt(|t + A centre|^2 + 80^2 / 5 * trace(A^T A)).
     *
     * \param a      (in) 4x4 row-major transform
     * \param b      (in) 4x4 row-major transform
     * \param centre (in) centre of the sphere, in the coordinates of a and b
     * \return RMS deviation (mm)
     */
    double CAT_VolRealignRms(const double a[16], const double b[16],
                             const double centre[3]);

#ifdef __cplusplus
}
#endif

#endif /* CAT_REALIGN_H */
