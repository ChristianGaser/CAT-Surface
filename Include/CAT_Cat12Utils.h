/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * C ports of the CAT12/CAT26 MATLAB helpers that the safeguards of CAT's
 * geodesic-shooting registration (run_Shooting in cat_main_registration.m:
 * separation of very large ventricles, filtering of Jacobian hot spots) call:
 * cat_vol_smooth3X, cat_vol_div, cat_vol_morph ('l', 'd', 'o', 'ldc', 'ldo'),
 * cat_vol_localstat (mean), cat_vol_laplace3R and cat_vbdist (Christian Gaser
 * and Robert Dahnke).
 *
 * All volumes are stored in MATLAB order (x fastest) with dims {nx, ny, nz}.
 * Binary volumes are unsigned char; any non-zero input voxel is foreground
 * (callers threshold float data as CAT does, vol > 0.5) and the output is
 * 0/1.  The morphological operations and the distance transform reproduce the
 * MATLAB originals voxel for voxel, including the crop to the bounding box and
 * the zero padding that decide what happens at the volume border; the
 * smoothing and the divergence agree to float rounding.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef _CAT_CAT12UTILS_H_
#define _CAT_CAT12UTILS_H_

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * \brief Gaussian smoothing of CAT12's cat_vol_smooth3X(S, s) (filter 's').
     *
     * NaN and +-Inf become 0, before and after the filter.  Depending on s
     * (in voxels):
     *  - 0 < s < 0.5: smooth3(S,'gaussian',3,0.5)*s + S*(1-s)
     *  - 0.5 <= s <= 1: smooth3(S,'gaussian',3,s)
     *  - s > 1 and every dimension > 6: the volume is halved (odd sizes padded
     *    by repeating the last plane, 2x2x2 means by cat_vol_interp3f), smoothed
     *    with s/2 recursively and interpolated back linearly to the original
     *    size (cat_vol_interp3f, borders clamped)
     *  - s >= 1 and a dimension <= 6: smooth3(S,'gaussian',5,s)
     *  - s == 0: unchanged; s < 0: unchanged; s NaN: all 0 (as in MATLAB)
     *
     * smooth3 is MATLAB's: a normalised Gaussian of size sz and standard
     * deviation sd voxels, the image padded by replicating its border, 'valid'
     * convolution.  The kernel is rounded to single as MATLAB's convn does, the
     * sum is accumulated in double; MATLAB accumulates in single in an
     * unpublished order, so the results differ by float rounding only (a few
     * 1e-7 of the data range).  The resampling steps are bit-identical to the
     * cat_vol_interp3f mex.  A volume with nz == 1 is filtered as a volume;
     * MATLAB takes its 2-D path (and fails for a non-integer 2*s+1).
     *
     * \param vol  (in/out) float volume, smoothed in place
     * \param dims (in)     {nx, ny, nz}
     * \param s    (in)     filter size in voxels (the original takes the mean of
     *                      a vector)
     */
    void cat12_vol_smooth3X(float *vol, const int dims[3], double s);

    /**
     * \brief Divergence of CAT12's cat_vol_div(Ym) with all defaults.
     *
     * The defaults are vx_vol = vx_volr = 1.5, norm = 0, sm = 0, the way
     * run_Shooting calls it on a velocity component:
     *  - NaN voxels within 2 voxels of a finite one (cat_vbdist inside
     *    cat_vol_morph(~Ynan,'d',2)) take the value of the nearest finite voxel
     *  - crop to the bounding box of (Ym ~= 0 & ~Ynan) plus one voxel
     *    (cat_vol_resize 'reduceBrain' with d = 2 and vx_vol = 1.5); an empty
     *    box keeps the whole volume
     *  - 'reduceV' to 1.5 mm is a no-op (vx_vol == vx_volr)
     *  - gradient of max(1/3, Ym) (NaN -> 1/3), each component divided by 1.5,
     *    and the derivative of each along its own axis, both by
     *    cat_vol_gradient3 (central differences / 2, one-sided differences / 2
     *    at the crop border)
     *  - 'dereduceV' with 'cubic' at the grid points: NaN -> 0, otherwise the
     *    identity, except that a +-Inf within the 4x4x4 cubic window makes the
     *    voxel NaN, as in cat_vol_interp3f
     *  - outside the crop and at the original NaN voxels the output is 0.
     *
     * The clamp at 1/3 makes this the divergence of the positive part of the
     * input only.  The arithmetic follows the original step by step (single
     * precision, no fused multiply-add), so the result is bit-identical to
     * MATLAB's.
     *
     * \param vol  (in)  float volume (one velocity component in voxels)
     * \param dims (in)  {nx, ny, nz}
     * \param div  (out) float volume of the same size; may be the same array as
     *                   vol
     */
    void cat12_vol_div(const float *vol, const int dims[3], float *div);

    /**
     * \brief Largest component of CAT12's cat_vol_morph(vol, 'l').
     *
     * Keeps the largest 6-connected component (spm_bwlabel(vol,6)).  Of
     * equally large components the one whose first voxel comes first in
     * memory order wins (spm_bwlabel numbers the components in that order and
     * MATLAB's sort is stable).  An empty volume stays empty.
     *
     * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
     * \param dims (in)     {nx, ny, nz}
     */
    void cat12_vol_morph_lab(unsigned char *vol, const int dims[3]);

    /**
     * \brief Several large components of CAT12's cat_vol_morph(vol, 'l', [nmax th]).
     *
     * Keeps the largest 6-connected component and, of the next ones in
     * decreasing size up to the nmax-th, those larger than th times the whole
     * foreground (0 < th < 1) or larger than |th| voxels (th < 0 or th > 1).
     * MATLAB returns the rank of each kept component as a single volume; this
     * returns their union, which is what the callers use after the next
     * threshold at 0.5.  cat12_vol_morph_lab() is the case nmax <= 1.  Used by
     * run_Shooting with [10 0.1] when more than two classes are registered.
     *
     * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
     * \param dims (in)     {nx, ny, nz}
     * \param nmax (in)     largest number of components to keep (as n(1))
     * \param th   (in)     relative or absolute size limit (as n(2))
     */
    void cat12_vol_morph_lab_multi(unsigned char *vol, const int dims[3], double nmax,
                                   double th);

    /**
     * \brief Dilation of CAT12's cat_vol_morph(vol, 'd', n, vx) ('cdilate').
     *
     * Dilation with a box of 2*round(n/vx)+1 voxels per axis (convn with a box
     * kernel, > 0); outside the volume is background, so nothing grows in from
     * the border.  Nothing happens if round(n) <= 0.  The crop to the bounding
     * box plus round((round(n)+1)/vx) voxels that the original uses for speed
     * does not change the result.
     *
     * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
     * \param dims (in)     {nx, ny, nz}
     * \param n    (in)     size in mm (need not be an integer)
     * \param vx   (in)     voxel size in mm; NULL is {1,1,1}, the default of the
     *                      original
     */
    void cat12_vol_morph_dilate(unsigned char *vol, const int dims[3], double n,
                                const double vx[3]);

    /**
     * \brief Opening of CAT12's cat_vol_morph(vol, 'o', n, vx) ('copen').
     *
     * ~cclose(~vol) with a box of 2*round(round(n)/vx)+1 voxels.  The closing
     * of the background is computed in its bounding box plus one voxel, padded
     * by ceil(round(n)/vx) voxels of zeros, so the region outside the volume
     * counts as foreground: objects touching the border are not eroded from
     * outside.  Nothing happens if round(n) <= 0.
     *
     * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
     * \param dims (in)     {nx, ny, nz}
     * \param n    (in)     size in mm
     * \param vx   (in)     voxel size in mm; NULL is {1,1,1}
     */
    void cat12_vol_morph_open(unsigned char *vol, const int dims[3], double n,
                              const double vx[3]);

    /**
     * \brief Closing and hole filling of CAT12's cat_vol_morph(vol, 'ldc', n, vx).
     *
     * 'distclose' followed by ~lab(~vol), which fills every background
     * component except the largest.  'distclose' crops to the bounding box
     * plus round((round(n)+1)/vx) voxels, pads by ceil(n/vx) voxels of zeros
     * and ORs the closing with the input.  For round(n) > 5 the closing is
     * cat_vbdist(vol) > n followed by cat_vbdist from that background, inside
     * the dilated region, >= n; otherwise a dilation and an erosion with the
     * kernel max(0, n - max(0, d - 0.5)) (d = distance in mm), each keeping
     * voxels whose convolution with the binary volume is >= 0.5.
     *
     * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
     * \param dims (in)     {nx, ny, nz}
     * \param n    (in)     distance in mm (run_Shooting: 12/tempres2)
     * \param vx   (in)     voxel size in mm; NULL is {1,1,1} (run_Shooting)
     */
    void cat12_vol_morph_labdistclose(unsigned char *vol, const int dims[3], double n,
                                      const double vx[3]);

    /**
     * \brief Opening and largest component of CAT12's cat_vol_morph(vol, 'ldo', n, vx).
     *
     * 'distopen' (= ~distclose(~vol), see cat12_vol_morph_labdistclose(), so
     * the outside of the volume counts as foreground) followed by the largest
     * 6-connected component.
     *
     * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
     * \param dims (in)     {nx, ny, nz}
     * \param n    (in)     distance in mm (run_Shooting: 3)
     * \param vx   (in)     voxel size in mm; NULL is {1,1,1} (run_Shooting:
     *                      tempres2 on every axis)
     */
    void cat12_vol_morph_labdistopen(unsigned char *vol, const int dims[3], double n,
                                     const double vx[3]);

    /**
     * \brief Local mean of CAT12's cat_vol_localstat(V, M, nb, 1).
     *
     * For every voxel of the mask with a finite, non-zero value: the mean of
     * the finite, non-zero values of the mask within the sphere of radius nb
     * voxels (nb = 1: the voxel and its six face neighbours) that lie inside
     * the volume.  Every other voxel is 0.  Zeros are dropped because the
     * original's filter0 = 0 removes them from the mask.  Accumulated in single
     * precision in the original's order, so the result is bit-identical.
     * libCAT's localstat3() differs: it replicates the border, keeps zeros,
     * leaves masked voxels unchanged and sums in double.
     *
     * \param vol  (in)  float volume
     * \param mask (in)  mask (non-zero = in); NULL is all voxels
     * \param dims (in)  {nx, ny, nz}
     * \param nb   (in)  radius in voxels, 1..10
     * \param out  (out) float volume of the same size; must not alias vol
     */
    void cat12_vol_localstat_mean(const float *vol, const unsigned char *mask,
                                  const int dims[3], int nb, float *out);

    /**
     * \brief Laplace filter of CAT12's cat_vol_laplace3R(V, M, TH).
     *
     * libCAT's laplace3R() is the same algorithm (Jacobi sweeps over the
     * face neighbours of the voxels of M, until the largest change is <= TH or
     * after 2000 sweeps; NaN becomes FLT_MAX and is not used as a neighbour).
     * The only difference is that the mex rounds TH to single precision, which
     * this wrapper does before calling laplace3R().
     *
     * \param vol  (in/out) float volume, filtered in place
     * \param mask (in)     filter region (non-zero = filter)
     * \param dims (in)     {nx, ny, nz}
     * \param th   (in)     largest change of a sweep that ends the iteration
     *                      (0.000001 <= th < 0.5 in the original)
     */
    void cat12_vol_laplace3R(float *vol, const unsigned char *mask, const int dims[3],
                             double th);

    /**
     * \brief Euclidean distance transform of CAT12's [D, I] = cat_vbdist(V, M, vx).
     *
     * Distance from every voxel of M to the nearest object voxel (V >= 0.5,
     * i.e. round(V) >= 1 for 0 <= V < 255.5) with a forward and a backward
     * sweep over a 13-voxel half neighbourhood, propagating the index of the
     * nearest object voxel and recomputing the exact distance to it.  Object
     * voxels get 0; voxels outside M and voxels the sweeps do not reach keep
     * FLT_MAX.  libCAT's euclidean_distance() runs the same sweeps but returns
     * 0 for unreached and masked voxels, uses V > 1e-6 as object, computes the
     * diagonal step lengths and the distances in single precision and does not
     * return the index, so it is not reused here.
     *
     * \param vol  (in)  float volume, object where >= 0.5
     * \param mask (in)  voxels to compute (non-zero); NULL is all voxels
     * \param dims (in)  {nx, ny, nz}
     * \param vx   (in)  voxel size; NULL is {1,1,1}
     * \param dist (out) float distances (FLT_MAX: not computed)
     * \param idx  (out) 0-based linear index of the nearest object voxel (the
     *                   voxel itself if none was found); NULL to skip
     */
    void cat12_vbdist(const float *vol, const unsigned char *mask, const int dims[3],
                      const double vx[3], float *dist, unsigned int *idx);

#ifdef __cplusplus
}
#endif

#endif
