# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrappers for libCAT volume-processing operations.

Functions that operate on 3-D volumes accept numpy arrays plus
shape/voxelsize metadata.  Functions that need full NIfTI header
information (e.g. marching cubes, vol2surf) accept a file path
and read the volume at the C level.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from libc.float cimport FLT_MAX

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Object_types, POLYGONS,
    get_polygons_ptr, get_object_type, delete_object,
)
from cat_surf._nifti_types cimport (
    nifti_image, nifti_image_read, nifti_image_free,
    NIFTI_TYPE_FLOAT32, DT_FLOAT32,
)
from cat_surf._convert cimport PolygonsMesh, _wrap_object
from cat_surf._convert import polygons_to_arrays
from cat_surf._volume cimport VolumeHandle, open_volume

cimport cat_surf._cat_funcs as C

cnp.import_array()


# ===================================================================
# SA-NLM denoising  (mirrors CAT_VolSanlm)
# ===================================================================
def vol_sanlm(volume, bint is_rician=False, double strength=1.0):
    """
    Apply spatial-adaptive non-local means (SA-NLM) denoising.

    Mirrors ``CAT_VolSanlm``.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Input volume (Fortran order, as used by NIfTI).
    is_rician : bool
        Use Rician noise model (MRI magnitude images).
    strength : float
        Filter strength scale (default 1.0).  Values >1 increase
        smoothing; values <1 reduce it.

    Returns
    -------
    denoised : ndarray, 3-D, float32
        Denoised volume (same shape as input).
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    cdef cnp.ndarray[cnp.float32_t, ndim=3] out = vol.copy(order='F')
    cdef int dims[3]
    dims[0] = out.shape[0]
    dims[1] = out.shape[1]
    dims[2] = out.shape[2]

    C.sanlm(<float *>out.data, 3, 1,
            1 if is_rician else 0, strength, dims)
    return out


# ===================================================================
# Blood-vessel correction  (mirrors CAT_VolBloodVesselCorrection)
# ===================================================================
def vol_blood_vessel_correction(volume, voxelsize=None):
    """
    Apply blood-vessel correction to a PVE label volume.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Input PVE label volume (Fortran order).
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.

    Returns
    -------
    corrected : ndarray, 3-D, float32
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    cdef cnp.ndarray[cnp.float32_t, ndim=3] out = vol.copy(order='F')
    cdef int dims[3]
    dims[0] = out.shape[0]
    dims[1] = out.shape[1]
    dims[2] = out.shape[2]

    cdef double vx[3]
    if voxelsize is not None:
        vs = np.asarray(voxelsize, dtype=np.float64).ravel()
        vx[0] = vs[0]; vx[1] = vs[1]; vx[2] = vs[2]
    else:
        vx[0] = 1.0; vx[1] = 1.0; vx[2] = 1.0

    C.blood_vessel_correction_pve_float(<float *>out.data, dims, vx)
    return out


# ===================================================================
# Projection-based thickness (PBT)  (mirrors CAT_VolThicknessPbt)
# ===================================================================
def vol_thickness_pbt(volume, voxelsize=None,
                      int n_avgs=2, int n_median_filter=2,
                      int median_subsample=4, double range_val=0.45,
                      double fill_thresh=0.5,
                      double correct_voxelsize=0.5,
                      double sulcal_width=2.5,
                      bint fast=False, bint verbose=False):
    """
    Compute projection-based cortical thickness (PBT).

    Mirrors ``CAT_VolThicknessPbt``.  Defaults match
    ``CAT_PbtOptionsInit`` from the library.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Tissue-label volume (PVE-style).
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    n_avgs : int
        Number of averages for distance estimation (default 2).
    n_median_filter : int
        Iterations of weighted local median filtering of the PPM
        (default 2).  Set to 0 to disable.
    median_subsample : int
        Subsampling size for the median filter (default 4).
    range_val : float
        Range extension for Euclidean distance masking (default 0.45).
    fill_thresh : float
        Hole-fill threshold for the PPM (default 0.5).  Set to 0 to
        disable filling.
    correct_voxelsize : float
        Half-voxel thickness-correction strength (default 0.5).
    sulcal_width : float
        Max distance from CSF boundary for sulcal PPM correction
        (default 2.5 mm).  Set to 0 to disable.
    fast : bool
        Fast/coarse thickness estimate only (default False).
    verbose : bool

    Returns
    -------
    GMT : ndarray, 3-D, float32
        Grey-matter thickness map.
    PPM : ndarray, 3-D, float32
        Percentage-position map.
    dist_CSF : ndarray, 3-D, float32
        Distance to CSF boundary.
    dist_WM : ndarray, 3-D, float32
        Distance to WM boundary.
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    cdef int dims[3]
    dims[0] = vol.shape[0]
    dims[1] = vol.shape[1]
    dims[2] = vol.shape[2]

    cdef double vx[3]
    if voxelsize is not None:
        vs = np.asarray(voxelsize, dtype=np.float64).ravel()
        vx[0] = vs[0]; vx[1] = vs[1]; vx[2] = vs[2]
    else:
        vx[0] = 1.0; vx[1] = 1.0; vx[2] = 1.0

    cdef int nvox = dims[0] * dims[1] * dims[2]

    cdef cnp.ndarray[cnp.float32_t, ndim=3] gmt = np.zeros_like(vol, order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] ppm = np.zeros_like(vol, order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] dcsf = np.zeros_like(vol, order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] dwm = np.zeros_like(vol, order='F')

    cdef C.CAT_PbtOptions opts
    C.CAT_PbtOptionsInit(&opts)
    opts.n_avgs = n_avgs
    opts.n_median_filter = n_median_filter
    opts.median_subsample = median_subsample
    opts.range = range_val
    opts.fill_thresh = fill_thresh
    opts.correct_voxelsize = correct_voxelsize
    opts.sulcal_width = sulcal_width
    opts.fast = 1 if fast else 0
    opts.verbose = 1 if verbose else 0

    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = np.asfortranarray(vol, dtype=np.float32)

    cdef int rc = C.CAT_VolComputePbt(
        <const float *>src.data,
        <float *>gmt.data,
        <float *>ppm.data,
        <float *>dcsf.data,
        <float *>dwm.data,
        dims, vx, &opts)

    if rc != 0:
        raise RuntimeError(f"CAT_VolComputePbt returned error code {rc}")

    return gmt, ppm, dcsf, dwm


# ===================================================================
# AMAP tissue segmentation  (mirrors CAT_VolAmap)
# ===================================================================
def vol_amap(volume, labels, voxelsize=None,
             int n_pure_classes=3, int n_iters=50, int sub=96,
             bint pve=True, double weight_mrf=0.0,
             int n_iters_icm=50, bint verbose=False,
             bint use_median=False, bint use_multistep=False,
             mrf_class_weights=None):
    """
    Adaptive MAP (AMAP) brain tissue segmentation.

    Direct binding to the libCAT ``Amap`` routine (called by
    ``CAT_VolAmap``).  This wrapper does NOT include the surrounding
    pipeline (bias correction, ORNLM, cleanup) the CLI performs - run
    those steps yourself or call the CLI for an end-to-end result.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Input intensity volume.
    labels : array_like, 3-D, uint8
        Initial tissue labels (3-class: CSF/GM/WM).
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    n_pure_classes : int
        Number of pure tissue classes (default 3: CSF, GM, WM).  The
        probability map will have ``5`` channels if ``pve`` is True,
        otherwise ``n_pure_classes`` channels.
    n_iters : int
        Number of Amap iterations (default 50).
    sub : int
        Sub-sampling factor for bias/mean estimation (default 96).
    pve : bool
        Use 5-class Partial Volume Estimation (default True).
    weight_mrf : float
        MRF prior weight in [0, 1] (default 0.0).
    n_iters_icm : int
        ICM iterations (default 50).
    verbose : bool
    use_median : bool
        Use local median instead of mean for peak estimation.
    use_multistep : bool
        Enable two-step subsampling schedule.
    mrf_class_weights : array_like, optional
        Per-class MRF weights.  Length must match the number of output
        classes (3 if ``pve`` is False, 5 otherwise).

    Returns
    -------
    prob : ndarray, 4-D, uint8, shape (X, Y, Z, n_out_classes)
        Tissue probability maps.  ``n_out_classes`` is 5 if ``pve``
        is True, else ``n_pure_classes``.
    label : ndarray, 3-D, uint8
        Final tissue labels.
    mean : ndarray, shape (n_pure_classes,), float64
        Estimated pure-class means.
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    lab = np.asfortranarray(labels, dtype=np.uint8)
    if lab.shape != vol.shape:
        raise ValueError("labels shape must match volume shape")

    cdef int dims[3]
    dims[0] = vol.shape[0]
    dims[1] = vol.shape[1]
    dims[2] = vol.shape[2]

    cdef double vx[3]
    if voxelsize is not None:
        vs = np.asarray(voxelsize, dtype=np.float64).ravel()
        vx[0] = vs[0]; vx[1] = vs[1]; vx[2] = vs[2]
    else:
        vx[0] = 1.0; vx[1] = 1.0; vx[2] = 1.0

    cdef int n_out = 5 if pve else n_pure_classes
    cdef int pve_flag = 1 if pve else 0

    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = vol.copy(order='F')
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] lab_out = lab.copy(order='F')

    prob_shape = (dims[0], dims[1], dims[2], n_out)
    cdef cnp.ndarray[cnp.uint8_t, ndim=4] prob = np.zeros(
        prob_shape, dtype=np.uint8, order='F')

    cdef cnp.ndarray[cnp.float64_t, ndim=1] mean = np.zeros(
        n_pure_classes, dtype=np.float64)

    cdef cnp.ndarray[cnp.float64_t, ndim=1] mrf_w
    cdef double *mrf_w_ptr = NULL
    if mrf_class_weights is not None:
        mrf_w = np.ascontiguousarray(mrf_class_weights, dtype=np.float64)
        if mrf_w.shape[0] != n_out:
            raise ValueError(
                f"mrf_class_weights length ({mrf_w.shape[0]}) != "
                f"output classes ({n_out})")
        mrf_w_ptr = <double *>mrf_w.data

    C.Amap(<float *>src.data,
           <unsigned char *>lab_out.data,
           <unsigned char *>prob.data,
           <double *>mean.data,
           n_pure_classes, n_iters, sub, dims,
           pve_flag, weight_mrf, vx, n_iters_icm,
           1 if verbose else 0,
           1 if use_median else 0,
           mrf_w_ptr,
           1 if use_multistep else 0)

    return prob, lab_out, mean


# ===================================================================
# Marching cubes  (mirrors CAT_VolMarchingCubes)
# ===================================================================
def vol_marching_cubes(volume, double threshold=0.5,
                       double pre_fwhm=2.0, int iter_laplacian=50,
                       dist_morph=None, int n_median_filter=2,
                       int n_iter=5, double strength_gyri_mask=0.1,
                       bint fast=False, label=None,
                       bint verbose=False):
    """
    Generate a surface mesh from a volume using marching cubes
    with topology correction (genus-0).

    Mirrors ``CAT_VolMarchingCubes``.

    Parameters
    ----------
    volume : str | (ndarray, affine) | nibabel-image-like
        Input volume.  Accepts a NIfTI file path, an ``(array, affine)``
        tuple, or a nibabel-image-like object.  Auto-converted to float32.
    threshold : float
        Isosurface threshold (default 0.5).
    pre_fwhm : float
        Pre-smoothing FWHM (default 2.0).
    iter_laplacian : int
        Iterations for final Laplacian smoothing (default 50).
    dist_morph : float, optional
        Morphological opening/closing distance.  None (default) lets the
        C routine auto-estimate it (FLT_MAX sentinel, matches the CLI).
        Positive values close, negative values open.
    n_median_filter : int
        Iterations of median filtering for artefact regions (default 2).
    n_iter : int
        Maximum number of topology-correction iterations (default 5).
    strength_gyri_mask : float
        Isovalue-correction strength using the gyri/sulci mask
        (default 0.1).  Only effective with ``label``.
    fast : bool
        Use fast variant (no preprocessing/topology correction).
    label : str | (ndarray, affine) | nibabel-image-like, optional
        Label volume for gyri/sulci masking.  Must match ``volume`` shape.
    verbose : bool

    Returns
    -------
    vertices : ndarray, shape (V, 3), float64
    faces    : ndarray, shape (F, 3), int32
    """
    cdef VolumeHandle vh = open_volume(volume)
    cdef VolumeHandle vh_label = None
    cdef float *label_data = NULL
    cdef nifti_image *nii_label = NULL

    if label is not None:
        vh_label = open_volume(label)
        if (vh_label.dims[0] != vh.dims[0]
                or vh_label.dims[1] != vh.dims[1]
                or vh_label.dims[2] != vh.dims[2]):
            vh.close()
            vh_label.close()
            raise ValueError("label volume must match input volume dimensions")
        label_data = vh_label.data

    cdef double dist_morph_val
    if dist_morph is None:
        # Sentinel matching the CLI default: apply_marching_cubes compares
        # the (double) dist_morph against FLT_MAX exactly, so the wrapper
        # must pass FLT_MAX (~3.4e38), NOT DBL_MAX (sys.float_info.max,
        # ~1.8e308) — otherwise the auto-estimate branch is skipped and
        # dist_close runs with an infinite radius, wiping the mask.
        dist_morph_val = FLT_MAX
    else:
        dist_morph_val = float(dist_morph)

    cdef object_struct *result = NULL
    try:
        if fast:
            result = C.apply_marching_cubes_fast(
                vh.data, vh.nii, threshold,
                iter_laplacian, 1 if verbose else 0)
        else:
            result = C.apply_marching_cubes(
                vh.data, vh.nii, label_data,
                threshold, pre_fwhm, iter_laplacian,
                dist_morph_val, n_median_filter, n_iter,
                strength_gyri_mask, 1 if verbose else 0)
    finally:
        vh.close()
        if vh_label is not None:
            vh_label.close()

    if result == NULL:
        raise RuntimeError("apply_marching_cubes returned NULL")

    if get_object_type(result) != POLYGONS:
        delete_object(result)
        raise RuntimeError("marching cubes did not return a polygon mesh")

    cdef PolygonsMesh mesh = _wrap_object(result, owns=True)
    return polygons_to_arrays(mesh)


# ===================================================================
# Volume smoothing  (mirrors CAT_VolSmooth)
# ===================================================================
def vol_smooth(volume, voxelsize=None, double fwhm=8.0,
               bint use_mask=False):
    """
    Smooth a 3-D volume with an isotropic Gaussian kernel.

    Mirrors ``CAT_VolSmooth``.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Input volume (Fortran memory order, as used by NIfTI).
    voxelsize : array_like, shape (3,), float64, optional
        Voxel size in mm.  Default ``[1, 1, 1]``.
    fwhm : float
        Full-width at half-maximum of the Gaussian kernel in mm.
        Default 8.0.
    use_mask : bool
        If True, use masked smoothing (zero voxels are excluded and
        corrected for).  Default False.

    Returns
    -------
    smoothed : ndarray, 3-D, float32
        Smoothed volume (same shape as input).
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    cdef cnp.ndarray[cnp.float32_t, ndim=3] out = vol.copy(order='F')
    cdef int dims[3]
    dims[0] = out.shape[0]
    dims[1] = out.shape[1]
    dims[2] = out.shape[2]

    cdef double vx[3]
    if voxelsize is not None:
        vs = np.asarray(voxelsize, dtype=np.float64).ravel()
        vx[0] = vs[0]; vx[1] = vs[1]; vx[2] = vs[2]
    else:
        vx[0] = 1.0; vx[1] = 1.0; vx[2] = 1.0

    cdef double s[3]
    s[0] = fwhm; s[1] = fwhm; s[2] = fwhm

    C.smooth3(<void *>out.data, dims, vx, s,
              1 if use_mask else 0, DT_FLOAT32)
    return out


