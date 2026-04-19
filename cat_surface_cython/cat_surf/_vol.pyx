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

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Object_types, POLYGONS,
    get_polygons_ptr, get_object_type, delete_object,
)
from cat_surf._nifti_types cimport (
    nifti_image, nifti_image_read, nifti_image_free,
    NIFTI_TYPE_FLOAT32,
)
from cat_surf._convert cimport PolygonsMesh, _wrap_object
from cat_surf._convert import polygons_to_arrays

cimport cat_surf._cat_funcs as C

cnp.import_array()


# ===================================================================
# SA-NLM denoising  (mirrors CAT_VolSanlm)
# ===================================================================
def vol_sanlm(volume, bint is_rician=False, double strength=3.0):
    """
    Apply spatial-adaptive non-local means (SA-NLM) denoising.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Input volume (Fortran order, as used by NIfTI).
    is_rician : bool
        Use Rician noise model (MRI magnitude images).
    strength : float
        Denoising strength (default 3.0).

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
                      int n_avgs=2, int n_median_filter=1,
                      int median_subsample=0, double range_val=100.0,
                      double fill_thresh=0.5,
                      double correct_voxelsize=0.0,
                      double sulcal_width=0.0,
                      bint fast=False, bint verbose=False):
    """
    Compute projection-based cortical thickness (PBT).

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Tissue-label volume (PVE-style).
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    n_avgs : int
    n_median_filter : int
    median_subsample : int
    range_val : float
    fill_thresh : float
    correct_voxelsize : float
    sulcal_width : float
    fast : bool
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
             int n_classes=3, int n_iters=200, int sub=16,
             int pve=5, double weight_mrf=0.15,
             int n_iters_icm=50, bint verbose=False,
             bint use_median=False, bint use_multistep=False):
    """
    Adaptive MAP (AMAP) brain tissue segmentation.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Input intensity volume.
    labels : array_like, 3-D, uint8
        Initial tissue labels.
    voxelsize : array_like, shape (3,), float64, optional
    n_classes : int
        Number of tissue classes (default 3: CSF, GM, WM).
    n_iters : int
    sub : int
        Sub-sampling factor for bias estimation.
    pve : int
        PVE model (5 = Pve5).
    weight_mrf : float
    n_iters_icm : int
    verbose : bool
    use_median : bool
    use_multistep : bool

    Returns
    -------
    prob : ndarray, 4-D, uint8, shape (X, Y, Z, n_classes)
        Tissue probability maps.
    label : ndarray, 3-D, uint8
        Final tissue labels.
    mean : ndarray, shape (n_classes,), float64
        Estimated class means.
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

    cdef int nvox = dims[0] * dims[1] * dims[2]

    # Amap modifies src in-place
    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = vol.copy(order='F')
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] lab_out = lab.copy(order='F')

    # prob: nvox * n_classes
    prob_shape = (dims[0], dims[1], dims[2], n_classes)
    cdef cnp.ndarray[cnp.uint8_t, ndim=4] prob = np.zeros(
        prob_shape, dtype=np.uint8, order='F')

    cdef cnp.ndarray[cnp.float64_t, ndim=1] mean = np.zeros(
        n_classes, dtype=np.float64)

    C.Amap(<float *>src.data,
           <unsigned char *>lab_out.data,
           <unsigned char *>prob.data,
           <double *>mean.data,
           n_classes, n_iters, sub, dims,
           pve, weight_mrf, vx, n_iters_icm,
           1 if verbose else 0,
           1 if use_median else 0,
           NULL,   # mrf_class_weights
           1 if use_multistep else 0)

    return prob, lab_out, mean


# ===================================================================
# Marching cubes  (mirrors CAT_VolMarchingCubes)
# ===================================================================
def vol_marching_cubes(str volume_file, double threshold=0.5,
                       double pre_fwhm=1.0, int iter_laplacian=100,
                       double dist_morph=1.0, int n_median_filter=3,
                       int n_iter=5, double strength_gyri_mask=0.0,
                       bint fast=False, str label_file=None,
                       bint verbose=False):
    """
    Generate a surface mesh from a volume using marching cubes
    with topology correction (genus-0).

    Parameters
    ----------
    volume_file : str
        Path to input NIfTI volume.
    threshold : float
        Isosurface threshold.
    pre_fwhm : float
        Pre-smoothing FWHM.
    iter_laplacian : int
    dist_morph : float
    n_median_filter : int
    n_iter : int
    strength_gyri_mask : float
    fast : bool
        Use fast variant (fewer parameters).
    label_file : str, optional
        Path to a label NIfTI volume.
    verbose : bool

    Returns
    -------
    vertices : ndarray, shape (V, 3), float64
    faces    : ndarray, shape (F, 3), int32
    """
    cdef bytes bvol = volume_file.encode("utf-8")
    cdef nifti_image *nii = nifti_image_read(bvol, 1)
    if nii == NULL:
        raise IOError(f"Failed to read NIfTI volume '{volume_file}'")

    cdef float *label_data = NULL
    cdef nifti_image *nii_label = NULL

    if label_file is not None:
        blabel = label_file.encode("utf-8")
        nii_label = nifti_image_read(blabel, 1)
        if nii_label != NULL:
            label_data = <float *>nii_label.data

    cdef object_struct *result

    if fast:
        result = C.apply_marching_cubes_fast(
            <float *>nii.data, nii, threshold,
            iter_laplacian, 1 if verbose else 0)
    else:
        result = C.apply_marching_cubes(
            <float *>nii.data, nii, label_data,
            threshold, pre_fwhm, iter_laplacian,
            dist_morph, n_median_filter, n_iter,
            strength_gyri_mask, 1 if verbose else 0)

    if nii_label != NULL:
        nifti_image_free(nii_label)
    nifti_image_free(nii)

    if result == NULL:
        raise RuntimeError("apply_marching_cubes returned NULL")

    if get_object_type(result) != POLYGONS:
        delete_object(result)
        raise RuntimeError("marching cubes did not return a polygon mesh")

    cdef PolygonsMesh mesh = _wrap_object(result, owns=True)
    return polygons_to_arrays(mesh)



