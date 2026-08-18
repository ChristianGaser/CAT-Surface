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
# ===================================================================
# Hessian sheetness and the oriented filters  (mirrors CAT_VolSheetness,
# CAT_VolLocalStat -oriented and CAT_VolSmooth -oriented)
# ===================================================================
cdef _pack_normal(normal, int dims[3]):
    """Interleave an (X, Y, Z, 3) normal field into the layout libCAT reads.

    The C side indexes the field as ``normal[3*i + c]`` with ``i`` the
    Fortran-order voxel index, i.e. the component varies fastest.  Reading
    an (X, Y, Z, 3) array in Fortran order and reshaping to ``(nvox, 3)``
    puts voxel ``i`` on row ``i``; making that C-contiguous and ravelling
    gives exactly ``3*i + c``.
    """
    arr = np.asarray(normal, dtype=np.float32)
    if arr.ndim != 4 or arr.shape[3] != 3:
        raise ValueError("normal must have shape (X, Y, Z, 3)")
    if (arr.shape[0] != dims[0] or arr.shape[1] != dims[1]
            or arr.shape[2] != dims[2]):
        raise ValueError("normal shape must match volume shape")
    return np.ascontiguousarray(
        arr.reshape(-1, 3, order='F'), dtype=np.float32).ravel()


cdef _unpack_normal(normal_flat, int dims[3]):
    """Inverse of :func:`_pack_normal`: ``3*i + c`` back to (X, Y, Z, 3)."""
    return np.reshape(np.asarray(normal_flat).reshape(-1, 3),
                      (dims[0], dims[1], dims[2], 3), order='F')


def vol_sheetness(volume, voxelsize=None,
                  double sigma_min=0.3, double sigma_max=1.0,
                  int n_scales=3, double alpha=0.5, double beta=0.5,
                  double c=-1.0, int polarity=0,
                  bint return_normal=False, bint verbose=False):
    """
    Multi-scale Hessian sheetness (plate) filter.

    Mirrors ``CAT_VolSheetness``.  Defaults match
    ``CAT_SheetnessOptionsInit`` from the library.

    Isotropic regularizers -- a local median, a Potts MRF, total
    variation -- penalize boundary area.  A thin structure has an extreme
    area-to-volume ratio, so deleting it is always the cheaper labelling
    ("shrinking bias"), which is why one and the same median filter opens
    glued sulci in one place and closes cerebellar fissures in another.
    A sheetness filter escapes this by using a *shape* prior instead: from
    the Hessian eigenvalues ``|l1| <= |l2| <= |l3|`` a sheet has
    ``|l1| ~ |l2| ~ 0`` with ``|l3|`` large, a tube has ``|l2| ~ |l3|``
    large, and a blob has all three large.  It therefore keeps thin sheets
    and ignores blobs, and it shrinks nothing, because it makes no
    statement about boundary length.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Scalar volume, e.g. a bias-corrected T1.
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    sigma_min, sigma_max : float
        Smallest and largest Gaussian scale in mm (defaults 0.3 and 1.0).
        The range should bracket the structures being looked for: a sulcal
        CSF sheet or a gyral WM blade that survives at all is one to three
        voxels thick at 0.5 mm.  Larger scales start responding to the
        cortical ribbon itself.
    n_scales : int
        Number of log-spaced scales in the range (default 3).
    alpha : float
        Plate-vs-tube sensitivity (default 0.5).
    beta : float
        Blob-vs-plate sensitivity (default 0.5).
    c : float
        Structure-vs-noise sensitivity.  Negative (default) selects the
        automatic estimate, which makes the filter independent of the
        intensity units of the input.
    polarity : int
        ``1`` accepts only bright sheets (ridges, e.g. gyral WM blades),
        ``-1`` only dark sheets (valleys, e.g. sulcal CSF), ``0`` either
        (default).  This guard is what stops a WM-blade detector from
        responding to sulci and vice versa.
    return_normal : bool
        Also return the unit sheet normals (default False).
    verbose : bool

    Returns
    -------
    sheetness : ndarray, 3-D, float32
        Response in [0, 1], the maximum over scales.
    normal : ndarray, 4-D, float32, shape (X, Y, Z, 3)
        Unit sheet normal at the winning scale.  Only when
        ``return_normal`` is True.

    References
    ----------
    Sato et al., Med Image Anal 2(2):143-168, 1998.
    Frangi et al., MICCAI 1998, LNCS 1496:130-137.
    Descoteaux et al., Med Image Anal 10(4):638-651, 2006.
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

    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = vol
    cdef cnp.ndarray[cnp.float32_t, ndim=3] sheet = np.zeros_like(vol, order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=1] nrm
    cdef float *nrm_ptr = NULL

    if return_normal:
        nrm = np.zeros(3 * nvox, dtype=np.float32)
        nrm_ptr = <float *>nrm.data

    cdef C.CAT_SheetnessOpts opts
    C.CAT_SheetnessOptionsInit(&opts)
    opts.sigma_min = sigma_min
    opts.sigma_max = sigma_max
    opts.n_scales = n_scales
    opts.alpha = alpha
    opts.beta = beta
    opts.c = c
    opts.polarity = polarity
    opts.verbose = 1 if verbose else 0

    cdef int rc = C.CAT_VolSheetness(<const float *>src.data,
                                     <float *>sheet.data, nrm_ptr, NULL,
                                     dims, vx, &opts)
    if rc != 0:
        raise RuntimeError(f"CAT_VolSheetness returned error code {rc}")

    if return_normal:
        return sheet, _unpack_normal(nrm, dims)
    return sheet


cdef _sheet_field(volume, guide, sheetness, normal, voxelsize,
                  double sigma_min, double sigma_max, int n_scales,
                  int polarity, double strength, bint verbose):
    """Resolve the (sheetness, normal) pair the oriented filters need.

    Estimates the field from ``guide`` (or from the volume itself) unless
    both parts were supplied.  Passing the intensity image as ``guide`` is
    what makes these filters usable on a label map: a T1 still shows the
    dip through a tight sulcus long after the label map committed to pure
    grey matter.
    """
    if sheetness is None or normal is None:
        sheetness, normal = vol_sheetness(
            volume if guide is None else guide, voxelsize=voxelsize,
            sigma_min=sigma_min, sigma_max=sigma_max, n_scales=n_scales,
            polarity=polarity, return_normal=True, verbose=verbose)
    sheetness = np.asfortranarray(sheetness, dtype=np.float32)
    if 0.0 <= strength < 1.0:
        sheetness = np.asfortranarray(sheetness * strength, dtype=np.float32)
    return sheetness, normal


def vol_oriented_median(volume, guide=None, sheetness=None, normal=None,
                        voxelsize=None, int iters=1,
                        double sigma_min=0.3, double sigma_max=1.0,
                        int n_scales=3, int polarity=0,
                        double strength=1.0, bint verbose=False):
    """
    Median filter over a sheetness-oriented neighbourhood.

    Mirrors ``CAT_VolLocalStat -stat 7 -oriented``.

    Drop-in replacement for an isotropic 3x3x3 median that cannot close a
    thin sheet.  A neighbour at offset ``d`` is admitted only when
    ``1 - s*(dhat.n)**2 > 0.5``, with ``s`` the local sheetness and ``n``
    the local sheet normal.  At ``s = 0`` every neighbour is admitted and
    the operator is bit-identical to the plain isotropic median, so
    behaviour away from thin structures does not change at all; at
    ``s = 1`` only offsets within 45 degrees of the sheet plane survive,
    so the filter averages *along* the sheet and never across it.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Volume to filter.
    guide : array_like, 3-D, float32, optional
        Volume the orientation is estimated from.  Defaults to ``volume``.
        Pass the intensity image when filtering a label map.
    sheetness : array_like, 3-D, float32, optional
        Precomputed sheetness; estimated from ``guide`` when omitted.
    normal : array_like, 4-D, float32, shape (X, Y, Z, 3), optional
        Precomputed unit sheet normals; estimated alongside ``sheetness``.
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    iters : int
        Number of successive passes (default 1).
    sigma_min, sigma_max, n_scales, polarity : see :func:`vol_sheetness`
    strength : float
        How far the filter may deviate from isotropic, in [0, 1]
        (default 1.0).  0 reproduces the isotropic median exactly.
    verbose : bool

    Returns
    -------
    out : ndarray, 3-D, float32
        Filtered volume.
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    cdef int dims[3]
    dims[0] = vol.shape[0]
    dims[1] = vol.shape[1]
    dims[2] = vol.shape[2]

    sheetness, normal = _sheet_field(vol, guide, sheetness, normal, voxelsize,
                                     sigma_min, sigma_max, n_scales, polarity,
                                     strength, verbose)

    cdef cnp.ndarray[cnp.float32_t, ndim=3] out = vol.copy(order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] sheet = sheetness
    cdef cnp.ndarray[cnp.float32_t, ndim=1] nrm = _pack_normal(normal, dims)

    cdef int rc = C.CAT_VolOrientedMedian(<float *>out.data,
                                          <const float *>sheet.data,
                                          <const float *>nrm.data,
                                          NULL, dims, iters)
    if rc != 0:
        raise RuntimeError(f"CAT_VolOrientedMedian returned error code {rc}")
    return out


def vol_oriented_smooth(volume, guide=None, sheetness=None, normal=None,
                        voxelsize=None, int iters=1, double sigma=0.5,
                        double sigma_min=0.3, double sigma_max=1.0,
                        int n_scales=3, int polarity=0,
                        double strength=1.0, bint verbose=False):
    """
    Sheetness-oriented (coherence-enhancing) smoothing.

    Mirrors ``CAT_VolSmooth -oriented``.

    Diffuses along the sheet and not across it: each 3x3x3 neighbour is
    weighted by ``(1 - s) + s*exp(-(dhat.n)**2 / (2*sigma**2))``, so
    tangential neighbours keep full weight while neighbours across the
    sheet are suppressed in proportion to the local sheetness.  With
    ``s = 0`` the weights are the plain distance weights and the operator
    reduces to ordinary local averaging.  This is an iterated local
    kernel, not a Gaussian: the amount of smoothing is set by ``iters``.

    Parameters
    ----------
    volume : array_like, 3-D, float32
        Volume to filter.
    guide : array_like, 3-D, float32, optional
        Volume the orientation is estimated from.  Defaults to ``volume``.
    sheetness : array_like, 3-D, float32, optional
        Precomputed sheetness; estimated from ``guide`` when omitted.
    normal : array_like, 4-D, float32, shape (X, Y, Z, 3), optional
        Precomputed unit sheet normals; estimated alongside ``sheetness``.
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    iters : int
        Number of successive passes (default 1).
    sigma : float
        Angular width of the anisotropy (default 0.5).  Smaller values
        confine the diffusion more tightly to the plane of the sheet.
    sigma_min, sigma_max, n_scales, polarity : see :func:`vol_sheetness`
    strength : float
        How far the filter may deviate from isotropic, in [0, 1]
        (default 1.0).  0 reproduces isotropic local averaging exactly.
    verbose : bool

    Returns
    -------
    out : ndarray, 3-D, float32
        Filtered volume.
    """
    vol = np.asfortranarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("volume must be 3-D")

    cdef int dims[3]
    dims[0] = vol.shape[0]
    dims[1] = vol.shape[1]
    dims[2] = vol.shape[2]

    sheetness, normal = _sheet_field(vol, guide, sheetness, normal, voxelsize,
                                     sigma_min, sigma_max, n_scales, polarity,
                                     strength, verbose)

    cdef cnp.ndarray[cnp.float32_t, ndim=3] out = vol.copy(order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] sheet = sheetness
    cdef cnp.ndarray[cnp.float32_t, ndim=1] nrm = _pack_normal(normal, dims)

    cdef int rc = C.CAT_VolOrientedSmooth(<float *>out.data,
                                          <const float *>sheet.data,
                                          <const float *>nrm.data,
                                          NULL, dims, sigma, iters)
    if rc != 0:
        raise RuntimeError(f"CAT_VolOrientedSmooth returned error code {rc}")
    return out


# ===================================================================
# Pre-PBT repair of a PVE label map  (mirrors CAT_VolSulcusRepair)
# ===================================================================
def vol_sulcus_repair(t1, label, voxelsize=None,
                      bint recover_csf=True, bint reconnect_gyri=True,
                      bint refine_pve=False,
                      double sheet_sigma_min=0.3, double sheet_sigma_max=1.0,
                      int sheet_n_scales=3,
                      double csf_min_dist=1.5, double csf_min_wmdist=0.75,
                      double csf_thresh=0.3, double csf_strength=0.8,
                      double wm_thresh=0.3, double wm_strength=0.8,
                      int wm_max_gap=2,
                      double band_min_dist=1.5, int band_window=4,
                      double band_strength=0.7,
                      bint return_sheetness=False, bint verbose=False):
    """
    Anatomy-aware repair of a PVE label map, to be run *before* PBT.

    Mirrors ``CAT_VolSulcusRepair``.  Defaults match
    ``CAT_SulcusRepairOptionsInit`` from the library.

    Three failure modes of a tissue classifier break central-surface
    extraction, and all three are failures of *evidence*, not of
    smoothness:

    1. **Glued sulci.**  Two banks of a tight sulcus are labelled as one
       thick GM band because no CSF was detected between them.  Typical in
       the occipital midline, where cortex is thin and contrast poorest.
    2. **Broken gyri.**  A thin gyral WM blade is interrupted by a small
       missegmentation, which corrupts the WM distance map and with it the
       thickness and the central surface along the whole blade.
    3. **Residual partial-volume error** in the band where (1) happens:
       the label map reports no CSF nearby while the T1 still shows an
       intensity dip across the sulcus.

    No amount of regularization fixes these, because regularization cannot
    create evidence -- it can only redistribute what the classifier already
    committed to.  Each step here goes back to the *intensity* image and
    recovers evidence the classifier discarded, using a Hessian sheetness
    filter as the shape prior.

    Every operation is one-sided and gated: the CSF recovery can only lower
    a label, never below what the intensity supports and never within
    ``csf_min_wmdist`` of the WM boundary; the gyral reconnection only
    fires where WM is present on two *opposite* sides, so it cannot fire
    inside a sulcus; the narrow-band refit only touches voxels with a
    genuine local minimum across the sheet.

    Parameters
    ----------
    t1 : array_like, 3-D, float32
        Bias-corrected intensity image, arbitrary units (rescaled
        internally onto the label axis).
    label : array_like, 3-D, float32
        PVE label image in [0, 3] with CSF = 1, GM = 2, WM = 3 -- the same
        convention :func:`vol_thickness_pbt` expects.
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    recover_csf : bool
        Run step 1, the sulcal CSF recovery (default True).
    reconnect_gyri : bool
        Run step 2, the gyral WM reconnection (default True).
    refine_pve : bool
        Run step 3, the narrow-band PVE refit (default False; it is the
        most aggressive of the three).
    sheet_sigma_min, sheet_sigma_max, sheet_n_scales : float, float, int
        Scale range of the sheetness filter; see :func:`vol_sheetness`.
    csf_min_dist : float
        Act on sulcal CSF only where the label map finds no CSF within
        this many mm (default 1.5).
    csf_min_wmdist : float
        Never carve closer than this many mm to the WM boundary
        (default 0.75).
    csf_thresh, csf_strength : float
        Sheetness threshold and blend weight for step 1.
    wm_thresh, wm_strength : float
        Sheetness threshold and blend weight for step 2.
    wm_max_gap : int
        Largest gap half-width in voxels the two-sided test will bridge
        (default 2).
    band_min_dist, band_window, band_strength : float, int, float
        Parameters of step 3.
    return_sheetness : bool
        Also return the sheetness map of the last executed step.
    verbose : bool

    Returns
    -------
    label : ndarray, 3-D, float32
        Repaired PVE label image.
    sheetness : ndarray, 3-D, float32
        Only when ``return_sheetness`` is True.

    References
    ----------
    Han et al., Proc SPIE Med Imag 4322:194-203, 2001 (ACE).
    Han et al., NeuroImage 23(3):997-1012, 2004 (CRUISE).
    Kim et al., NeuroImage 27(1):210-221, 2005 (CLASP, CSF skeleton).
    """
    vol_t1 = np.asfortranarray(t1, dtype=np.float32)
    vol_lab = np.asfortranarray(label, dtype=np.float32)
    if vol_t1.ndim != 3:
        raise ValueError("t1 must be 3-D")
    if vol_lab.shape != vol_t1.shape:
        raise ValueError("label shape must match t1 shape")

    cdef int dims[3]
    dims[0] = vol_t1.shape[0]
    dims[1] = vol_t1.shape[1]
    dims[2] = vol_t1.shape[2]

    cdef double vx[3]
    if voxelsize is not None:
        vs = np.asarray(voxelsize, dtype=np.float64).ravel()
        vx[0] = vs[0]; vx[1] = vs[1]; vx[2] = vs[2]
    else:
        vx[0] = 1.0; vx[1] = 1.0; vx[2] = 1.0

    cdef C.CAT_SulcusRepairOpts opts
    C.CAT_SulcusRepairOptionsInit(&opts)
    opts.sheet_sigma_min = sheet_sigma_min
    opts.sheet_sigma_max = sheet_sigma_max
    opts.sheet_n_scales = sheet_n_scales
    opts.csf_min_dist = csf_min_dist
    opts.csf_min_wmdist = csf_min_wmdist
    opts.csf_thresh = csf_thresh
    opts.csf_strength = csf_strength
    opts.wm_thresh = wm_thresh
    opts.wm_strength = wm_strength
    opts.wm_max_gap = wm_max_gap
    opts.band_min_dist = band_min_dist
    opts.band_window = band_window
    opts.band_strength = band_strength
    opts.verbose = 1 if verbose else 0

    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = vol_t1
    cdef cnp.ndarray[cnp.float32_t, ndim=3] lab = vol_lab.copy(order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] sheet
    cdef float *sheet_ptr = NULL
    cdef int rc

    if return_sheetness:
        sheet = np.zeros_like(lab, order='F')
        sheet_ptr = <float *>sheet.data

    # Order matters: the CSF recovery changes the label map the gyral
    # reconnection reads, and both change the CSF distance the narrow-band
    # refit is gated on.
    if recover_csf:
        rc = C.CAT_VolRecoverSulcalCSF(<const float *>src.data,
                                       <float *>lab.data, sheet_ptr,
                                       dims, vx, &opts)
        if rc != 0:
            raise RuntimeError(f"CAT_VolRecoverSulcalCSF returned error code {rc}")

    if reconnect_gyri:
        rc = C.CAT_VolReconnectGyri(<const float *>src.data,
                                    <float *>lab.data, sheet_ptr,
                                    dims, vx, &opts)
        if rc != 0:
            raise RuntimeError(f"CAT_VolReconnectGyri returned error code {rc}")

    if refine_pve:
        rc = C.CAT_VolRefinePveNarrowBand(<const float *>src.data,
                                          <float *>lab.data, dims, vx, &opts)
        if rc != 0:
            raise RuntimeError(f"CAT_VolRefinePveNarrowBand returned error code {rc}")

    if return_sheetness:
        return lab, sheet
    return lab


def vol_thickness_pbt(volume, voxelsize=None,
                      int n_avgs=2, int n_median_filter=2,
                      int median_subsample=4, double range_val=0.45,
                      double fill_thresh=0.5,
                      double correct_thickness=0.25,
                      double sulcal_width=2.5,
                      bint pve_distance=False,
                      bint oriented_filter=False,
                      double oriented_strength=1.0,
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
    correct_thickness : float
        Additive thickness correction in mm (default 0.25), compensating
        the systematic border shift of the segmentation.  Replaces the
        former ``correct_voxelsize``, which was given in voxels and so
        scaled with ``voxelsize``; to reproduce a previous call, pass
        ``correct_voxelsize * mean(voxelsize)``.
    sulcal_width : float
        Max distance from CSF boundary for sulcal PPM correction
        (default 2.5 mm).  Set to 0 to disable.
    pve_distance : bool
        EXPERIMENTAL, default False.  Correct the WM/CSF distance maps
        for partial volume: the distance transform measures to the
        nearest source voxel *centre*, and this uses that voxel's
        partial volume to measure to the tissue boundary instead.
        Mirrors the correction in CAT12's ``cat_vol_pbtsimpleCS4``.
        Shifts the thickness calibration, so ``correct_thickness``
        has to be re-derived when enabling it.
    oriented_filter : bool
        Replace the isotropic 3x3x3 median filters by sheetness-oriented
        ones (default False).  An isotropic median penalizes boundary
        area, so it removes thin structures whichever side of the label
        boundary they lie on: the same filter that opens a glued sulcus
        closes a cerebellar fissure.  The oriented variant admits only
        neighbours lying in the plane of the locally detected sheet, so
        it averages along a thin structure and never across it.  Where
        no sheet is detected it is identical to the isotropic median.
    oriented_strength : float
        How far the oriented filters may deviate from isotropic, in
        [0, 1] (default 1.0).  0 reproduces the isotropic filters
        exactly.
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
    opts.correct_thickness = correct_thickness
    opts.sulcal_width = sulcal_width
    opts.pve_distance = 1 if pve_distance else 0
    opts.oriented_filter = 1 if oriented_filter else 0
    opts.oriented_strength = oriented_strength
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
             mrf_class_weights=None, int mrf_aniso=0,
             double mrf_aniso_strength=1.0, double mrf_aniso_sigma=0.5,
             sheetness=None, normal=None):
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
    mrf_aniso : int
        Anisotropic MRF prior: ``0`` off (default), ``1`` locally
        varying beta, ``2`` direction-weighted Potts.  The isotropic
        Potts prior penalizes boundary area, so it always finds it
        cheaper to delete a thin structure than to keep it -- the reason
        a stronger MRF closes cerebellar fissures and glues tight sulci
        while it removes noise.  Mode 1 scales beta by
        ``1 - strength * sheetness``, so the prior stops pulling on a
        sheet.  Mode 2 keeps beta and instead down-weights neighbours
        lying *across* the sheet, so the prior still denoises within a
        sulcal bank but can no longer merge two banks facing each other.
        Both are exact no-ops where the sheetness is zero.  Requires
        ``weight_mrf > 0`` and ``n_iters_icm > 0``.
    mrf_aniso_strength : float
        Strength of the anisotropic relaxation in [0, 1] (default 1.0).
        0 reproduces the isotropic prior exactly.
    mrf_aniso_sigma : float
        Angular width of the direction weighting for ``mrf_aniso=2``
        (default 0.5).  Smaller values confine the prior more tightly to
        the plane of the sheet.
    sheetness : array_like, 3-D, float32, optional
        Precomputed sheetness field in [0, 1].  Estimated from
        ``volume`` with :func:`vol_sheetness` defaults when omitted and
        ``mrf_aniso`` is non-zero.
    normal : array_like, 4-D, float32, shape (X, Y, Z, 3), optional
        Precomputed unit sheet normals, required by ``mrf_aniso=2`` and
        estimated alongside ``sheetness`` when omitted.

    Returns
    -------
    prob : ndarray, 4-D, uint8, shape (X, Y, Z, n_out_classes)
        Tissue probability maps.  ``n_out_classes`` is 5 if ``pve``
        is True, else ``n_pure_classes``.
    label : ndarray, 3-D, uint8
        Final tissue labels.
    mean : ndarray, shape (n_out_classes,), float64
        Estimated class means, in the class order of ``CAT_Amap.h``.  With
        ``pve`` that is ``[CSF, CSF/GM, GM, GM/WM, WM]``, so the pure-class
        means are at indices 0, 2 and 4; without it, ``[CSF, GM, WM]``.
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

    # Amap() raises its class count to 5 in the PVE branch and rewrites the
    # caller's mean[0..4] there (CAT_Amap.c), so the buffer has to hold the
    # mixed classes too -- CAT_VolAmap.c sizes it the same way.  Allocating
    # only n_pure_classes overruns it by 16 bytes on every pve=True call.
    cdef cnp.ndarray[cnp.float64_t, ndim=1] mean = np.zeros(
        max(n_out, n_pure_classes), dtype=np.float64)

    cdef cnp.ndarray[cnp.float64_t, ndim=1] mrf_w
    cdef double *mrf_w_ptr = NULL
    if mrf_class_weights is not None:
        mrf_w = np.ascontiguousarray(mrf_class_weights, dtype=np.float64)
        if mrf_w.shape[0] != n_out:
            raise ValueError(
                f"mrf_class_weights length ({mrf_w.shape[0]}) != "
                f"output classes ({n_out})")
        mrf_w_ptr = <double *>mrf_w.data

    cdef C.CAT_MrfAnisoField aniso
    cdef C.CAT_MrfAnisoField *aniso_ptr = NULL
    cdef cnp.ndarray[cnp.float32_t, ndim=3] sheet_arr
    cdef cnp.ndarray[cnp.float32_t, ndim=1] normal_flat

    if mrf_aniso != 0 and weight_mrf > 0.0 and n_iters_icm > 0:
        if sheetness is None or (mrf_aniso == 2 and normal is None):
            sheet_est, normal_est = vol_sheetness(
                vol, voxelsize=voxelsize, polarity=0, return_normal=True,
                verbose=verbose)
            if sheetness is None:
                sheetness = sheet_est
            if normal is None:
                normal = normal_est

        sheet_arr = np.asfortranarray(sheetness, dtype=np.float32)
        if (sheet_arr.shape[0] != dims[0] or sheet_arr.shape[1] != dims[1]
                or sheet_arr.shape[2] != dims[2]):
            raise ValueError("sheetness shape must match volume shape")

        aniso.mode = mrf_aniso
        aniso.sheetness = <const float *>sheet_arr.data
        aniso.normal = NULL
        aniso.strength = mrf_aniso_strength
        aniso.sigma = mrf_aniso_sigma

        if mrf_aniso == 2:
            normal_flat = _pack_normal(normal, dims)
            aniso.normal = <const float *>normal_flat.data

        aniso_ptr = &aniso

    C.AmapAniso(<float *>src.data,
                <unsigned char *>lab_out.data,
                <unsigned char *>prob.data,
                <double *>mean.data,
                n_pure_classes, n_iters, sub, dims,
                pve_flag, weight_mrf, vx, n_iters_icm,
                1 if verbose else 0,
                1 if use_median else 0,
                mrf_w_ptr,
                1 if use_multistep else 0,
                aniso_ptr)

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


