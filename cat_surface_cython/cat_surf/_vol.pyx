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
# Blood-vessel correction
#
# The standalone CAT_VolBloodVesselCorrection tool was removed; this is
# still the correction CAT_VolThicknessPbt applies to its input PVE map
# (its -no-bvc switch turns it off), so the array API is kept for anyone
# who wants to inspect or apply it on its own.
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
# Hessian sheetness and the oriented median  (mirrors CAT_VolSheetness
# and CAT_VolLocalStat -oriented)
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
                  double c=-1.0, double gain=1.0, int polarity=0,
                  double normalize=1.0, bint skeletonize=False,
                  bint signed_response=False,
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
        intensity units of the input.  That estimate is half the largest
        Hessian norm in the volume, so on a whole head it is set by the
        scalp/air step and the cortical response collapses towards zero;
        lowering it is the principled way to bring the response back up.
    gain : float
        Overall gain on the response (default 1.0).  The blunt alternative
        to ``c`` when the intensity units are unknown: the map is
        multiplied by this and clamped to ``[0, 1]``.  Values above 1
        amplify a response too weak to reach the thresholds the consumers
        gate on -- notably the hard 0.5 of :func:`vol_oriented_median`,
        below which it is exactly the isotropic median.  Because the map
        is linear and fixes zero, a zero response stays zero at any gain,
        so the oriented operators remain exact no-ops where no sheet was
        found.
    polarity : int
        ``1`` accepts only bright sheets (ridges, e.g. gyral WM blades),
        ``-1`` only dark sheets (valleys, e.g. sulcal CSF), ``0`` either
        (default).  This guard is what stops a WM-blade detector from
        responding to sulci and vice versa.
    signed_response : bool
        Return the polarity as a sign (default False).  With ``polarity=0``
        the filter answers on both kinds of sheet but returns a magnitude,
        so a consumer cannot tell a sulcus from a gyral blade.  This makes
        a valley (dark sheet, a dip in a PPM) come out negative and a ridge
        (bright sheet, a white-matter blade) positive, in [-1, 1] — so the
        two move in opposite directions under a single operation.  Adding
        the map to an image lowers it along sulci and raises it along blades
        at once, which a global isovalue shift cannot do.  Do not hand a
        signed map to the oriented filters: they clamp to [0, 1].
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
    opts.gain = gain
    opts.polarity = polarity
    opts.normalize = normalize
    opts.skeletonize = 1 if skeletonize else 0
    opts.signed_response = 1 if signed_response else 0
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
    if strength != 1.0:
        # Same gain the C side applies, done here so that it also reaches a
        # precomputed field, which never passes through CAT_VolSheetness.
        sheetness = np.asfortranarray(
            np.clip(sheetness * strength, 0.0, 1.0), dtype=np.float32)
    return sheetness, normal


def vol_oriented_median(volume, guide=None, sheetness=None, normal=None,
                        voxelsize=None, int iters=1, double cutoff=0.0,
                        double sigma_min=0.3, double sigma_max=1.0,
                        int n_scales=3, int polarity=0,
                        double strength=1.0, bint verbose=False):
    """
    Median filter over a sheetness-oriented neighbourhood.

    Mirrors ``CAT_VolLocalStat -stat 7 -oriented``.

    Drop-in replacement for an isotropic 3x3x3 median that cannot close a
    thin sheet.  A neighbour at offset ``d`` is admitted only when
    ``s*(dhat.n)**2 < cutoff``, with ``s`` the local sheetness and ``n``
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
    cutoff : float
        Admission cutoff (default 0, which selects 0.10).  The 9 offsets
        lying in the sheet plane are always admitted; the 6 face
        neighbours drop out at ``s = cutoff``, the 12 edge neighbours at
        ``2*cutoff`` and the 8 corners at ``3*cutoff``.  A one-voxel-thick
        sheet is preserved from ``s = 2*cutoff`` upwards, where the
        in-plane offsets first make up half the admitted set, so the
        cutoff should be about half the sheetness the data reaches.
    sigma_min, sigma_max, n_scales, polarity : see :func:`vol_sheetness`
    strength : float
        Overall gain on the sheetness before the filter uses it
        (default 1.0).  0 reproduces the isotropic median exactly.
        Values above 1 amplify a response too weak to matter: this filter
        admits every neighbour while the sheetness stays below ``cutoff``,
        and only preserves a one-voxel sheet from ``2*cutoff`` upwards.
        See ``gain`` in :func:`vol_sheetness`.
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
                                          NULL, dims, cutoff, iters)
    if rc != 0:
        raise RuntimeError(f"CAT_VolOrientedMedian returned error code {rc}")
    return out


# ===================================================================
# ===================================================================
def vol_open_ppm_sulci(ppm, voxelsize=None, double isovalue=0.5,
                       sheetness=None, normal=None,
                       double sigma_min=0.3, double sigma_max=1.0,
                       int n_scales=3, double sheet_normalize=1.0,
                       bint sheet_skeleton=False,
                       double sheet_strength=1.0,
                       double thresh=0.1, double band=0.25,
                       double margin=0.05, double strength=1.0,
                       bint verbose=False):
    """
    Push buried sulcal valleys in a PPM below the isovalue.

    Mirrors the correction ``CAT_VolMarchingCubes -strength-sulci`` applies
    internally; defaults match ``CAT_PpmSulciOptionsInit``.

    A buried sulcus is a valley in the percentage position map whose floor
    never drops below the isovalue, so an isosurface at that value bridges
    the two banks.  No intensity image is needed, which matters because by
    the time the surface is extracted the T1 is gone: the PPM carries the
    geometry itself.  Crossing a sulcus it runs 1 (WM) -> 0.5 -> ~0 (pial)
    -> 0.5 -> 1, and crossing a gyral blade it runs 0 -> 0.5 -> ~1 -> 0.5
    -> 0, so a sulcus is a *valley* and a blade a *ridge*, and the polarity
    guard of the sheetness filter separates them exactly.

    Three conditions must hold together, and it is the conjunction that
    makes this safe rather than just another filter: the value lies in
    ``(isovalue, isovalue + band)`` so solid white matter is never touched;
    the dark-sheet response exceeds ``thresh``; and the value is a genuine
    local minimum *along the sheet normal*, which a gyral crown can never
    satisfy.  The correction only ever lowers a value.

    Parameters
    ----------
    ppm : array_like, 3-D, float32
        Percentage position map in [0, 1].
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    isovalue : float
        Threshold the surface will be extracted at (default 0.5).
    sheetness, normal : array_like, optional
        Precomputed dark-sheet field; estimated from ``ppm`` when omitted.
    sigma_min, sigma_max, n_scales : float, float, int
        Scale range of the sheetness filter; see :func:`vol_sheetness`.
    sheet_strength : float
        Gain on the sheetness response (default 1.0).  The filter's
        automatic noise scale is half the largest Hessian norm in the
        volume, so where a strong structure dominates it the raw response
        can fall below ``thresh`` and nothing happens.  Raise this when
        that is the case; ``verbose`` reports where the response landed.
    thresh : float
        Sheetness below this is ignored (default 0.1).
    band : float
        Only values in ``(isovalue, isovalue + band)`` are touched
        (default 0.25).
    margin : float
        How far below the isovalue an opened valley is pushed
        (default 0.05).
    strength : float
        Blend towards that target, in [0, 1] (default 1.0).  0 is an
        exact no-op.
    verbose : bool

    Returns
    -------
    out : ndarray, 3-D, float32
        Corrected PPM.
    """
    vol = np.asfortranarray(ppm, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("ppm must be 3-D")

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

    cdef cnp.ndarray[cnp.float32_t, ndim=3] out = vol.copy(order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] sheet_arr
    cdef cnp.ndarray[cnp.float32_t, ndim=1] normal_flat
    cdef const float *sheet_ptr = NULL
    cdef const float *normal_ptr = NULL

    if sheetness is not None and normal is not None:
        sheet_arr = np.asfortranarray(sheetness, dtype=np.float32)
        if (sheet_arr.shape[0] != dims[0] or sheet_arr.shape[1] != dims[1]
                or sheet_arr.shape[2] != dims[2]):
            raise ValueError("sheetness shape must match ppm shape")
        normal_flat = _pack_normal(normal, dims)
        sheet_ptr = <const float *>sheet_arr.data
        normal_ptr = <const float *>normal_flat.data

    cdef C.CAT_PpmSulciOpts opts
    C.CAT_PpmSulciOptionsInit(&opts)
    opts.sigma_min = sigma_min
    opts.sigma_max = sigma_max
    opts.n_scales = n_scales
    opts.sheet_normalize = sheet_normalize
    opts.sheet_skeleton = 1 if sheet_skeleton else 0
    opts.sheet_strength = sheet_strength
    opts.thresh = thresh
    opts.band = band
    opts.margin = margin
    opts.strength = strength
    opts.verbose = 1 if verbose else 0

    cdef int rc = C.CAT_VolOpenPpmSulci(<float *>out.data, sheet_ptr, normal_ptr,
                                        dims, vx, isovalue, &opts)
    if rc != 0:
        raise RuntimeError(f"CAT_VolOpenPpmSulci returned error code {rc}")
    return out


# ===================================================================
# Boundary offset  (mirrors CAT_VolBoundaryOffset)
# ===================================================================
def vol_boundary_offset(label, t1w, voxelsize=None,
                        double ref_fwhm=10.0,
                        double erosion_mm=3.0,
                        double t_lo=0.25,
                        double t_hi=0.75,
                        double search_mm=4.0,
                        double step_mm=0.25,
                        double width_pct=88.0,
                        double gain=0.5,
                        double max_offset_mm=1.5,
                        double smooth_fwhm=8.0,
                        bint return_width=False,
                        bint verbose=False):
    """
    Myelination-induced displacement of the GM/WM boundary, in mm.

    Mirrors ``CAT_VolBoundaryOffset``.  Defaults match
    ``CAT_BoundaryOffsetOptionsInit`` from the library -- the single
    source of truth.

    In the primary motor and somatosensory strip, and along the line of
    Gennari in V1, the deep cortical layers carry enough myelin to
    approach white-matter intensity on T1w.  The classifier follows the
    intensity, so the GM/WM boundary is placed too far out and the cortex
    comes back too thin.

    The displacement is a fraction of a voxel, which is why relabelling
    cannot express it -- a label can only move in whole isovalue
    crossings.  What *can* carry it is PBT's distance field, which already
    holds the boundary to sub-voxel precision, so this returns the
    displacement in millimetres to be added to ``dist_WM``.

    **The observable is the width of the intensity transition, not its
    position.**  The label map is derived from the intensity, so the label
    boundary and the intensity boundary agree by construction and their
    difference measures nothing.  What distinguishes myelinated cortex is
    the shape of the profile across the boundary: healthy cortex steps
    from grey to white over about the partial-volume width, while a
    myelinated deep layer turns that step into a ramp one to three
    millimetres wide.  The excess over this brain's own healthy majority
    is the signal.

    Parameters
    ----------
    label : array_like, 3-D, float32
        PVE label volume, values in [0, 3] (1=CSF, 2=GM, 3=WM).
    t1w : array_like, 3-D, float32
        T1w intensity volume on the same grid.  Intensity units are
        irrelevant: the profile is read in a contrast-normalized
        coordinate built from local tissue references.
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    ref_fwhm : float
        FWHM in mm of the local GM and WM intensity references
        (default 10.0).  Because both ends of the normalized coordinate
        are local, a multiplicative bias cancels out of it.
    erosion_mm : float
        Erosion in mm of the WM intensity control set (default 3.0), so
        the myelinated band does not calibrate the reference it is being
        measured against.
    t_lo, t_hi : float
        Level crossings of the normalized profile bounding the transition
        (defaults 0.25 and 0.75).  Their separation is the width.
    search_mm : float
        Half-length in mm of the profile searched along the normal
        (default 4.0).
    step_mm : float
        Sampling step in mm along the profile (default 0.25).
    width_pct : float
        Upper percentile of the measured widths above which a transition
        is read as myelinated (default 88.0) -- the share of the boundary
        the correction may touch.

        A *central* location is no use here: whatever one is chosen the
        healthy population straddles it, every voxel above acquires a
        positive offset, and the whole cortex ends up displaced.  At p25
        that was 99.8% of the boundary on a real 0.75 mm subject, at
        0.05-0.1 mm, because p25 sits below the mode of a distribution
        whose healthy peak is one voxel wide.  Nor does a location plus a
        multiple of the spread transfer between subjects, because the
        spread follows resolution and noise rather than anatomy: one
        setting selected 12.0% of the boundary on a 0.75 mm scan and 1.4%
        on a 1x1x1.25 mm one.

        Fixing the share does not fix the correction.  The displacement is
        the excess *beyond* the threshold, so a brain whose widths are
        tightly grouped is barely corrected however large the share.
        Measured at p88 on two subjects the mean displacement came out at
        0.19 and 0.13 mm, and it moves by under 10% across p85-p92.
    gain : float
        Displacement per unit excess width (default 0.5).  For a symmetric
        ramp the cytoarchitectonic boundary sits about half the excess
        beyond the intensity midpoint.
    max_offset_mm : float
        Clamp on the displacement (default 1.5).
    smooth_fwhm : float
        FWHM in mm of smoothing applied to the width *within* the boundary
        sheet (default 8.0), before any threshold.  Myelination is
        centimetre-scale and stereotyped while the per-voxel profile is
        noisy, and the band is a thin surface, so a normalized convolution
        restricted to it averages along the boundary and not across it.
        Doing this before the threshold is what makes the threshold
        meaningful: on a 0.75 mm subject it takes the healthy spread from
        0.121 to 0.080 mm and the p99 of the width from 4.06 to 2.01 mm --
        the extreme per-voxel values were noise.
    return_width : bool
        Also return the raw transition width (default False).  This is the
        map to inspect first -- the displacement is only a rescaling of it.
    verbose : bool

    Returns
    -------
    offset : ndarray, 3-D, float32
        Displacement in mm, 0 outside the transition band.
    width : ndarray, 3-D, float32
        Transition width in mm, smoothed within the sheet and before the
        threshold -- the quantity the decision is made on, and the map to
        inspect first.  Only when ``return_width`` is True.
    """
    lab = np.asfortranarray(label, dtype=np.float32)
    t1 = np.asfortranarray(t1w, dtype=np.float32)
    if lab.ndim != 3:
        raise ValueError("label must be 3-D")
    if t1.shape != lab.shape:
        raise ValueError("label and t1w must have the same shape")

    cdef int dims[3]
    dims[0] = lab.shape[0]
    dims[1] = lab.shape[1]
    dims[2] = lab.shape[2]

    cdef double vx[3]
    if voxelsize is not None:
        vs = np.asarray(voxelsize, dtype=np.float64).ravel()
        vx[0] = vs[0]; vx[1] = vs[1]; vx[2] = vs[2]
    else:
        vx[0] = 1.0; vx[1] = 1.0; vx[2] = 1.0

    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = lab
    cdef cnp.ndarray[cnp.float32_t, ndim=3] img = t1
    cdef cnp.ndarray[cnp.float32_t, ndim=3] off = np.zeros_like(lab, order='F')
    cdef cnp.ndarray[cnp.float32_t, ndim=3] wid
    cdef float *wid_ptr = NULL

    if return_width:
        wid = np.zeros_like(lab, order='F')
        wid_ptr = <float *>wid.data

    cdef C.CAT_BoundaryOffsetOpts opts
    C.CAT_BoundaryOffsetOptionsInit(&opts)
    opts.ref_fwhm = ref_fwhm
    opts.erosion_mm = erosion_mm
    opts.t_lo = t_lo
    opts.t_hi = t_hi
    opts.search_mm = search_mm
    opts.step_mm = step_mm
    opts.width_pct = width_pct
    opts.gain = gain
    opts.max_offset_mm = max_offset_mm
    opts.smooth_fwhm = smooth_fwhm
    opts.verbose = 1 if verbose else 0

    cdef long rc = C.CAT_VolBoundaryOffset(<const float *>src.data,
                                           <const float *>img.data,
                                           <float *>off.data, wid_ptr,
                                           dims, vx, &opts)
    if rc < 0:
        raise RuntimeError(f"CAT_VolBoundaryOffset returned error code {rc}")

    if return_width:
        return off, wid
    return off


def vol_thickness_pbt(volume, voxelsize=None,
                      int n_avgs=-1, int n_median_filter=-1,
                      int median_subsample=-1, double range_val=-1.0,
                      double fill_thresh=-1.0,
                      correct_thickness=None,
                      double sulcal_width=-1.0,
                      bint pve_distance=False,
                      bint sulcal_barrier=False,
                      double barrier_q=-1.0, double barrier_tmin=-1.0,
                      double barrier_gmtfactor=-1.0, double barrier_gmtpct=-1.0, double barrier_ramp=-1.0, double barrier_local=-1.0,
                      double barrier_gmtmax=-1.0,
                      double barrier_dmin=-1.0,
                      double barrier_halfwidth=-1.0,
                      bint oriented_filter=False,
                      double oriented_strength=-1.0, double oriented_cutoff=-1.0,
                      bint fast=False, bint verbose=False):
    """
    Compute projection-based cortical thickness (PBT).

    Mirrors ``CAT_VolThicknessPbt``.  Every argument defaults to a
    sentinel meaning "not asked for", so the value actually used comes
    from ``CAT_PbtOptionsInit`` in the library -- the single source of
    truth.  ``correct_thickness`` takes ``None`` rather than a negative
    sentinel, because a negative correction is a legitimate value.

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
    sulcal_barrier : bool
        Stop the CSF distance at the sulcal medial surface (default
        False).  Where the classifier lost the CSF in a sulcus there is
        no boundary for the CSF distance to stop at, so the front from
        one bank runs through the fused grey matter into the other; the
        thickness follows it, the PPM never drops below the isovalue, and
        the buried sulcus is created there — in the distance map, long
        before marching cubes is asked to draw it.  The midline the front
        should have stopped at is recovered from geometry rather than
        intensity: it is where the fronts from the two banks collide, so
        no intensity image and no per-subject threshold are involved.
        Applied as a minimum, which makes it self-limiting — where CSF
        *was* segmented it is nearer than the midline and the result is
        identical to leaving this off, and the distance can only shrink,
        so an overestimated thickness can be fixed but a correct one can
        never be inflated.
    barrier_q : float
        Shock threshold of the medial set (default 0.0, which selects
        ``CAT_MEDIAL_SET_Q``).  Lower is stricter and gives a thinner
        midline; the result is flat over a wide band and then falls off,
        so the default sits in the middle of the flat region.
    barrier_gmtfactor : float
        Multiple of the median thickness at which the gate sits (default
        2.0; <= 0 disables the gate).  The criterion in its natural form —
        a glued sulcus is *two* cortices back to back, so the threshold
        belongs at twice the typical thickness of the brain being
        processed rather than at a fixed millimetre value that is only
        right for the cortex it was tuned on.  The median is taken over
        ``dist_WM + dist_CSF`` in the GM band, which for a band of locally
        constant thickness is that thickness exactly, and a median is
        unmoved by the glued minority the gate exists to catch.  It is the
        *pre-projection* proxy though, running ~15% high against the
        reported GMT, so 2.0 here is nearer 2.3x the final thickness.
    barrier_gmtmax : float
        Only bound where the implied thickness is impossible for cortex,
        in mm (default 5.0; 0 disables).  This is the gate that matters.
        A glued sulcus does not merely look thick, it looks like *two*
        cortices back to back — 5-6 mm where 2-3 mm is normal — so the
        implied thickness at the voxel, ``dist_WM + dist_CSF``, separates
        the two populations cleanly.  Gating on the CSF distance alone
        does not: plenty of ordinary voxels sit far from CSF simply by
        being near the white matter.  Measured on an ADNI subject, the
        gate takes the capped fraction of the cortex from 20.5% to 2.4%
        and the mean thickness from 1.81 mm back to 2.28 mm against an
        unbarriered 2.31 mm.
    barrier_dmin : float
        Only bound a CSF distance already implausible for real cortex, in
        mm (default 2.0; 0 disables).  Bounding always shrinks the
        thickness, so without this the mean thickness of the brain becomes
        a readout of how often the barrier fired and moves with
        ``barrier_q``.
    barrier_tmin : float
        Ignore collisions closer than this to the WM boundary, in mm
        (default 0.5).  Stops the barrier carving into thin cortex from
        the inside.
    barrier_halfwidth : float
        Half the width the barrier stands for, in mm (default 0.0, which
        selects one voxel).  The distance transform measures to the
        nearest medial voxel centre while the grey matter ends at the
        surface of that set.  It does not move the 0.5 crossing — it
        trades against thickness accuracy only.
    oriented_filter : bool
        Replace the isotropic 3x3x3 median filters by sheetness-oriented
        ones (default False).  An isotropic median penalizes boundary
        area, so it removes thin structures whichever side of the label
        boundary they lie on: the same filter that opens a glued sulcus
        closes a cerebellar fissure.  The oriented variant admits only
        neighbours lying in the plane of the locally detected sheet, so
        it averages along a thin structure and never across it.  Where
        no sheet is detected it is identical to the isotropic median.
    oriented_cutoff : float
        Admission cutoff of the oriented medians (default 0, which
        selects 0.10).  A one-voxel-thick sheet is preserved from a
        sheetness of ``2*cutoff`` upwards; see
        :func:`vol_oriented_median`.
    oriented_strength : float
        Overall gain on the sheetness before the filter uses it
        (default 1.0).  0 reproduces the isotropic filters exactly.
        Values above 1 amplify a response too weak to matter: the oriented
        median admits every neighbour while the sheetness stays below
        ``oriented_cutoff``, and only preserves a one-voxel sheet from
        ``2*oriented_cutoff`` upwards.  See ``gain`` in
        :func:`vol_sheetness`.
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
    # CAT_PbtOptionsInit() is the single source of truth for the defaults; a
    # sentinel here means "not asked for" and the library default stands.
    # correct_thickness takes None rather than a negative sentinel because a
    # negative correction is a legitimate value.
    if n_avgs           >= 0:   opts.n_avgs = n_avgs
    if n_median_filter  >= 0:   opts.n_median_filter = n_median_filter
    if median_subsample >= 0:   opts.median_subsample = median_subsample
    if range_val        >= 0.0: opts.range = range_val
    if fill_thresh      >= 0.0: opts.fill_thresh = fill_thresh
    if correct_thickness is not None: opts.correct_thickness = correct_thickness
    if sulcal_width     >= 0.0: opts.sulcal_width = sulcal_width
    opts.pve_distance = 1 if pve_distance else 0
    opts.sulcal_barrier = 1 if sulcal_barrier else 0
    if barrier_q         >= 0.0: opts.barrier_q = barrier_q
    if barrier_gmtfactor >= 0.0: opts.barrier_gmtfactor = barrier_gmtfactor
    if barrier_gmtpct    >= 0.0: opts.barrier_gmtpct = barrier_gmtpct
    if barrier_ramp      >= 0.0: opts.barrier_ramp = barrier_ramp
    if barrier_local     >= 0.0: opts.barrier_local = barrier_local
    if barrier_gmtmax    >= 0.0: opts.barrier_gmtmax = barrier_gmtmax
    if barrier_dmin      >= 0.0: opts.barrier_dmin = barrier_dmin
    if barrier_tmin      >= 0.0: opts.barrier_tmin = barrier_tmin
    if barrier_halfwidth >= 0.0: opts.barrier_halfwidth = barrier_halfwidth
    opts.oriented_filter = 1 if oriented_filter else 0
    if oriented_strength >= 0.0: opts.oriented_strength = oriented_strength
    if oriented_cutoff   >= 0.0: opts.oriented_cutoff = oriented_cutoff
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
                       double strength_sulci=-1.0, double sulci_cutoff=-1.0,
                       double sulci_sheet_strength=-1.0,
                       double sulci_thresh=-1.0, double sulci_band=-1.0,
                       double sulci_normalize=-1.0, int sulci_skeleton=-1,
                       double sulci_sigma_factor=-1.0,
                       double sulci_sigma_min=-1.0, double sulci_sigma_max=-1.0,
                       int sulci_scales=-1, double sheet_offset=-1.0,
                       double sheet_offset_gyri=-1.0,
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
    strength_sulci : float
        Strength of the buried-sulcus correction on the PPM, and its
        on/off gate (default 0 = off).  Raising it enables the
        correction; the other ``sulci_*`` arguments are the values tuned
        on real data and take effect only once this is set.  A buried sulcus is a valley in the PPM
        whose floor never drops below the isovalue, so the two banks fuse
        when the isosurface is extracted.  No intensity image is needed:
        crossing a sulcus the PPM runs 1 → 0.5 → ~0 → 0.5 → 1 and
        crossing a gyral blade it runs 0 → 0.5 → ~1 → 0.5 → 0, so a
        sulcus is a valley and a blade a ridge, and a Hessian sheetness
        filter separates them by polarity alone.  The field is used three
        times: floors just above the isovalue are pushed below it, the
        gyral boost of ``strength_gyri_mask`` is damped there (otherwise
        strengthening a thin WM finger lifts the neighbouring sulcal
        floor back over the isovalue), and the median filter is oriented
        so it cannot re-close what was opened.
    sulci_cutoff : float
        Admission cutoff of that oriented median (default 0.0, which
        selects ``CAT_ORIENTED_MEDIAN_CUTOFF``).  See
        :func:`vol_oriented_median`.
    sulci_sheet_strength : float
        Gain on the sheetness response (default 1.0).  The filter's
        automatic noise scale is half the largest Hessian norm in the
        volume, so where a strong structure dominates it the raw response
        can sit well below ``sulci_thresh`` and the whole correction goes
        inert.  Measure the response first with :func:`vol_sheetness` on
        the PPM and pick a gain that puts its p99 near twice
        ``sulci_thresh``; ``verbose=True`` reports both numbers.
    sulci_skeleton : bool
        Thin the valley field to its medial sheet before any threshold is
        applied (default False).  The plate response is as wide as the
        Gaussian that produced it, so a large ``sulci_sigma_max`` locates
        the valleys well but answers several voxels into the banks on
        either side, and the per-voxel gates then reach into that tissue
        too.  Suppressing everything that is not a maximum along its own
        normal collapses the band onto its ridge line -- one voxel at any
        scale -- and leaves the value on the ridge unchanged.  Runs before
        the percentile anchor, so p99.9 is then taken over ridge values.
    sheet_offset : float
        Signed sheetness offset added to the PPM before extraction, in map
        units (library default 0.6; 0 disables it).  A sulcus is a valley
        and a gyral blade a ridge, so the signed map is negative on one and
        positive on the other and one addition lowers sulci while raising
        blades -- which a global isovalue shift cannot do.
    sheet_offset_gyri : float
        The same offset for the raising half alone.  Negative (default)
        means "use ``sheet_offset``", i.e. the signed map is applied whole.
        Only this half can re-glue banks: raising a ridge protects a thin
        blade but lifts the sulcal floor beside it too, and on a real PPM
        it is the larger of the two -- at ``sheet_offset=0.6`` on an OASIS
        subject it lifted 134214 voxels over the isovalue against 85450
        pushed under.  Lowering it keeps the blade protection without the
        net inflation; it does not reduce the number of sulci opened.
        ``verbose=True`` reports both crossing counts.
    sulci_normalize : float
        Value the p99.9 of the valley response is scaled to (default 1.0);
        0 keeps the raw response.  This is what makes ``sulci_thresh``
        mean the same thing on every dataset.
    sulci_thresh : float
        Sheetness below this is ignored (default 0.1).
    sulci_band : float
        Only values in ``(threshold, threshold + band)`` are opened
        (default 0.25).
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

    # apply_marching_cubes() writes the topology change map back over its input
    # -- that is the optional third output of the CLI, which writes it to a file
    # and exits.  A binding cannot do that: it would hand the caller back a
    # destroyed array, and a second call on the same array would run on the
    # change map instead of the volume and silently return the same surface
    # whatever the options were.  Work on a private copy.
    cdef cnp.ndarray[cnp.float32_t, ndim=1] _mc_scratch = np.empty(
        vh.dims[0] * vh.dims[1] * vh.dims[2], dtype=np.float32)
    memcpy(<void *>_mc_scratch.data, <const void *>vh.data,
           sizeof(float) * <size_t>(vh.dims[0] * vh.dims[1] * vh.dims[2]))

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
    cdef C.CAT_PpmSulciOpts sulci_opts
    cdef C.CAT_PpmSulciOpts *sulci_ptr = NULL
    if strength_sulci != 0.0:
        # CAT_PpmSulciOptionsInit() is the single source of truth for the
        # defaults; a negative value here means "not asked for" and the library
        # default stands.  Duplicating the numbers is how this binding and the
        # CLI drifted apart -- the binding was running the sheetness at a tenth
        # of the library's gain, which left the whole correction inert.
        C.CAT_PpmSulciOptionsInit(&sulci_opts)
        if sulci_sheet_strength >= 0.0: sulci_opts.sheet_strength = sulci_sheet_strength
        if sulci_normalize      >= 0.0: sulci_opts.sheet_normalize = sulci_normalize
        if sulci_skeleton       >= 0:   sulci_opts.sheet_skeleton = sulci_skeleton
        if sulci_thresh         >= 0.0: sulci_opts.thresh = sulci_thresh
        if sulci_band           >= 0.0: sulci_opts.band = sulci_band
        if sulci_sigma_factor   >= 0.0: sulci_opts.sigma_factor = sulci_sigma_factor
        if sulci_sigma_min      >= 0.0: sulci_opts.sigma_min = sulci_sigma_min
        if sulci_sigma_max      >= 0.0:
            sulci_opts.sigma_max = sulci_sigma_max
            # An explicit sigma_max means an explicit sigma_max: the derivation
            # from the PPM's own thickness runs whenever sigma_factor is positive
            # and would overwrite it.  Passing both keeps the factor.
            if sulci_sigma_factor < 0.0:
                sulci_opts.sigma_factor = 0.0
        if sulci_scales         >= 1:   sulci_opts.n_scales = sulci_scales
        if sheet_offset    >= 0.0: sulci_opts.offset = sheet_offset
        if sheet_offset_gyri >= 0.0: sulci_opts.offset_gyri = sheet_offset_gyri
        if strength_sulci  >= 0.0: sulci_opts.strength = strength_sulci
        if sulci_cutoff    >= 0.0: sulci_opts.cutoff = sulci_cutoff
        sulci_opts.verbose = 1 if verbose else 0
        sulci_ptr = &sulci_opts

    try:
        if fast:
            result = C.apply_marching_cubes_fast(
                <float *>_mc_scratch.data, vh.nii, threshold,
                iter_laplacian, 1 if verbose else 0)
        else:
            result = C.apply_marching_cubes(
                <float *>_mc_scratch.data, vh.nii, label_data,
                threshold, pre_fwhm, iter_laplacian,
                dist_morph_val, n_median_filter, n_iter,
                strength_gyri_mask, sulci_ptr, 1 if verbose else 0)
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


