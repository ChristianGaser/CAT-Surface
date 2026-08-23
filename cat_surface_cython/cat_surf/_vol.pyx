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
# Shape triage of implausibly thick cortex  (mirrors CAT_VolThicknessQC)
# ===================================================================
def vol_thickness_qc(gmt, label=None, voxelsize=None,
                     double thresh=4.5, double plate_radius=2.5,
                     double min_volume=20.0, int conn=26,
                     bint return_classmap=False, bint verbose=False):
    """
    Triage implausibly thick cortex by shape.

    Mirrors ``CAT_VolThicknessQC``.  A thickness map normally carries a
    population of voxels above any defensible value -- cortex does not
    exceed roughly 4.5 mm -- and two different faults produce them:

    ``plate``
        two sulcal banks with no CSF between them, measured as one band.
        There is a sulcus to recover; :func:`vol_sulcus_repair` is the
        tool.
    ``solid``
        cortex merged with subcortical grey matter, or a genuinely thick
        region such as a temporal or orbitofrontal pole.  There is no
        sulcus, and carving one invents anatomy.

    Thickness alone cannot separate them; shape can.  A glued sulcus is a
    band at most about 5 mm across however far it runs along the sulcus,
    so the largest sphere fitting inside it has a radius of about 2.5 mm,
    while a solid mass has no such bound.  The discriminator is therefore
    the maximum inscribed radius of each connected component, which
    depends neither on where it sits nor on how large it is.

    Parameters
    ----------
    gmt : array_like, 3-D, float32
        Cortical thickness map in mm, e.g. from
        :func:`vol_thickness_pbt`.
    label : array_like, 3-D, float32, optional
        PVE label map in [0..3] confining the search to the cortical
        band.  Recommended; without it the thickness map is used as it
        stands.
    voxelsize : array_like, shape (3,), float64, optional
        Voxel dimensions in mm.  Default ``[1, 1, 1]``.
    thresh : float
        Thickness above which a voxel is implausible, in mm (default
        4.5).
    plate_radius : float
        Largest inscribed radius still reported as a plate, in mm
        (default 2.5).
    min_volume : float
        Discard components below this volume, in mm^3 (default 20).
    conn : int
        Voxel connectivity: 6, 18 or 26 (default 26).
    return_classmap : bool
        Also return a map tagging each flagged voxel 1 (plate) or 2
        (solid), 0 elsewhere.
    verbose : bool

    Returns
    -------
    components : list of dict
        One entry per component, with keys ``n_voxels``, ``volume_mm3``,
        ``max_radius``, ``centroid``, ``gmt_mean``, ``gmt_max`` and
        ``shape`` (``"plate"`` or ``"solid"``), ordered largest first.
    classmap : ndarray, 3-D, float32
        Only when ``return_classmap`` is True.
    """
    vol = np.asfortranarray(gmt, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError("gmt must be 3-D")

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

    cdef cnp.ndarray[cnp.float32_t, ndim=3] src = vol
    cdef cnp.ndarray[cnp.float32_t, ndim=3] lab
    cdef const float *lab_ptr = NULL
    if label is not None:
        lab = np.asfortranarray(label, dtype=np.float32)
        if (lab.shape[0] != dims[0] or lab.shape[1] != dims[1]
                or lab.shape[2] != dims[2]):
            raise ValueError("label shape must match gmt shape")
        lab_ptr = <const float *>lab.data

    cdef cnp.ndarray[cnp.float32_t, ndim=3] cls
    cdef float *cls_ptr = NULL
    if return_classmap:
        cls = np.zeros_like(vol, order='F')
        cls_ptr = <float *>cls.data

    cdef C.CAT_ThicknessQCOpts opts
    C.CAT_ThicknessQCOptionsInit(&opts)
    opts.thresh = thresh
    opts.plate_radius = plate_radius
    opts.min_volume = min_volume
    opts.conn = conn
    opts.verbose = 1 if verbose else 0

    cdef C.CAT_ThicknessComponent *comps = NULL
    cdef int n_comps = 0
    cdef int rc = C.CAT_VolThicknessQC(<const float *>src.data, lab_ptr, dims,
                                       vx, &opts, &comps, &n_comps, cls_ptr)
    if rc != 0:
        raise RuntimeError(f"CAT_VolThicknessQC returned error code {rc}")

    out = []
    cdef int i
    try:
        for i in range(n_comps):
            out.append(dict(
                n_voxels=int(comps[i].n_voxels),
                volume_mm3=float(comps[i].volume_mm3),
                max_radius=float(comps[i].max_radius),
                centroid=(float(comps[i].centroid[0]),
                          float(comps[i].centroid[1]),
                          float(comps[i].centroid[2])),
                gmt_mean=float(comps[i].gmt_mean),
                gmt_max=float(comps[i].gmt_max),
                shape="plate" if comps[i].shape == C.CAT_QC_PLATE else "solid",
            ))
    finally:
        free(comps)

    out.sort(key=lambda c: c["volume_mm3"], reverse=True)

    if return_classmap:
        return out, cls
    return out


# ===================================================================
# Pre-PBT repair of a PVE label map  (mirrors CAT_VolSulcusRepair)
# ===================================================================
def vol_sulcus_repair(t1, label, voxelsize=None,
                      bint recover_csf=True, strengthen_wm=None,
                      bint refine_pve=False,
                      double sheet_sigma_min=0.3, double sheet_sigma_max=1.0,
                      int sheet_n_scales=3, double sheet_strength=1.0,
                      double csf_min_dist=1.5, double csf_min_wmdist=0.75,
                      double csf_thresh=0.1, double csf_strength=0.8,
                      double wm_thresh=0.1, double wm_strength=0.8,
                      double wm_min_int=2.1, int wm_max_gap=3,
                      double wm_sulcus_guard=1.0,
                      double sheet_normalize=1.0,
                      bint sheet_skeleton=False,
                      double band_min_dist=1.5, int band_window=4,
                      double band_strength=0.7,
                      bint return_sheetness=False, bint verbose=False,
                      reconnect_gyri=None):
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
    strengthen_wm : bool
        Run step 2, the WM blade strengthening (default True).  It
        rescues the fine white-matter fingers reaching into the gyral
        crowns, which partial volume drags towards GM so that the
        classifier drops their last millimetre.
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
    sheet_strength : float
        Overall gain on the sheetness before the thresholds see it
        (default 1.0).  Every step is gated on the response clearing
        ``csf_thresh`` or ``wm_thresh`` and then ramps its blend weight
        from zero at that threshold, so a uniformly weak response makes
        the whole call a no-op.  See ``gain`` in :func:`vol_sheetness`.
    csf_thresh, csf_strength : float
        Sheetness threshold and blend weight for step 1 (defaults 0.1
        and 0.8).  The weight ramps from 0 at ``csf_thresh`` to
        ``csf_strength`` at ``2*csf_thresh``, so the threshold is half
        the sheetness at which the correction acts at full strength --
        the same relation :func:`vol_oriented_median` has to its cutoff.
        Match it to the response the data produces.
    wm_thresh, wm_strength : float
        Sheetness threshold and blend weight for step 2 (defaults 0.1
        and 0.8), ramping exactly as the CSF pair does.
    wm_min_int : float
        Intensity floor for strengthening a blade, on the 1..3 label
        axis where GM = 2 and WM = 3 (default 2.1).  A blade tip is
        dragged towards GM by partial volume, so this sits just above
        pure GM rather than half way to WM.
    reconnect_gyri : bool, optional
        Deprecated alias of ``strengthen_wm``.
    wm_sulcus_guard : float
        How strongly a neighbouring sulcus vetoes the blade strengthening,
        in [0, 1] (default 1.0; 0 disables the guard).  A blade tip and
        the sulcal floor behind it are one voxel apart where the cortex
        is thin, so raising the tip towards WM closes the sulcus — the
        failure is common in the occipital lobe, where the banks are
        already almost touching.  The polarity guard alone does not catch
        it: the bright ridge really is there, it is the *neighbouring*
        dark sheet that must not be filled.  A second dark-sheet pass is
        run, kept only where the intensity confirms a genuine sulcal
        floor, dilated by one voxel, and used to damp the blend weight.
    wm_max_gap : int
        How far from existing WM, in voxels, a blade may still be
        strengthened (default 3).
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
    # step 2 was called "reconnect_gyri" while it bridged two-sided gaps; it now
    # strengthens blade tips, which is where blades actually break.  The old
    # keyword still works so that existing callers keep running.
    if reconnect_gyri is not None:
        if strengthen_wm is not None and bool(strengthen_wm) != bool(reconnect_gyri):
            raise ValueError(
                "strengthen_wm and its deprecated alias reconnect_gyri disagree")
        strengthen_wm = reconnect_gyri
    cdef bint do_strengthen_wm = True if strengthen_wm is None else bool(strengthen_wm)

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
    opts.sheet_strength = sheet_strength
    opts.csf_min_dist = csf_min_dist
    opts.csf_min_wmdist = csf_min_wmdist
    opts.csf_thresh = csf_thresh
    opts.csf_strength = csf_strength
    opts.wm_thresh = wm_thresh
    opts.wm_strength = wm_strength
    opts.wm_min_int = wm_min_int
    opts.wm_max_gap = wm_max_gap
    opts.wm_sulcus_guard = wm_sulcus_guard
    opts.sheet_normalize = sheet_normalize
    opts.sheet_skeleton = 1 if sheet_skeleton else 0
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

    if do_strengthen_wm:
        rc = C.CAT_VolStrengthenWmBlades(<const float *>src.data,
                                         <float *>lab.data, sheet_ptr,
                                         dims, vx, &opts)
        if rc != 0:
            raise RuntimeError(
                f"CAT_VolStrengthenWmBlades returned error code {rc}")

    if refine_pve:
        rc = C.CAT_VolRefinePveNarrowBand(<const float *>src.data,
                                          <float *>lab.data, dims, vx, &opts)
        if rc != 0:
            raise RuntimeError(f"CAT_VolRefinePveNarrowBand returned error code {rc}")

    if return_sheetness:
        return lab, sheet
    return lab


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


def vol_thickness_pbt(volume, voxelsize=None,
                      int n_avgs=2, int n_median_filter=2,
                      int median_subsample=4, double range_val=0.45,
                      double fill_thresh=0.5,
                      double correct_thickness=0.25,
                      double sulcal_width=2.5,
                      bint pve_distance=False,
                      bint sulcal_barrier=False,
                      double barrier_q=0.0, double barrier_tmin=0.5,
                      double barrier_gmtfactor=2.0, double barrier_gmtmax=0.0,
                      double barrier_dmin=2.0,
                      double barrier_halfwidth=0.0,
                      bint oriented_filter=False,
                      double oriented_strength=1.0, double oriented_cutoff=0.0,
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
    opts.n_avgs = n_avgs
    opts.n_median_filter = n_median_filter
    opts.median_subsample = median_subsample
    opts.range = range_val
    opts.fill_thresh = fill_thresh
    opts.correct_thickness = correct_thickness
    opts.sulcal_width = sulcal_width
    opts.pve_distance = 1 if pve_distance else 0
    opts.sulcal_barrier = 1 if sulcal_barrier else 0
    opts.barrier_q = barrier_q
    opts.barrier_gmtfactor = barrier_gmtfactor
    opts.barrier_gmtmax = barrier_gmtmax
    opts.barrier_dmin = barrier_dmin
    opts.barrier_tmin = barrier_tmin
    opts.barrier_halfwidth = barrier_halfwidth
    opts.oriented_filter = 1 if oriented_filter else 0
    opts.oriented_strength = oriented_strength
    opts.oriented_cutoff = oriented_cutoff
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
                       double strength_sulci=0.0, double sulci_cutoff=0.0,
                       double sulci_sheet_strength=1.0,
                       double sulci_thresh=0.3, double sulci_band=0.25,
                       double sulci_normalize=1.0,
                       bint sulci_skeleton=False, double sheet_offset=0.0,
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
        Strength of the buried-sulcus correction on the PPM, in [0, 1]
        (default 0.0 = off).  A buried sulcus is a valley in the PPM
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
        inert — the same reason :func:`vol_sulcus_repair` has its own
        ``sheet_strength``.  The two are *not* interchangeable: that one
        is measured on the intensity image, this one on the PPM, which is
        the better-conditioned of the two.  Pass ``verbose=True`` to see
        the p99 and maximum of the response next to the threshold.
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
    if strength_sulci > 0.0:
        C.CAT_PpmSulciOptionsInit(&sulci_opts)
        sulci_opts.sheet_normalize = sulci_normalize
        sulci_opts.sheet_skeleton = 1 if sulci_skeleton else 0
        sulci_opts.offset = sheet_offset
        sulci_opts.sheet_strength = sulci_sheet_strength
        sulci_opts.thresh = sulci_thresh
        sulci_opts.band = sulci_band
        sulci_opts.strength = strength_sulci
        sulci_opts.cutoff = sulci_cutoff
        sulci_opts.verbose = 1 if verbose else 0
        sulci_ptr = &sulci_opts

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


