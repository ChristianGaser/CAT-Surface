# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrappers for boundary-based registration (BBR) and volume-based
rigid co-registration — mirrors the CAT_SurfBBReg command-line tool.

All volume inputs are given as NIfTI file paths so the full header
(voxel-to-world matrix) is available for spatial sampling.  Surfaces are
passed as ``(vertices, faces)`` numpy array pairs.

For 4-D EPI series the middle frame is automatically selected, matching
the FreeSurfer ``mri_segreg`` convention.

Smoothing is applied
- *symmetrically* to both images for the volume-registration initialisation
  step (``fwhm``), and
- to the *moving* (EPI) image only before the BBR optimiser starts
  (``fwhm_bbr``) — the BBR cost function then sees the same pre-smoothed
  frame at every iteration, matching ``mri_segreg --fwhm`` behaviour.

The defaults of :func:`bbreg` are those of ``CAT_SurfBBReg``, set in the
library by ``CAT_BBRegOptionsInit()``.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

from cat_surf._nifti_types cimport (
    nifti_image, nifti_image_free, DT_FLOAT32,
)
from cat_surf._volume cimport _build_nii
from cat_surf._bic_types cimport polygons_struct
from cat_surf._convert cimport PolygonsMesh
from cat_surf._convert import arrays_to_polygons

cimport cat_surf._cat_funcs as C

cnp.import_array()


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

cdef object _m16_to_numpy(double m[16]):
    """Return a (4, 4) float64 ndarray from a row-major C double[16]."""
    out = np.empty((4, 4), dtype=np.float64)
    cdef double[:, :] v = out
    cdef int r, c
    for r in range(4):
        for c in range(4):
            v[r, c] = m[r * 4 + c]
    return out


cdef void _smooth_inplace(float *buf, int dims[3],
                           double dx, double dy, double dz,
                           double fwhm_mm) noexcept:
    """Gaussian smooth *buf* (float32, Fortran order) in-place."""
    cdef double vs[3]
    vs[0] = dx if dx > 0.0 else 1.0
    vs[1] = dy if dy > 0.0 else 1.0
    vs[2] = dz if dz > 0.0 else 1.0
    cdef double fv[3]
    fv[0] = fwhm_mm / vs[0]
    fv[1] = fwhm_mm / vs[1]
    fv[2] = fwhm_mm / vs[2]
    C.smooth3(<void *>buf, dims, vs, fv, 0, DT_FLOAT32)


def _bbreg_defaults():
    """Return the defaults of ``CAT_BBRegOptionsInit()`` as a dict.

    The source of truth for the defaults of :func:`bbreg` and
    ``CAT_SurfBBReg``; the smoke test checks the signature against it.
    """
    cdef C.CAT_BBRegOpts opts
    C.CAT_BBRegOptionsInit(&opts)
    return opts


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def bbreg(str volume_file,
          lh_surface=None,
          rh_surface=None,
          str ref_file=None,
          lh_mask=None,
          rh_mask=None,
          lh_thickness=None,
          rh_thickness=None,
          double gm_proj_frac=0.0,
          double wm_dist=0.5,
          double gm_dist=0.5,
          double slope=0.5,
          int invert_contrast=-1,
          double fwhm=5.0,
          double grid_range_mm=4.0,
          double grid_range_rad=0.07,
          int grid_steps=4,
          int max_iter=200,
          double tol=1e-6,
          int vol_reg_levels=4,
          int vol_reg_bins=64,
          int vol_reg_iter=30,
          bint verbose=False,
          int dof=6,
          double fwhm_bbr=2.0):
    """
    Boundary-Based Registration (BBR) of a volume to cortical surfaces.

    Mirrors the ``CAT_SurfBBReg`` command-line tool.  When ``ref_file`` is
    provided an NMI volume-registration step initialises the transform before
    BBR.  When no surfaces are given the function returns the volume-
    registration result directly (standalone mode).  With ``dof`` 9 or 12
    the rigid BBR result is refined with scales (and shears), as
    ``CAT_SurfBBReg -dof``.

    The defaults are those of ``CAT_SurfBBReg`` (``CAT_BBRegOptionsInit()``
    in the library), so the two give the same transform for the same inputs.
    Before cat-surf 1.1.3 the defaults here were others -- a translation
    grid of +-15 mm in 5 mm steps, ``tol`` 1e-5 and no smoothing -- and
    ``fwhm`` smoothed for BBR as well, which ``fwhm_bbr`` now does; in
    1.1.3 ``wm_dist`` was 2.0, ``CAT_SurfBBReg``'s default until then.

    Parameters
    ----------
    volume_file : str
        Moving volume (EPI / BOLD), any NIfTI format.  For 4-D series the
        middle frame is used automatically.
    lh_surface, rh_surface : tuple (vertices, faces) or None
        Left / right hemisphere white-matter surfaces as ``(V×3 float64,
        F×3 int32)`` numpy arrays.
    ref_file : str, optional
        T1w reference volume for NMI initialisation.  Required for standalone
        volume-registration mode (no surfaces).
    lh_mask, rh_mask : array_like, shape (V,), float32, optional
        Per-vertex cortex mask (1 = included, 0 = excluded).
    lh_thickness, rh_thickness : array_like, shape (V,), float32, optional
        Per-vertex cortical thickness (mm).  Used with ``gm_proj_frac > 0``.
    gm_proj_frac : float
        GM sampling offset as fraction of local thickness (0 = use gm_dist).
    wm_dist : float
        WM-side sampling offset in mm (default 0.5).  0.5 rather than
        bbregister's 2 mm: on BOLD-like images with a known transform 2 mm
        doubled the error with T1Prep's surfaces and changed nothing with
        FreeSurfer's.
    gm_dist : float
        GM-side absolute sampling offset in mm (default 0.5).
    slope : float
        BBR cost saturation slope (default 0.5).
    invert_contrast : int
        ``-1`` auto-detect (default), ``0`` force T1/FLAIR, ``1`` force
        T2/BOLD.
    fwhm : float
        Gaussian pre-smoothing FWHM in mm of both volumes for the NMI
        initialisation (default 5.0, 0 = off).
    grid_range_mm : float
        Stage-1 brute-force translation grid half-width (mm, default 4.0).
    grid_range_rad : float
        Stage-1 rotation grid half-width (radians, default 0.07); unused,
        the grid searches translations only.
    grid_steps : int
        Stage-1 steps on either side of the start (default 4).
    max_iter : int
        Powell maximum iterations for BBR (default 200).
    tol : float
        Powell convergence tolerance (default 1e-6).
    vol_reg_levels : int
        Gaussian pyramid levels for NMI initialisation (default 4).
    vol_reg_bins : int
        NMI histogram bins per axis (default 64).
    vol_reg_iter : int
        Powell iterations per pyramid level (default 30).
    verbose : bool
    dof : int
        Degrees of freedom of the BBR: 6 (rigid, default), 9 (+ scales) or
        12 (+ shears).  9 and 12 refine the rigid result and need surfaces;
        the NMI initialisation is always rigid.  The BBR contrast peaks where
        the boundary lies midway between the WM and GM samples, which scales
        can follow: with ``wm_dist`` larger than the GM offset they come out
        larger by about (wm_dist - GM offset) / 2 mm at the cortex.
    fwhm_bbr : float
        Gaussian pre-smoothing FWHM in mm of the moving volume before BBR
        (default 2.0, 0 = off).

    Returns
    -------
    matrix : ndarray, shape (4, 4), float64
        Transform in the standard output convention: moving (EPI) RAS
        → fixed (T1) RAS, rigid unless ``dof`` is 9 or 12.  Equivalent to
        the matrix written by ``CAT_SurfBBReg`` / ``bbregister``.
    cost : float
        Final BBR cost (lower = better).  Returns 0.0 when only volume
        registration is performed (no surfaces).
    """
    # --- C-level variables (all declared at function scope) ----------------
    cdef float        *vol_flat       = NULL
    cdef nifti_image  *nii            = NULL
    cdef float        *ref_flat       = NULL
    cdef nifti_image  *ref_nii        = NULL
    cdef float        *vol_smooth_buf = NULL   # smoothed moving for vol-reg
    cdef float        *ref_smooth_buf = NULL   # smoothed fixed  for vol-reg
    cdef float        *vol_bbr_buf    = NULL   # smoothed moving for BBR
    cdef float        *vol3d                   # pointer into vol_flat
    cdef float        *vol_for_reg             # vol3d or vol_smooth_buf
    cdef float        *ref_for_reg             # ref_flat or ref_smooth_buf
    cdef float        *vol_bbr                 # vol3d or vol_bbr_buf
    cdef C.CAT_RigidParams p_init, p_best
    cdef C.CAT_AffineParams a_init, a_best
    cdef C.CAT_SurfData    surfs[2]
    cdef int  n_surfs = 0
    cdef int  dims[3], ref_dims[3]
    cdef int  nvox3, ref_nvox3, nt, ic
    cdef double m_best[16], m_inv[16]
    cdef double final_cost = 0.0

    # Python-level variables that keep optional arrays alive
    lh_mask_np = None
    rh_mask_np = None
    lh_thick_np = None
    rh_thick_np = None
    cdef PolygonsMesh lh_mesh = None
    cdef PolygonsMesh rh_mesh = None

    if dof not in (6, 9, 12):
        raise ValueError(f"dof must be 6, 9 or 12, not {dof}")
    if dof != 6 and lh_surface is None and rh_surface is None:
        raise ValueError(
            f"dof {dof} needs surfaces; the volume registration alone is rigid.")

    try:
        # --- 1. Load moving volume -----------------------------------------
        bvol = volume_file.encode("utf-8")
        nii = C.read_nifti_float(bvol, &vol_flat, 0)
        if nii == NULL or vol_flat == NULL:
            raise IOError(f"Failed to read volume '{volume_file}'")

        dims[0] = nii.nx; dims[1] = nii.ny; dims[2] = nii.nz
        nvox3 = dims[0] * dims[1] * dims[2]
        nt = nii.nt if nii.nt > 1 else 1
        vol3d = vol_flat + <size_t>(nt // 2) * nvox3
        if verbose and nt > 1:
            print(f"4D volume: {nt} frames — using mid-frame {nt // 2} "
                  f"for registration")

        # --- 2. Initialise rigid parameters --------------------------------
        p_init.tx = p_init.ty = p_init.tz = 0.0
        p_init.rx = p_init.ry = p_init.rz = 0.0

        # --- 3. Build surface array ----------------------------------------
        if lh_surface is not None:
            lh_verts, lh_faces = lh_surface
            lh_mesh = arrays_to_polygons(lh_verts, lh_faces)
            surfs[n_surfs].surface = lh_mesh.ptr()
            surfs[n_surfs].gm_proj_frac = gm_proj_frac
            if lh_mask is not None:
                lh_mask_np = np.ascontiguousarray(lh_mask, dtype=np.float32)
                surfs[n_surfs].cortex_mask = \
                    <float *>cnp.PyArray_DATA(lh_mask_np)
            else:
                surfs[n_surfs].cortex_mask = NULL
            if lh_thickness is not None:
                lh_thick_np = np.ascontiguousarray(
                    lh_thickness, dtype=np.float32)
                surfs[n_surfs].thickness = \
                    <float *>cnp.PyArray_DATA(lh_thick_np)
            else:
                surfs[n_surfs].thickness = NULL
            n_surfs += 1

        if rh_surface is not None:
            rh_verts, rh_faces = rh_surface
            rh_mesh = arrays_to_polygons(rh_verts, rh_faces)
            surfs[n_surfs].surface = rh_mesh.ptr()
            surfs[n_surfs].gm_proj_frac = gm_proj_frac
            if rh_mask is not None:
                rh_mask_np = np.ascontiguousarray(rh_mask, dtype=np.float32)
                surfs[n_surfs].cortex_mask = \
                    <float *>cnp.PyArray_DATA(rh_mask_np)
            else:
                surfs[n_surfs].cortex_mask = NULL
            if rh_thickness is not None:
                rh_thick_np = np.ascontiguousarray(
                    rh_thickness, dtype=np.float32)
                surfs[n_surfs].thickness = \
                    <float *>cnp.PyArray_DATA(rh_thick_np)
            else:
                surfs[n_surfs].thickness = NULL
            n_surfs += 1

        if n_surfs == 0 and ref_file is None:
            raise ValueError(
                "Provide at least one surface (-lh/-rh) or a reference "
                "volume (ref_file) for standalone volume registration.")

        # --- 4. NMI volume-registration initialisation --------------------
        if ref_file is not None:
            bref = ref_file.encode("utf-8")
            ref_nii = C.read_nifti_float(bref, &ref_flat, 0)
            if ref_nii == NULL or ref_flat == NULL:
                raise IOError(f"Failed to read reference volume '{ref_file}'")
            ref_dims[0] = ref_nii.nx
            ref_dims[1] = ref_nii.ny
            ref_dims[2] = ref_nii.nz
            ref_nvox3   = ref_dims[0] * ref_dims[1] * ref_dims[2]

            vol_for_reg = vol3d
            ref_for_reg = ref_flat

            if fwhm > 0.0:
                # Smoothed copy of mid-frame (moving)
                vol_smooth_buf = <float *>malloc(
                    <size_t>nvox3 * sizeof(float))
                if vol_smooth_buf == NULL:
                    raise MemoryError()
                memcpy(vol_smooth_buf, vol3d,
                       <size_t>nvox3 * sizeof(float))
                _smooth_inplace(vol_smooth_buf, dims,
                                nii.dx, nii.dy, nii.dz, fwhm)
                vol_for_reg = vol_smooth_buf

                # Smoothed copy of reference (fixed)
                ref_smooth_buf = <float *>malloc(
                    <size_t>ref_nvox3 * sizeof(float))
                if ref_smooth_buf == NULL:
                    raise MemoryError()
                memcpy(ref_smooth_buf, ref_flat,
                       <size_t>ref_nvox3 * sizeof(float))
                _smooth_inplace(ref_smooth_buf, ref_dims,
                                ref_nii.dx, ref_nii.dy, ref_nii.dz, fwhm)
                ref_for_reg = ref_smooth_buf

            C.CAT_VolumeReg_register_NMI(
                ref_for_reg, ref_nii, ref_dims,
                vol_for_reg, nii,     dims,
                &p_init,
                vol_reg_levels, vol_reg_bins, vol_reg_iter,
                1 if verbose else 0)

            free(vol_smooth_buf); vol_smooth_buf = NULL
            free(ref_smooth_buf); ref_smooth_buf = NULL

        # --- 5. Pre-smooth for BBR (once, before optimiser starts) ---------
        vol_bbr = vol3d
        if n_surfs > 0 and fwhm_bbr > 0.0:
            vol_bbr_buf = <float *>malloc(<size_t>nvox3 * sizeof(float))
            if vol_bbr_buf == NULL:
                raise MemoryError()
            memcpy(vol_bbr_buf, vol3d, <size_t>nvox3 * sizeof(float))
            _smooth_inplace(vol_bbr_buf, dims,
                            nii.dx, nii.dy, nii.dz, fwhm_bbr)
            vol_bbr = vol_bbr_buf

        # --- 6. Standalone vol-reg exit (no surfaces) ----------------------
        if n_surfs == 0:
            C.CAT_BBReg_params_to_matrix(&p_init, m_best)
            C.CAT_BBReg_invert_matrix(m_best, m_inv)
            return _m16_to_numpy(m_inv), 0.0

        # --- 7. Auto-detect contrast ---------------------------------------
        ic = invert_contrast
        if ic < 0:
            ic = C.CAT_BBReg_detect_contrast(
                &p_init, surfs, n_surfs,
                vol_bbr, nii, dims,
                wm_dist, gm_dist, 1 if verbose else 0)
            if ic < 0:
                if verbose:
                    print("Contrast auto-detection inconclusive "
                          "— defaulting to T1/FLAIR.")
                ic = 0

        # --- 8. BBR optimisation ------------------------------------------
        final_cost = C.CAT_BBReg_optimise(
            &p_init, &p_best, surfs, n_surfs,
            vol_bbr, nii, dims,
            wm_dist, gm_dist, slope, ic,
            grid_range_mm, grid_range_rad, grid_steps,
            max_iter, tol, 1 if verbose else 0)

        # --- 9. Scales (and shears) refine the rigid result ----------------
        if dof > 6:
            C.CAT_BBReg_affine_init(&p_best, surfs, n_surfs, &a_init)
            final_cost = C.CAT_BBReg_optimise_affine(
                &a_init, &a_best, dof, surfs, n_surfs,
                vol_bbr, nii, dims,
                wm_dist, gm_dist, slope, ic,
                max_iter, tol, 1 if verbose else 0)
            C.CAT_BBReg_affine_to_matrix(&a_best, m_best)
            C.CAT_BBReg_invert_affine(m_best, m_inv)
        else:
            C.CAT_BBReg_params_to_matrix(&p_best, m_best)
            C.CAT_BBReg_invert_matrix(m_best, m_inv)
        return _m16_to_numpy(m_inv), final_cost

    finally:
        free(vol_bbr_buf)
        free(vol_smooth_buf)
        free(ref_smooth_buf)
        if ref_nii != NULL:
            free(ref_flat)
            ref_nii.data = NULL
            nifti_image_free(ref_nii)
        if nii != NULL:
            free(vol_flat)
            nii.data = NULL
            nifti_image_free(nii)


def bbreg_detect_contrast(str volume_file,
                           lh_surface=None,
                           rh_surface=None,
                           lh_mask=None,
                           rh_mask=None,
                           double wm_dist=0.5,
                           double gm_dist=0.5,
                           bint verbose=False):
    """
    Detect image contrast type (T1/FLAIR vs T2/BOLD) from WM/GM intensities.

    Samples the volume at WM-side and GM-side positions along the surface
    normals using the identity transform.  Requires at least one surface and
    a roughly pre-aligned volume.

    Parameters
    ----------
    volume_file : str
        NIfTI volume to test (EPI / BOLD).  Mid-frame used for 4-D inputs.
    lh_surface, rh_surface : tuple (vertices, faces) or None
    lh_mask, rh_mask : array_like, shape (V,), float32, optional
    wm_dist : float
        WM-side sampling offset in mm (default 0.5, as for :func:`bbreg`).
    gm_dist : float
        GM-side sampling offset in mm (default 0.5, as for :func:`bbreg`).
    verbose : bool

    Returns
    -------
    int
        ``0`` — T1 / FLAIR (WM brighter than GM).
        ``1`` — T2 / BOLD  (WM darker  than GM).
        ``-1`` — undetermined (fewer than 100 valid samples or
        ``|contrast| < 1 %``).
    """
    cdef float       *vol_flat = NULL
    cdef nifti_image *nii      = NULL
    cdef float       *vol3d
    cdef C.CAT_RigidParams p_zero
    cdef C.CAT_SurfData    surfs[2]
    cdef int  n_surfs = 0
    cdef int  dims[3], nvox3, nt

    lh_mask_np = None
    rh_mask_np = None
    cdef PolygonsMesh lh_mesh = None
    cdef PolygonsMesh rh_mesh = None

    try:
        bvol = volume_file.encode("utf-8")
        nii = C.read_nifti_float(bvol, &vol_flat, 0)
        if nii == NULL or vol_flat == NULL:
            raise IOError(f"Failed to read volume '{volume_file}'")

        dims[0] = nii.nx; dims[1] = nii.ny; dims[2] = nii.nz
        nvox3 = dims[0] * dims[1] * dims[2]
        nt = nii.nt if nii.nt > 1 else 1
        vol3d = vol_flat + <size_t>(nt // 2) * nvox3

        p_zero.tx = p_zero.ty = p_zero.tz = 0.0
        p_zero.rx = p_zero.ry = p_zero.rz = 0.0

        if lh_surface is not None:
            lh_verts, lh_faces = lh_surface
            lh_mesh = arrays_to_polygons(lh_verts, lh_faces)
            surfs[n_surfs].surface     = lh_mesh.ptr()
            surfs[n_surfs].gm_proj_frac = 0.0
            if lh_mask is not None:
                lh_mask_np = np.ascontiguousarray(lh_mask, dtype=np.float32)
                surfs[n_surfs].cortex_mask = \
                    <float *>cnp.PyArray_DATA(lh_mask_np)
            else:
                surfs[n_surfs].cortex_mask = NULL
            surfs[n_surfs].thickness = NULL
            n_surfs += 1

        if rh_surface is not None:
            rh_verts, rh_faces = rh_surface
            rh_mesh = arrays_to_polygons(rh_verts, rh_faces)
            surfs[n_surfs].surface      = rh_mesh.ptr()
            surfs[n_surfs].gm_proj_frac = 0.0
            if rh_mask is not None:
                rh_mask_np = np.ascontiguousarray(rh_mask, dtype=np.float32)
                surfs[n_surfs].cortex_mask = \
                    <float *>cnp.PyArray_DATA(rh_mask_np)
            else:
                surfs[n_surfs].cortex_mask = NULL
            surfs[n_surfs].thickness = NULL
            n_surfs += 1

        if n_surfs == 0:
            raise ValueError(
                "At least one surface (lh_surface or rh_surface) is required.")

        return C.CAT_BBReg_detect_contrast(
            &p_zero, surfs, n_surfs,
            vol3d, nii, dims,
            wm_dist, gm_dist, 1 if verbose else 0)

    finally:
        if nii != NULL:
            free(vol_flat)
            nii.data = NULL
            nifti_image_free(nii)


def volume_register_nmi(str fixed_file,
                         str moving_file,
                         int n_levels=4,
                         int n_bins=64,
                         int max_iter=30,
                         bint verbose=False):
    """
    Cross-modal rigid registration using Normalised Mutual Information (NMI).

    Equivalent to FreeSurfer ``mri_coreg``.  Appropriate for EPI ↔ T1w
    pairs.  For same-modality registration use :func:`volume_register_robust`.

    A Gaussian image pyramid is built for both volumes; NMI is maximised
    at each level (coarse → fine) using Powell's direction-set method.

    Parameters
    ----------
    fixed_file : str
        Reference volume (T1w), any NIfTI format.
    moving_file : str
        Moving volume (EPI / BOLD).  Mid-frame used for 4-D inputs.
    n_levels : int
        Gaussian pyramid levels (default 4).
    n_bins : int
        Histogram bins per axis (default 64).
    max_iter : int
        Powell iterations per pyramid level (default 30).
    verbose : bool

    Returns
    -------
    matrix : ndarray, shape (4, 4), float64
        Rigid transform moving (EPI) RAS → fixed (T1w) RAS.
    nmi : float
        Final NMI value (higher = better alignment).
    """
    cdef float       *fixed_flat  = NULL
    cdef nifti_image *fixed_nii   = NULL
    cdef float       *moving_flat = NULL
    cdef nifti_image *moving_nii  = NULL
    cdef float       *moving3d
    cdef int  fixed_dims[3], moving_dims[3]
    cdef int  nvox3_moving, nt
    cdef C.CAT_RigidParams p
    cdef double m[16], m_inv[16]
    cdef double nmi

    p.tx = p.ty = p.tz = p.rx = p.ry = p.rz = 0.0

    try:
        bfixed = fixed_file.encode("utf-8")
        fixed_nii = C.read_nifti_float(bfixed, &fixed_flat, 0)
        if fixed_nii == NULL or fixed_flat == NULL:
            raise IOError(f"Failed to read fixed volume '{fixed_file}'")
        fixed_dims[0] = fixed_nii.nx
        fixed_dims[1] = fixed_nii.ny
        fixed_dims[2] = fixed_nii.nz

        bmoving = moving_file.encode("utf-8")
        moving_nii = C.read_nifti_float(bmoving, &moving_flat, 0)
        if moving_nii == NULL or moving_flat == NULL:
            raise IOError(f"Failed to read moving volume '{moving_file}'")
        moving_dims[0] = moving_nii.nx
        moving_dims[1] = moving_nii.ny
        moving_dims[2] = moving_nii.nz
        nvox3_moving = moving_dims[0] * moving_dims[1] * moving_dims[2]
        nt = moving_nii.nt if moving_nii.nt > 1 else 1
        moving3d = moving_flat + <size_t>(nt // 2) * nvox3_moving

        nmi = C.CAT_VolumeReg_register_NMI(
            fixed_flat, fixed_nii, fixed_dims,
            moving3d,   moving_nii, moving_dims,
            &p, n_levels, n_bins, max_iter,
            1 if verbose else 0)

        C.CAT_BBReg_params_to_matrix(&p, m)
        C.CAT_BBReg_invert_matrix(m, m_inv)
        return _m16_to_numpy(m_inv), nmi

    finally:
        if moving_nii != NULL:
            free(moving_flat)
            moving_nii.data = NULL
            nifti_image_free(moving_nii)
        if fixed_nii != NULL:
            free(fixed_flat)
            fixed_nii.data = NULL
            nifti_image_free(fixed_nii)


def volume_register_robust(str fixed_file,
                            str moving_file,
                            int n_levels=4,
                            double sat_k=4.685,
                            int max_iter=20,
                            bint iscale=True,
                            bint verbose=False):
    """
    Same-modality rigid registration with Tukey biweight M-estimation.

    Implements the linearised intensity-matching criterion with Tukey
    biweight robust estimation (FreeSurfer ``mri_robust_register --cost
    ROB``).  Appropriate for same-modality pairs (T1w↔T1w, EPI↔EPI).
    For cross-modal pairs use :func:`volume_register_nmi`.

    Parameters
    ----------
    fixed_file : str
        Reference volume, any NIfTI format.
    moving_file : str
        Moving volume.  Mid-frame used for 4-D inputs.
    n_levels : int
        Gaussian pyramid levels (default 4).
    sat_k : float
        Tukey saturation multiplier (default 4.685 ≈ Gaussian 95 %).
    max_iter : int
        IRLS iterations per pyramid level (default 20).
    iscale : bool
        Estimate a global intensity scale between the two volumes along with
        the motion (default True), like ``mri_robust_register --iscale``.
        Without it a reference brighter than the moving volume, such as the
        median of the non-steady-state volumes of a BOLD series, makes most
        of the brain a Tukey outlier.
    verbose : bool

    Returns
    -------
    matrix : ndarray, shape (4, 4), float64
        Rigid transform moving RAS → fixed RAS.
    residual : float
        Final normalised residual (lower = better alignment).
    """
    cdef float       *fixed_flat  = NULL
    cdef nifti_image *fixed_nii   = NULL
    cdef float       *moving_flat = NULL
    cdef nifti_image *moving_nii  = NULL
    cdef float       *moving3d
    cdef int  fixed_dims[3], moving_dims[3]
    cdef int  nvox3_moving, nt
    cdef C.CAT_RigidParams p
    cdef C.CAT_VolumeRegOpts opts
    cdef double m[16], m_inv[16]
    cdef double residual
    cdef double scale = 0.0

    p.tx = p.ty = p.tz = p.rx = p.ry = p.rz = 0.0
    C.CAT_VolumeRegOptionsInit(&opts)
    opts.n_levels = n_levels
    opts.sat_k = sat_k
    opts.max_iter = max_iter
    opts.iscale = 1 if iscale else 0
    opts.verbose = 1 if verbose else 0

    try:
        bfixed = fixed_file.encode("utf-8")
        fixed_nii = C.read_nifti_float(bfixed, &fixed_flat, 0)
        if fixed_nii == NULL or fixed_flat == NULL:
            raise IOError(f"Failed to read fixed volume '{fixed_file}'")
        fixed_dims[0] = fixed_nii.nx
        fixed_dims[1] = fixed_nii.ny
        fixed_dims[2] = fixed_nii.nz

        bmoving = moving_file.encode("utf-8")
        moving_nii = C.read_nifti_float(bmoving, &moving_flat, 0)
        if moving_nii == NULL or moving_flat == NULL:
            raise IOError(f"Failed to read moving volume '{moving_file}'")
        moving_dims[0] = moving_nii.nx
        moving_dims[1] = moving_nii.ny
        moving_dims[2] = moving_nii.nz
        nvox3_moving = moving_dims[0] * moving_dims[1] * moving_dims[2]
        nt = moving_nii.nt if moving_nii.nt > 1 else 1
        moving3d = moving_flat + <size_t>(nt // 2) * nvox3_moving

        residual = C.CAT_VolumeReg_register(
            fixed_flat, fixed_nii, fixed_dims,
            moving3d,   moving_nii, moving_dims,
            &p, &scale, &opts)

        C.CAT_BBReg_params_to_matrix(&p, m)
        C.CAT_BBReg_invert_matrix(m, m_inv)
        return _m16_to_numpy(m_inv), residual

    finally:
        if moving_nii != NULL:
            free(moving_flat)
            moving_nii.data = NULL
            nifti_image_free(moving_nii)
        if fixed_nii != NULL:
            free(fixed_flat)
            fixed_nii.data = NULL
            nifti_image_free(fixed_nii)


# ---------------------------------------------------------------------------
# Head-motion correction (CAT_VolRealign, an mcflirt replacement)
# ---------------------------------------------------------------------------

def vol_realign(series, affine, ref=None, ref_affine=None, int refvol=-1,
                n_levels=None, max_iter=None, sat_k=None, iscale=None,
                bint resample=True, bint verbose=False):
    """
    Head-motion correction of a 4-D series, in mcflirt's conventions.

    Every volume is registered rigidly to the reference with the robust
    registration of :func:`volume_register_robust` (global intensity scale
    included), outward from the middle volume, each starting from its
    neighbour's estimate.  The results are what ``CAT_VolRealign`` -- and
    mcflirt with ``-mats -plots -rmsabs -rmsrel`` -- write.

    Parameters
    ----------
    series : ndarray, shape (nx, ny, nz, nt) or (nx, ny, nz)
        The series; cast to float32.
    affine : ndarray, shape (4, 4)
        Voxel-to-world matrix of the series.
    ref : ndarray, shape (mx, my, mz), optional
        Separate reference volume (mcflirt's ``-reffile``); the realigned
        series is resampled into its grid.  Default: volume ``refvol``.
    ref_affine : ndarray, shape (4, 4), optional
        Voxel-to-world matrix of ``ref``; default ``affine``.
    refvol : int
        Reference volume of the series without ``ref`` (0-based); -1, the
        default, is ``nt // 2`` as in mcflirt.
    n_levels, max_iter, sat_k, iscale : optional
        Options of each registration; ``None`` keeps the default of
        ``CAT_VolRealignOptionsInit()`` (3 levels, 20 iterations, 4.685,
        intensity scale on).
    resample : bool
        Also return the realigned series (default True).
    verbose : bool

    Returns
    -------
    dict
        ``realigned`` (float32, reference grid, with ``resample``),
        ``ras`` and ``fsl`` (nt, 4, 4) volume -> reference transforms in
        world space and in FSL's convention (mcflirt's MAT files), ``par``
        (nt, 6) rx ry rz (rad) tx ty tz (mm), ``rms_abs`` (nt,) and
        ``rms_rel`` (nt-1,) RMS displacements (mm) with ``rms_abs_mean`` and
        ``rms_rel_mean``, ``scale`` (nt,) intensity scales and
        ``ref_index`` (-1 for a separate reference).
    """
    cdef cnp.ndarray s = np.asfortranarray(series, dtype=np.float32)
    if s.ndim == 3:
        s = np.asfortranarray(s[..., None])
    if s.ndim != 4:
        raise ValueError("series must be a 3-D or 4-D array")
    cdef cnp.ndarray aff = np.ascontiguousarray(affine, dtype=np.float64)
    cdef cnp.ndarray vol0 = np.asfortranarray(s[..., 0])
    cdef cnp.ndarray r = None, raff = None
    if ref is not None:
        r = np.asfortranarray(ref, dtype=np.float32)
        if r.ndim == 4:
            r = np.asfortranarray(r[..., 0])
        if r.ndim != 3:
            raise ValueError("ref must be a 3-D array")
        raff = np.ascontiguousarray(aff if ref_affine is None else ref_affine,
                                    dtype=np.float64)

    cdef int nt = s.shape[3]
    cdef int dims[3]
    cdef int rdims[3]
    dims[0] = s.shape[0]; dims[1] = s.shape[1]; dims[2] = s.shape[2]
    cdef nifti_image *snii = NULL
    cdef nifti_image *rnii = NULL
    cdef C.CAT_VolRealignOpts opts
    cdef C.CAT_VolRealignResult res
    cdef cnp.ndarray out = None
    cdef int rc, t

    C.CAT_VolRealignOptionsInit(&opts)
    opts.refvol = refvol
    if n_levels is not None:
        opts.reg.n_levels = n_levels
    if max_iter is not None:
        opts.reg.max_iter = max_iter
    if sat_k is not None:
        opts.reg.sat_k = sat_k
    if iscale is not None:
        opts.reg.iscale = 1 if iscale else 0
    opts.verbose = 1 if verbose else 0
    res.nt = 0
    res.ras = res.fsl = res.par = NULL
    res.rms_abs = res.rms_rel = res.scale = NULL

    try:
        snii = _build_nii(vol0, aff)
        if r is not None:
            rnii = _build_nii(r, raff)
            rdims[0] = r.shape[0]; rdims[1] = r.shape[1]; rdims[2] = r.shape[2]
            rc = C.CAT_VolRealignEstimate(<float *>cnp.PyArray_DATA(s), snii, dims, nt,
                                          <float *>cnp.PyArray_DATA(r), rnii, rdims,
                                          &opts, &res)
        else:
            rdims[0] = dims[0]; rdims[1] = dims[1]; rdims[2] = dims[2]
            rc = C.CAT_VolRealignEstimate(<float *>cnp.PyArray_DATA(s), snii, dims, nt,
                                          NULL, NULL, NULL, &opts, &res)
        if rc != 0:
            raise RuntimeError("CAT_VolRealignEstimate failed")

        result = {
            "ras": np.array(<double[:16 * nt]>res.ras).reshape(nt, 4, 4),
            "fsl": np.array(<double[:16 * nt]>res.fsl).reshape(nt, 4, 4),
            "par": np.array(<double[:6 * nt]>res.par).reshape(nt, 6),
            "rms_abs": np.array(<double[:nt]>res.rms_abs),
            "rms_rel": (np.array(<double[:nt - 1]>res.rms_rel) if nt > 1
                        else np.zeros(0)),
            "rms_abs_mean": res.rms_abs_mean,
            "rms_rel_mean": res.rms_rel_mean,
            "scale": np.array(<double[:nt]>res.scale),
            "ref_index": res.ref_index,
        }
        if resample:
            out = np.empty((rdims[0], rdims[1], rdims[2], nt), dtype=np.float32,
                           order="F")
            rc = C.CAT_VolRealignResample(<float *>cnp.PyArray_DATA(s), snii, dims, nt,
                                          rnii if rnii != NULL else snii, rdims,
                                          &res, <float *>cnp.PyArray_DATA(out))
            if rc != 0:
                raise MemoryError("CAT_VolRealignResample failed")
            result["realigned"] = out
        return result
    finally:
        C.CAT_VolRealignFree(&res)
        if snii != NULL:
            snii.data = NULL
            nifti_image_free(snii)
        if rnii != NULL:
            rnii.data = NULL
            nifti_image_free(rnii)


def _realign_write(str base, dict result, bint mats, bint plots,
                   bint rmsabs, bint rmsrel):
    """Write mcflirt's text outputs of a :func:`vol_realign` result.

    Calls ``CAT_VolRealignWrite()``, so the files are those ``CAT_VolRealign``
    writes: ``<base>.mat/MAT_%04d``, ``<base>.par``, ``<base>_abs.rms``,
    ``<base>_rel.rms`` and the two ``_mean.rms`` files.
    """
    cdef cnp.ndarray fsl = np.ascontiguousarray(result["fsl"], dtype=np.float64).ravel()
    cdef cnp.ndarray par = np.ascontiguousarray(result["par"], dtype=np.float64).ravel()
    cdef cnp.ndarray ra = np.ascontiguousarray(result["rms_abs"], dtype=np.float64)
    cdef cnp.ndarray rr = np.ascontiguousarray(
        np.append(result["rms_rel"], 0.0), dtype=np.float64)
    cdef C.CAT_VolRealignResult res
    res.nt = ra.shape[0]
    res.ref_index = result.get("ref_index", -1)
    res.ras = NULL
    res.fsl = <double *>cnp.PyArray_DATA(fsl)
    res.par = <double *>cnp.PyArray_DATA(par)
    res.rms_abs = <double *>cnp.PyArray_DATA(ra)
    res.rms_rel = <double *>cnp.PyArray_DATA(rr)
    res.scale = NULL
    res.rms_abs_mean = result["rms_abs_mean"]
    res.rms_rel_mean = result["rms_rel_mean"]
    if C.CAT_VolRealignWrite(base.encode("utf-8"), &res, mats, plots,
                             rmsabs, rmsrel) != 0:
        raise IOError(f"could not write the outputs of '{base}'")


def _realign_image_name(str base):
    """File name of the realigned series for an output base (as mcflirt)."""
    cdef char buf[4096]
    C.CAT_VolRealignImageName(base.encode("utf-8"), buf, sizeof(buf))
    return buf.decode("utf-8")
