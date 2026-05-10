# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrappers for boundary-based registration (BBR) and volume-based
rigid co-registration — mirrors the CAT_SurfBBReg command-line tool.

All volume inputs are given as NIfTI file paths so the full header
(voxel-to-world matrix) is available for spatial sampling.  Surfaces are
passed as ``(vertices, faces)`` numpy array pairs.

For 4-D EPI series the middle frame is automatically selected, matching
the FreeSurfer ``mri_segreg`` convention.

Smoothing (``fwhm > 0``) is applied
- *symmetrically* to both images for the volume-registration initialisation
  step, and
- to the *moving* (EPI) image only before the BBR optimiser starts —
  the BBR cost function then sees the same pre-smoothed frame at every
  iteration, matching ``mri_segreg --fwhm`` behaviour.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

from cat_surf._nifti_types cimport (
    nifti_image, nifti_image_free, DT_FLOAT32,
)
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
          double fwhm=0.0,
          double grid_range_mm=15.0,
          double grid_range_rad=0.3,
          int grid_steps=3,
          int max_iter=200,
          double tol=1e-5,
          int vol_reg_levels=4,
          int vol_reg_bins=64,
          int vol_reg_iter=30,
          bint verbose=False):
    """
    Boundary-Based Registration (BBR) of a volume to cortical surfaces.

    Mirrors the ``CAT_SurfBBReg`` command-line tool.  When ``ref_file`` is
    provided an NMI volume-registration step initialises the transform before
    BBR.  When no surfaces are given the function returns the volume-
    registration result directly (standalone mode).

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
        WM-side sampling offset in mm (default 0.5).
    gm_dist : float
        GM-side absolute sampling offset in mm (default 0.5).
    slope : float
        BBR cost saturation slope (default 0.5).
    invert_contrast : int
        ``-1`` auto-detect (default), ``0`` force T1/FLAIR, ``1`` force
        T2/BOLD.
    fwhm : float
        Gaussian pre-smoothing FWHM in mm.  Applied symmetrically to both
        volumes for vol-reg, and to the moving volume before BBR.  0 = off.
    grid_range_mm : float
        Stage-1 brute-force translation grid half-width (mm).
    grid_range_rad : float
        Stage-1 brute-force rotation grid half-width (radians).
    grid_steps : int
        Stage-1 steps per DOF.
    max_iter : int
        Powell maximum iterations for BBR.
    tol : float
        Powell convergence tolerance.
    vol_reg_levels : int
        Gaussian pyramid levels for NMI initialisation.
    vol_reg_bins : int
        NMI histogram bins per axis.
    vol_reg_iter : int
        Powell iterations per pyramid level.
    verbose : bool

    Returns
    -------
    matrix : ndarray, shape (4, 4), float64
        Rigid transform in the standard output convention: moving (EPI) RAS
        → fixed (T1) RAS.  Equivalent to the matrix written by
        ``CAT_SurfBBReg`` / ``bbregister``.
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
        if fwhm > 0.0:
            vol_bbr_buf = <float *>malloc(<size_t>nvox3 * sizeof(float))
            if vol_bbr_buf == NULL:
                raise MemoryError()
            memcpy(vol_bbr_buf, vol3d, <size_t>nvox3 * sizeof(float))
            _smooth_inplace(vol_bbr_buf, dims,
                            nii.dx, nii.dy, nii.dz, fwhm)
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
        WM-side sampling offset in mm.
    gm_dist : float
        GM-side sampling offset in mm.
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
    cdef double m[16], m_inv[16]
    cdef double residual

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

        residual = C.CAT_VolumeReg_register(
            fixed_flat, fixed_nii, fixed_dims,
            moving3d,   moving_nii, moving_dims,
            &p, n_levels, sat_k, max_iter,
            1 if verbose else 0)

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
