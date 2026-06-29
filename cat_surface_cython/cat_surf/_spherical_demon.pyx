# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrapper for Spherical Demons spherical surface registration — mirrors
the ``CAT_SurfSphericalDemon`` command-line tool.

Warps a source sphere onto a template sphere by matching curvature features of
the corresponding cortical surfaces with a diffeomorphic, multi-resolution
Spherical Demons flow (Yeo et al., IEEE TMI 2010). All registration logic lives
in ``Lib/CAT_WarpDemons.c``.

The wrapper returns the warped source sphere as ``(vertices, faces)``, at the
resolution of the input source sphere.
"""

import numpy as np
cimport numpy as cnp

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Status, OK,
    create_object, get_polygons_ptr, delete_object, POLYGONS,
)
from cat_surf._convert cimport PolygonsMesh, _wrap_object
from cat_surf._convert import arrays_to_polygons, polygons_to_arrays

cimport cat_surf._warp_demons as WD

cnp.import_array()


def spherical_demon(source_surface,
                    source_sphere,
                    target_surface,
                    target_sphere,
                    int n_points=20480,
                    int n_steps=2,
                    int iters=100,
                    int curvtype0=5,
                    int curvtype1=5,
                    int curvtype2=5,
                    int curvtype3=0,
                    double fwhm_flow=5.0,
                    double fwhm_curv=1.0,
                    double fwhm_disp=10.0,
                    double max_step_deg=15.0,
                    double sigma_x=2.0,
                    double rate=1.0,
                    double step_factor=1.0,
                    bint rotate=True,
                    bint smooth_velocity=True,
                    bint smooth_displacement=True,
                    bint use_hessian=True,
                    bint use_line_search=True,
                    bint use_expmap=True,
                    bint use_tangent=False,
                    bint use_geodesic=False,
                    std_map=None,
                    double std_exp=1.0,
                    cortex_mask=None,
                    double l_dist=0.0,
                    bint verbose=False,
                    bint debug=False):
    """
    Warp a source sphere onto a template sphere via Spherical Demons.

    Mirrors ``CAT_SurfSphericalDemon``.  The source/template surfaces drive the
    feature matching (curvature), the spheres carry the deformation.  Each
    surface must share its vertex count and topology with its sphere.

    Parameters
    ----------
    source_surface, source_sphere : (vertices, faces)
        Subject's cortical surface and its spherical mapping.
    target_surface, target_sphere : (vertices, faces)
        Template surface and its sphere.
    n_points : int
        Finest pyramid resolution (default 20480); each coarser level uses
        1/4 the points.
    n_steps : int
        Number of multi-resolution pyramid levels, coarse to fine
        (1..CAT_WARP_DEMONS_MAX_STEPS, default 2).
    iters : int
        Maximum iterations per level (default 100).
    curvtype0, curvtype1, curvtype2, curvtype3 : int
        Curvature feature per level (coarse to fine).  0 mean curvature
        (3 mm, deg), 1 gaussian, 2 curvedness, 3 shape index, 4 mean
        curvature (rad), 5 sulcal-depth-like.  Defaults 5, 5, 5, 0.
    fwhm_flow : float
        Velocity-update smoothing FWHM (default 5.0).
    fwhm_curv : float
        Curvature pre-smoothing FWHM applied to both surfaces (default 1.0).
    fwhm_disp : float
        Displacement-field smoothing FWHM, the elastic prior (default 10.0).
    max_step_deg : float
        Clamp per-iteration |dtheta,dphi| to this many degrees; <=0 disables
        (default 15.0).
    sigma_x : float
        Spherical Demons Tikhonov regularization weight (default 2.0).
    rate : float
        Per-iteration multiplier for ``fwhm_flow`` (default 1.0 = constant).
    step_factor : float
        Global step-size factor (default 1.0).
    rotate : bool
        Rigid rotation pre-alignment on the coarsest level (default True).
    smooth_velocity : bool
        Smooth the velocity update, a fluid prior (default True).
    smooth_displacement : bool
        Smooth the displacement field, an elastic prior (default True).
    use_hessian : bool
        Per-vertex Gauss-Newton 2x2 Hessian update (default True).
    use_line_search : bool
        Adaptive step backtracking on stalled correlation (default True).
    use_expmap : bool
        Diffeomorphic scaling-and-squaring exponential map (default True).
        If False, falls back to an additive theta/phi flow.
    use_tangent : bool
        Prototype: compute the update in a per-vertex tangent-plane frame (as
        in Spherical Demons) instead of the global lat-lon chart.  Requires
        ``use_expmap``.  Default False.
    use_geodesic : bool
        Compose the diffeomorphic exp-map warp with geodesic (slerp)
        barycentric interpolation on the sphere instead of
        linear-then-renormalize.  Slightly more accurate, slightly slower.
        Default False.
    std_map : array_like or None
        Optional per-vertex standard deviation of the (mean-curvature) feature,
        defined on the TEMPLATE mesh (one value per ``target_sphere`` vertex).
        When given, the Gauss-Newton data term is locally weighted by
        1/variance (atlas-style template registration, as in Spherical Demons);
        the map is resampled internally to each pyramid level.  Default None.
    std_exp : float
        Exponent on the precision weight, ``w = (1/variance) ** std_exp``
        (default 1.0 = SD's 1/variance).  Raise above 1 to sharpen a
        low-contrast std map; 0 disables local weighting.  Only used when
        ``std_map`` is given.
    cortex_mask : array_like or None
        Optional per-vertex cortex mask on the TEMPLATE mesh (one value per
        ``target_sphere`` vertex; 0 excludes a vertex, e.g. the medial wall,
        >0 includes it).  Excludes non-cortex from the data term
        (FreeSurfer-style).  Independent of ``std_map`` and may be combined
        with it.  Default None.
    l_dist : float
        Weight of the metric-distortion regularizer (FreeSurfer-style distance
        term): each iteration takes a gradient step pulling warped neighbour
        distances back toward the original sphere metric, resisting local
        stretch/fold while still allowing large smooth warps.  0 disables it
        (default).  Try small values (e.g. 0.05-0.2).
    verbose : bool
        Print per-iteration progress (default False).
    debug : bool
        Write intermediate debug files (default False).

    Returns
    -------
    warped_sphere_vertices : ndarray, shape (V, 3), float64
    warped_sphere_faces    : ndarray, shape (F, 3), int32
        The deformed source sphere at the input source-sphere resolution.
    """
    if not (1 <= n_steps <= WD.CAT_WARP_DEMONS_MAX_STEPS):
        raise ValueError(
            "n_steps must be in [1, %d]" % WD.CAT_WARP_DEMONS_MAX_STEPS)
    if n_points <= 0:
        raise ValueError("n_points must be positive")
    if iters <= 0:
        raise ValueError("iters must be positive")

    src_v, src_f = source_surface
    sph_v, sph_f = source_sphere
    trg_v, trg_f = target_surface
    tsph_v, tsph_f = target_sphere

    cdef PolygonsMesh src_mesh = arrays_to_polygons(
        np.ascontiguousarray(src_v, dtype=np.float64),
        np.ascontiguousarray(src_f, dtype=np.int32))
    cdef PolygonsMesh sph_mesh = arrays_to_polygons(
        np.ascontiguousarray(sph_v, dtype=np.float64),
        np.ascontiguousarray(sph_f, dtype=np.int32))
    cdef PolygonsMesh trg_mesh = arrays_to_polygons(
        np.ascontiguousarray(trg_v, dtype=np.float64),
        np.ascontiguousarray(trg_f, dtype=np.int32))
    cdef PolygonsMesh tsph_mesh = arrays_to_polygons(
        np.ascontiguousarray(tsph_v, dtype=np.float64),
        np.ascontiguousarray(tsph_f, dtype=np.int32))

    # Fill options from the library defaults, then apply overrides (mirrors the
    # CLI in Progs/CAT_SurfSphericalDemon.c).
    cdef WD.CAT_WarpDemonsOptions opt
    WD.CAT_WarpDemonsDefaults(&opt)
    opt.n_points            = n_points
    opt.n_steps             = n_steps
    opt.iters               = iters
    opt.curvtype[0]         = curvtype0
    opt.curvtype[1]         = curvtype1
    opt.curvtype[2]         = curvtype2
    opt.curvtype[3]         = curvtype3
    opt.rotate              = 1 if rotate else 0
    opt.smooth_velocity     = 1 if smooth_velocity else 0
    opt.smooth_displacement = 1 if smooth_displacement else 0
    opt.use_hessian         = 1 if use_hessian else 0
    opt.use_line_search     = 1 if use_line_search else 0
    opt.use_expmap          = 1 if use_expmap else 0
    opt.use_tangent         = 1 if use_tangent else 0
    opt.geodesic            = 1 if use_geodesic else 0
    opt.l_dist              = l_dist
    opt.fwhm_flow           = fwhm_flow
    opt.fwhm_curv           = fwhm_curv
    opt.fwhm_disp           = fwhm_disp
    opt.rate                = rate
    opt.max_step_deg        = max_step_deg
    opt.sigma_x             = sigma_x
    opt.step_factor         = step_factor
    opt.std_exp             = std_exp
    opt.verbose             = 1 if verbose else 0
    opt.debug               = 1 if debug else 0

    # Optional template std map / cortex mask -> local data-term weighting.
    # Held in std_arr / mask_arr for the call so the opt pointers stay valid.
    cdef cnp.ndarray[cnp.float64_t, ndim=1] std_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] mask_arr
    cdef int n_trg = tsph_mesh.ptr().n_points
    if std_map is not None:
        std_arr = np.ascontiguousarray(std_map, dtype=np.float64).ravel()
        if std_arr.shape[0] != n_trg:
            raise ValueError(
                "std_map must have one value per template vertex (%d), got %d"
                % (n_trg, std_arr.shape[0]))
        opt.std_map = &std_arr[0]
    if cortex_mask is not None:
        mask_arr = np.ascontiguousarray(cortex_mask, dtype=np.float64).ravel()
        if mask_arr.shape[0] != n_trg:
            raise ValueError(
                "cortex_mask must have one value per template vertex (%d), got %d"
                % (n_trg, mask_arr.shape[0]))
        opt.cortex_mask = &mask_arr[0]

    # Build the coarse-to-fine pyramid: finest level = n_points, each coarser
    # level uses 1/4 of the points (mirrors the CLI).
    cdef int lvl, npl = n_points
    for lvl in range(opt.n_steps - 1, -1, -1):
        opt.level_points[lvl] = npl
        npl //= 4

    cdef object_struct *out_obj = create_object(POLYGONS)
    if out_obj == NULL:
        raise MemoryError()
    cdef polygons_struct *warped = get_polygons_ptr(out_obj)

    cdef Status st = WD.CAT_WarpDemonsRegister(
        src_mesh.ptr(), sph_mesh.ptr(),
        trg_mesh.ptr(), tsph_mesh.ptr(),
        warped, &opt)
    if st != OK:
        delete_object(out_obj)
        raise RuntimeError("CAT_WarpDemonsRegister failed")

    cdef PolygonsMesh warped_mesh = _wrap_object(out_obj, True)
    return polygons_to_arrays(warped_mesh)
