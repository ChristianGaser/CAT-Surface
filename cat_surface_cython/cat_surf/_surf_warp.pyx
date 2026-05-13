# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrapper for DARTEL-based spherical surface registration — mirrors
the ``CAT_SurfWarp`` command-line tool.

Warps a source sphere onto a template sphere by matching curvature maps
of the corresponding cortical surfaces with a multi-step DARTEL flow.

The wrapper returns the warped source sphere as ``(vertices, faces)``.
For the Jacobian-determinant or flow-field outputs available in the CLI,
fall back to invoking the binary directly.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memset

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Status, OK,
    create_object, get_polygons_ptr, delete_object,
    copy_polygons, Object_types, POLYGONS,
)

from libc.math cimport M_PI as PI_C
from cat_surf._convert cimport PolygonsMesh
from cat_surf._convert import arrays_to_polygons, polygons_to_arrays

cimport cat_surf._cat_funcs as C
cimport cat_surf._dartel as D

cnp.import_array()


def surf_warp(source_surface,
              source_sphere,
              target_surface,
              target_sphere,
              int loop=6,
              int n_steps=2,
              int n_runs=2,
              int rtype=1,
              int code=1,
              double mu=0.125,
              double lambda_=0.0,
              double lmreg=1e-3,
              int muchange=4,
              double murate=1.25,
              int curvtype0=5,
              int curvtype1=5,
              int curvtype2=2,
              double fwhm=5.0,
              double fwhm_surf=0.0,
              tuple size=(512, 256),
              int multires_levels=0,
              int n_triangles=81920,
              bint rotate=True,
              bint avg=False,
              bint verbose=False):
    """
    Warp a source sphere onto a template sphere via DARTEL.

    Mirrors ``CAT_SurfWarp``.  Returns the warped source sphere - the
    Jacobian-determinant, flow-field, and curvature-PGM outputs of the
    CLI are not exposed.

    Parameters
    ----------
    source_surface, source_sphere : (vertices, faces)
        Subject's surface and its spherical mapping.
    target_surface, target_sphere : (vertices, faces)
        Template surface and its sphere.
    loop : int
        Outer DARTEL loops (max 6, default 6).
    n_steps : int
        Number of DARTEL steps (1–3, default 2).
    n_runs : int
        Total runs / repetitions of the DARTEL pipeline (default 2).
    rtype : int
        Regularization type: 0 elastic, 1 membrane (default), 2 bending.
    code : int
        Objective: 0 sum-of-squares, 1 symmetric (default), 2 multinomial.
    mu : float
        Initial regularization parameter mu (default 0.125).
    lambda_ : float
        Regularization parameter lambda (default 0.0).
    lmreg : float
        Levenberg–Marquardt regularization (default 1e-3).
    muchange : int
        Divide mu by ``murate`` after every ``muchange`` loops (default 4).
    murate : float
        Mu divisor (default 1.25).
    curvtype0, curvtype1, curvtype2 : int
        Curvature type per step (see CAT_SurfWarp docstring; default
        5, 5, 2).
    fwhm : float
        Curvature-map smoothing FWHM (default 5.0).  Decreases per step.
    fwhm_surf : float
        Surface smoothing FWHM (default 0.0).
    size : (int, int)
        Curvature-map dimensions (default ``(512, 256)``).
    multires_levels : int
        Multi-resolution coarse levels in [0, 3] (0 = disabled).
    n_triangles : int
        Spherical resampling triangle count (default 81920).
    rotate : bool
        Run rotation-only initialization before DARTEL (default True).
    avg : bool
        Mirrors the ``-avg`` CLI flag.  On the LAST DARTEL run, rotate
        the source/template pole by 90 degrees, solve a second flow,
        un-rotate, and merge the two warped spheres via
        ``average_xz_surf``.  Reduces pole-distortion at modest extra
        cost (default False).
    verbose : bool

    Returns
    -------
    warped_sphere_vertices : ndarray, shape (V, 3), float64
    warped_sphere_faces    : ndarray, shape (F, 3), int32
    """
    if not (1 <= n_steps <= 3):
        raise ValueError("n_steps must be in [1, 3]")
    if loop <= 0 or loop > 100:
        raise ValueError("loop must be in [1, 100]")
    if not (0 <= multires_levels <= 3):
        raise ValueError("multires_levels must be in [0, 3]")

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

    cdef polygons_struct *src_sph_ptr = sph_mesh.ptr()
    cdef polygons_struct *trg_sph_ptr = tsph_mesh.ptr()

    # Normalise spheres: centre at origin, project to unit sphere.
    C.translate_to_center_of_mass(src_sph_ptr)
    cdef int i
    for i in range(src_sph_ptr.n_points):
        C.set_vector_length(&src_sph_ptr.points[i], 1.0)
    C.translate_to_center_of_mass(trg_sph_ptr)
    for i in range(trg_sph_ptr.n_points):
        C.set_vector_length(&trg_sph_ptr.points[i], 1.0)

    # Build DARTEL parameter array (one entry per loop).
    cdef D.dartel_prm *prm = <D.dartel_prm *>malloc(
        sizeof(D.dartel_prm) * 100)
    if prm == NULL:
        raise MemoryError()
    memset(prm, 0, sizeof(D.dartel_prm) * 100)

    cdef double mu_local = mu
    cdef double lam_local = lambda_
    cdef int j
    cdef int run
    cdef Status st
    for j in range(loop):
        prm[j].rparam[0] = 1.0
        prm[j].rparam[1] = 1.0
        prm[j].rtype = rtype
        prm[j].cycles = 3
        prm[j].its = 3
        prm[j].code = code
        prm[j].lmreg = lmreg
        prm[j].rparam[2] = mu_local
        prm[j].rparam[3] = lam_local
        prm[j].rparam[4] = lam_local / 2.0
        prm[j].k = j
        if (j + 1) % muchange == 0:
            mu_local /= murate
        lam_local /= 5.0

    cdef int dm[3]
    dm[0] = size[0]
    dm[1] = size[1]
    dm[2] = 1

    cdef size_t xy_size = <size_t>dm[0] * <size_t>dm[1]
    cdef double *flow = <double *>malloc(sizeof(double) * xy_size * 2)
    if flow == NULL:
        free(prm)
        raise MemoryError()
    memset(flow, 0, sizeof(double) * xy_size * 2)

    cdef double fwhm_local = fwhm
    cdef double fwhm_surf_local = fwhm_surf
    cdef D.CAT_SurfWarpDartelOptions opt
    opt.multires_levels = multires_levels
    opt.n_triangles = n_triangles
    opt.verbose = 1 if verbose else 0
    opt.debug = 0
    opt.rotate = 1 if rotate else 0
    opt.curvtype0 = curvtype0
    opt.curvtype1 = curvtype1
    opt.curvtype2 = curvtype2
    opt.fwhm = &fwhm_local
    opt.fwhm_surf = &fwhm_surf_local
    opt.jacdet_file = NULL

    cdef double rot[3]
    rot[0] = 0.0; rot[1] = 0.0; rot[2] = 0.0

    # apply_warp with the !INVERSE_WARPING (i.e. 1) flag matches the CLI.
    cdef int forward_warp = 0 if C.INVERSE_WARPING else 1

    # ------------------------------------------------------------ avg buffers
    cdef double *flow2 = NULL
    cdef double rotmat[9]
    # The CLI uses raw `polygons_struct *` (malloc'd, struct field copied in
    # via rotate_polygons / copy_polygons).  We mirror that here -- using
    # bicpl object_struct wrappers so initialize_polygons / delete_object
    # do the right thing on cleanup.
    cdef object_struct *rsrc_obj    = NULL
    cdef object_struct *rs_sph_obj  = NULL
    cdef object_struct *rtrg_obj    = NULL
    cdef object_struct *rt_sph_obj  = NULL
    cdef object_struct *as_sph_obj  = NULL
    cdef polygons_struct *rsrc_p    = NULL
    cdef polygons_struct *rs_sph_p  = NULL
    cdef polygons_struct *rtrg_p    = NULL
    cdef polygons_struct *rt_sph_p  = NULL
    cdef polygons_struct *as_sph_p  = NULL

    if avg:
        flow2 = <double *>malloc(sizeof(double) * xy_size * 2)
        if flow2 == NULL:
            free(prm); free(flow)
            raise MemoryError()
        memset(flow2, 0, sizeof(double) * xy_size * 2)

        rsrc_obj   = create_object(POLYGONS)
        rs_sph_obj = create_object(POLYGONS)
        rtrg_obj   = create_object(POLYGONS)
        rt_sph_obj = create_object(POLYGONS)
        as_sph_obj = create_object(POLYGONS)
        rsrc_p   = get_polygons_ptr(rsrc_obj)
        rs_sph_p = get_polygons_ptr(rs_sph_obj)
        rtrg_p   = get_polygons_ptr(rtrg_obj)
        rt_sph_p = get_polygons_ptr(rt_sph_obj)
        as_sph_p = get_polygons_ptr(as_sph_obj)

        # Pre-rotate the TEMPLATE pair by +90 deg about Y (same matrix the
        # CLI uses).  These rotated targets are reused across iterations.
        C.rotation_to_matrix(rotmat, 0.0, PI_C / 2.0, 0.0)
        C.rotate_polygons(trg_mesh.ptr(), rtrg_p, rotmat)
        C.rotate_polygons(trg_sph_ptr,    rt_sph_p, rotmat)

    try:
        if rotate:
            st = D.CAT_SurfWarpSolveDartelFlow(
                src_mesh.ptr(), src_sph_ptr,
                trg_mesh.ptr(), trg_sph_ptr,
                prm, dm, n_steps, rot, flow, -1, &opt)
            if st != OK:
                raise RuntimeError(
                    "CAT_SurfWarpSolveDartelFlow (rotate-init) failed")

        for run in range(n_runs):
            st = D.CAT_SurfWarpSolveDartelFlow(
                src_mesh.ptr(), src_sph_ptr,
                trg_mesh.ptr(), trg_sph_ptr,
                prm, dm, n_steps, rot, flow, loop, &opt)
            if st != OK:
                raise RuntimeError("CAT_SurfWarpSolveDartelFlow failed")

            if avg and run == (n_runs - 1):
                # Rotate the SOURCE pair to put the pole somewhere else,
                # solve a second DARTEL flow into flow2, un-rotate the
                # warped rotated sphere, then average with the main one.
                C.rotation_to_matrix(rotmat, 0.0, PI_C / 2.0, 0.0)
                C.rotate_polygons(src_mesh.ptr(), rsrc_p,   rotmat)
                C.rotate_polygons(src_sph_ptr,    rs_sph_p, rotmat)

                st = D.CAT_SurfWarpSolveDartelFlow(
                    rsrc_p, rs_sph_p,
                    rtrg_p, rt_sph_p,
                    prm, dm, n_steps, rot, flow2, loop, &opt)
                if st != OK:
                    raise RuntimeError(
                        "CAT_SurfWarpSolveDartelFlow (avg-rotated) failed")

                C.apply_warp(src_sph_ptr, src_sph_ptr, flow,  dm, forward_warp)
                C.apply_warp(rs_sph_p,    rs_sph_p,    flow2, dm, forward_warp)

                # Un-rotate by -90 deg about Y.
                C.rotation_to_matrix(rotmat, 0.0, -PI_C / 2.0, 0.0)
                C.rotate_polygons(rs_sph_p, NULL, rotmat)

                C.average_xz_surf(rs_sph_p, src_sph_ptr, as_sph_p)
                copy_polygons(as_sph_p, src_sph_ptr)
            else:
                C.apply_warp(src_sph_ptr, src_sph_ptr, flow, dm,
                             forward_warp)
    finally:
        if flow2 != NULL:
            free(flow2)
        if rsrc_obj   != NULL: delete_object(rsrc_obj)
        if rs_sph_obj != NULL: delete_object(rs_sph_obj)
        if rtrg_obj   != NULL: delete_object(rtrg_obj)
        if rt_sph_obj != NULL: delete_object(rt_sph_obj)
        if as_sph_obj != NULL: delete_object(as_sph_obj)
        free(flow)
        free(prm)

    return polygons_to_arrays(sph_mesh)
