# cython: language_level=3
"""
Isolated Spherical Demons declarations.

Mirrors the public interface in ``Include/CAT_WarpDemons.h``: the options
struct, its defaults initializer, and the registration entry point. Kept in a
dedicated pxd (like ``_dartel.pxd``) so only ``_spherical_demon.pyx`` pulls in
the header.
"""

from cat_surf._bic_types cimport polygons_struct, Status


cdef extern from "CAT_WarpDemons.h":
    int CAT_WARP_DEMONS_MAX_STEPS

    ctypedef struct CAT_WarpDemonsOptions:
        int    n_points
        int    level_points[4]
        int    n_steps
        int    curvtype[4]
        int    iters
        int    rotate
        int    smooth_velocity
        int    smooth_displacement
        int    use_hessian
        int    use_line_search
        int    use_expmap
        int    use_tangent
        double fwhm_flow
        double fwhm_curv
        double fwhm_disp
        double rate
        double max_step_deg
        double sigma_x
        double step_factor
        double *std_map
        double std_exp
        double *cortex_mask
        double l_dist
        int    verbose
        int    debug

    void CAT_WarpDemonsDefaults(CAT_WarpDemonsOptions *opt)

    Status CAT_WarpDemonsRegister(polygons_struct *src,
                                  polygons_struct *src_sphere,
                                  polygons_struct *trg,
                                  polygons_struct *trg_sphere,
                                  polygons_struct *warped_src_sphere,
                                  const CAT_WarpDemonsOptions *opt)
