# cython: language_level=3
"""
Isolated DARTEL declarations.

The 3rdparty ``dartel.h`` header has no include guards AND defines a
macro ``#define S 1.0000000001``.  Including ``CAT_SurfWarpDartel.h``
in a file that also uses Cython's stdint-based buffer helpers triggers
preprocessor pollution and compilation failures.  By isolating these
declarations in their own pxd, only ``_surf_warp.pyx`` pulls them in.
"""

from cat_surf._bic_types cimport polygons_struct, Status


cdef extern from "CAT_SurfWarpDartel.h":
    cdef struct dartel_prm:
        int     rtype
        double  rparam[5]
        double  lmreg
        int     cycles
        int     its
        int     k
        int     code

    ctypedef struct CAT_SurfWarpDartelOptions:
        int         multires_levels
        int         n_triangles
        int         verbose
        int         debug
        int         rotate
        int         curvtype0
        int         curvtype1
        int         curvtype2
        double     *fwhm
        double     *fwhm_surf
        const char *jacdet_file

    Status CAT_SurfWarpSolveDartelFlow(
        polygons_struct *src,
        polygons_struct *src_sphere,
        polygons_struct *trg,
        polygons_struct *trg_sphere,
        dartel_prm *prm,
        int dm[3],
        int n_steps,
        double rot[3],
        double *flow,
        int n_loops,
        const CAT_SurfWarpDartelOptions *opt)
