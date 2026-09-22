/*
 * CAT_SurfWarpDartel.h
 *
 * DARTEL-based spherical registration driver used by CAT_SurfWarp.
 */

#ifndef CAT_SURFWARP_DARTEL_H
#define CAT_SURFWARP_DARTEL_H

#include <bicpl.h>

#include "dartel.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int multires_levels;   /* 0 disables multires, otherwise #coarse levels (1-3) */
    int n_triangles;       /* resampling resolution for spherical surfaces */

    int verbose;
    int debug;
    int rotate;            /* used for warnings only */

    int curvtype0;
    int curvtype1;
    int curvtype2;

    /* In/out: modified across steps (historical CAT_SurfWarp behaviour). */
    double *fwhm;
    double *fwhm_surf;

    /* Optional: if set, compute and write Jacobian det. values. */
    const char *jacdet_file;
} CAT_SurfWarpDartelOptions;

/**
 * \brief Solve a multi-resolution DARTEL flow for spherical registration.
 *
 * \param src        (in)  source surface mesh
 * \param src_sphere (in)  spherical source mesh
 * \param trg        (in)  target surface mesh
 * \param trg_sphere (in)  spherical target mesh
 * \param prm        (in)  DARTEL parameter array
 * \param dm         (in)  sheet dimensions for 2D mapping
 * \param n_steps    (in)  number of smoothing/registration steps
 * \param rot        (in/out) rotation vector updated for initial alignment
 * \param flow       (out) output flow field (2 * dm[0] * dm[1])
 * \param n_loops    (in)  number of DARTEL loops per step
 * \param opt        (in)  solver options
 * \return OK on success, ERROR on failure
 */
Status CAT_SurfWarpSolveDartelFlow(polygons_struct *src,
                                   polygons_struct *src_sphere,
                                   polygons_struct *trg,
                                   polygons_struct *trg_sphere,
                                   struct dartel_prm *prm, int dm[3],
                                   int n_steps, double rot[3], double *flow,
                                   int n_loops,
                                   const CAT_SurfWarpDartelOptions *opt);

#ifdef __cplusplus
}
#endif

#endif
