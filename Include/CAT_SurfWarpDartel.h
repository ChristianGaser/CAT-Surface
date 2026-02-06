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

Status CAT_SurfWarpSolveDartelFlow(
    polygons_struct *src,
    polygons_struct *src_sphere,
    polygons_struct *trg,
    polygons_struct *trg_sphere,
    struct dartel_prm *prm,
    int dm[3],
    int n_steps,
    double rot[3],
    double *flow,
    int n_loops,
    const CAT_SurfWarpDartelOptions *opt);

#ifdef __cplusplus
}
#endif

#endif
