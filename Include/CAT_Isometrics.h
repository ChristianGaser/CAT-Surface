/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_ISOMETRICS_H_
#define _CAT_ISOMETRICS_H_

#include <bicpl.h>

#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

#define SELECT_OFF     0
#define SELECT_ON      1
#define LARGE_ONLY     1
#define QUIET_OFF      0
#define QUIET_ON       1

struct metricdata {
        polygons_struct *polygons;
        int *n_neigh;
        int **neigh;
        struct pointdata **ptdata;
};

struct pointdata {
        double *lengths;
        Vector *norm;
        double *areas;
};

/**
 * \brief Build per-vertex metric data for isometric mapping.
 *
 * \param polygons (in)  input surface mesh
 * \return Allocated metricdata structure (caller must free)
 */
struct metricdata * getmetricdata(polygons_struct *polygons);
/**
 * \brief Smooth the map to reduce area distortion.
 *
 * \param brain      (in)  metric data for the original surface
 * \param map        (in/out) map surface to optimize
 * \param maxiters   (in)  maximum iterations
 * \param selectflag (in)  SELECT_ON to accept only improvements
 * \param tolerance  (in)  stopping tolerance for improvement
 * \return Number of iterations performed
 */
int smooth(struct metricdata *brain, polygons_struct *map, int maxiters,
           int selectflag, double tolerance);
/**
 * \brief Correct area distortion using weighted triangle centers.
 *
 * \param brain      (in)  metric data for the original surface
 * \param map        (in/out) map surface to optimize
 * \param maxiters   (in)  maximum iterations
 * \param selectflag (in)  SELECT_ON to accept only improvements
 * \param tolerance  (in)  stopping tolerance for improvement
 * \return Number of iterations performed, or -1 on mismatch
 */
int distortcorrect(struct metricdata *brain, polygons_struct *map, int maxiters,
                   int selectflag, double tolerance);
/**
 * \brief Optimize stretch to preserve edge lengths and avoid flips.
 *
 * \param brain      (in)  metric data for the original surface
 * \param map        (in/out) map surface to optimize
 * \param maxiters   (in)  maximum iterations
 * \param selectflag (in)  SELECT_ON to accept only improvements
 * \param largeonly  (in)  LARGE_ONLY to restrict updates to large errors
 * \param tolerance  (in)  stopping tolerance for improvement
 * \return Number of iterations performed, or -1 on mismatch
 */
int stretch(struct metricdata *brain, polygons_struct *map, int maxiters,
            int selectflag, int largeonly, double tolerance);

#endif
