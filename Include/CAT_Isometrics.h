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
#include <ParseArgv.h>

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

struct metricdata * getmetricdata(polygons_struct *);
int smooth(struct metricdata *, polygons_struct *, int, int, double);
int distortcorrect(struct metricdata *, polygons_struct *, int, int, double);
int stretch(struct metricdata *, polygons_struct *, int, int, int, double);

#endif
