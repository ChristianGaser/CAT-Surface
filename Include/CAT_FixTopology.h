/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <bicpl.h>
#include <ParseArgv.h>

/* argument defaults */
int bw = 512;
int lim = 128;
int n_triangles = 81920;
double max_refine_length = 2;
double laplace_thresh = 0.01;
char *reparam_file = NULL;
int holes = 0;
int handles = 0;

object_struct ** fix_topology_sph(polygons_struct *, polygons_struct *, int,
                int, int, char *, double, int, double);
