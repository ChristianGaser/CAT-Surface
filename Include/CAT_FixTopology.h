/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <bicpl.h>
#include <ParseArgv.h>

Volume volume;

/* argument defaults */
int bw = 1024;
int lim = 64;
int n_triangles = 81920;
double max_refine_length = 1.5;
char *reparam_file = NULL;
int do_surface_deform = 0;

object_struct **
fix_topology_sph(polygons_struct *, polygons_struct *, int, int, int, char *, double, int);
