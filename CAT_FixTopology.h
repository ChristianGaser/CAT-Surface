/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_FixTopology.c 238 2012-03-14 23:35:41Z gaser $
 */

#include <bicpl.h>
#include <ParseArgv.h>

Volume volume;

/* argument defaults */
int bw = 1024;
int lim = 64;
int n_triangles = 327680;
char *t1_file = NULL;
char *reparam_file = NULL;

object_struct **
fix_topology_sph(polygons_struct *, polygons_struct *, int, Volume, char *, int, int, char *);
