/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#ifndef _CAT_FIXTOPOLOGY_H_
#define _CAT_FIXTOPOLOGY_H_

#include <bicpl.h>

/**
 * \brief Fix topological defects on a sphere using spherical harmonics.
 *
 * \param surface          (in)  original surface mesh
 * \param sphere           (in)  spherical parameterization
 * \param n_triangles      (in)  target triangle count for resampling
 * \param bw               (in)  spherical harmonic bandwidth
 * \param lim              (in)  Butterworth filter limit
 * \param reparam_file     (in)  optional reparameterization sphere file
 * \param max_refine_length (in) max edge length for refinement (<=0 disables)
 * \param force            (in)  force label for holes/handles (0 = auto)
 * \param laplace_thresh   (in)  Laplace filtering threshold (0 disables)
 * \return object list containing corrected surface (POLYGONS)
 */
object_struct ** fix_topology_sph(polygons_struct *surface,
                                  polygons_struct *sphere, int n_triangles,
                                  int bw, int lim, char *reparam_file,
                                  double max_refine_length, int force,
                                  double laplace_thresh);

#endif
