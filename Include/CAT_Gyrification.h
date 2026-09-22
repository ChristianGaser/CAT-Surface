/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_GYRIFICATION_H_
#define _CAT_GYRIFICATION_H_

#include <bicpl.h>

#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define BW_SPH 1024

/**
 * \brief Compute global and local gyrification index using SPH-based resampling.
 *
 * \param surface     (in)  input cortical surface mesh
 * \param sphere      (in)  spherical parameterization of the surface
 * \param file        (in)  output filename for local gyrification values
 * \param n_triangles (in)  triangle count for resampling (currently unused)
 * \param reparam     (in)  optional reparameterized sphere (currently unused)
 * \return Global gyrification index (area ratio)
 */
double gyrification_index_sph(polygons_struct *surface, polygons_struct *sphere,
                              char *file, int n_triangles,
                              polygons_struct *reparam);

#endif
