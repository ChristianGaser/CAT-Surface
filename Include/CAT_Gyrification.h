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

#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define BW_SPH 1024

double gyrification_index_sph(polygons_struct *, polygons_struct *, char *, int,
                             polygons_struct *);
void find_conformal_map(polygons_struct *polygons);

#endif
