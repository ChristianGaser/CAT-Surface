/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_RESAMPLE_H_
#define _CAT_RESAMPLE_H_

#include <bicpl.h>
#include <limits.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

void correct_shift_scale_sphere(polygons_struct *, polygons_struct *,
            polygons_struct **, polygons_struct **);
object_struct ** resample_tetrahedron(polygons_struct *, polygons_struct *,
            int, double *, double *);
object_struct ** resample_surface(polygons_struct *, polygons_struct *,
            int, double *, double *);
object_struct ** resample_surface_to_target_sphere(polygons_struct *, 
            polygons_struct *, polygons_struct *, double *, double *, int);
void resample_values_sphere_noscale(polygons_struct *, polygons_struct *, 
            double *, double *);
void resample_values_sphere(polygons_struct *, polygons_struct *, double *, 
            double *, int);
void resample_spherical_surface(polygons_struct *, polygons_struct *, polygons_struct *,
               double *, double *, int );

#endif
