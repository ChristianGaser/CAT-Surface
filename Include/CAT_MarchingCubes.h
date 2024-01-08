/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl/marching.h>

 void extract_isosurface(
    float            *vol,
    int              sizes[3],
    double           min_label,
    double           max_label,
    mat44            nii_mat,
    Marching_cubes_methods method,
    BOOLEAN          binary_flag,
    double           min_threshold,
    double           max_threshold,
    double           valid_low,
    double           valid_high,
    polygons_struct  *polygons);
