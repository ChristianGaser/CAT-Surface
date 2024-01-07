/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/curvature.c,v 1.22 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_vertex_curvatures
@INPUT      : polygons
              smoothing_distance
              low_threshold
@OUTPUT     : curvatures
@RETURNS    : 
@DESCRIPTION: Computes the curvatures at each vertex of the polygon, using
              1 of two methods.  If smoothing distance is zero, computes
              instantaneous curvature in terms of a fractional relative
              curvature.  If non-zero, returns +/- angle in degrees of the
              smoothed curvature.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1994    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_polygon_vertex_curvatures(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    Real              smoothing_distance,
    Real              low_threshold,
    Real              curvatures[] )
{
    int              size, point_index, vertex_index, poly;
    Real             curvature, base_length;
    Smallest_int     *point_done;
    Point            centroid;
    Vector           normal;
    float            *distances;
    BOOLEAN          initialized;
    progress_struct  progress;

    compute_polygon_normals( polygons );

    ALLOC( point_done, polygons->n_points );

    for_less( point_index, 0, polygons->n_points )
        point_done[point_index] = FALSE;

    if( smoothing_distance > 0.0 )
    {
        ALLOC( distances, polygons->n_points );
        initialized = FALSE;
    }

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Computing Curvatures" );

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( vertex_index, 0, size )
        {
            point_index = polygons->indices[
                POINT_INDEX(polygons->end_indices,poly,vertex_index)];

            if( !point_done[point_index] )
            {
                point_done[point_index] = TRUE;

                if( smoothing_distance <= 0.0 )
                {
                    compute_points_centroid_and_normal( polygons, point_index,
                                        n_neighbours[point_index],
                                        neighbours[point_index],
                                        &centroid, &normal, &base_length,
                                        &curvature );

                }
                else
                {
		    /* There is some state being carried around in
		       the "distances" list.  Negative numbers are being
		       used as markers in some manner.
		    */
                    curvature = get_smooth_surface_curvature( polygons,
                                  n_neighbours, neighbours,
                                  poly, vertex_index,
                                  initialized, distances, smoothing_distance );

                    initialized = TRUE;
                }

                if( FABS( curvature ) < low_threshold )
                    curvature = 0.0;

                curvatures[point_index] = curvature;
            }
        }

        update_progress_report( &progress, poly + 1 );
    }

    terminate_progress_report( &progress );

    if( smoothing_distance > 0.0 )
        FREE( distances );
}
