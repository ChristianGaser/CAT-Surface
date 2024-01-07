/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1994,1995 David MacDonald,
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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/smooth_curvature.c,v 1.17 2005/08/17 22:30:25 bert Exp $";
#endif

static  int  get_smoothing_points(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    int               n_found,
    int               list[],
    Real              smoothing_distance,
    float             distances[],
    Point             *smoothing_points[] );

static  Real  get_average_curvature(
    Point   *point,
    Vector  *normal,
    int     n_smoothing_points,
    Point   smoothing_points[] );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_smooth_surface_curvature
@INPUT      : polygons
              poly
              vertex
              smoothing_distance
@OUTPUT     : 
@RETURNS    : curvature between -180 and 180
@DESCRIPTION: Computes the smooth surface curvature by finding a set of
              points the appropriate distance from the vertex and computing
              an average curvature of these points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1994    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Compute a measure of surface curvature.
 *
 * Returns a curvature measure for surface at \a vertex.
 *
 * Imagine drawing a circle on the surface of radius \a smoothing_distance,
 * centred on \a vertex.  The radius of the circle should be measured 
 * geodesically on the surface.  In this function it is approximated by
 * the link distance on the polyhedron vertex graph.  The disc will cut
 * through a number of polyhedron edges.  Call these points of intersection
 * the \em smoothing_points.  Imagine the cone with apex at \a vertex
 * that best fits the smoothing_points.  Denote the cone angle by
 * \em angle, measured in degrees.  The curvature magnitude is
 * 180 - angle.  The sign is negative if the cone points away from
 * the surface normal at \a vertex, otherwise positive.
 *
 * The precise calculation of the cone angle is as follows.  Let
 * \em centroid denote the mean location of the smoothing_points.
 * The cone angle is twice the mean angle formed by the centroid,
 * \a vertex, and the i'th smoothing point.
 *
 */
BICAPI  Real  get_smooth_surface_curvature(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    int               poly,
    int               vertex,
    BOOLEAN           distances_initialized,
    float             distances[],
    Real              smoothing_distance )
{
    int      n_smoothing_points, point_index, n_found, *list, p;
    Point    *smoothing_points;
    Real     curvature;
    BOOLEAN  alloced_distances;

    if( distances == NULL )
    {
        alloced_distances = TRUE;
        distances_initialized = FALSE;
        ALLOC( distances, polygons->n_points );
    }
    else
        alloced_distances = FALSE;

    point_index = polygons->indices[
                  POINT_INDEX(polygons->end_indices,poly,vertex)];

    /* If !distances_initialized, this routine sets them all to -1.
       Then the mesh distances from point_index are computed, up to 
       smoothing_distance.  The list of vertices visited is returned
       and the distances for the visited vertices are re-set to -1
       below.
    */
    n_found = compute_distances_from_point( polygons, n_neighbours, neighbours,
                                       &polygons->points[point_index], poly,
                                       smoothing_distance,
                                       distances_initialized, distances,
                                       &list );

    n_smoothing_points = get_smoothing_points( polygons,
                                               n_neighbours, neighbours,
                                               n_found, list,
                                               smoothing_distance,
                                               distances, &smoothing_points );

    if( alloced_distances )
        FREE( distances );
    else
    {
        for_less( p, 0, n_found )
            distances[list[p]] = -1.0f;
    }

    if( n_found > 0 )
        FREE( list );

    if( n_smoothing_points > 0 )
    {
        curvature = get_average_curvature( &polygons->points[point_index],
                                           &polygons->normals[point_index],
                                           n_smoothing_points,
                                           smoothing_points );
    }
    else
	/* I don't think this should ever be reached.  (?)
	 */
        curvature = 0.0;

    if( n_smoothing_points > 0 )
        FREE( smoothing_points );

    return( curvature );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_smoothing_points
@INPUT      : polygons
              smoothing_distance
@OUTPUT     : distances
              smoothing_points
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1994    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
/*! \brief Return intersection of disc with mesh edges.
 *
 * Returns the points at distance \a smoothing_distance from a given
 * start point.  Distances are measured in the graph, so they
 * are approximately geodesic.  The array \a distances must be
 * filled in with distances from the start point, or -1 to indicate
 * that the vertex lies at distance greater than \a smoothing_distance.
 */
static  int  get_smoothing_points(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    int               n_found,
    int               list[],
    Real              smoothing_distance,
    float             distances[],
    Point             *smoothing_points[] )
{
    int      p, neigh, point_index, prev_index, inside, outside;
    int      n_smoothing_points;
    Real     dist, ratio;
    Point    point;

    n_smoothing_points = 0;

    for_less( p, 0, n_found )
    {
        point_index = list[p];

        if( distances[point_index] < 0.0f )
            handle_internal_error( "get_smoothing_points" );

        for_less( neigh, 0, n_neighbours[point_index] )
        {
            prev_index = neighbours[point_index][neigh];

	    /* The array distances[] is expected to be initialized to -1
	       before starting the graph search in
	       compute_distance_from_point().  Upon exit from the graph
	       search, distance -1 indicates that the vertex was not visited;
	       i.e. it is further from the initial point than the maximum
	       distance.
	    */
            if( distances[prev_index] < 0.0f )
            {
                inside = point_index;
                outside = prev_index;

                dist = (Real) distances[inside] +
                                   distance_between_points(
                                         &polygons->points[inside],
                                         &polygons->points[outside] );

                if( dist != (Real) distances[inside] )
                {
                    ratio = (smoothing_distance - (Real) distances[inside]) /
                            (dist - (Real) distances[inside]);
                    INTERPOLATE_POINTS( point,
                                        polygons->points[inside],
                                        polygons->points[outside],
                                        ratio );
                    ADD_ELEMENT_TO_ARRAY( *smoothing_points, n_smoothing_points,
                                          point, DEFAULT_CHUNK_SIZE );
                }
            }
        }
    }

    return( n_smoothing_points );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_average_curvature
@INPUT      : point
              n_smoothing_points
              smoothing_points
@OUTPUT     : 
@RETURNS    : angle
@DESCRIPTION: Gets average angle computed from centroid of points, point,
              and each of the smoothing points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1994    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  get_average_curvature(
    Point   *point,
    Vector  *normal,
    int     n_smoothing_points,
    Point   smoothing_points[] )
{
    int      i;
    Real     sum_curvature, curvature, angle, sign_curvature;
    Point    centroid;
    Vector   offset;

    get_points_centroid( n_smoothing_points, smoothing_points, &centroid );

    SUB_POINTS( offset, *point, centroid );
    if( DOT_VECTORS( offset, *normal ) < 0.0 )
        sign_curvature = -1.0;
    else
        sign_curvature = 1.0;

    sum_curvature = 0.0;
    for_less( i, 0, n_smoothing_points )
    {
        angle = RAD_TO_DEG * get_angle_between_points( &smoothing_points[i],
                                                       point, &centroid );

        curvature = 180.0 - 2.0 * angle;
        sum_curvature += curvature;
    }

    curvature = sum_curvature * sign_curvature / (Real) n_smoothing_points;

#ifdef DEBUG
{
static  BOOLEAN  first = TRUE;
Real  **tags;
if( first )
{
    first = FALSE;
    ALLOC2D( tags, n_smoothing_points+1, 3 );
    tags[0][X] = RPoint_x(*point);
    tags[0][Y] = RPoint_y(*point);
    tags[0][Z] = RPoint_z(*point);
    for_less( i, 0, n_smoothing_points )
    {
        tags[i+1][X] = RPoint_x(smoothing_points[i]);
        tags[i+1][Y] = RPoint_y(smoothing_points[i]);
        tags[i+1][Z] = RPoint_z(smoothing_points[i]);
    }
    (void) output_tag_file( "smooth_curvature.tag", "", 1, n_smoothing_points+1,
                            tags, NULL, NULL, NULL, NULL, NULL );
    FREE2D( tags );
}
}
#endif

    return( curvature );
}
