/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

/* Some of the code is used from caret 5.3 (BrainModelSurfaceCurvature.cxx)  */

#include <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include "CAT_Curvature.h"

#define PI2 0.6366197724    /* 2/pi */

Vector projectToPlane(Vector projected, Vector basis[2])
{
   Vector xyz;
   Vector_x(xyz) = DOT_VECTORS(projected,basis[0]);
   Vector_y(xyz) = DOT_VECTORS(projected,basis[1]);
   Vector_z(xyz) = 0.0;
   return(xyz);
}

Vector projection(Vector vector, Vector normal)
{
   Vector xyz;
   const Real t2 = DOT_VECTORS(vector, normal);
   Vector_x(xyz) = Vector_x(vector) - (t2*Vector_x(normal));
   Vector_y(xyz) = Vector_y(vector) - (t2*Vector_y(normal));
   Vector_z(xyz) = Vector_z(vector) - (t2*Vector_z(normal));
   return(xyz);
}

void leastSquares(const int num, Vector dc[],
                                         Vector dn[],
                                         Real *k1,
                                         Real *k2)
{
   int i;
   Real sum1 = 0.0;
   Real sum2 = 0.0;
   Real sum3 = 0.0;
   Real wx   = 0.0;
   Real wy   = 0.0;
   Real wxy  = 0.0;
   
   for (i = 0; i < num; i++) {
      sum1 += (Vector_x(dc[i]) * Vector_x(dn[i]));
      sum2 += ((Vector_x(dc[i]) * Vector_y(dn[i])) + (Vector_y(dc[i]) * Vector_x(dn[i])));
      sum3 += (Vector_y(dc[i]) * Vector_y(dn[i]));
      wx   += (Vector_x(dc[i]) * Vector_x(dc[i]));
      wy   += (Vector_y(dc[i]) * Vector_y(dc[i]));
      wxy  += (Vector_x(dc[i]) * Vector_y(dc[i]));
   }
   
   const Real wx2  = wx  * wx;
   const Real wy2  = wy  * wy;
   const Real wxy2 = wxy * wxy;
   
   Real a = 0.0;
   Real b = 0.0;
   Real c = 0.0;
   const Real t1 = (wx + wy) * (-wxy2 + wx * wy);
   if (t1 > 0.0) {
      a = (sum3 * wxy2 - sum2 * wxy * wy +
           sum1 * (-wxy2 + wx * wy + wy2)) / t1;
      b = (-(sum3 * wx * wxy) + sum2 * wx * wy - sum1 * wxy * wy) / t1;
      c = (-(sum2 * wx * wxy) + sum1 * wxy2 +
           sum3 * (wx2 - wxy2 + wx * wy)) / t1;
   }
   
   const Real trC = a + c;
   const Real detC = a * c - b * b;
   const Real temp = trC * trC - 4 * detC;
   *k1 = 0.0;
   *k2 = 0.0;
   if (temp > 0.0) {
      const Real deltaPlus = sqrt(temp);
      const Real deltaMinus = -deltaPlus;
      *k1 = (trC + deltaPlus)  / 2.0;
      *k2 = (trC + deltaMinus) / 2.0;
   }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_points_centroid_and_normal_cg
@INPUT      : polygons
              point_index
              n_neighbours
              neighbours
@OUTPUT     : centroid
              normal
              base_length
              curvature
@RETURNS    : 
@DESCRIPTION: Computes the centroid and normal of the neighbours of a vertex,
              as well as a measure of the size of the polygon defined by the
              neighbours, and the relative curvature of the surface at the
              vertex.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  compute_points_centroid_and_normal_cg(
    polygons_struct  *polygons,
    int              point_index,
    int              n_neighbours,
    int              neighbours[],
    Point            *centroid,
    Vector           *normal,
    Real             *base_length,
    int              curvtype,
    Real             *curvparameter)
{
#define  MAX_NEIGHBOURS   1000
    int              i, j, k;
    Point            neigh_points[MAX_NEIGHBOURS];
    Vector           deltaNormal[MAX_NEIGHBOURS],deltaCoord[MAX_NEIGHBOURS];
    Vector           basis[2], t1, dn[MAX_NEIGHBOURS], dc[MAX_NEIGHBOURS];
    Vector           projected;
    Real             k1, k2;

    if( n_neighbours > 2 )
    {
        for_less( i, 0, n_neighbours )
            neigh_points[i] = polygons->points[neighbours[i]];

        get_points_centroid( n_neighbours, neigh_points, centroid );

        find_polygon_normal( n_neighbours, neigh_points, normal );
        
        for_less( i, 0, n_neighbours ) {
            SUB_VECTORS(deltaNormal[i],polygons->normals[neighbours[i]],polygons->normals[point_index]);
            SUB_VECTORS(deltaCoord[i], polygons->points[neighbours[i]], polygons->points[point_index]);
        }
        basis[0] = projection(deltaCoord[0], polygons->normals[point_index]);
        NORMALIZE_VECTOR(basis[0],basis[0]);
        fill_Vector( t1, 0.0, 0.0, 0.0 );
        
        SUB_VECTORS(t1, t1, basis[0]);
        CROSS_VECTORS(basis[1], t1, polygons->normals[point_index]);
        NORMALIZE_VECTOR(basis[1],basis[1]);

        for_less( i, 0, n_neighbours ) {
            projected = projection(deltaNormal[i], polygons->normals[point_index]);
            dn[i] = projectToPlane(projected, basis);
            
            projected = projection(deltaCoord[i], polygons->normals[point_index]);
            dc[i] = projectToPlane(projected, basis);
        }

        leastSquares(n_neighbours, dc, dn, &k1, &k2);

        switch( curvtype )
        {
        case 1: /* gauss curvature */
            *curvparameter = k1 * k2;
            break;
        case 2: /* curvedness */
            *curvparameter = sqrt((k1*k1 + k2*k2) / 2.0);
            break;
        case 3: /* shape index */
            *curvparameter = -PI2 * atan((k1 + k2) / (k1 - k2));
            break;        
        case 4: /* mean curvature */
            *curvparameter = (k1 + k2) / 2.0;
            break;
        }
    }
    else
    {
        *centroid = polygons->points[point_index];
        fill_Vector( *normal, 0.0, 0.0, 0.0 );
        *base_length = 1.0;
        *curvparameter = 0.0;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_vertex_curvatures_cg
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

public  void  get_polygon_vertex_curvatures_cg(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    Real              smoothing_distance,
    int               curvtype,    
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
                /* if smoothing_distance is > 0 mean curvature will be calculated
                and averaged over smoothing_distance */
                if( smoothing_distance <= 0.0 )
                {
                    compute_points_centroid_and_normal_cg( polygons, point_index,
                                        n_neighbours[point_index],
                                        neighbours[point_index],
                                        &centroid, &normal, &base_length,
                                        curvtype, &curvature);

                } else {
                    curvature = get_smooth_surface_curvature( polygons,
                                  n_neighbours, neighbours,
                                  poly, vertex_index,
                                  initialized, distances, smoothing_distance );

                    initialized = TRUE;
                }

                curvatures[point_index] = curvature;
            }
        }

        update_progress_report( &progress, poly + 1 );
    }

    terminate_progress_report( &progress );

    if( smoothing_distance > 0.0 )
        FREE( distances );
}
