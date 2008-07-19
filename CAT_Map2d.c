/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <float.h>
#include  <CAT_Blur2d.h>
#include  <CAT_Surf.h>

#define  BINTREE_FACTOR   0.5


double compute_clockwise_rotation2( double x, double y )
{
    double radians;

    if( x == 0.0 )
    {
        if( y < 0.0 )
            return( PI / 2.0 );
        else if( y > 0.0 )
            return( 3.0 * PI / 2.0 );
        else
            return( 0.0 );
    }
    else if( y == 0.0 )
    {
        if( x > 0.0 )
            return( 0.0 );
        else
            return( PI );
    }
    else
    {
        radians = - (double) atan2( (double) y, (double) x );

        if( radians < 0.0 )
            radians += 2.0 * PI;

        return( radians );
    }
}

public  double point_to_uv_radius(
    Point            *point,
    double           *u,
    double           *v )
{
    double x, y, z, theta, phi, scale;

    scale = 100.0;
    
    x = (double) Point_x(*point)/scale;
    y = (double) Point_y(*point)/scale;
    z = (double) Point_z(*point)/scale;

    phi = acos( z );
    theta = compute_clockwise_rotation( y, x );
    *u = theta / (PI * 2.0);
    *v = phi / PI;
    return(sqrt(x*x + y*y + z*z));
}

public  void  point_to_uv(
    Point            *point,
    double           *u,
    double           *v )
{
    double x, y, z, theta, phi;

    x = (double) Point_x(*point);
    y = (double) Point_y(*point);
    z = (double) Point_z(*point);

    phi = acos( z );
    theta = compute_clockwise_rotation2( y, x );
    *u = theta / (PI * 2.0);
    *v = phi / PI;
}

public void  uv_to_point(
    double           u,
    double           v,
    Point            *point )
{
    double x, y, z, theta, phi, cos_u;
    double sin_u, cos_v, sin_v;
    
    /* shift theta by 90 to obtain correct position of midline */    
    theta = u * PI * 2.0 + PI/2.0;
    phi =   v * PI;
 
    cos_u = cos( theta );
    sin_u = sin( theta );
    cos_v = cos( phi );
    sin_v = sin( phi );

    z = cos_v;
    x = sin_v * cos_u;
    y = sin_v * sin_u;

    fill_Point( *point, x, y, z );
}

void map_smoothed_curvature_to_sphere(
    polygons_struct     *polygons,
    double              *values,
    double              *data,
    double              fwhm,
    int                 *size_map
)
{
    int                 *n_neighbours, **neighbours;
    int                 i, j, n_iter;
    double              sigma, value, *smooth_values;
    polygons_struct     unit_sphere;
    int                 x, y, point_index, value_index;
    int                 poly, size, ind;
    double              u, v;
    double              weights[1000], mn, mx;
    Point               unit_point, on_sphere_point, centre;
    Point               poly_points[1000];

    get_all_polygon_point_neighbours( polygons, &n_neighbours, &neighbours);

    // if values is empty calculate curvature
    if (values == (double *)0) {
        values = (double *)malloc(sizeof(double)*polygons->n_points);
        get_polygon_vertex_curvatures( polygons, n_neighbours, neighbours,
                                     3.0, 0.0, values );
    }
    smooth_values = (double *)malloc(sizeof(double)*polygons->n_points);

    /* calculate n_iter for sigma = 1.0 */
    n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482);
    if( n_iter == 0 )
        n_iter = 1;

    /* select sigma according fwhm */
    if( fwhm > 50.0 )
            sigma = 8.0;
    else if( fwhm > 30.0)
            sigma = 3.0;
    else if( fwhm > 20.0)
            sigma = 2.0;
    else sigma = 1.0;
				
    for ( j=0; j<n_iter; j++ )
    {
        for ( i=0; i<polygons->n_points; i++ )
        {
            heatkernel_blur_points( polygons->n_points, polygons->points, values,
                            n_neighbours[i], neighbours[i], i, sigma, NULL, &value );
            smooth_values[i] = value;
        }
        for ( i=0; i<polygons->n_points; i++ )
            values[i] = smooth_values[i];
    }

    /* create a unit sphere with same number of triangles as skin surface */
    fill_Point( centre, 0.0, 0.0, 0.0 );

    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0,
                                 polygons->n_items, &unit_sphere );

    create_polygons_bintree( &unit_sphere,
                             ROUND( (double) unit_sphere.n_items *
                                    BINTREE_FACTOR ) );
        
    for( x = 0; x < size_map[0]; x++ ) {
        for( y = 0; y < size_map[1]; y++ ) {
        
            u = ((double) x)/(double) (size_map[0] - 1);
            v = ((double) y)/(double) (size_map[1] - 1);
    
            uv_to_point(u, v, &unit_point );
            
            poly = find_closest_polygon_point( &unit_point, &unit_sphere,
                                                 &on_sphere_point );
            
            size = get_polygon_points( &unit_sphere, poly, poly_points );

            get_polygon_interpolation_weights( &on_sphere_point, size,
                                                     poly_points, weights );

            value = 0.0;
            for_less( i, 0, size ) {
                ind = unit_sphere.indices[
                         POINT_INDEX(unit_sphere.end_indices,poly,i)];
                value += weights[i] * values[ind];
            }
            value_index = x + (size_map[0]*y);
            data[value_index] = value;    
        }

    }

    // scale data to uint8 range
    mn = FLT_MAX; mx = -FLT_MAX;
    for( i=0; i<size_map[0]*size_map[1]; i++) {
	if(data[i] > mx) mx = data[i];
	if(data[i] < mn) mn = data[i];
    }
	
    for( i=0; i<size_map[0]*size_map[1]; i++) 
        data[i] = 1.0*(data[i]-mn)/(mx - mn);
		
    delete_polygons( &unit_sphere );
    free(smooth_values);
}

void map_sheet2d_to_sphere(
    double *sheet2d,
    double *values,
    polygons_struct *polygons,
    int    interpolate,
    int    *size_map
)
{
    double               tmp_x, tmp_y;
    double               u, v;
    Point                unit_point, on_sphere_point, centre;
    polygons_struct      unit_sphere;
    int                  i, x, y, point_index;

    /*--- create a unit sphere with same number of triangles as skin surface */
    fill_Point( centre, 0.0, 0.0, 0.0 );

    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0,
                               polygons->n_items, &unit_sphere );

    create_polygons_bintree( &unit_sphere,
                             ROUND( (double) unit_sphere.n_items *
                                    BINTREE_FACTOR ) );

	tmp_x = (double)size_map[0] - 1.0;
	tmp_y = (double)size_map[1] - 1.0;

    for ( i=0; i<polygons->n_points; i++ ) {

        point_to_uv(&unit_sphere.points[i], &u, &v);

		x = (int) round(u*tmp_x - 1.0);
		y = (int) round(v*tmp_y - 1.0);
		if (interpolate) {
			double xp = u*tmp_x - 1.0 - x;
			double yp = v*tmp_y - 1.0 - y;
			double xm = 1.0 - xp;
			double ym = 1.0 - yp;
			double H00 = sheet2d[bound(x,  y,  size_map)];
			double H01 = sheet2d[bound(x,  y+1,size_map)];
			double H10 = sheet2d[bound(x+1,y,  size_map)];
			double H11 = sheet2d[bound(x+1,y+1,size_map)];
			
			values[i] = (ym * ( xm * H00 + xp * H10) + 
		 		yp * ( xm * H01 + xp * H11));
		} else values[i] = sheet2d[x + y*size_map[0]];
		
    }
    delete_polygons( &unit_sphere );
}