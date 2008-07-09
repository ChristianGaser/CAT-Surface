/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

/* Some of the code is used from caret 5.3 (BrainModelSurface.cxx)  */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <float.h>
#include "CAT_Surf.h"

#define pi 3.14159265358979323846264338327510

int bound(int i, int j, int dm[])
{
	int i1,j1;
	
	i1 = (((i)>=0)     ? (i)%(dm[0]) : ((dm[0])+(i)%(dm[0]))%dm[0]);
	if(j>=dm[1]) {
		i1 = dm[0]   - i1 - 1;
		j1 = 2*dm[1] - j  - 1; 	
	} else if(j<0) {
		i1 = - i1 - 1;
		j1 = - j; 	
	} else j1 = j;
	return(i1 + dm[0]*j1);
}

float* get_surface_ratio(float r, polygons_struct *polygons)
{
	int		i,j,x,y,z,a,b,c,nan=0;
	float	*lf;
	int		size, poly_size;
	float	area,asum,*avol=(float*)calloc(256*256*256,sizeof(float));
	char	str[512];
    Point   points[MAX_POINTS_PER_POLYGON];
	
	lf=(float*)calloc(polygons->n_points,sizeof(float));
		
	for(i=0;i<256*256*256;i++)
	  avol[i] = 0.0;

	for(i=0;i<polygons->n_items;i++)
	{
        size = get_polygon_points( polygons, i, points );
	    area = get_polygon_surface_area( size, points );

        poly_size = GET_OBJECT_SIZE( *polygons, i );

		for(j=0;j<poly_size;j++)
		{
		  *points = polygons->points[polygons->indices[POINT_INDEX(polygons->end_indices,i,j)]];
          x = 128 + (int) Point_x(*points);
		  y = 128 + (int) Point_y(*points);
		  z = 128 + (int) Point_z(*points);
		  if(x>=0&&x<256 && y>=0&&y<256 && z>=0&&z<256)
				avol[z*256*256+y*256+x]+=area;
		}
	}
	
	for(i=0;i<polygons->n_points;i++)
	{
		if(i%100==0)
		{
		  sprintf(str,"%i/%i",i,polygons->n_points);
		  printf(str);
		  for(j=0;j<strlen(str);j++)
			printf("\b");
		  fflush(stdout);
		}

		asum=0;
		for(x=-r;x<=r;x++)
		  for(y=-r;y<=r;y++)
		    for(z=-r;z<=r;z++)
			  if(x*x + y*y + z*z < r*r)
			  {
		        *points = polygons->points[i];
		        a = x + 128 + (int) Point_x(*points);
		        b = y + 128 + (int) Point_y(*points);
		        c = z + 128 + (int) Point_z(*points);
			    if(a>=0&&a<256 && b>=0&&b<256 && c>=0&&c<256)
					asum+=avol[c*65536+b*256+a];
			  }
		// the area of a triangle completely inside the sphere
		// will be added 3 times (one per vertex)
		// the area of a triangle partially inside the sphere		
		// will be accounted proportionally to the number
		// of vertices inside the sphere
		lf[i]=asum/(pi*r*r)/3.0;
		
		if(!(lf[i]==lf[i]))
			nan++;
	}
	free(avol);
	
	if(nan)
		printf("ERROR: there are %i NaN\n",nan);
	
	return lf;
}

float get_area_of_points(
    polygons_struct     *polygons,
    float               *area_values)
{
    Smallest_int    *point_done;
    float           poly_size, area;
    Point           points[MAX_POINTS_PER_POLYGON];
    int             point_index, poly, vertex_index, size;
    float           surface_area = 0.0;

    ALLOC( point_done, polygons->n_points );

    for_less( point_index, 0, polygons->n_points )
        point_done[point_index] = FALSE;

    for_less( poly, 0, polygons->n_items ) {    
        size = get_polygon_points( polygons, poly, points );
	    area = get_polygon_surface_area( size, points );
        surface_area += area;

        poly_size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( vertex_index, 0, poly_size ) {
            point_index = polygons->indices[
                POINT_INDEX(polygons->end_indices,poly,vertex_index)];
            if( !point_done[point_index] ) {
                point_done[point_index] = TRUE;
                area_values[point_index] = area;
            }
        }
    }
    FREE( point_done );
    return ( surface_area );
}

void translate_to_center_of_mass(
    polygons_struct     *polygons)
{
    int i, j;
    float c[3] = {0.0, 0.0, 0.0};

    for_less( i, 0, polygons->n_points )
        for_less( j, 0, 3 )
            c[j] += Point_coord(polygons->points[i],j);

    for_less( j, 0, 3 )
        c[j] /= polygons->n_points;

    for_less( i, 0, polygons->n_points )
        for_less( j, 0, 3 )
            Point_coord(polygons->points[i],j) -= c[j];
}

float get_largest_dist(
    polygons_struct     *polygons)
{
    int i, j;
    float radius = 0.0, dist = 0.0;

    for_less( i, 0, polygons->n_points ) {
        dist = 0.0;
        for_less( j, 0, 3 )
            dist += Point_coord(polygons->points[i],j)*Point_coord(polygons->points[i],j);
        radius = MAX(radius, dist);        
    }
    return( sqrt(radius) );
}


void set_vector_length(
    Point   *p,
    float   newLength)
{
    int   j;
    const float len = sqrt(Point_coord(*p,0)*Point_coord(*p,0) + 
        Point_coord(*p,1)*Point_coord(*p,1) + Point_coord(*p,2)*Point_coord(*p,2));
    if (len > 0) {
        const float scale = newLength / len;
        for_less( j, 0, 3 )
            Point_coord(*p,j) *= scale;
    }
}

void get_radius_of_points(
    polygons_struct     *polygons,
    float               *radius)
{
    int i;
    Vector xyz;
    
    for_less( i, 0, polygons->n_points ) {
        fill_Vector( xyz, Point_coord(polygons->points[i],0), Point_coord(polygons->points[i],1), Point_coord(polygons->points[i],2));
        radius[i] = MAGNITUDE(xyz); 
    }

}

void get_bounds(
    polygons_struct     *polygons,
    float bounds[6])
{
    int i;

    bounds[0] = FLT_MAX;
    bounds[1] = -FLT_MAX;
    bounds[2] = FLT_MAX;
    bounds[3] = -FLT_MAX;
    bounds[4] = FLT_MAX;
    bounds[5] = -FLT_MAX;

    for_less( i, 0, polygons->n_points ) {
        bounds[0] = MIN(bounds[0], Point_coord(polygons->points[i],0));        
        bounds[1] = MAX(bounds[1], Point_coord(polygons->points[i],0));        
        bounds[2] = MIN(bounds[2], Point_coord(polygons->points[i],1));        
        bounds[3] = MAX(bounds[3], Point_coord(polygons->points[i],1));        
        bounds[4] = MIN(bounds[4], Point_coord(polygons->points[i],2));        
        bounds[5] = MAX(bounds[5], Point_coord(polygons->points[i],2));        
    }
}

int  count_edges(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[] )
{
    int    point, n, nn, n_edges, n_duplicate_edges, neigh, n_points;

    n_points = polygons->n_points;

    n_edges = 0;
    n_duplicate_edges = 0;
    for_less( point, 0, n_points )
    {
        for_less( n, 0, n_neighbours[point] )
        {
            neigh = neighbours[point][n];
            for_less( nn, n+1, n_neighbours[point] )
            {
                if( neighbours[point][nn] == neigh )
                    break;
            }
            if( nn < n_neighbours[point] )
                ++n_duplicate_edges;
            else if( point < neigh )
                ++n_edges;
        }
    }

    if( n_duplicate_edges > 0 )
        print( "N duplicate edges: %d\n", n_duplicate_edges );

    return( n_edges );
}

void apply_warp(
  polygons_struct  *polygons,
  double           *flow,
  int              *size_map,
  int              *shift
)
{

  Point             centre, unit_point, *new_points, trans_point;
  polygons_struct   unit_sphere;
  double            inflow_x, inflow_y, u, v, x, y, z, ux, vy;
  double            indx, indy;
  int               p, ind;

  fill_Point( centre, 0.0, 0.0, 0.0 );

  create_polygons_bintree( polygons,
               round( (double) polygons->n_items * BINTREE_FACTOR ) );

  create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0, polygons->n_items,
                 &unit_sphere );

  create_polygons_bintree( &unit_sphere,
               round( (double) unit_sphere.n_items * BINTREE_FACTOR ) );

  ALLOC( new_points, polygons->n_points );
  
  inflow_x = (double)size_map[0] - 1.0;
  inflow_y = (double)size_map[1] - 1.0;

  for( p = 0; p < polygons->n_points; p++ )
  {
    map_point_to_unit_sphere( polygons, &polygons->points[p],
                  &unit_sphere, &unit_point );

    point_to_uv(&unit_point, &u, &v);

    indx = u*inflow_x;
    indy = v*inflow_y;    
    ind  = (int)round(indx) + size_map[0]*(int)round(indy);

    ux = (flow[ind] - 1.0 - indx + shift[0])/inflow_x;
    vy = (flow[ind + size_map[0]*size_map[1]] - 1.0 - indy + shift[1])/inflow_y;
    
    u += ux;
    v += vy;

    // wrap borders
	while( u < 0.0 )  u += 1.0;
	while( u >= 1.0 ) u -= 1.0;
	if( v < 0.0 )     v = 0.0;
	if( v > 1.0 )     v = 1.0;

    uv_to_point(u, v, &unit_point );

    x = Point_x( unit_point );
    y = Point_y( unit_point );
    z = Point_z( unit_point );

    fill_Point( trans_point, x, y, z );

    map_unit_sphere_to_point( &unit_sphere, &trans_point,
                  polygons, &new_points[p] );

  }

  for( p = 0; p < polygons->n_points; p++ )
    polygons->points[p] = new_points[p];

  compute_polygon_normals( polygons );

}

int euler_characteristic(
    polygons_struct   *polygons)
{
    int n_edges, *n_neighbours, **neighbours;

    create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL );
    n_edges = count_edges( polygons, n_neighbours, neighbours );
    
    delete_polygon_point_neighbours( polygons, n_neighbours,
                                     neighbours, NULL, NULL );

    create_polygon_point_neighbours( polygons, FALSE, &n_neighbours,
                                     &neighbours, NULL, NULL );

    n_edges = count_edges( polygons, n_neighbours, neighbours );

    return( polygons->n_items + polygons->n_points - n_edges );
}

void convert_ellipsoid_to_sphere_with_surface_area(
    polygons_struct     *polygons,
    float         desiredSurfaceArea)
{
    int     i, j;
   
    // Determine radius of the output sphere.
    // Note: Sphere surface area = 4 * PI * Radius * Radius
    float sphereRadius = sqrt(desiredSurfaceArea / (4.0 * PI));
   
    // Determine lengths of the axes
    float bounds[6];
    get_bounds(polygons, bounds);
    const float A = (fabs(bounds[0]) + fabs(bounds[1])) * 0.5;
    const float B = (fabs(bounds[2]) + fabs(bounds[3])) * 0.5;
    const float C = (fabs(bounds[4]) + fabs(bounds[5])) * 0.5;

    // Convert the coordinates from ellipsoid to sphere
    const float aSquared = A * A;
    const float bSquared = B * B;
    const float cSquared = C * C;
   
    float xyz[3] = { 0.0, 0.0, 0.0 };
    for_less( i, 0, polygons->n_points ) {
        for_less( j, 0, 3 )
            xyz[j] = Point_coord(polygons->points[i],j);
         
        //  ellipsoidal coordinates
        //
        //  x*x   y*y   z*z
        //  --- + --- + --- = 1.0
        //  A*A   B*B   C*C
        const float t1 = (xyz[0]*xyz[0]) / (aSquared);
        const float t2 = (xyz[1]*xyz[1]) / (bSquared);
        const float t3 = (xyz[2]*xyz[2]) / (cSquared);
        const float f = sqrt(t1 + t2 + t3);
        if (f != 0.0) {
            xyz[0] /= f;
            xyz[1] /= f;
            xyz[2] /= f;
        }
         
        // Push coordinate onto the sphere
        xyz[0] = (sphereRadius * xyz[0]) / A;
        xyz[1] = (sphereRadius * xyz[1]) / B;
        xyz[2] = (sphereRadius * xyz[2]) / C;
        for_less( j, 0, 3 )
            Point_coord(polygons->points[i],j) = xyz[j];
    }
    // Determine lengths of the axes
    get_bounds(polygons, bounds);
}

void  linear_smoothing(
    polygons_struct     *polygons,
    float                strength,
    int                 iterations,
    int                 smoothEdgesEveryXIterations,
    int                 *smoothOnlyTheseNodes,
    int                 projectToSphereEveryXIterations)
{
    int     i, j, k, l;
    int     *n_neighbours, **neighbours;

    BOOLEAN smoothSubsetOfNodes = 0;
    const float invStrength = 1.0 - strength;
    float sphereRadius = get_largest_dist( polygons );
    
    create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL );
    
    if (smoothOnlyTheseNodes != NULL) {
         smoothSubsetOfNodes = 1;
    }

    for_less( k, 1, iterations )
    {
        // see if edges should be smoothed
        BOOLEAN smoothEdges = 0;
        if (smoothEdgesEveryXIterations > 0)
            if ((k % smoothEdgesEveryXIterations) == 0)
                smoothEdges = 1;

        for_less( i, 0, polygons->n_points ) {
            BOOLEAN smoothIt = smoothEdges;
            if (smoothIt && smoothSubsetOfNodes)
                smoothIt = (smoothOnlyTheseNodes)[i];
        
            if (smoothIt) {

                if (n_neighbours[i] > 0) { 
                
                    float neighXYZ[3] = {0.0, 0.0, 0.0};    
                    for_less( j, 0, n_neighbours[i] )
                            for_less( l, 0, 3 )
                                neighXYZ[l] += (float) Point_coord(polygons->points[neighbours[i][j]],l);
                    // Update the nodes position
                    for_less( l, 0, 3 )
                        Point_coord(polygons->points[i],l) = (Point_coord(polygons->points[i],l) * invStrength) + 
                            (neighXYZ[l]/(float) n_neighbours[i] * strength);
                }
            }
        }
        
        // If the surface should be projected to a sphere
        if (projectToSphereEveryXIterations > 0)
            if ((k % projectToSphereEveryXIterations) == 0) 
                for_less( i, 0, polygons->n_points ) 
                    set_vector_length(&polygons->points[i], sphereRadius);
    }
}

void  areal_smoothing(
    polygons_struct     *polygons,
    float                strength,
    int                 iterations,
    int                 smoothEdgesEveryXIterations,
    int                 *smoothOnlyTheseNodes,
    int                 projectToSphereEveryXIterations)
{
    int     i, j, k, l;
    int     *n_neighbours, **neighbours;
    float   *area_values;
    Point   points[1000];

    BOOLEAN smoothSubsetOfNodes = 0;
    const float invStrength = 1.0 - strength;
    float sphereRadius = get_largest_dist( polygons );
    
    create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL );
    
    if (smoothOnlyTheseNodes != NULL)
         smoothSubsetOfNodes = 1;

    for_less( k, 1, iterations )
    {
        // see if edges should be smoothed
        BOOLEAN smoothEdges = 0;
        if (smoothEdgesEveryXIterations > 0)
            if ((k % smoothEdgesEveryXIterations) == 0)
                smoothEdges = 1;

        for_less( i, 0, polygons->n_points ) {
            BOOLEAN smoothIt = smoothEdges;
            if (smoothIt && smoothSubsetOfNodes)
                smoothIt = (smoothOnlyTheseNodes)[i];
        
            if (smoothIt) {

                if (n_neighbours[i] > 1) { 
                
                    float tileAreas[n_neighbours[i]];
                    float tileCenters[n_neighbours[i]*3];
                    float totalArea = 0.0;    
                    
                    // Get 2 consecutive neighbors of this node
                    for_less( j, 0, n_neighbours[i] ) {        
                    
                        const int n1 = neighbours[i][j];
                        int next = j + 1;
                        if (next >= n_neighbours[i])
                            next = 0;
                        const n2 = neighbours[i][next];
                    
                        // Area of the triangle
                        points[0] = polygons->points[i];
                        points[1] = polygons->points[n1];
                        points[2] = polygons->points[n2];
                        tileAreas[j] = get_polygon_surface_area(3, points);
                        totalArea += tileAreas[j];
                        
                        // Save center of this tile
                        for_less( l, 0, 3 ) {
                            tileCenters[j*3+l] = (Point_coord(polygons->points[i],l) + 
                                                  Point_coord(polygons->points[n1],l) +
                                                  Point_coord(polygons->points[n2],l)) / 3.0;
                        }
                    }
                    
                    // Compute the influence of the neighboring nodes
                    float xyz[3] = {0.0, 0.0, 0.0};
                    for_less( j, 0, n_neighbours[i] ) {
                        if (tileAreas[j] > 0.0) {
                            const float weight = tileAreas[j] / totalArea;
                            for_less( l, 0, 3 ) {
                                xyz[l] += weight * tileCenters[j*3+l];
                            }
                        }
                    }
                    // Update the nodes position
                    for_less( l, 0, 3 )
                        Point_coord(polygons->points[i],l) = (Point_coord(polygons->points[i],l) * invStrength) + 
                            (xyz[l] * strength);
                }
            }
        }
        
        // If the surface should be projected to a sphere
        if (projectToSphereEveryXIterations > 0)
            if ((k % projectToSphereEveryXIterations) == 0)
                for_less( i, 0, polygons->n_points )
                    set_vector_length(&polygons->points[i], sphereRadius);
    }
}

void  distance_smoothing(
    polygons_struct     *polygons,
    float               strength,
    int                 iterations,
    int                 smoothEdgesEveryXIterations,
    int                 *smoothOnlyTheseNodes,
    int                 projectToSphereEveryXIterations)
{
    int     i, j, k, l;
    int     *n_neighbours, **neighbours;

    BOOLEAN smoothSubsetOfNodes = 0;
    const float invStrength = 1.0 - strength;
    float sphereRadius = get_largest_dist( polygons );
    
    create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL );

    if (smoothOnlyTheseNodes != NULL)
         smoothSubsetOfNodes = 1;

    for_less( k, 1, iterations )
    {
        // see if edges should be smoothed
        BOOLEAN smoothEdges = 0;
        if (smoothEdgesEveryXIterations > 0)
            if ((k % smoothEdgesEveryXIterations) == 0)
                smoothEdges = 1;

        for_less( i, 0, polygons->n_points ) {
            BOOLEAN smoothIt = smoothEdges;
            if (smoothIt && smoothSubsetOfNodes)
                smoothIt = (smoothOnlyTheseNodes)[i];
        
            if (smoothIt) {

                if (n_neighbours[i] > 1) { 
                
	                // use distance weighting only for a limited number of neighbours
    	            // (usually the tetras have 5-7 neighbours)
	                if (n_neighbours[i] < 20) { 
    	                float totalDistance = 0.0;    
        	            float tileDistance[MAX_POINTS_PER_POLYGON];
            	        for_less( j, 0, n_neighbours[i] ) {  
		        	        tileDistance[j] = distance_between_points( &polygons->points[i],
    		        	        &polygons->points[neighbours[i][j]] );
                        	totalDistance += tileDistance[j];                      
	                    }

	                    // Compute the influence of the neighboring nodes
    	                float xyz[3] = {0.0, 0.0, 0.0};
        	            for_less( j, 0, n_neighbours[i] ) {
            	            if (tileDistance[j] > 0.0) {
                	            const float weight = tileDistance[j] / totalDistance;
                    	        for_less( l, 0, 3 )
                        	        xyz[l] += weight * (float) Point_coord(polygons->points[neighbours[i][j]],l);
	                        }
    	                }
    	                // Update the nodes position
        	            for_less( l, 0, 3 )
            	            Point_coord(polygons->points[i],l) = (Point_coord(polygons->points[i],l) * invStrength) + 
                	            (xyz[l] * strength);
                	            
					// use linear smoothing if number of neighbours is too large
					} else {
    	                float neighXYZ[3] = {0.0, 0.0, 0.0};    
        	            for_less( j, 0, n_neighbours[i] )
            	                for_less( l, 0, 3 )
                	                neighXYZ[l] += (float) Point_coord(polygons->points[neighbours[i][j]],l);
                    	// Update the nodes position
                    	for_less( l, 0, 3 )
                        	Point_coord(polygons->points[i],l) = (Point_coord(polygons->points[i],l) * invStrength) + 
                            	(neighXYZ[l]/(float) n_neighbours[i] * strength);
					}
                }
            }
        }
        
        // If the surface should be projected to a sphere
        if (projectToSphereEveryXIterations > 0)
            if ((k % projectToSphereEveryXIterations) == 0)
                for_less( i, 0, polygons->n_points )
                    set_vector_length(&polygons->points[i], sphereRadius);
    }
}


void inflate_surface_and_smooth_fingers(
    polygons_struct     *polygonsIn,
    const int numberSmoothingCycles,
    const float regularSmoothingStrength,
    const int regularSmoothingIterations,
    const float inflationFactorIn,
    const float compressStretchThreshold,
    const float fingerSmoothingStrength,
    const int fingerSmoothingIterations)
{

    polygons_struct     *polygons;
    int i, j, k, cycles;
    const float inflationFactor = inflationFactorIn - 1.0;
    float   *averageCompressedStretched, *maximumLinearDistortion, *averageArealCompression;
    float   *compressedStretched, *stretching, *area_values, *area_valuesIn;
    int     *n_neighbours, **neighbours, *needSmoothing;
    object_struct    *out_object;
    
    // Copy the fiducial surface since it will be modified (translated to center of mass)
    out_object = create_object( POLYGONS );
    polygons = get_polygons_ptr(out_object);
    copy_polygons(polygonsIn, polygons);
    
    // Translate the fiducial to center of mass
    translate_to_center_of_mass( polygons );
    translate_to_center_of_mass( polygonsIn );

    // Get bounds of fiducial surface
    float fiducialBounds[6];
    get_bounds(polygons, fiducialBounds);

    const float xdiff = fiducialBounds[1] - fiducialBounds[0];
    const float ydiff = fiducialBounds[3] - fiducialBounds[2];
    const float zdiff = fiducialBounds[5] - fiducialBounds[4];

    float fiducialSurfaceArea = get_polygons_surface_area( polygons );

    ALLOC( averageCompressedStretched, polygons->n_points );
    ALLOC( maximumLinearDistortion, polygons->n_points );
    ALLOC( averageArealCompression, polygons->n_points );
    ALLOC( compressedStretched, polygons->n_points );
    ALLOC( stretching, polygons->n_points );

    float surfaceAreaRatio = 0.0;

    for (cycles = 0; cycles < (numberSmoothingCycles + 1); cycles++) {
   
        if (cycles < numberSmoothingCycles) {
            // Step 6a: Apply Smoothing to AUX coord
            areal_smoothing(polygonsIn,
                        regularSmoothingStrength, 
                        regularSmoothingIterations,
                        1, NULL, 0);

            // Step 6b: Incrementally Inflate AUX Surface by Ellipsoidal Projection
            for_less( i, 0, polygons->n_points ) {
                float xyz[3];
                for_less( j, 0, 3 )
                    xyz[j] = Point_coord(polygonsIn->points[i],j);
                
                const float x = xyz[0] / xdiff;
                const float y = xyz[1] / ydiff;
                const float z = xyz[2] / zdiff;

                const float r = sqrt(x*x + y*y + z*z);

                const float k = 1.0 + inflationFactor * (1.0 - r);
                xyz[0] *= k;
                xyz[1] *= k;
                xyz[2] *= k;
                for_less( j, 0, 3 )
                    Point_coord(polygonsIn->points[i],j) = xyz[j];
            }
        }
      
        // Step 6c: Calculate surface area of this surface
        const float inflatedSurfaceArea = get_polygons_surface_area( polygonsIn );
      
        // Ratio of inflated and spherical surfaces
        surfaceAreaRatio = inflatedSurfaceArea / fiducialSurfaceArea;
      
        create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL );

        ALLOC( area_values, polygons->n_points);
        get_area_of_points( polygons, area_values);
        ALLOC( area_valuesIn, polygonsIn->n_points);
        get_area_of_points( polygonsIn, area_valuesIn);

        // Step 6d: Calculate compress/stretched value for each node
        for_less( i, 0, polygons->n_points ) {
      
        // Get position of node in both this and fiducial surface
        float nodePos[3];
        for_less( k, 0, 3 )
            nodePos[k] = Point_coord(polygonsIn->points[i],k);
        float nodeFiducialPos[3];
        for_less( k, 0, 3 )
            nodeFiducialPos[k] = Point_coord(polygons->points[i],k);
         
        maximumLinearDistortion[i] = 0.0;
        averageArealCompression[i] = 0.0;
         
        float numValidNeighbors = 0;
         
        // Loop through the neighbors
        for_less( j, 0, n_neighbours[i] ) {                        
         
            const int neighbor = neighbours[i][j];
            
            float neighPos[3];
            for_less( k, 0, 3 )
                neighPos[k] = Point_coord(polygonsIn->points[neighbor],k);
            float neighFiducialPos[3];
            for_less( k, 0, 3 )
                neighFiducialPos[k] = Point_coord(polygons->points[neighbor],k);
            
            // calculate maximum linear distortion on aux and ref surface
            float dx = neighPos[0] - nodePos[0];
            float dy = neighPos[1] - nodePos[1];
            float dz = neighPos[2] - nodePos[2];
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            dx = neighFiducialPos[0] - nodeFiducialPos[0];
            dy = neighFiducialPos[1] - nodeFiducialPos[1];
            dz = neighFiducialPos[2] - nodeFiducialPos[2];

            float fiducialDist = sqrt(dx*dx + dy*dy + dz*dz);
            if (fiducialDist > 0.0) {
               const float ratio = dist / fiducialDist;
               if (ratio > maximumLinearDistortion[i]) {
                  maximumLinearDistortion[i] = ratio;
               }
            }
            
            // Next neighbour
            int jNext = j + 1;
            if (jNext >= n_neighbours[i]) {
               jNext = 0;
            }
            
            // compute area of tiles on aux and ref surfaces
            const float tileArea = area_valuesIn[i];
            const float fiducialTileArea = area_values[i];
                                 
            //
            // average areal compression of tiles associated with node
            //
            float distort = 0.0;
            if (tileArea > 0.0) {
               distort = fiducialTileArea / tileArea;
            }
            else {
               if (fiducialTileArea != 0.0) {
                  distort = 10000.0;  // big dist since denominator zero
               }
               else {
                  distort = 1.0;  // if both zero then use assume same area
               }
            }
            
            // Zero will cause -inf
            if (distort < 0.00000001) {
               distort = 0.00000001;
            }
            
            averageArealCompression[i] += distort; //arealCompression[i] += distort;  //log(distort) / log2;
            numValidNeighbors += 1.0;
         }
                  
         if (numValidNeighbors > 0) {
            averageArealCompression[i] /= numValidNeighbors;
         }
      
         // compressed/stretched for node
         compressedStretched[i] = maximumLinearDistortion[i]
                                * averageArealCompression[i]  //arealCompression[i] 
                                * surfaceAreaRatio;
         
         // stretching for node
         stretching[i] = maximumLinearDistortion[i]
                                * sqrt(averageArealCompression[i]  //arealCompression[i] 
                                * surfaceAreaRatio);
    }

    // average compressed/stretched for all nodes by averaging with neighbors
    for_less( i, 0, polygons->n_points ) {
         averageCompressedStretched[i] = compressedStretched[i];
         
         if (n_neighbours[i] > 0) {
            for (j = 0; j < n_neighbours[i]; j++) {
               const int n = neighbours[i][j];
               averageCompressedStretched[i] += compressedStretched[n];
            }
            averageCompressedStretched[i] /= ((float)(n_neighbours[i] + 1));
         }
    }

    // Step 6e: Flag highly compressed/stretched nodes for targeted smoothing
    int numDistortionAboveThreshold = 0;
    float maxDistortion = -FLT_MAX;
    float minDistortion =  FLT_MAX;
    ALLOC( needSmoothing, polygons->n_points);
      
    for_less( i, 0, polygons->n_points ) {
         if (averageCompressedStretched[i] > compressStretchThreshold) {
            numDistortionAboveThreshold++;
            needSmoothing[i] = 1;
         }
         else needSmoothing[i] = 0;
         if (averageCompressedStretched[i] > maxDistortion) 
            maxDistortion = averageCompressedStretched[i];
         if (averageCompressedStretched[i] < minDistortion) 
            minDistortion = averageCompressedStretched[i];
    }  
             
    if (cycles < numberSmoothingCycles) {
         // Step 6f: Targeted smoothing
         areal_smoothing(polygonsIn,
                        fingerSmoothingStrength,
                        fingerSmoothingIterations,
                        1,
                        needSmoothing, 0);
      }
    }
      
    compute_polygon_normals( polygonsIn );

    FREE( area_values );
    FREE( area_valuesIn );
    FREE( averageCompressedStretched );
    FREE( maximumLinearDistortion );
    FREE( averageArealCompression );
    FREE( compressedStretched );
    FREE( stretching );
    FREE( needSmoothing );
}

