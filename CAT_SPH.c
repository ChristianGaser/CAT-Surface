/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  "Cat_SPH.h"

// not yet working (pointers!!!)
int read_SPHxyz(
    char *filename,
    int *bandwidth,
    double *rcoeffsx,
    double *rcoeffsy,
    double *rcoeffsz,
    double *icoeffsx,
    double *icoeffsy,
    double *icoeffsz
)
{
    FILE  *fp;
    char  line[256];
    int   n_dims, i;
    int   dataformat;

    // read coefficients
    if((fp = fopen(filename, "rb")) == NULL) {
        fprintf(stderr, "Error opening file %s.\n", filename);
        return(1);
    }

    fgets(line, 256, fp);
    if (strncmp(line,"SPH", 3)) {
        fprintf(stderr, "Wrong file format.\n");
        return(1);
    }
    fgets(line, 256, fp);
    sscanf(line, "%d %d %d", bandwidth, &n_dims, &dataformat);

    if (dataformat !=1) {
        fprintf(stderr, "Wrong dataformat: Data should be real and not complex.\n");
        return(1);
    }

    if (n_dims !=3) {
        fprintf(stderr, "Data dimension should be 3.\n");
        return(1);
    }

    int size = (*bandwidth) * (*bandwidth);
    rcoeffsx = (double *) malloc(sizeof(double)*size);
    rcoeffsy = (double *) malloc(sizeof(double)*size);
    rcoeffsz = (double *) malloc(sizeof(double)*size);
    icoeffsx = (double *) malloc(sizeof(double)*size);
    icoeffsy = (double *) malloc(sizeof(double)*size);
    icoeffsz = (double *) malloc(sizeof(double)*size);

    if (fread(rcoeffsx, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(icoeffsx, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(rcoeffsy, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(icoeffsy, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(rcoeffsz, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(icoeffsz, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
        
    return(fclose(fp));
}

int write_SPHxyz(
    char *filename,
    int bandwidth,
    double *rcoeffsx,
    double *rcoeffsy,
    double *rcoeffsz,
    double *icoeffsx,
    double *icoeffsy,
    double *icoeffsz
)
{
    FILE  *fp;

    // output coefficients
    fp = fopen( filename, "w" );
    fprintf(fp,"SPH\n%d %d %d\n",bandwidth, 3, 1);
    if (fwrite(rcoeffsx, sizeof(double), bandwidth*bandwidth, fp) != bandwidth*bandwidth) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    if (fwrite(icoeffsx, sizeof(double), bandwidth*bandwidth, fp) != bandwidth*bandwidth) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    if (fwrite(rcoeffsy, sizeof(double), bandwidth*bandwidth, fp) != bandwidth*bandwidth) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    if (fwrite(icoeffsy, sizeof(double), bandwidth*bandwidth, fp) != bandwidth*bandwidth) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    if (fwrite(rcoeffsz, sizeof(double), bandwidth*bandwidth, fp) != bandwidth*bandwidth) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    if (fwrite(icoeffsz, sizeof(double), bandwidth*bandwidth, fp) != bandwidth*bandwidth) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    return(fclose(fp));
}

void sample_sphere_from_sph(
    double *rdatax,
    double *rdatay,
    double *rdataz,
    polygons_struct *polygons_sphere,
    int n_triangles,
    int bandwidth
)
{
    Point                centre, new_point;    
    double               H00, H01, H10, H11, valuex, valuey, valuez;
    double               u, v;
    int                  i, size_map[2], bandwidth2, n_polygons;

    bandwidth2 = bandwidth*2;
    
	// check tetrahedral topology
	// best areal distribution of triangles is achieved for 20 edges
	n_polygons = n_triangles;
	
    while( n_polygons != 20 && n_polygons > 8 && n_polygons % 4 == 0 )
         n_polygons /= 4;

	if( n_polygons != 20 ) {
		fprintf(stderr,"Warning: Number of triangles %d is not recommend because\n", n_triangles);
		fprintf(stderr,"tetrahedral topology is not optimal.\n");
		fprintf(stderr,"Please try 20*(4*x) triangles (e.g. 81920).\n");
	}
    
    fill_Point( centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0,
                               n_triangles, polygons_sphere );

    create_polygons_bintree( polygons_sphere,
                             round( (double) polygons_sphere->n_items * BINTREE_FACTOR ) );

    size_map[0] = bandwidth2;
    size_map[1] = bandwidth2;

    for ( i=0; i<polygons_sphere->n_points; i++ ) {

        point_to_uv(&polygons_sphere->points[i], &u, &v);

        // interpolate points
        double x0 = (double)(bandwidth2 - 1) * u;
        double y0 = (double)(bandwidth2 - 1) * v;
        int x = (int) x0;
        int y = (int) y0;
        double xp = x0 - x;
        double yp = y0 - y;
        double xm = 1.0 - xp;
        double ym = 1.0 - yp;
        
        H00 = rdatax[bound(x,  y,  size_map)];
        H01 = rdatax[bound(x,  y+1,size_map)];
        H10 = rdatax[bound(x+1,y,  size_map)];
        H11 = rdatax[bound(x+1,y+1,size_map)];

        valuex = (ym * ( xm * H00 + xp * H10) + 
		    yp * ( xm * H01 + xp * H11));
 
        H00 = rdatay[bound(x,  y,  size_map)];
        H01 = rdatay[bound(x,  y+1,size_map)];
        H10 = rdatay[bound(x+1,y,  size_map)];
        H11 = rdatay[bound(x+1,y+1,size_map)];

        valuey = (ym * ( xm * H00 + xp * H10) + 
		    yp * ( xm * H01 + xp * H11));

        H00 = rdataz[bound(x,  y,  size_map)];
        H01 = rdataz[bound(x,  y+1,size_map)];
        H10 = rdataz[bound(x+1,y,  size_map)];
        H11 = rdataz[bound(x+1,y+1,size_map)];

        valuez = (ym * ( xm * H00 + xp * H10) + 
		    yp * ( xm * H01 + xp * H11));

        fill_Point( new_point,valuex,valuey,valuez);
        polygons_sphere->points[i] = new_point;
    }

    compute_polygon_normals( polygons_sphere );
}

void  replaceSPH(
    int bandwidth,
    int bandwidth_limited,
    double *coeffs,
    double *coeffs_filter)
{
    int l, m, i;

    for(l=0; l<bandwidth; l++) {
	    for(m=-l; m<l+1; m++) {
            i = seanindex(m,l,bandwidth);
            if(l<bandwidth_limited) 
                coeffs[i] = coeffs_filter[i];
        }
    }
}


void shape_description(
    int bandwidth, 
    double *rcoeffsx, 
    double *rcoeffsy,
    double *rcoeffsz,
    double *icoeffsx, 
    double *icoeffsy, 
    double *icoeffsz, 
    double *shape_descriptor)
{

    int k,l,m;

    for(l=0; l<bandwidth; l++)
    {
        shape_descriptor[l] = 0;
        for(m=-l; m<l+1; m++)
        {
            k = seanindex(m,l,bandwidth);
            shape_descriptor[l] += rcoeffsx[k]*rcoeffsx[k]+rcoeffsy[k]*rcoeffsy[k]+rcoeffsz[k]*rcoeffsz[k];
            shape_descriptor[l] += icoeffsx[k]*icoeffsx[k]+icoeffsy[k]*icoeffsy[k]+icoeffsz[k]*icoeffsz[k];
        }
    }
}

void  limit_bandwidth(
    int bandwidth,
    int bandwidth_limited,
    double *coeffs_old,
    double *coeffs_new)
{
    int l, m, i;

    for(l=0; l<bandwidth; l++) {
	    for(m=-l; m<l+1; m++) {
            i = seanindex(m,l,bandwidth);
            if(l<bandwidth_limited && l<bandwidth) {
                coeffs_new[i] = coeffs_old[i];
            } else coeffs_new[i] = 0;
        }
    }
}

void get_sph_coeffs_of_realdata(
    double     *rdata,
    int        bandwidth,
    int        dataformat,
    double     *rcoeffs,
    double     *icoeffs)
{
    double               *workspace, *weights, *idata;
    int                  rank, howmany_rank, n_objects, bandwidth2;
    fftw_plan            dctPlan ;
    fftw_plan            fftPlan ;
    fftw_iodim           dims[1], howmany_dims[1];

    bandwidth2 = bandwidth*2;
    
    idata = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    weights = (double *) malloc(sizeof(double) * 4 * bandwidth);
    workspace = (double *) malloc(sizeof(double) * 
				  ((10 * (bandwidth*bandwidth)) + (24 * bandwidth)));

    // set imag data to zero
    memset( idata, 0, sizeof(double) * bandwidth2*bandwidth2 );

    // forward DCT
    dctPlan = fftw_plan_r2r_1d( bandwidth2, weights, rdata,
			      FFTW_REDFT10, FFTW_ESTIMATE ) ;

    rank = 1 ;
    dims[0].n = bandwidth2 ;
    dims[0].is = 1 ;
    dims[0].os = bandwidth2 ;
    howmany_rank = 1 ;
    howmany_dims[0].n = bandwidth2 ;
    howmany_dims[0].is = bandwidth2 ;
    howmany_dims[0].os = 1 ;

    /* forward fft */
    fftPlan = fftw_plan_guru_split_dft( rank, dims,
				      howmany_rank, howmany_dims,
				      rdata, idata,
				      workspace, workspace+(4*bandwidth*bandwidth),
				      FFTW_ESTIMATE );

    /* now make the weights */
    makeweights( bandwidth, weights );

    FST_semi_fly(rdata, idata,
	       rcoeffs, icoeffs,
	       bandwidth,
	       workspace,
	       dataformat,
	       bandwidth, /* Use seminaive for all orders */
	       &dctPlan,
	       &fftPlan,
	       weights ) ;

    fftw_destroy_plan( fftPlan );
    fftw_destroy_plan( dctPlan );
    free( idata );
    free( workspace );
    free( weights );
}

void get_realdata_from_sph_coeffs(
    double     *rdata,
    int        bandwidth,
    int        dataformat,
    double     *rcoeffs,
    double     *icoeffs
)
{
    double               *workspace, *weights, *idata;
    int                  rank, howmany_rank, n_objects, bandwidth2;
    fftw_plan            idctPlan ;
    fftw_plan            ifftPlan ;
    fftw_iodim           dims[1], howmany_dims[1];

    bandwidth2 = bandwidth*2;
    
    idata =     (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    weights =   (double *) malloc(sizeof(double) * 4 * bandwidth);
    workspace = (double *) malloc(sizeof(double) * 
				  ((10 * (bandwidth*bandwidth)) + (24 * bandwidth)));

    // inverse DCT
    idctPlan = fftw_plan_r2r_1d( bandwidth2, weights, rcoeffs,
			      FFTW_REDFT01, FFTW_ESTIMATE ) ;

    rank = 1 ;
    dims[0].n = bandwidth2 ;
    dims[0].is = bandwidth2 ;
    dims[0].os = 1 ;
    howmany_rank = 1 ;
    howmany_dims[0].n = bandwidth2 ;
    howmany_dims[0].is = 1 ;
    howmany_dims[0].os = bandwidth2 ;

    /* forward fft */
    ifftPlan = fftw_plan_guru_split_dft( rank, dims,
				      howmany_rank, howmany_dims,
				      rcoeffs, icoeffs,
				      workspace, workspace+(4*bandwidth*bandwidth),
				      FFTW_ESTIMATE );

    /* now make the weights */
    makeweights( bandwidth, weights );

    InvFST_semi_fly(rcoeffs, icoeffs,
	       rdata, idata,
	       bandwidth,
	       workspace,
	       dataformat,
	       bandwidth, /* Use seminaive for all orders */
	       &idctPlan,
	       &ifftPlan ) ;

    fftw_destroy_plan( ifftPlan );
    fftw_destroy_plan( idctPlan );
    free( idata );
    free( workspace );
    free( weights );
}

void get_equally_sampled_coords_of_polygon(
    polygons_struct *polygons,
    polygons_struct *polygons_sphere,
    int             bandwidth,
    double          xcoord[],
    double          ycoord[],
    double          zcoord[]
)
{
    int                  i, j, x, y;
    double               value, u, v;
    Point                unit_point, on_sphere_point, new_point;
    Point                poly_points[1000], poly_points_src[1000], scaled_point;
    int                  poly, size, ind, bandwidth2;
    double               weights[1000], centre[3];

    // Determine centre of sphere based on bounds of input (to correct for shiftings)
    float bounds[6];
    get_bounds(polygons_sphere, bounds);
    for( j=0; j<3; j++ )
        centre[j] = bounds[2*j]+bounds[2*j+1];
        
    // set centre and radius
    // make radius slightly smaller to get sure that the inner side of handles will be found as
    // nearest point on the surface
    for( i=0; i<polygons_sphere->n_points; i++ ) {
        set_vector_length(&polygons_sphere->points[i], 0.99);
        for( j=0; j<3; j++ )
            Point_coord(polygons_sphere->points[i],j) -= centre[j];        

    }

    create_polygons_bintree( polygons_sphere,
                             round( (double) polygons_sphere->n_items *
                                    BINTREE_FACTOR ) );

    bandwidth2 = bandwidth*2;
    
    for( x=0; x < bandwidth2; x++ ) {
        for( y=0; y < bandwidth2; y++ ) {
        
            u = ((double) x)/(double) (bandwidth2 - 1);
            v = ((double) y)/(double) (bandwidth2 - 1);
  
            uv_to_point(u, v, &unit_point );
            
            poly = find_closest_polygon_point( &unit_point, polygons_sphere,
                                               &on_sphere_point );
            
            size = get_polygon_points( polygons_sphere, poly, poly_points_src );

            get_polygon_interpolation_weights( &on_sphere_point, size,
                                                   poly_points_src, weights );

            if( get_polygon_points( polygons, poly, poly_points ) != size )
                handle_internal_error( "map_point_between_polygons" );
                
            fill_Point( new_point, 0.0, 0.0, 0.0 );
            for( i=0; i<size; i++ ) {
                SCALE_POINT( scaled_point, poly_points[i], weights[i] );
                ADD_POINTS( new_point, new_point, scaled_point );
            }
            xcoord[x + (bandwidth2*y)] = Point_x(new_point);
            ycoord[x + (bandwidth2*y)] = Point_y(new_point);
            zcoord[x + (bandwidth2*y)] = Point_z(new_point);
        }

    }
}