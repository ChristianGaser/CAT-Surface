/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include "MarchingCubes.h"
#include "CAT_Separate.h"
#include "CAT_NiftiIO.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

static  STRING    dimension_names_3D[] = { MIzspace, MIyspace, MIxspace };
static  STRING    dimension_names[] = { MIyspace, MIxspace };

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: marching_cubes  input.nii  output_surface_file  threshold\n\
\n\
     Creates a polygonal surface of either the thresholded volume, or the\n\
     boundary of the region of values between min and max threshold.\n\
     and extracts the largest component\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               input_filename, output_filename;
    Volume               volume;
    double               min_threshold, val, xw, yw, zw;
    BOOLEAN              binary_flag;
    int                  i, j, k, c, spatial_axes[N_DIMENSIONS];
    int                  n_out, sizes[MAX_DIMENSIONS];
    object_struct        *object, **object2, *object3;
    General_transform    voxel_to_world_transform;
    polygons_struct		 *polygons;
    Point                point;
    MCB                  *mcb ;
    int                  obj_type=4;
    int                  Resx, Resy, Resz;
    double               real_voxel[N_DIMENSIONS], voxel_pos[N_DIMENSIONS];

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_real_argument( 0.5, &min_threshold );

    if (input_volume_all(input_filename, 3, dimension_names_3D,
                                     NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                     TRUE, &volume, NULL) != OK)
        return( 0 );

    get_volume_sizes( volume, sizes );
    
    copy_general_transform( get_voxel_to_world_transform(volume),
                            &voxel_to_world_transform );

    for_less( c, 0, N_DIMENSIONS )
        spatial_axes[c] = volume->spatial_axes[c];

   /* It is really weird, but only this combination worked */
    spatial_axes[0] = 2;
    spatial_axes[1] = 1;
    spatial_axes[2] = 0;

    mcb = MarchingCubes(-1, -1, -1);
    Resx = sizes[0]; Resy = sizes[1]; Resz = sizes[2];
    set_resolution( mcb, Resx, Resy, Resz) ;

    init_all(mcb) ;
   
    for (i = 0; i < sizes[0]; i++) {
        for (j = 0; j < sizes[1]; j++) {
            for (k = 0; k < sizes[2]; k++) {
                val = get_volume_real_value( volume, i, j, k, 0, 0);
                mcb->data[i + j*sizes[0] + k*sizes[0]*sizes[1]] = val - min_threshold;
            }
        }
    }
    
    run(mcb) ;
    clean_temps(mcb) ;

    object  = create_object( POLYGONS );
    object3 = create_object( POLYGONS );

    polygons = get_polygons_ptr( object );
    
    /* convert mcb structure to BIC polygon data*/
    polygons->n_items = mcb->ntrigs;
    polygons->n_points = mcb->nverts;
    ALLOC(polygons->points, polygons->n_points);
    ALLOC(polygons->normals, polygons->n_points);
    ALLOC(polygons->end_indices, polygons->n_items);
    polygons->bintree = (bintree_struct_ptr) NULL;
    for (i = 0; i < polygons->n_items; i++)
        polygons->end_indices[i] = (i + 1) * 3;
    ALLOC(polygons->indices,
        polygons->end_indices[polygons->n_items-1]);
    for (i = 0; i < polygons->n_points; i++) {
        real_voxel[0] = mcb->vertices[i].x;
        real_voxel[1] = mcb->vertices[i].y;
        real_voxel[2] = mcb->vertices[i].z;

        for(j= 0; j <N_DIMENSIONS; j++ )
        {
            if( spatial_axes[j] >= 0 )
                voxel_pos[j] = real_voxel[spatial_axes[j]];
            else
                voxel_pos[j] = 0.0;
        }
        general_transform_point( &voxel_to_world_transform,
                             voxel_pos[X], voxel_pos[Y], voxel_pos[Z],
                             &xw, &yw, &zw );

        Point_x(point) = xw;
        Point_y(point) = yw;
        Point_z(point) = zw;
        polygons->points[i] = point;
    }
    for (i = 0; i < polygons->n_items; i++) {
        polygons->indices[POINT_INDEX(polygons->end_indices, i, 0)] = mcb->triangles[i].v3;
        polygons->indices[POINT_INDEX(polygons->end_indices, i, 1)] = mcb->triangles[i].v2;
        polygons->indices[POINT_INDEX(polygons->end_indices, i, 2)] = mcb->triangles[i].v1;
    }

    compute_polygon_normals(polygons);
        
    check_polygons_neighbours_computed( polygons );
    n_out = separate_polygons( polygons, -1, &object2 );

	if( n_out > 2) printf("Extract largest of %d components.\n",n_out);
	
    triangulate_polygons( get_polygons_ptr(object2[0]), get_polygons_ptr(object3) );
    
    printf( "Euler characteristics is %d...\n", euler_characteristic(get_polygons_ptr(object3)));

    (void) output_graphics_any_format( output_filename, ASCII_FORMAT, 1, &object3, NULL);

    delete_volume( volume );
    delete_general_transform( &voxel_to_world_transform );

    return( 0 );
}



