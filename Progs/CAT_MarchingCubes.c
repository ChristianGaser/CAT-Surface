/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl/marching.h>
#include "CAT_Separate.h"
#include "CAT_NiftiIO.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"

#define  CHUNK_SIZE   1000000

private  void  extract_isosurface(
    Volume            volume,
    double              min_label,
    double              max_label,
    int               spatial_axes[],
    General_transform *voxel_to_world_transform,
    Marching_cubes_methods  method,
    BOOLEAN           binary_flag,
    double              min_threshold,
    double              max_threshold,
    double              valid_low,
    double              valid_high,
    polygons_struct   *polygons );
    
private  void  extract_surface(
    Marching_cubes_methods  method,
    BOOLEAN           binary_flag,
    double              min_threshold,
    double              max_threshold,
    double              valid_low,
    double              valid_high,
    int               x_size,
    int               y_size,
    double              ***slices,
    double              min_label,
    double              max_label,
    double              ***label_slices,
    int               slice_index,
    BOOLEAN           right_handed,
    int               spatial_axes[],
    General_transform *voxel_to_world_transform,
    int               ***point_ids[],
    polygons_struct   *polygons );

static  STRING    dimension_names_3D[] = { MIzspace, MIyspace, MIxspace };
static  STRING    dimension_names[] = { MIyspace, MIxspace };

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: marching_cubes  input.nii  output.obj  threshold\n\
       marching_cubes  input.nii  output.obj  min_threshold max_threshold\n\
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
    double                 min_threshold, max_threshold;
    double                 min_label, max_label;
    double                 valid_low, valid_high;
    BOOLEAN              binary_flag;
    int                  i, c, spatial_axes[N_DIMENSIONS];
    int                  int_method, n_out;
    Marching_cubes_methods  method;
    object_struct        *object, **object2, *object3;
    General_transform    voxel_to_world_transform;
    polygons_struct		 *polygons;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_real_argument( 0.5, &min_threshold );
    (void) get_real_argument( min_threshold-1.0, &max_threshold );
    (void) get_int_argument( (int) MARCHING_NO_HOLES,  &int_method );
    method = (Marching_cubes_methods) int_method;

    if( max_threshold < min_threshold )
    {
        binary_flag = FALSE;
        max_threshold = min_threshold;
    }
    else
        binary_flag = TRUE;

    (void) get_real_argument( 0.0, &valid_low );
    (void) get_real_argument( -1.0, &valid_high );

    min_label = 0.0;
    max_label = -1.0;

    if (input_volume_all(input_filename, 3, dimension_names_3D,
                                     NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                     TRUE, &volume, NULL) != OK)
        return( 0 );

    copy_general_transform( get_voxel_to_world_transform(volume),
                            &voxel_to_world_transform );

    for_less( c, 0, N_DIMENSIONS )
        spatial_axes[c] = volume->spatial_axes[c];

   /* It is really weird, but only this combination worked */
    spatial_axes[0] = 0;
    spatial_axes[1] = 2;
    spatial_axes[2] = 1;

    object  = create_object( POLYGONS );
    object3 = create_object( POLYGONS );

    extract_isosurface( volume,
                        min_label, max_label,
                        spatial_axes,
                        &voxel_to_world_transform,
                        method, binary_flag,
                        min_threshold, max_threshold,
                        valid_low, valid_high, get_polygons_ptr(object) );

    polygons = get_polygons_ptr( object );
    check_polygons_neighbours_computed( polygons );
    n_out = separate_polygons( polygons, -1, &object2 );

	if( n_out > 2) fprintf(stderr,"Extract largest of %d components.\n",n_out);
	
    triangulate_polygons( get_polygons_ptr(object2[0]), get_polygons_ptr(object3) );
    
    fprintf(stderr, "Euler characteristics is %d...\n", euler_characteristic(get_polygons_ptr(object3)));

    (void) output_graphics_any_format( output_filename, ASCII_FORMAT, 1, &object3, NULL);

    delete_volume( volume );
    delete_marching_cubes_table();
    delete_general_transform( &voxel_to_world_transform );

    return( 0 );
}

private  void  clear_slice(
    Volume            volume,
    double              **slice )
{
    int    x, y, sizes[MAX_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    for_less( x, 0, sizes[0] )
    for_less( y, 0, sizes[1] )
        slice[x][y] = 0.0;
}

private  void  input_slice(
    Volume            volume,
    double            **slice,
    int               z )
{
    int    x, y, sizes[MAX_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    for_less( x, 0, sizes[0] )
    for_less( y, 0, sizes[1] )
        slice[x][y] = get_volume_real_value( volume, x, y, z, 0, 0);
}

private  double  get_slice_value(
    double      ***slices,
    int       x_size,
    int       y_size,
    int       z,
    int       x,
    int       y )
{
    if( x < 0 || x >= x_size || y < 0 || y >= y_size )
        return( 0.0 );
    else
        return( slices[z][x][y] );
}

private  void  clear_points(
    int         x_size,
    int         y_size,
    int         max_edges,
    int         ***point_ids )
{
    int    x, y, edge;

    for_less( x, 0, x_size+2 )
    {
        for_less( y, 0, y_size+2 )
        {
            for_less( edge, 0, max_edges )
                point_ids[x][y][edge] = -1;
        }
    }
}

private  void   get_world_point(
    double                slice,
    double                x,
    double                y,
    int                 spatial_axes[],
    General_transform   *voxel_to_world_transform,
    Point               *point )
{
    int            c;
    double           xw, yw, zw;
    double           real_voxel[N_DIMENSIONS], voxel_pos[N_DIMENSIONS];

    real_voxel[0] = slice;
    real_voxel[1] = x;
    real_voxel[2] = y;

    for_less( c, 0, N_DIMENSIONS )
    {
        if( spatial_axes[c] >= 0 )
            voxel_pos[c] = real_voxel[spatial_axes[c]];
        else
            voxel_pos[c] = 0.0;
    }

    general_transform_point( voxel_to_world_transform,
                             voxel_pos[X], voxel_pos[Y], voxel_pos[Z],
                             &xw, &yw, &zw );

    fill_Point( *point, xw, yw, zw );
}

private  void  extract_isosurface(
    Volume            volume,
    double              min_label,
    double              max_label,
    int               spatial_axes[],
    General_transform *voxel_to_world_transform,
    Marching_cubes_methods  method,
    BOOLEAN           binary_flag,
    double              min_threshold,
    double              max_threshold,
    double              valid_low,
    double              valid_high,
    polygons_struct   *polygons )
{
    int             n_slices, sizes[MAX_DIMENSIONS], x_size, y_size, slice;
    int             ***point_ids[2], ***tmp_point_ids;
    int             max_edges;
    double            **slices[2], **tmp_slices;
    double            **label_slices[2];
    progress_struct progress;
    Surfprop        spr;
    Point           point000, point100, point010, point001;
    Vector          v100, v010, v001, perp;
    BOOLEAN         right_handed;

    get_world_point( 0.0, 0.0, 0.0, spatial_axes, voxel_to_world_transform,
                     &point000 );
    get_world_point( 1.0, 0.0, 0.0, spatial_axes, voxel_to_world_transform,
                     &point100 );
    get_world_point( 0.0, 1.0, 0.0, spatial_axes, voxel_to_world_transform,
                     &point010 );
    get_world_point( 0.0, 0.0, 1.0, spatial_axes, voxel_to_world_transform,
                     &point001 );

    SUB_POINTS( v100, point100, point000 );
    SUB_POINTS( v010, point010, point000 );
    SUB_POINTS( v001, point001, point000 );
    CROSS_VECTORS( perp, v100, v010 );

    right_handed = DOT_VECTORS( perp, v001 ) >= 0.0;

    get_volume_sizes( volume, sizes );
    x_size = sizes[X];
    y_size = sizes[Y];
    n_slices = sizes[Z];

    ALLOC2D( slices[0], x_size, y_size );
    ALLOC2D( slices[1], x_size, y_size );

    max_edges = get_max_marching_edges( method );

    ALLOC3D( point_ids[0], x_size+2, y_size+2, max_edges );
    ALLOC3D( point_ids[1], x_size+2, y_size+2, max_edges );

    clear_slice( volume, slices[1] );

    clear_points( x_size, y_size, max_edges, point_ids[0] );
    clear_points( x_size, y_size, max_edges, point_ids[1] );

    Surfprop_a(spr) = 0.3f;
    Surfprop_d(spr) = 0.6f;
    Surfprop_s(spr) = 0.6f;
    Surfprop_se(spr) = 30.0f;
    Surfprop_t(spr) = 1.0f;
    initialize_polygons( polygons, WHITE, &spr );

    initialize_progress_report( &progress, FALSE, n_slices+1,
                                "Extracting Surface" );

    for_less( slice, -1, n_slices )
    {
        tmp_slices = slices[0];
        slices[0] = slices[1];
        slices[1] = tmp_slices;
        if( slice < n_slices - 1 )
            input_slice( volume, slices[1], slice );
        else
            clear_slice( volume, slices[1] );

        tmp_point_ids = point_ids[0];
        point_ids[0] = point_ids[1];
        point_ids[1] = tmp_point_ids;
        clear_points( x_size, y_size, max_edges, point_ids[1] );

        extract_surface( method, binary_flag, min_threshold, max_threshold,
                         valid_low, valid_high,
                         x_size, y_size, slices,
                         min_label, max_label, label_slices, slice,
                         right_handed, spatial_axes, voxel_to_world_transform,
                         point_ids, polygons );

        update_progress_report( &progress, slice+2 );
    }

    terminate_progress_report( &progress );

    if( polygons->n_points > 0 )
    {
        ALLOC( polygons->normals, polygons->n_points );
        compute_polygon_normals( polygons );
    }

    FREE2D( slices[0] );
    FREE2D( slices[1] );

    FREE3D( point_ids[0] );
    FREE3D( point_ids[1] );
}

private  int   get_point_index(
    int                 x,
    int                 y,
    int                 slice_index,
    int                 x_size,
    int                 y_size,
    voxel_point_type    *point,
    double                corners[2][2][2],
    int                 spatial_axes[],
    General_transform   *voxel_to_world_transform,
    BOOLEAN             binary_flag,
    double                min_threshold,
    double                max_threshold,
    int                 ***point_ids[],
    polygons_struct     *polygons )
{
    int            voxel[N_DIMENSIONS], edge, point_index;
    int            edge_voxel[N_DIMENSIONS];
    double           v[N_DIMENSIONS];
    Point          world_point;
    Point_classes  point_class;

    voxel[X] = x + point->coord[X];
    voxel[Y] = y + point->coord[Y];
    voxel[Z] = point->coord[Z];
    edge = point->edge_intersected;

    point_index = point_ids[voxel[Z]][voxel[X]+1][voxel[Y]+1][edge];
    if( point_index < 0 )
    {
        edge_voxel[X] = point->coord[X];
        edge_voxel[Y] = point->coord[Y];
        edge_voxel[Z] = point->coord[Z];
        point_class = get_isosurface_point( corners, edge_voxel, edge,
                                            binary_flag,
                                            min_threshold, max_threshold, v );

        get_world_point( v[Z] + (Real) slice_index,
                         v[X] + (Real) x, v[Y] + (Real) y,
                         spatial_axes, voxel_to_world_transform, &world_point );

        point_index = polygons->n_points;
        ADD_ELEMENT_TO_ARRAY( polygons->points, polygons->n_points,
                              world_point, CHUNK_SIZE );

        point_ids[voxel[Z]][voxel[X]+1][voxel[Y]+1][edge] = point_index;
    }

    return( point_index );
}

private  void  extract_surface(
    Marching_cubes_methods  method,
    BOOLEAN           binary_flag,
    double              min_threshold,
    double              max_threshold,
    double              valid_low,
    double              valid_high,
    int               x_size,
    int               y_size,
    double              ***slices,
    double              min_label,
    double              max_label,
    double              ***label_slices,
    int               slice_index,
    BOOLEAN           right_handed,
    int               spatial_axes[],
    General_transform *voxel_to_world_transform,
    int               ***point_ids[],
    polygons_struct   *polygons )
{
    int                x, y, *sizes, tx, ty, tz, n_polys, ind;
    int                p, point_index, poly, size, start_points, dir;
    voxel_point_type   *points;
    double               corners[2][2][2], label;
    BOOLEAN            valid;

    for_less( x, -1, x_size )
    {
        for_less( y, -1, y_size )
        {
            valid = TRUE;
            for_less( tx, 0, 2 )
            for_less( ty, 0, 2 )
            for_less( tz, 0, 2 )
            {
                corners[tx][ty][tz] = get_slice_value( slices, x_size, y_size,
                                                       tz, x + tx, y + ty );
                if( valid_low <= valid_high &&
                    (corners[tx][ty][tz] < min_threshold ||
                    corners[tx][ty][tz] > max_threshold) &&
                    (corners[tx][ty][tz] < valid_low ||
                     corners[tx][ty][tz] > valid_high) )
                    valid = FALSE;

                if( min_label <= max_label )
                {
                    label = get_slice_value( label_slices, x_size, y_size,
                                             tz, x + tx, y + ty );
                    if( label < min_label || label > max_label )
                        corners[tx][ty][tz] = 0.0;
                }
            }

            if( !valid )
                continue;

            n_polys = compute_isosurface_in_voxel( method, x, y, slice_index,
                                  corners, binary_flag, min_threshold,
                                  max_threshold, &sizes, &points );

            if( n_polys == 0 )
                continue;

#ifdef DEBUG
{
int   ind, i;
print( "Okay %d %d: \n", x, y );
ind = 0;
for_less( poly, 0, n_polys )
{
    for_less( i, 0, sizes[poly] )
    {
        print( "%d %d %d : %d\n", points[ind].coord[0],
               points[ind].coord[1],
               points[ind].coord[2],
               points[ind].edge_intersected );
        ++ind;
    }

    print( "\n" );
}
}
#endif

            if( right_handed )
            {
                start_points = 0;
                dir = 1;
            }
            else
            {
                start_points = sizes[0]-1;
                dir = -1;
            }

            for_less( poly, 0, n_polys )
            {
                size = sizes[poly];

                start_new_polygon( polygons );

                /*--- orient polygons properly */

                for_less( p, 0, size )
                {
                    ind = start_points + p * dir;
                    point_index = get_point_index( x, y, slice_index,
                                   x_size, y_size, &points[ind], corners,
                                   spatial_axes, voxel_to_world_transform,
                                   binary_flag, min_threshold, max_threshold,
                                   point_ids, polygons );

                    ADD_ELEMENT_TO_ARRAY( polygons->indices,
                             polygons->end_indices[polygons->n_items-1],
                             point_index, CHUNK_SIZE );
                }

                if( right_handed )
                    start_points += size;
                else if( poly < n_polys-1 )
                    start_points += sizes[poly+1];
            }
        }
    }
}

