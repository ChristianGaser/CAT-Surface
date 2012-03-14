#include  <volume_io/internal_volume_io.h>
#include  <bicpl/marching.h>

#define  CHUNK_SIZE   1000000

private  int   separate_polygons(
    polygons_struct    *polygons,
    int                desired_index,
    object_struct      **out[] );

private  void  triangulate_polygons(
    polygons_struct  *polygons,
    polygons_struct  *triangles );

private  void  extract_isosurface(
    Volume            volume,
    Volume            label_volume,
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
    Volume               volume, label_volume;
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

//    delete_volume( volume );

//    volume = create_volume( 2, dimension_names, NC_UNSPECIFIED, FALSE,
 //                           0.0, 0.0 );

    label_volume = NULL;

    object  = create_object( POLYGONS );
    object3 = create_object( POLYGONS );

    extract_isosurface( volume,
                        label_volume, min_label, max_label,
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
    (void) output_graphics_any_format( output_filename, ASCII_FORMAT, 1, &object3 );

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
    double              **slice )
{
    int    sizes[MAX_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    get_volume_value_hyperslab_2d( volume, 0, 0, sizes[0], sizes[1],
                                   &slice[0][0] );
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
    Volume            label_volume,
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

    if( label_volume != NULL )
    {
        ALLOC2D( label_slices[0], x_size, y_size );
        ALLOC2D( label_slices[1], x_size, y_size );
    }

    max_edges = get_max_marching_edges( method );

    ALLOC3D( point_ids[0], x_size+2, y_size+2, max_edges );
    ALLOC3D( point_ids[1], x_size+2, y_size+2, max_edges );

    clear_slice( volume, slices[1] );
    if( label_volume != NULL )
        clear_slice( volume, label_slices[1] );

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
            input_slice( volume, slices[1] );
        else
            clear_slice( volume, slices[1] );

        if( label_volume != NULL )
        {
            tmp_slices = label_slices[0];
            label_slices[0] = label_slices[1];
            label_slices[1] = tmp_slices;
            if( slice < n_slices - 1 )
                input_slice( label_volume, label_slices[1] );
            else
                clear_slice( volume, label_slices[1] );
        }

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

    if( label_volume != NULL )
    {
        FREE2D( label_slices[0] );
        FREE2D( label_slices[1] );
    }

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

private  BOOLEAN  recursive_triangulate_one_polygon(
    int     size,
    int     poly[],
    int     n_neighbours[],
    int     *neighbours[],
    int     indices[] )
{
    int       *left, *right, n_left, n_right, p1, p2;
    int       start_index, end_index, count, i, n;
    BOOLEAN   found;

    if( size < 3 )
        handle_internal_error( "recursive_triangulate_one_polygon" );
    
    if( size == 3 )
    {
        indices[0] = poly[0];
        indices[1] = poly[1];
        indices[2] = poly[2];
        return( TRUE );
    }

    found = FALSE;

    ALLOC( left, size );
    ALLOC( right, size );

    for_less( start_index, 0, size-2 )
    {
        p1 = poly[start_index];

        for_less( end_index, start_index+2, size )
        {
            if( start_index == 0 && end_index == size-1 )
                continue;

            p2 = poly[end_index];

            count = 0;
            for_less( n, 0, n_neighbours[p1] )
            {
                if( neighbours[p1][n] == p2 )
                    ++count;
            }

            if( count != 1 )
                continue;

            n_left = 0;
            for_inclusive( i, 0, start_index )
                left[n_left++] = poly[i];
            for_less( i, end_index, size )
                left[n_left++] = poly[i];

            n_right = 0;
            for_inclusive( i, start_index, end_index )
                right[n_right++] = poly[i];

            if( n_left + n_right != size + 2 )
            {
                handle_internal_error( "n_left" );
            }

            if( recursive_triangulate_one_polygon( n_left, left,
                                       n_neighbours, neighbours, indices ) &&
                recursive_triangulate_one_polygon( n_right, right,
                                       n_neighbours, neighbours,
                                       &indices[3*(n_left-2)] ) )
            {
                found = TRUE;
                break;
            }
        }

        if( found )
            break;
    }

    FREE( left );
    FREE( right );

    return( found );
}

private  BOOLEAN  triangulate_one_polygon(
    int     size,
    int     poly[],
    int     n_neighbours[],
    int     *neighbours[],
    int     indices[] )
{
    return( recursive_triangulate_one_polygon( size, poly,
                                n_neighbours, neighbours, indices ) );
}

private  void  triangulate_polygons(
    polygons_struct  *polygons,
    polygons_struct  *triangles )
{
    int                poly, size, index, ind, n_matches, n;
    int                *n_neighbours, **neighbours, *indices, max_size;
    progress_struct    progress;
    BOOLEAN            done;

    create_polygon_point_neighbours( polygons, TRUE,
                                     &n_neighbours, &neighbours,
                                     NULL, NULL );

    *triangles = *polygons;

    triangles->colour_flag = ONE_COLOUR;
    ALLOC( triangles->colours, 1 );
    triangles->colours[0] = polygons->colours[0];

    triangles->points = polygons->points;
    triangles->normals = polygons->normals;

    triangles->n_items = 0;
    for_less( poly, 0, polygons->n_items )
        triangles->n_items += GET_OBJECT_SIZE( *polygons, poly ) - 2;

    ALLOC( triangles->indices, 3 * triangles->n_items );

    max_size = 0;
    for_less( poly, 0, polygons->n_items )
        max_size = MAX( max_size, GET_OBJECT_SIZE( *polygons, poly ) );

    ALLOC( indices, max_size );

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Triangulating" );

    ind = 0;
    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );
        for_less( index, 0, size )
        {
            indices[index] = polygons->indices[POINT_INDEX(
                                            polygons->end_indices,poly,index)];
        }

        for_less( index, 2, size-1 )
        {
            n_matches = 0;
            for_less( n, 0, n_neighbours[indices[0]] )
            {
                if( neighbours[indices[0]][n] == indices[index] )
                    ++n_matches;
            }

            if( n_matches != 1 )
                break;
        }

        if( size > 3 && index < size-1 )
        {
            done = triangulate_one_polygon( size, indices,
                                              n_neighbours, neighbours,
                                              &triangles->indices[ind] );
            if( !done )
                print( "Could not find good triangulation: %d\n", poly );
            else
                ind += 3 * (size-2);
        }
        else
            done = FALSE;

        if( !done )
        {
            for_less( index, 1, size-1 )
            {
                triangles->indices[ind] = indices[0];
                ++ind;
                triangles->indices[ind] = indices[index];
                ++ind;
                triangles->indices[ind] = indices[index+1];
                ++ind;
            }
        }

        update_progress_report( &progress, poly+1 );
    }

    terminate_progress_report( &progress );

    if( ind != 3 * triangles->n_items )
        handle_internal_error( "Summation of ind" );

    FREE( polygons->end_indices );
    FREE( polygons->indices );
    FREE( indices );

    delete_polygon_point_neighbours( polygons, n_neighbours, neighbours,
                                     NULL, NULL );

    ALLOC( triangles->end_indices, triangles->n_items );
    for_less( poly, 0, triangles->n_items )
        triangles->end_indices[poly] = 3 * (poly + 1);
}

private  int   make_connected_components(
    polygons_struct    *polygons,
    Smallest_int       polygon_classes[],
    int                n_in_class[] )
{
    int                poly, current_poly, edge, size;
    int                neigh;
    int                n_components;
    Smallest_int       not_done;
    QUEUE_STRUCT(int)  queue;

    n_components = 0;

    not_done = (Smallest_int) 255;

    for_less( poly, 0, polygons->n_items )
        polygon_classes[poly] = not_done;

    for_less( poly, 0, polygons->n_items )
    {
        if( polygon_classes[poly] != not_done )
            continue;

        if( n_components == 255 )
        {
            ++n_components;
            break;
        }

        INITIALIZE_QUEUE( queue );
        INSERT_IN_QUEUE( queue, poly );
        polygon_classes[poly] = (Smallest_int) n_components;
        n_in_class[n_components] = 1;

        while( !IS_QUEUE_EMPTY(queue) )
        {
            REMOVE_FROM_QUEUE( queue, current_poly );
            size = GET_OBJECT_SIZE( *polygons, current_poly );

            for_less( edge, 0, size )
            {
                neigh = polygons->neighbours[
                    POINT_INDEX(polygons->end_indices,current_poly,edge)];
                if( neigh >= 0 &&
                    polygon_classes[neigh] == not_done )
                {
                    polygon_classes[neigh] = (Smallest_int) n_components;
                    ++n_in_class[n_components];
                    INSERT_IN_QUEUE( queue, neigh );
                }
            }
        }

        DELETE_QUEUE( queue );

        ++n_components;
    }

    return( n_components );
}

private  int   separate_polygons(
    polygons_struct    *polygons,
    int                desired_index,
    object_struct      **out[] )
{
    int                point, ind, p_ind, poly, vertex, size, i, j, tmp;
    int                point_index, *new_point_ids, n_objects, comp, c;
    int                biggest;
    Smallest_int       *poly_classes;
    int                n_components, *n_in_class, *ordered;
    polygons_struct    *new_poly;

    ALLOC( poly_classes, polygons->n_items );
    ALLOC( n_in_class, 256 );
    ALLOC( ordered, 256 );

    n_components = make_connected_components( polygons, poly_classes,
                                              n_in_class );

    for_less( i, 0, n_components )
        ordered[i] = i;

    for_less( i, 0, n_components-1 )
    {
        biggest = i;
        for_less( j, i+1, n_components )
        {
            if( n_in_class[ordered[j]] > n_in_class[ordered[biggest]] )
                biggest = j;
        }

        tmp = ordered[i];
        ordered[i] = ordered[biggest];
        ordered[biggest] = tmp;
    }

    ALLOC( new_point_ids, polygons->n_points );

    n_objects = 0;

    for_less( c, 0, n_components )
    {
        if( desired_index >= 0 && c != desired_index )
            continue;

        comp = ordered[c];

        for_less( point, 0, polygons->n_points )
            new_point_ids[point] = -1;

        SET_ARRAY_SIZE( *out, n_objects, n_objects+1,
                        DEFAULT_CHUNK_SIZE);
        (*out)[n_objects] = create_object( POLYGONS );
        new_poly = get_polygons_ptr( (*out)[n_objects] );
        ++n_objects;
        initialize_polygons( new_poly, WHITE, NULL );
        if( desired_index >= 0 )
        {
            new_poly->points = polygons->points;
            new_poly->normals = polygons->normals;
            new_poly->indices = polygons->indices;
            new_poly->n_items = 0;

            for_less( poly, 0, polygons->n_items )
            {
                if( poly_classes[poly] != (Smallest_int) comp )
                    continue;
                size = GET_OBJECT_SIZE( *polygons, poly );
                ++new_poly->n_items;
                for_less( vertex, 0, size )
                {
                    point_index = polygons->indices[
                              POINT_INDEX(polygons->end_indices,poly,vertex)];
                    if( new_point_ids[point_index] < 0 )
                        new_point_ids[point_index] = 0;
                }
            }

            ALLOC( new_poly->end_indices, new_poly->n_items );

            ind = 0;
            for_less( point, 0, polygons->n_points )
            {
                if( new_point_ids[point] >= 0 )
                {
                    new_point_ids[point] = ind;
                    new_poly->points[ind] = new_poly->points[point];
                    new_poly->normals[ind] = new_poly->normals[point];
                    ++ind;
                }
            }

            new_poly->n_points = ind;

            p_ind = 0;
            ind = 0;
            for_less( poly, 0, polygons->n_items )
            {
                if( poly_classes[poly] != (Smallest_int) comp )
                    continue;

                size = GET_OBJECT_SIZE( *polygons, poly );

                for_less( vertex, 0, size )
                {
                    point_index = polygons->indices[
                         POINT_INDEX(polygons->end_indices,poly,vertex)];
                    new_poly->indices[ind] = new_point_ids[point_index];
                    ++ind;
                }

                new_poly->end_indices[p_ind] = ind;
                ++p_ind;
            }
        }
        else
        {
            ind = 0;
            for_less( poly, 0, polygons->n_items )
            {
                if( poly_classes[poly] != (Smallest_int) comp )
                    continue;

                size = GET_OBJECT_SIZE( *polygons, poly );
                for_less( vertex, 0, size )
                {
                    point_index = polygons->indices[
                              POINT_INDEX(polygons->end_indices,poly,vertex)];

                    if( new_point_ids[point_index] < 0 )
                    {
                        new_point_ids[point_index] = new_poly->n_points;
                        ADD_ELEMENT_TO_ARRAY( new_poly->points,
                                              new_poly->n_points,
                                              polygons->points[point_index],
                                              DEFAULT_CHUNK_SIZE );
                        --new_poly->n_points;
                        ADD_ELEMENT_TO_ARRAY( new_poly->normals,
                                              new_poly->n_points,
                                              polygons->normals[point_index],
                                              DEFAULT_CHUNK_SIZE );
                    }

                    ADD_ELEMENT_TO_ARRAY( new_poly->indices, ind,
                                          new_point_ids[point_index],
                                          DEFAULT_CHUNK_SIZE );
                }

                ADD_ELEMENT_TO_ARRAY( new_poly->end_indices, new_poly->n_items,
                                      ind, DEFAULT_CHUNK_SIZE );
            }

            REALLOC( new_poly->points, new_poly->n_points );
            REALLOC( new_poly->normals, new_poly->n_points );
            REALLOC( new_poly->end_indices, new_poly->n_items );
            REALLOC( new_poly->indices, ind );
        }
    }

    FREE( poly_classes );
    FREE( new_point_ids );
    FREE( n_in_class );
    FREE( ordered );

    return( n_objects );
}

