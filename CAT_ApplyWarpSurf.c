/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Map2d.h>

#define  BINTREE_FACTOR   0.5

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s  surface.obj ShiftField output.obj\n\n";

    print_error( usage_str, executable );
}

int  main(
    int    argc,
    char   *argv[] )
{
    STRING              surface_filename, output_filename, vector_filename;
    FILE		        *inputfile;
    File_formats        format;
    int                 n_objects, p, ind, i;
    Real                u, v, x, y, z;
    Point               centre, unit_point, *new_points, trans_point;
    polygons_struct     *polygons, unit_sphere;
    object_struct       **objects;
    Header              raw_header;
    Real                tmp_x, tmp_y;
    Vector2D		    *image;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &surface_filename ) ||
        !get_string_argument( NULL, &vector_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_graphics_any_format( surface_filename,
                             &format, &n_objects, &objects ) != OK ||
        n_objects != 1 || get_object_type(objects[0]) != POLYGONS )
    {
        print_error( "Surface file must contain one polygons struct\n" );
        return( 1 );
    }

    polygons = get_polygons_ptr( objects[0] );

    if((inputfile = fopen(vector_filename, "rb")) == NULL)
    {
        fprintf(stderr, "Error: Couldn't read file %s.\n", vector_filename);
        return(1);
    }

    fread(&raw_header, 1, sizeof(raw_header), inputfile);
    ALLOC(image, raw_header.x*raw_header.y*2);
    fread(image, raw_header.x*raw_header.y*2, sizeof(float), inputfile);
    fclose(inputfile);

    /* don't ask me why, but sometimes there are artifacts in the vectorfile indicated
    by very large/small values. We set these values to zero */
    for_less( i, 0, raw_header.x*raw_header.y )
    {
        if(( image[i].x > 30) | ( image[i].x < -30))
            image[i].x = 0;
        if(( image[i].y > 30) | ( image[i].y < -30))
            image[i].y = 0;
    }

    create_polygons_bintree( polygons,
                             ROUND( (Real) polygons->n_items * BINTREE_FACTOR ) );

    fill_Point( centre, 0.0, 0.0, 0.0 );

    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0, polygons->n_items,
                               &unit_sphere );

    create_polygons_bintree( &unit_sphere,
                             ROUND( (Real) unit_sphere.n_items * BINTREE_FACTOR ) );

    ALLOC( new_points, polygons->n_points );
    
    tmp_x = raw_header.x - 1;
	tmp_y = raw_header.y - 1;

    for_less( p, 0, polygons->n_points )
    {

        map_point_to_unit_sphere( polygons, &polygons->points[p],
                                  &unit_sphere, &unit_point );

        point_to_uv(&unit_point, &u, &v);

        ind = ROUND(u*tmp_x) + raw_header.x*ROUND(v*tmp_y);
        
        u -= image[ind].x/(Real) tmp_x;
        v -= image[ind].y/(Real) tmp_y;

        uv_to_point(u, v, &unit_point );

        x = Point_x( unit_point );
        y = Point_y( unit_point );
        z = Point_z( unit_point );

        fill_Point( trans_point, x, y, z );

        map_unit_sphere_to_point( &unit_sphere, &trans_point,
                                  polygons, &new_points[p] );

    }

    for_less( p, 0, polygons->n_points )
        polygons->points[p] = new_points[p];

    if( output_graphics_any_format( output_filename, format, n_objects,
                              objects ) != OK )
        return( 1 );

    return( 0 );
}
