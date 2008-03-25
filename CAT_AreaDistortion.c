/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Surf.h>

int  main(
   int   argc,
   char  *argv[] )
{
   STRING               object_filename, object2_filename, output_filename;
   FILE                 *file;
   File_formats         format;
   int                  poly, n_objects, point_index, vertex_index, size;
   object_struct        **objects, **objects2;
   polygons_struct      *polygons, *polygons2;
   float                 poly_size, area, surface_area, surface_area2;
   float                 *area_values, *area_values2, ratio, value;
   Point                points[MAX_POINTS_PER_POLYGON];
   Smallest_int         *point_done;

   initialize_argument_processing( argc, argv );

   if( !get_string_argument( NULL, &object_filename ) ||
       !get_string_argument( NULL, &object2_filename ) ||
       !get_string_argument( NULL, &output_filename ))
   {
       print_error(
             "Usage: %s  object_file object_file2 output_file\n",
             argv[0] );
       return( 1 );
   }

   if( input_graphics_file( object_filename,
                            &format, &n_objects, &objects ) != OK )
       return( 1 );

   if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
   {
       print( "File must contain 1 polygons object.\n" );
       return( 1 );
   }

   if( input_graphics_file( object2_filename,
                            &format, &n_objects, &objects2 ) != OK )
       return( 1 );

   if( n_objects != 1 || get_object_type( objects2[0] ) != POLYGONS )
   {
       print( "File must contain 1 polygons object.\n" );
       return( 1 );
   }

       if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
               return( 1 );

   polygons = get_polygons_ptr( objects[0] );
   polygons2 = get_polygons_ptr( objects2[0] );

   ALLOC( area_values, polygons->n_points );
   surface_area = get_area_of_points( polygons, area_values );

   ALLOC( area_values2, polygons2->n_points );
   surface_area2 = get_area_of_points( polygons2, area_values2 );

       ratio = surface_area/surface_area2;

   for_less( point_index, 0, polygons->n_points )
   {
       value = ratio*area_values2[point_index]/area_values[point_index];
       if( output_real( file, value ) != OK || output_newline( file ) != OK )
           break;
   }
       (void) close_file( file );

   delete_object_list( n_objects, objects );
   delete_object_list( n_objects, objects2 );

   return( 0 );
}

