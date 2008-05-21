/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Surf.h>

#include  "ParseArgv.h"

typedef  enum { RATIO, LOG_RATIO, PERCENTAGE }
              Distortion_method;

/* argument defaults */
Distortion_method     method = RATIO; /* area disortion method - default: ratio */

/* the argument table */
ArgvInfo argTable[] = {
  { "-ratio", ARGV_CONSTANT, (char *) 0, 
    (char *) &method,
    "Use area ratio (Default): ratio = a2/a1." },
  { "-log", ARGV_CONSTANT, (char *) 1, 
    (char *) &method,
    "Use log of area ratio: ratio = log10(a2/a1)." },
  { "-percent", ARGV_CONSTANT, (char *) 2, 
    (char *) &method,
    "Use percent change of area: ratio = 100*(a2-a1)/a1." },
  
  { NULL, ARGV_END, NULL, NULL, NULL }
};


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

    /* Call ParseArgv */
    if ( ParseArgv( &argc, argv, argTable, 0 ) || ( argc < 3 ) ) {
      (void) fprintf( stderr,"\nUsage: %s [options] object_file object_file2 output_file\n", argv[0] );
      (void) fprintf( stderr,"\nCalculate area distortion between two surfaces.\n" );
      (void) fprintf(stderr, "       %s -help\n\n", argv[0]);
      return( 1 );
    }

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

   if( input_graphics_any_format( object_filename,
                            &format, &n_objects, &objects ) != OK )
       return( 1 );

   if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
   {
       print( "File must contain 1 polygons object.\n" );
       return( 1 );
   }

   if( input_graphics_any_format( object2_filename,
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
       switch( method )
       {
       case RATIO:
         value = ratio*area_values2[point_index]/area_values[point_index];
         break;
       case LOG_RATIO:
         value = log10(ratio*area_values2[point_index]/area_values[point_index]);
         break;
       case PERCENTAGE:
         value = 100*(ratio*area_values2[point_index]-area_values[point_index])/area_values[point_index];
         break;
       }

       if( output_real( file, value ) != OK || output_newline( file ) != OK )
           break;
   }
       (void) close_file( file );

   delete_object_list( n_objects, objects );
   delete_object_list( n_objects, objects2 );

   return( 0 );
}

