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
Usage: %s  surface.obj FlowField output.obj\n\n";

  print_error( usage_str, executable );
}

int  main(
  int  argc,
  char   *argv[] )
{
  STRING            surface_filename, output_filename, flow_filename;
  FILE		        *inputfile;
  File_formats      format;
  int               n_objects, i, xy_size, shift[2];
  polygons_struct   *polygons;
  object_struct     **objects;
  double            *inflow, *flow, *flow1;
  int               size_map[2];

  initialize_argument_processing( argc, argv );

  if( !get_string_argument( NULL, &surface_filename ) ||
    !get_string_argument( NULL, &flow_filename ) ||
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

  if((inputfile = fopen(flow_filename, "rb")) == NULL)
  {
    fprintf(stderr, "Error: Couldn't read file %s.\n", flow_filename);
    return(1);
  }

  fread(&size_map, 2, sizeof(int), inputfile);
  fread(&shift, 2, sizeof(int), inputfile);
  xy_size = size_map[0]*size_map[1];

  flow  = (double *)malloc(sizeof(double)*xy_size*2);
  flow1 = (double *)malloc(sizeof(double)*xy_size*2);
  inflow   = (double *)malloc(sizeof(double)*xy_size*2);
  fread(inflow, xy_size*2, sizeof(double), inputfile);
  fclose(inputfile);

  expdef(size_map, 10, inflow, flow, flow1, (double *)0, (double *)0);  
  free( flow1 );
  free( inflow );

  apply_warp(polygons, flow, size_map, shift);  

  if( output_graphics_any_format( output_filename, format, n_objects,
                objects ) != OK )
    return( 1 );

  free( flow );

  return( 0 );
}
