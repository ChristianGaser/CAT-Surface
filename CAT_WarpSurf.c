/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <float.h>
#include  <ParseArgv.h>
#include  "CAT_SheetIO.h"
#include  "CAT_Map2d.h"
#include  "CAT_Blur2d.h"
#include  "CAT_Surf.h"
  
struct dartel_prm {
  int rform;         // regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy
  double rparam[5];  // regularization parameters
  double lmreg;      // LM regularization
  int cycles;        // number of cycles for full multi grid (FMG)
  int its;           // number of relaxation iterations in each multigrid cycle
  int k;             // time steps for solving the PDE
  int code;          // objective function: 0 - sum of squares; 1 - symmetric sum of squares
};

// defaults
char *param_filename   = NULL;
char *source_filename  = NULL;
char *target_filename  = NULL;
char *jacdet_filename  = NULL;
char *output_filename  = NULL;
char *pgm_filename     = NULL;
char *outflow_filename = NULL;
char *inflow_filename  = NULL;
int translate = 0;
int code = 1;
int loop = 6;
int verbose = 0;
int rform = 1;
double reg   = 0;
double lmreg = 0.0001;
double fwhm  = 10.0;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &source_filename, 
     "Input file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &target_filename, 
     "Template file."},
  {"-j", ARGV_STRING, (char *) 1, (char *) &jacdet_filename, 
     "Jacobian determinant (-1)."},
  {"-w", ARGV_STRING, (char *) 1, (char *) &output_filename, 
     "Warped input."},
  {"-o", ARGV_STRING, (char *) 1, (char *) &pgm_filename, 
     "Warped map as pgm file."},
  {"-if", ARGV_STRING, (char *) 1, (char *) &inflow_filename, 
     "Warped input."},
  {"-of", ARGV_STRING, (char *) 1, (char *) &outflow_filename, 
     "Warped input."},
  {"-p", ARGV_STRING, (char *) 1, (char *) &param_filename, 
     "Parameter file."},
  {"-code", ARGV_INT, (char *) 1, (char *) &code,
     "Objective function (code): 0 - sum of squares; 1 - symmetric sum of squares."},
  {"-rform", ARGV_INT, (char *) 1, (char *) &rform,
     "Regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy."},
  {"-reg", ARGV_FLOAT, (char *) 1, (char *) &reg,
     "Regularization parameter."},
  {"-lmreg", ARGV_FLOAT, (char *) 1, (char *) &lmreg,
     "LM regularization."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm,
     "Filter size for curvature map in FWHM."},
  {"-loop", ARGV_INT, (char *) 1, (char *) &loop,
     "Number of outer loops for default parameters (max. 6)."},
  {"-shift", ARGV_CONSTANT, (char *) TRUE, (char *) &translate,
     "Shift map before warping."},
  {"-v", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Be verbose."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

private void translate_to_template(
  double          *source,
  double          *target,
  int             *size_map,
  int             *shift
)
{
  int mn = -10, mx = 10;
  int i, sx, sy, x, y;
  double *shifted, sdiff, diff, sdiff0 = FLT_MAX;
  
  shifted = (double *)malloc(sizeof(double)*size_map[0]*size_map[1]);

// restricted to shifts on x-scale only

  for(sx = mn; sx < mx; sx++) {
    for(sy = 0; sy < 1; sy++) {
      sdiff = 0.0;
      for(x = 0; x < size_map[0]; x++) {
        for(y = 0; y < size_map[1]; y++) {
          i = x + (size_map[0]*y);
          diff = target[i] - source[bound(x+sx, y+sy, size_map)];
          sdiff += diff*diff;
        }
      }
      if (sdiff < sdiff0) {
        shift[0] = sx;
        shift[1] = sy;
        sdiff0 = sdiff;
      }
    }
  }
  
  for(x = 0; x < size_map[0]; x++) {
    for(y = 0; y < size_map[1]; y++) {
      i = x + (size_map[0]*y);
      shifted[i] = source[bound(x+shift[0], y+shift[1], size_map)];
    }
  }

  for(i = 0; i < size_map[0]*size_map[1]; i++) source[i] = shifted[i];

  free(shifted);
  
}

private void map2d_smoothed_curvature(
  polygons_struct *polygons,
  double          *data,
  double          fwhm,
  int             *size_map
)
{
  int             *n_neighbours, **neighbours;
  int             i, j, n_iter;
  double          sigma, value, *values, *smooth_values;
  polygons_struct unit_sphere;
  int             x, y, point_index, value_index;
  int             poly, size, ind;
  double          u, v;
  double          weights[1000], mn, mx;
  Point           unit_point, on_sphere_point, centre;
  Point           poly_points[1000];

  ALLOC( values, polygons->n_points );
  ALLOC( smooth_values, polygons->n_points );

  get_all_polygon_point_neighbours( polygons, &n_neighbours, &neighbours);

  get_polygon_vertex_curvatures( polygons, n_neighbours, neighbours,
                   3.0, 0.0, values );

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
  FREE(values);
  FREE(smooth_values);
}

int  main(
  int   argc,
  char  *argv[] )
{
  File_formats     format;
  FILE             *fp, *fp_flow;
  char             line[1024];
  polygons_struct  *polygons_source, *polygons_target;
  int              x, y, i, j, it, it0, it1, n_objects, it_scratch, xy_size;
  double           value;
  double           *map_source, *map_target, *map_warp, *flow, *flow1, *inflow, *scratch, *jd, *jd1;
  object_struct    **objects, *object;
  double           ll[3];
  static double    param[3] = {1.0, 1.0, 1.0};
  int              size_map[3] = {360, 360, 1}, shift[2];

  /* get the arguments from the command line */

  // Get arguments
  if (ParseArgv(&argc, argv, argTable, 0) ||
    ( source_filename == NULL) ||
    ( target_filename == NULL) ||
    ((jacdet_filename == NULL) && (output_filename == NULL) && (pgm_filename == NULL)  && (outflow_filename == NULL))) {
      (void) fprintf(stderr, 
      "\nUsage: %s [options]\n",
           argv[0]);
      (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
  }

  if( input_graphics_any_format( target_filename, &format, &n_objects, &objects ) != OK )
    return( 1 );
    /* check that the surface file contains a polyhedron */
    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS ) {
      print( "Surface file must contain 1 polygons object.\n" );
      return( 1 );
  }
  /* get a pointer to the surface */
  polygons_target = get_polygons_ptr( objects[0] );

  if( input_graphics_any_format( source_filename, &format, &n_objects, &objects ) != OK )
    return( 1 );
    /* check that the surface file contains a polyhedron */
    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS ) {
      print( "Surface file must contain 1 polygons object.\n" );
      return( 1 );
  }
  
  struct dartel_prm* prm = (struct dartel_prm*)malloc(sizeof(struct dartel_prm)*100);

  // first three entrys of param are equal
  for (j = 0; j < loop; j++)
    for (i = 0; i < 3; i++)  prm[j].rparam[i] = param[i];

  // read values from parameter file
  if (param_filename != NULL) {
    if((fp = fopen(param_filename, "r")) == 0) {
          fprintf(stderr, "Couldn't open parameter file %s.\n", param_filename);
          return(0);
    }
    
    loop = 0;

    fprintf(stderr,"Read parameters from %s\n",param_filename);
    while (fgets(line, sizeof(line), fp)) {
      if (!isalnum(line[0]))  continue;           
      // check for 9 values in each line
      if (sscanf(line,"%d %lf %lf %lf %d %d %d %d", &prm[loop].rform, &prm[loop].rparam[3], &prm[loop].rparam[4],
           &prm[loop].lmreg, &prm[loop].cycles, &prm[loop].its, &prm[loop].k, &prm[loop].code) != 8)
        continue;
      loop++;
    }
    fclose(fp);

    if (loop == 0) {
      fprintf(stderr,"Could not read parameter file %s. Check that each line contains 9 values\n", param_filename);
      return(0);
    }

  // use predefined values
  } else {
      // some entry are equal
    for (j = 0; j < loop; j++) {
      prm[j].rform = rform;
      prm[j].cycles = 3;
      prm[j].its = 3;
      prm[j].code = code;
      prm[j].lmreg = lmreg;
    }
    for (i = 0; i < 24; i++) {
      prm[i].rparam[3] = reg;  prm[i].rparam[4] = reg/2;  prm[i].k = i; reg /= 5;
    }
  }

  if (verbose) {
    fprintf(stderr,"___________________________________________________________________________\n");
    fprintf(stderr,"Parameters\n");
    fprintf(stderr,"___________________________________________________________________________\n");
    fprintf(stderr,"Regularization (0 - elastic; 1 - membrane; 2 - bending):\t\t%d\n", prm[0].rform);
    fprintf(stderr,"Number of cycles for full multi grid (FMG):\t\t\t\t%d\n", prm[0].cycles);
    fprintf(stderr,"Number of relaxation iterations in each multigrid cycle:\t\t%d\n", prm[0].its);
    fprintf(stderr,"Objective function (0 - sum of squares; 1 - sym. sum of squares):\t%d\n", prm[0].code);
    fprintf(stderr,"Levenberg-Marquardt regularization:\t\t\t\t\t%g\n", prm[0].lmreg);
    fprintf(stderr,"\n%d Iterative loops\n", loop);
    fprintf(stderr,"\nRegularization parameters 1:\t"); for (i = 0; i < loop; i++) fprintf(stderr,"%8g\t",prm[i].rparam[3]); fprintf(stderr,"\n");
    fprintf(stderr,"Regularization parameters 2:\t"); for (i = 0; i < loop; i++) fprintf(stderr,"%8g\t",prm[i].rparam[4]); fprintf(stderr,"\n");
    fprintf(stderr,"Time steps for solving the PDE:\t"); for (i = 0; i < loop; i++) fprintf(stderr,"%8d\t",prm[i].k); fprintf(stderr,"\n\n");
  }
  
  /* get a pointer to the surface */
  polygons_source = get_polygons_ptr( objects[0] );
  
  xy_size = size_map[0]*size_map[1];

  flow  = (double *)malloc(sizeof(double)*xy_size*2);
  flow1 = (double *)malloc(sizeof(double)*xy_size*2);
  inflow   = (double *)malloc(sizeof(double)*xy_size*2);
  map_source = (double *)malloc(sizeof(double)*xy_size);
  map_target = (double *)malloc(sizeof(double)*xy_size);
  map_warp   = (double *)malloc(sizeof(double)*xy_size);

  map2d_smoothed_curvature(polygons_source, map_source, fwhm, size_map);
  map2d_smoothed_curvature(polygons_target, map_target, fwhm, size_map);
               
  if (translate) {
    translate_to_template(map_source, map_target, size_map, shift);
    fprintf(stderr,"%d %d\n", shift[0], shift[1]);
  }  
  else {
    shift[0] = 0; shift[1] = 0;
  }
  
  for (i = 0; i < xy_size; i++) map_warp[i]  = 0.0;

  if (inflow_filename != NULL) {  
    if((fp_flow = fopen(inflow_filename, "rb")) == NULL)
    {
        fprintf(stderr, "Error: Couldn't read file %s.\n", inflow_filename);
        return(1);
    }
    fprintf(stderr,"Shift is not considered!!!\n");
    fread(&size_map, 2, sizeof(int), fp_flow);
    fread(&shift, 2, sizeof(int), fp_flow);
    fread(inflow, xy_size*2, sizeof(double), fp_flow);
    fclose(fp_flow);
    size_map[2] = 1;
  } else {
    for (i = 0; i < xy_size*2; i++) {
      inflow[i]   = 0.0;
    }
  }

  it = 0;
  for (it0 = 0; it0 < loop; it0++) {
    int it_scratch = dartel_scratchsize((int *)size_map, prm[it0].code);
    scratch = (double *)malloc(sizeof(double)*it_scratch);
    for (it1 = 0; it1 < prm[it0].its; it1++) {
      it++;
      dartel(size_map, prm[it0].k, inflow, map_target, map_source, (double *)0, prm[it0].rform, 
             prm[it0].rparam, prm[it0].lmreg, prm[it0].cycles, prm[it0].its, prm[it0].code, flow, ll, scratch);
      fprintf(stderr, "%02d:\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n", it, ll[0], ll[1], ll[0]+ll[1], ll[2]);
    for (i = 0; i < xy_size*2; i++) inflow[i] = flow[i];
    }
    free(scratch);
  }

  if (outflow_filename != NULL) {  
    if((fp_flow = fopen(outflow_filename, "wb")) == NULL)
    {
        fprintf(stderr, "Error: Couldn't write file %s.\n", outflow_filename);
        return(1);
    }
    fwrite(&size_map, 2, sizeof(int), fp_flow);
    fwrite(&shift, 2, sizeof(int), fp_flow);
    fwrite(inflow, xy_size*2, sizeof(double), fp_flow);
    fclose(fp_flow);
  }

  // get deformations and jacobian det. from flow field
  if (jacdet_filename != NULL) {
    jd    = (double *)malloc(sizeof(double)*xy_size);
    jd1   = (double *)malloc(sizeof(double)*xy_size);
    expdefdet(size_map, 10, inflow, flow, flow1, jd, jd1);
    if( write_pgm(jacdet_filename, jd1, size_map[0], size_map[1]) != 0 )
      return( 1 );
    free(jd);
    free(jd1);
  }
  else expdef(size_map, 10, inflow, flow, flow1, (double *)0, (double *)0);  

  free( flow1 );
  free( inflow );

  if (pgm_filename != NULL) { 
    for (i = 0; i < xy_size; i++) {
      x = (int) flow[i] - 1.0;
      y = (int) flow[i + xy_size] - 1.0;
      double xp = flow[i] - 1.0 - x;
      double yp = flow[i + xy_size] - 1.0 - y;
      double xm = 1.0 - xp;
      double ym = 1.0 - yp;
      double H00 = map_source[bound(x,  y,  size_map)];
      double H01 = map_source[bound(x,  y+1,size_map)];
      double H10 = map_source[bound(x+1,y,  size_map)];
      double H11 = map_source[bound(x+1,y+1,size_map)];

      map_warp[i] = (ym * ( xm * H00 + xp * H10) + 
		 yp * ( xm * H01 + xp * H11));
    }
    if( write_pgm(pgm_filename, map_warp, size_map[0], size_map[1]) != 0 )
      return( 1 );
  }

  if (output_filename != NULL) { 
    apply_warp(polygons_source, flow, size_map, shift);  
  
    if( output_graphics_any_format( output_filename, format, n_objects,
                objects ) != OK )
      return( 1 );
  }
  
if( write_pgm("source.pgm", map_source, size_map[0], size_map[1]) != 0 )
return( 1 );

if( write_pgm("target.pgm", map_target, size_map[0], size_map[1]) != 0 )
return( 1 );


  delete_object_list( n_objects, objects );
  free( map_source );
  free( map_warp );
  free( map_target );
  free( flow );
  free(prm);
    
  return( 0 );
}
