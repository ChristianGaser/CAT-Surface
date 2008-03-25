/* Constants */
#define public
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

/* Definitions for accessing information on each dimension */
#define NUMBER_OF_DIMENSIONS 3      /* Number of dimensions of volume */
#define SLICE 0                     /* Index for slice information */
#define ROW 1                       /* Index for row information */
#define COLUMN 2                    /* Index for column information */

/* Default ncopts values for error handling */
#define NC_OPTS_VAL NC_VERBOSE | NC_FATAL

#define  MIN( x, y )  ( ((x) <= (y)) ? (x) : (y) )
#define  MAX( x, y )  ( ((x) >= (y)) ? (x) : (y) )

/* Types */
typedef struct {
   long nslices;             /* Number of slices */
   long nrows;               /* Number of rows in image */
   long ncolumns;            /* Number of columns in image */
   double maximum;           /* Volume maximum */
   double minimum;           /* Volume minimum */
   double step[NUMBER_OF_DIMENSIONS];     /* Step sizes for dimensions */
   double start[NUMBER_OF_DIMENSIONS];    /* Start positions for dimensions */
   char dimension_names[NUMBER_OF_DIMENSIONS][MAX_NC_NAME];
                             /* Dimension names */
} Volume_Info;

/* Function prototypes */
public int get_volume_info(char *infile, Volume_Info *volume_info, int nc_type);
public void setup_input_icv(int icvid, int nc_type);
public void get_dimension_info(char *infile, int icvid, 
                               Volume_Info *volume_info);
public void close_volume(int icvid);
public void get_volume_slice_double(int icvid, Volume_Info *volume_info, 
                             int slice_num, double *image);
public void get_volume_slice_float(int icvid, Volume_Info *volume_info, 
                             int slice_num, float *image);
public void get_volume_slice_byte(int icvid, Volume_Info *volume_info, 
                             int slice_num, unsigned char *image);
public int save_volume_info(int input_icvid, char *outfile, char *arg_string, 
                            Volume_Info *volume_info, int nc_type);
public void setup_output_icv(int icvid, int nc_type);
public void setup_variables(int inmincid, int mincid, 
                            Volume_Info *volume_info, 
                            char *arg_string, int nc_type);
public void setup_image_variables(int inmincid, int mincid, 
                                  int ndims, int dim[], int nc_type);
public void update_history(int mincid, char *arg_string);
public void save_volume_slice_double(int icvid, Volume_Info *volume_info, 
                              int slice_num, double *image);
public void save_volume_slice_float(int icvid, Volume_Info *volume_info, 
                              int slice_num, float *image);
public void save_volume_slice_byte(int icvid, Volume_Info *volume_info, 
                              int slice_num, unsigned char *image);
