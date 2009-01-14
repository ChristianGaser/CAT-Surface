/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

/* Definitions for accessing information on each dimension */
#define NUM_DIMENSIONS 3            /* Number of dimensions of volume */
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
   double step[NUM_DIMENSIONS];     /* Step sizes for dimensions */
   double start[NUM_DIMENSIONS];    /* Start positions for dimensions */
   char dimension_names[NUM_DIMENSIONS][MAX_NC_NAME]; /* Dimension names */
} Volume_Info;

/* Function prototypes */
int  get_volume_info(char *, Volume_Info *, int);
void setup_input_icv(int, int);
void get_dimension_info(char *, int, Volume_Info *);
void close_volume(int);
void get_volume_slice_double(int, Volume_Info *, int, double *);
void get_volume_slice_float(int, Volume_Info *, int, float *);
void get_volume_slice_byte(int, Volume_Info *, int, unsigned char *);
int  save_volume_info(int, char *, char *, Volume_Info *, int);
void setup_output_icv(int, int);
void setup_variables(int, int, Volume_Info *, char *, int);
void setup_image_variables(int, int, int, int [], int);
void update_history(int, char *);
void save_volume_slice_double(int, Volume_Info *, int, double *);
void save_volume_slice_float(int, Volume_Info *, int, float *);
void save_volume_slice_byte(int, Volume_Info *, int, unsigned char *);
