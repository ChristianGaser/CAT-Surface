/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <minc.h>
#include <time_stamp.h>
#include "CAT_SliceIO.h"


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_volume_info
@INPUT      : infile - input file name
@OUTPUT     : volume_info - input volume information
@RETURNS    : Id of icv created
@DESCRIPTION: Routine to read volume information for a file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public int get_volume_info(char *infile, Volume_Info *volume_info, int nc_type)
{
   int icvid, mincid;

   /* Create and set up icv for input */
   icvid = miicv_create();
   setup_input_icv(icvid, nc_type);

   /* Open the image file */
   mincid = miopen(infile, NC_NOWRITE);

   /* Attach the icv to the file */
   (void) miicv_attach(icvid, mincid, ncvarid(mincid, MIimage));

   /* Get dimension information */
   get_dimension_info(infile, icvid, volume_info);

   /* Get the volume min and max */
   (void) miicv_inqdbl(icvid, MI_ICV_NORM_MIN, &volume_info->minimum);
   (void) miicv_inqdbl(icvid, MI_ICV_NORM_MAX, &volume_info->maximum);

   return icvid;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : setup_input_icv
@INPUT      : icvid - id of icv to set up
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to set up an icv
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void setup_input_icv(int icvid, int nc_type)
{
   /* Set desired type */
   (void) miicv_setint(icvid, MI_ICV_TYPE, nc_type);
   if (nc_type == NC_BYTE)
       (void) miicv_setstr(icvid, MI_ICV_SIGN, MI_UNSIGNED);
   else
       (void) miicv_setstr(icvid, MI_ICV_SIGN, MI_SIGNED);   

   /* Set range of values */
   if (nc_type == NC_BYTE) {
       (void) miicv_setint(icvid, MI_ICV_VALID_MIN, 0);
       (void) miicv_setint(icvid, MI_ICV_VALID_MAX, 255);
   } else {
       (void) miicv_setint(icvid, MI_ICV_VALID_MIN, 0);
       (void) miicv_setint(icvid, MI_ICV_VALID_MAX, -1);   
   }

   /* Do normalization so that all pixels are on same scale */
   (void) miicv_setint(icvid, MI_ICV_DO_NORM, TRUE);

   /* Make sure that any out of range values are mapped to lowest value
      of type (for input only) */
   (void) miicv_setint(icvid, MI_ICV_DO_FILLVALUE, TRUE);

   /* We want to ensure that images have X, Y and Z dimensions in the
      positive direction, giving patient left on left and for drawing from
      bottom up. If we wanted patient right on left and drawing from
      top down, we would set to MI_ICV_NEGATIVE. */
   (void) miicv_setint(icvid, MI_ICV_DO_DIM_CONV, TRUE);
   (void) miicv_setint(icvid, MI_ICV_DO_SCALAR, TRUE);
   (void) miicv_setint(icvid, MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(icvid, MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(icvid, MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_dimension_info
@INPUT      : infile - name of input file
              icvid - id of the image conversion variable
@OUTPUT     : volume - input volume data
@RETURNS    : (nothing)
@DESCRIPTION: Routine to get dimension information from an input file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void get_dimension_info(char *infile, int icvid, 
                               Volume_Info *volume_info)
{
   int mincid, imgid, varid;
   int idim, ndims;
   int dim[MAX_VAR_DIMS];
   long *dimlength;
   char *dimname;
   int offset;
   int missing_one_dimension;

   /* Get the minc file id and the image variable id */
   (void) miicv_inqint(icvid, MI_ICV_CDFID, &mincid);
   (void) miicv_inqint(icvid, MI_ICV_VARID, &imgid);
   if ((mincid == MI_ERROR) || (imgid == MI_ERROR)) {
      (void) fprintf(stderr, "File %s is not attached to an icv!\n",
                     infile);
      exit(EXIT_FAILURE);
   }

   /* Get the list of dimensions subscripting the image variable */
   (void) ncvarinq(mincid, imgid, NULL, NULL, &ndims, dim, NULL);
   (void) miicv_inqint(icvid, MI_ICV_NUM_DIMS, &ndims);

   /* Check that we have two or three dimensions */
   if ((ndims != NUMBER_OF_DIMENSIONS) &&
       (ndims != NUMBER_OF_DIMENSIONS-1)) {
      (void) fprintf(stderr, 
                     "File %s does not have %d or %d dimensions\n",
                     infile, NUMBER_OF_DIMENSIONS,
                     NUMBER_OF_DIMENSIONS-1);
      exit(EXIT_FAILURE);
   }

   /* Pretend that we have three dimensions */
   offset = ndims - NUMBER_OF_DIMENSIONS;
   missing_one_dimension = (offset < 0);
   ndims = NUMBER_OF_DIMENSIONS;

   /* Loop through dimensions, checking them and getting their sizes */
   for (idim=0; idim<ndims; idim++) {

      /* Get pointers to the appropriate dimension size and name */
      switch (idim) {
      case 0: dimlength = &(volume_info->nslices) ; break;
      case 1: dimlength = &(volume_info->nrows) ; break;
      case 2: dimlength = &(volume_info->ncolumns) ; break;
      }
      dimname = volume_info->dimension_names[idim];

      /* Get dimension name and size */
      if (missing_one_dimension && (idim==0)) {
         (void) strcpy(dimname, "unknown");
         *dimlength = 1;
      }
      else {
         (void) ncdiminq(mincid, dim[idim+offset], dimname, dimlength);
      }

   }

   /* Get dimension step and start (defaults = 1 and 0). For slices,
      we read straight from the variable, but for image dimensions, we
      get the step and start from the icv. If we didn't have 
      MI_ICV_DO_DIM_CONV set to TRUE, then we would have to get the image 
      dimension step and start straight from the variables. */
   for (idim=0; idim<ndims; idim++) {
      volume_info->step[idim] = 1.0;
      volume_info->start[idim] = 0.0;
   }
   ncopts = 0;
   (void) miicv_inqdbl(icvid, MI_ICV_ADIM_STEP, &volume_info->step[COLUMN]);
   (void) miicv_inqdbl(icvid, MI_ICV_ADIM_START, &volume_info->start[COLUMN]);
   (void) miicv_inqdbl(icvid, MI_ICV_BDIM_STEP, &volume_info->step[ROW]);
   (void) miicv_inqdbl(icvid, MI_ICV_BDIM_START, &volume_info->start[ROW]);
   if ((varid=ncvarid(mincid, volume_info->dimension_names[SLICE])) 
                                != MI_ERROR) {
      (void) miattget1(mincid, varid, MIstep, NC_DOUBLE, 
                       &volume_info->step[SLICE]);
      (void) miattget1(mincid, varid, MIstart, NC_DOUBLE,
                       &volume_info->start[SLICE]);
   }
   ncopts = NC_OPTS_VAL;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : close_volume
@INPUT      : icvid - id of open icv
@OUTPUT     : (none)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to close a minc file and free the associated icv
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : March 16, 1994 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void close_volume(int icvid)
{
   int mincid;

   /* Get the minc file id and close the file */
   ncopts = 0;
   if (miicv_inqint(icvid, MI_ICV_CDFID, &mincid) != MI_ERROR) {
      (void) miclose(mincid);
   }
   ncopts = NC_OPTS_VAL;

   /* Free the icv */
   (void) miicv_free(icvid);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_volume_slice_double
@INPUT      : icvid - id of icv to read
              volume_info - info for volume
              slice_num - number of slice to read in (counting from zero)
@OUTPUT     : image - image that is read in
@RETURNS    : (nothing)
@DESCRIPTION: Routine to read in an image.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void get_volume_slice_double(int icvid, Volume_Info *volume_info, 
                             int slice_num, double *image)
{
   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int offset, ndims;

   /* Get number of dimensions */
   (void) miicv_inqint(icvid, MI_ICV_NUM_DIMS, &ndims);
   offset = ndims - NUMBER_OF_DIMENSIONS;

   /* Check slice_num */
   if (slice_num >= volume_info->nslices) {
      (void) fprintf(stderr, "Slice %d is not in the file.\n",
                     slice_num);
      exit(EXIT_FAILURE);
   }

   /* Set up the start and count variables for reading the volume */
   (void) miset_coords(3, 0, start);
   if (offset >= 0) {
      start[offset] = slice_num;
      count[offset] = 1;
   }
   count[1+offset] = volume_info->nrows;
   count[2+offset] = volume_info->ncolumns;

   /* Read in the volume */
   (void) miicv_get(icvid, start, count, image);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_volume_slice_float
@INPUT      : icvid - id of icv to read
              volume_info - info for volume
              slice_num - number of slice to read in (counting from zero)
@OUTPUT     : image - image that is read in
@RETURNS    : (nothing)
@DESCRIPTION: Routine to read in an image.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void get_volume_slice_float(int icvid, Volume_Info *volume_info, 
                             int slice_num, float *image)
{
   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int offset, ndims;

   /* Get number of dimensions */
   (void) miicv_inqint(icvid, MI_ICV_NUM_DIMS, &ndims);
   offset = ndims - NUMBER_OF_DIMENSIONS;

   /* Check slice_num */
   if (slice_num >= volume_info->nslices) {
      (void) fprintf(stderr, "Slice %d is not in the file.\n",
                     slice_num);
      exit(EXIT_FAILURE);
   }

   /* Set up the start and count variables for reading the volume */
   (void) miset_coords(3, 0, start);
   if (offset >= 0) {
      start[offset] = slice_num;
      count[offset] = 1;
   }
   count[1+offset] = volume_info->nrows;
   count[2+offset] = volume_info->ncolumns;

   /* Read in the volume */
   (void) miicv_get(icvid, start, count, image);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_volume_slice_byte
@INPUT      : icvid - id of icv to read
              volume_info - info for volume
              slice_num - number of slice to read in (counting from zero)
@OUTPUT     : image - image that is read in
@RETURNS    : (nothing)
@DESCRIPTION: Routine to read in an image.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void get_volume_slice_byte(int icvid, Volume_Info *volume_info, 
                             int slice_num, unsigned char *image)
{
   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int offset, ndims;

   /* Get number of dimensions */
   (void) miicv_inqint(icvid, MI_ICV_NUM_DIMS, &ndims);
   offset = ndims - NUMBER_OF_DIMENSIONS;

   /* Check slice_num */
   if (slice_num >= volume_info->nslices) {
      (void) fprintf(stderr, "Slice %d is not in the file.\n",
                     slice_num);
      exit(EXIT_FAILURE);
   }

   /* Set up the start and count variables for reading the volume */
   (void) miset_coords(3, 0, start);
   if (offset >= 0) {
      start[offset] = slice_num;
      count[offset] = 1;
   }
   count[1+offset] = volume_info->nrows;
   count[2+offset] = volume_info->ncolumns;

   /* Read in the volume */
   (void) miicv_get(icvid, start, count, image);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : save_volume_info
@INPUT      : input_icvid - input file icvid (MI_ERROR means no input volume)
              outfile - output file name
              arg_string - string giving argument list
              volume_info - volume information
@OUTPUT     : (nothing)
@RETURNS    : icv of output file
@DESCRIPTION: Routine to save a 3-D volume, copying information
              from an optional input file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public int save_volume_info(int input_icvid, char *outfile, char *arg_string, 
                            Volume_Info *volume_info, int nc_type)
{
   int mincid, icvid, inmincid;

   /* Create output file */
   mincid = micreate(outfile, NC_CLOBBER);

   /* Open the input file if it is provided */
   inmincid = MI_ERROR;
   if (input_icvid != MI_ERROR) {
      (void) miicv_inqint(input_icvid, MI_ICV_CDFID, &inmincid);
   }

   /* Set up variables and put output file in data mode */
   setup_variables(inmincid, mincid, volume_info, arg_string, nc_type);

   /* Create an icv and set it up */
   icvid = miicv_create();
   setup_output_icv(icvid, nc_type);

   /* Attach the icv to the file */
   (void) miicv_attach(icvid, mincid, ncvarid(mincid, MIimage));

   return icvid;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : setup_output_icv
@INPUT      : icvid - id of icv to set up
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to set up an icv
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void setup_output_icv(int icvid, int nc_type)
{
   /* Set desired type */
   (void) miicv_setint(icvid, MI_ICV_TYPE, nc_type);
   if (nc_type == NC_BYTE)
       (void) miicv_setstr(icvid, MI_ICV_SIGN, MI_UNSIGNED);
   else
       (void) miicv_setstr(icvid, MI_ICV_SIGN, MI_SIGNED);   

   if (nc_type == NC_BYTE) {
       (void) miicv_setint(icvid, MI_ICV_VALID_MIN, 0);
       (void) miicv_setint(icvid, MI_ICV_VALID_MAX, 255);
   } else {
       (void) miicv_setint(icvid, MI_ICV_VALID_MIN, 0);
       (void) miicv_setint(icvid, MI_ICV_VALID_MAX, -1);   
   }

   /* No normalization so that pixels are scaled to the slice */
   (void) miicv_setint(icvid, MI_ICV_DO_NORM, FALSE);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : setup_variables
@INPUT      : inmincid - id of input minc file (MI_ERROR if no file)
              mincid - id of output minc file
              volume_info - volume information
              arg_string - string giving argument list
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to set up variables in the output minc file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void setup_variables(int inmincid, int mincid, 
                            Volume_Info *volume_info, 
                            char *arg_string, int nc_type)
{
   int dim[MAX_VAR_DIMS], ndims, idim, varid;
   int excluded_vars[10], nexcluded;

   /* Create the dimensions */
   ndims = NUMBER_OF_DIMENSIONS;
   dim[0] = ncdimdef(mincid, volume_info->dimension_names[0], 
                     volume_info->nslices);
   dim[1] = ncdimdef(mincid, volume_info->dimension_names[1], 
                     volume_info->nrows);
   dim[2] = ncdimdef(mincid, volume_info->dimension_names[2], 
                     volume_info->ncolumns);

   /* If an input file is provided, copy all header info from that file except
      image, image-max, image-min */
   if (inmincid != MI_ERROR) {

      /* Look for the image variable and the image-max/min variables so that
         we can exclude them from the copy */
      nexcluded = 0;
      excluded_vars[nexcluded] = ncvarid(inmincid, MIimage);
      if (excluded_vars[nexcluded] != MI_ERROR) nexcluded++;
      excluded_vars[nexcluded] = ncvarid(inmincid, MIimagemax);
      if (excluded_vars[nexcluded] != MI_ERROR) nexcluded++;
      excluded_vars[nexcluded] = ncvarid(inmincid, MIimagemin);
      if (excluded_vars[nexcluded] != MI_ERROR) nexcluded++;

      /* Copy the variable definitions */
      (void) micopy_all_var_defs(inmincid, mincid, nexcluded, excluded_vars);

   }

   /* Set up the dimension variables. If the variable doesn't exist, create
      it (either no input file or variable did not exist in it). If the
      dimensions are not standard, then no variable is created. */

   for (idim=0; idim < ndims; idim++) {
      ncopts = 0;
      varid = ncvarid(mincid, volume_info->dimension_names[idim]);
      if (varid == MI_ERROR) {
         varid = micreate_std_variable(mincid, 
                                       volume_info->dimension_names[idim],
                                       NC_INT, 0, NULL);
      }
      ncopts = NC_OPTS_VAL;
      if (varid != MI_ERROR) {
         (void) miattputdbl(mincid, varid, MIstep, 
                            volume_info->step[idim]);
         (void) miattputdbl(mincid, varid, MIstart, 
                            volume_info->start[idim]);
      }
   }
   
   /* Create the image, image-max and image-min variables */
   setup_image_variables(inmincid, mincid, ndims, dim, nc_type);

   /* Add the time stamp to the history */
   update_history(mincid, arg_string);

   /* Put the file in data mode */
   (void) ncendef(mincid);

   /* Copy over variable values */
   if (inmincid != MI_ERROR) {
      (void) micopy_all_var_values(inmincid, mincid,
                                   nexcluded, excluded_vars);
   }

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : setup_image_variables
@INPUT      : inmincid - id of input minc file (MI_ERROR if no file)
              mincid - id of output minc file
              ndims - number of dimensions
              dim - list of dimension ids
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to set up image, image-max and image-min variables 
              in the output minc file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void setup_image_variables(int inmincid, int mincid, 
                                  int ndims, int dim[], int nc_type)
{
   int imgid, maxid, minid;
   double valid_range[2];

   if (nc_type == NC_BYTE) {
       valid_range[0] = 0;
       valid_range[1] = 255;
   } else if (nc_type == NC_FLOAT) {
       valid_range[0] = -1e15;
       valid_range[1] = 1e15;   
   } else if (nc_type == NC_DOUBLE) {
       valid_range[0] = -1e30;
       valid_range[1] = 1e30;   
   }
      
   /* Create the image max and min variables (varying over slices) */
   maxid = micreate_std_variable(mincid, MIimagemax, NC_DOUBLE, 1, dim);
   minid = micreate_std_variable(mincid, MIimagemin, NC_DOUBLE, 1, dim);
   if (inmincid != MI_ERROR) {
      (void) micopy_all_atts(inmincid, ncvarid(inmincid, MIimagemax),
                             mincid, maxid);
      (void) micopy_all_atts(inmincid, ncvarid(inmincid, MIimagemin),
                             mincid, minid);
   }

   /* Create the image variable, copy attributes, set the signtype attribute,
      set the valid range attribute and delete valid max/min attributes */
   imgid = micreate_std_variable(mincid, MIimage, nc_type, ndims, dim);
   if (inmincid != MI_ERROR) {
      (void) micopy_all_atts(inmincid, ncvarid(inmincid, MIimage),
                             mincid, imgid);
      ncopts = 0;
      (void) ncattdel(mincid, imgid, MIvalid_max);
      (void) ncattdel(mincid, imgid, MIvalid_min);
      ncopts = NC_OPTS_VAL;
   }
   if (nc_type == NC_BYTE)
       (void) miattputstr(mincid, imgid, MIsigntype, MI_UNSIGNED);
   else
       (void) miattputstr(mincid, imgid, MIsigntype, MI_SIGNED);
   (void) miset_valid_range(mincid, imgid, valid_range);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : update_history
@INPUT      : mincid - id of output minc file
              arg_string - string giving list of arguments
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to update the history global variable in the output 
              minc file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void update_history(int mincid, char *arg_string)
{
   nc_type datatype;
   int att_length;
   char *string;

   /* Get the history attribute length */
   ncopts=0;
   if ((ncattinq(mincid, NC_GLOBAL, MIhistory, &datatype,
                 &att_length) == MI_ERROR) ||
       (datatype != NC_CHAR))
      att_length = 0;
   att_length += strlen(arg_string) + 1;

   /* Allocate a string and get the old history */
   string = (void *) malloc(att_length);
   string[0] = '\0';
   (void) miattgetstr(mincid, NC_GLOBAL, MIhistory, att_length, 
                      string);
   ncopts = NC_OPTS_VAL;

   /* Add the new command and put the new history. */
   (void) strcat(string, arg_string);
   (void) miattputstr(mincid, NC_GLOBAL, MIhistory, string);
   free(string);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : save_volume_slice_double
@INPUT      : icvid - id of icv to write
              volume_info - volume information (minimum and maximum are
                 ignored)
              slice_num - number of slice to write
              image - image to write
              slice_min - minimum real value for slice
              slice_max - maximum real value for slice
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to write out a slice.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void save_volume_slice_double(int icvid, Volume_Info *volume_info, 
                              int slice_num, double *image)
{
   int mincid, i;
   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   double slice_min, slice_max;

   slice_min = 1e15; slice_max = -1e15;
   
   /* Get the minc file id */
   (void) miicv_inqint(icvid, MI_ICV_CDFID, &mincid);

   /* Set up the start and count variables for writinging the volume */
   (void) miset_coords(3, 0, start);
   start[0] = slice_num;
   count[0] = 1;
   count[1] = volume_info->nrows;
   count[2] = volume_info->ncolumns;
   
   for (i=0; i<(count[1]*count[2]); i++) {
       slice_min = MIN(slice_min,image[i]);
       slice_max = MAX(slice_max,image[i]);
    }
    
   /* Write out the slice min and max */
   (void) mivarput1(mincid, ncvarid(mincid, MIimagemin), start, NC_DOUBLE,
                    NULL, &slice_min);
   (void) mivarput1(mincid, ncvarid(mincid, MIimagemax), start, NC_DOUBLE,
                    NULL, &slice_max);

   /* Write out the volume */
   (void) miicv_put(icvid, start, count, image);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : save_volume_slice_float
@INPUT      : icvid - id of icv to write
              volume_info - volume information (minimum and maximum are
                 ignored)
              slice_num - number of slice to write
              image - image to write
              slice_min - minimum real value for slice
              slice_max - maximum real value for slice
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to write out a slice.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void save_volume_slice_float(int icvid, Volume_Info *volume_info, 
                              int slice_num, float *image)
{
   int mincid, i;
   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   double slice_min, slice_max;

   slice_min = 1e15; slice_max = -1e15;
   
   /* Get the minc file id */
   (void) miicv_inqint(icvid, MI_ICV_CDFID, &mincid);

   /* Set up the start and count variables for writinging the volume */
   (void) miset_coords(3, 0, start);
   start[0] = slice_num;
   count[0] = 1;
   count[1] = volume_info->nrows;
   count[2] = volume_info->ncolumns;
   
   for (i=0; i<(count[1]*count[2]); i++) {
       slice_min = MIN(slice_min,image[i]);
       slice_max = MAX(slice_max,image[i]);
    }

   /* Write out the slice min and max */
   (void) mivarput1(mincid, ncvarid(mincid, MIimagemin), start, NC_DOUBLE,
                    NULL, &slice_min);
   (void) mivarput1(mincid, ncvarid(mincid, MIimagemax), start, NC_DOUBLE,
                    NULL, &slice_max);

   /* Write out the volume */
   (void) miicv_put(icvid, start, count, image);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : save_volume_slice_byte
@INPUT      : icvid - id of icv to write
              volume_info - volume information (minimum and maximum are
                 ignored)
              slice_num - number of slice to write
              image - image to write
              slice_min - minimum real value for slice
              slice_max - maximum real value for slice
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Routine to write out a slice.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 26, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void save_volume_slice_byte(int icvid, Volume_Info *volume_info, 
                              int slice_num, unsigned char *image)
{
   int mincid, i;
   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   double slice_min, slice_max;

   slice_min = 1e15; slice_max = -1e15;
   
   /* Get the minc file id */
   (void) miicv_inqint(icvid, MI_ICV_CDFID, &mincid);

   /* Set up the start and count variables for writinging the volume */
   (void) miset_coords(3, 0, start);
   start[0] = slice_num;
   count[0] = 1;
   count[1] = volume_info->nrows;
   count[2] = volume_info->ncolumns;
   
   for (i=0; i<(count[1]*count[2]); i++) {
       slice_min = MIN(slice_min,image[i]);
       slice_max = MAX(slice_max,image[i]);
    }

   /* Write out the slice min and max */
   (void) mivarput1(mincid, ncvarid(mincid, MIimagemin), start, NC_DOUBLE,
                    NULL, &slice_min);
   (void) mivarput1(mincid, ncvarid(mincid, MIimagemax), start, NC_DOUBLE,
                    NULL, &slice_max);

   /* Write out the volume */
   (void) miicv_put(icvid, start, count, image);

}
