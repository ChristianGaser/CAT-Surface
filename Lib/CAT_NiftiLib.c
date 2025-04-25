/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_NiftiLib.h"

int
equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2) {

    if ((nii_ptr->nx != nii_ptr2->nx) ||
         (nii_ptr->ny != nii_ptr2->ny) ||
         (nii_ptr->nz != nii_ptr2->nz) ||
         (fabs(nii_ptr->dx - nii_ptr2->dx) > 0.001) ||
         (fabs(nii_ptr->dy - nii_ptr2->dy) > 0.001) ||
         (fabs(nii_ptr->dz - nii_ptr2->dz) > 0.001)) {
        fprintf(stderr,"Warning: Image %s and image %s differ.\n",nii_ptr->fname, nii_ptr2->fname);
        if (nii_ptr->nx != nii_ptr2->nx) fprintf(stderr,"nx: %d vs. %d.\n",nii_ptr->nx, nii_ptr2->nx);
        if (nii_ptr->ny != nii_ptr2->ny) fprintf(stderr,"ny: %d vs. %d.\n",nii_ptr->ny, nii_ptr2->ny);
        if (nii_ptr->nz != nii_ptr2->nz) fprintf(stderr,"nz: %d vs. %d.\n",nii_ptr->nz, nii_ptr2->nz);
        if (nii_ptr->dx != nii_ptr2->dx) fprintf(stderr,"dx: %g vs. %g.\n",nii_ptr->dx, nii_ptr2->dx);
        if (nii_ptr->dy != nii_ptr2->dy) fprintf(stderr,"dy: %g vs. %g.\n",nii_ptr->dy, nii_ptr2->dy);
        if (nii_ptr->dz != nii_ptr2->dz) fprintf(stderr,"dz: %g vs. %g.\n",nii_ptr->dz, nii_ptr2->dz);
        return(0);        
    }
    return(1);
}


/* Explicitly set all of the fields of the NIfTI I/O header structure to
 * something reasonable. Right now this is overkill since a simple memset() 
 * would do the same job, but I want this function to help me keep track
 * of all of the header fields and to allow me to easily override a default
 * if it becomes useful.
 */
void
init_nifti_header(nifti_image *nii_ptr)
{
    int i, j;
    double val;

    nii_ptr->ndim = 0;

    nii_ptr->nx = nii_ptr->ny = nii_ptr->nz = nii_ptr->nt = nii_ptr->nu = 
        nii_ptr->nv = nii_ptr->nw = 0;

    for (i = 0; i < MAX_NII_DIMS; i++) {
        /* Fix suggested by Hyun-Pil Kim (hpkim@ihanyang.ac.kr):
        Use 1 as the default, not zero */
        nii_ptr->dim[i] = 1;
    }

    nii_ptr->nvox = 0;
    nii_ptr->nbyper = 0;
    nii_ptr->datatype = DT_UNKNOWN;

    nii_ptr->dx = nii_ptr->dy = nii_ptr->dz = nii_ptr->dt = nii_ptr->du = 
        nii_ptr->dv = nii_ptr->dw = 0.0;
  
    for (i = 0; i < MAX_NII_DIMS; i++) 
        nii_ptr->pixdim[i] = 0.0;

    nii_ptr->scl_slope = 0.0;
    nii_ptr->scl_inter = 0.0;
    nii_ptr->cal_min = 0.0;
    nii_ptr->cal_max = 0.0;

    nii_ptr->qform_code = NIFTI_XFORM_UNKNOWN;
    nii_ptr->sform_code = NIFTI_XFORM_UNKNOWN;

    nii_ptr->freq_dim = 0;
    nii_ptr->phase_dim = 0;
    nii_ptr->slice_dim = 0;

    nii_ptr->slice_code = 0;
    nii_ptr->slice_start = 0;
    nii_ptr->slice_end = 0;
    nii_ptr->slice_duration = 0.0;

    nii_ptr->quatern_b = 0.0;
    nii_ptr->quatern_c = 0.0;
    nii_ptr->quatern_d = 0.0;
    nii_ptr->qoffset_x = 0.0;
    nii_ptr->qoffset_y = 0.0;
    nii_ptr->qoffset_z = 0.0;
    nii_ptr->qfac = 0.0;

    nii_ptr->toffset = 0.0;

    nii_ptr->xyz_units = NIFTI_UNITS_MM; /* Default spatial units */
    nii_ptr->time_units = NIFTI_UNITS_SEC; /* Default time units */

    nii_ptr->nifti_type = FT_ANALYZE;
    nii_ptr->intent_code = 0;
    nii_ptr->intent_p1 = 0.0;
    nii_ptr->intent_p2 = 0.0;
    nii_ptr->intent_p3 = 0.0;
    memset(nii_ptr->intent_name, 0, sizeof (nii_ptr->intent_name));

    memset(nii_ptr->descrip, 0, sizeof (nii_ptr->descrip));

    memset(nii_ptr->aux_file, 0, sizeof (nii_ptr->aux_file));
  
    nii_ptr->fname = NULL;
    nii_ptr->iname = NULL;
    nii_ptr->iname_offset = 0;
    nii_ptr->swapsize = 0;
    nii_ptr->byteorder = 1; /* Use LSB first by default. */
    nii_ptr->data = NULL;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            if (i==j) val = 1.0; else val = 0.0;
            nii_ptr->qto_xyz.m[i][j] = val;
            nii_ptr->qto_ijk.m[i][j] = val;
            nii_ptr->sto_xyz.m[i][j] = val;
            nii_ptr->sto_ijk.m[i][j] = val;
        }
    }
}  


int
write_nifti_double(const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr)
{
    nifti_image *nii_ptr, nii_rec;
    char *extension, buffer[1024];
    int i;
    double img_range[2], dt_range[2], val;
    
    if ((data_type != DT_UINT8) && 
        (data_type != DT_INT8) && 
        (data_type != DT_UINT16) && 
        (data_type != DT_INT16) && 
        (data_type != DT_UINT32) && 
        (data_type != DT_INT32) && 
        (data_type != DT_INT64) && 
        (data_type != DT_FLOAT32) && 
        (data_type != DT_FLOAT64)) {
        fprintf(stderr,"Datatype %d not supported to write data.\n",data_type);
        return(0);
    }
    
    if (!in_ptr) {
        nii_ptr = &nii_rec;
        init_nifti_header(nii_ptr);
    } else nii_ptr = in_ptr;
    
    extension = nifti_find_file_extension(output_filename);
    
    /* if no valid extension was found use .nii */
    if (!extension) {
        fprintf(stderr,"No valid extension found for output filename %s.\n",output_filename);
        return(0);
    }

    nii_ptr->nifti_type = 1;
    
    if (strcmp(extension,".img") == 0)
        nii_ptr->nifti_type = 2;
    
    if (strcmp(extension,".hdr") == 0) {
        nii_ptr->nifti_type = 2;
        strcpy(extension,".img");
    }

    output_filename = nifti_makebasename(output_filename);

    nii_ptr->nx = dim[0];
    nii_ptr->ny = dim[1];
    nii_ptr->nz = dim[2];
    nii_ptr->pixdim[0] = vox[0];
    nii_ptr->pixdim[1] = vox[1];
    nii_ptr->pixdim[2] = vox[2];
    nii_ptr->dx = vox[0];
    nii_ptr->dy = vox[1];
    nii_ptr->dz = vox[2];

    nii_ptr->data = NULL;
    nii_ptr->scl_slope = slope;
    nii_ptr->scl_inter = 0;
    nii_ptr->nvox = dim[0]*dim[1]*dim[2];
    nii_ptr->ndim = 3;

    /* for floating data no rescaling is necessary */
    if ((slope == 0.0) && ((data_type == DT_FLOAT32) || (data_type == DT_FLOAT64)))
        slope = 1.0;
         
    if (slope == 0.0) {
        /* find min and max of image */
        img_range[0] = FLT_MAX;
        img_range[1] = -FLT_MAX;
        for (i = 0; i < nii_ptr->nvox; i++) {
            img_range[0] = MIN(image[i], img_range[0]);
            img_range[1] = MAX(image[i], img_range[1]);
        }
    }

    switch (data_type) {
    case DT_UINT8:
        dt_range[0] = 0.0; dt_range[1] = 255.0;
        break;
    case DT_INT8:
        dt_range[0] = -128.0; dt_range[1] = 127.0;
        break;
    case DT_INT16:
        dt_range[0] = -32768.0; dt_range[1] = 32767.0;
        break;
    case DT_UINT16:
        dt_range[0] = 0.0; dt_range[1] = 65535.0;
        break;
    case DT_INT32:
        dt_range[0] = -2147483648.0; dt_range[1] = 2147483647.0;
        break;
    case DT_UINT32:
        dt_range[0] = 0.0; dt_range[1] = 4294967295.0;
        break;
    case DT_INT64:
        dt_range[0] = -9223372036854775808.0; dt_range[1] = 9223372036854775807.0;
        break;
    }

    /* rescale data to fit into data range */
    if (slope == 0.0) {
        nii_ptr->scl_slope = ((img_range[1] - img_range[0]) / 
                                                    (dt_range[1] - dt_range[0]));
        nii_ptr->scl_inter = img_range[0] - (dt_range[0] * nii_ptr->scl_slope);
    }

    nii_ptr->datatype = data_type;
    switch (data_type) {
    case DT_UINT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (unsigned char *)malloc(sizeof(unsigned char)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((unsigned char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (signed char *)malloc(sizeof(signed char)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((signed char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (unsigned short *)malloc(sizeof(unsigned short)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((unsigned short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (signed short *)malloc(sizeof(signed short)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((signed short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (unsigned int *)malloc(sizeof(unsigned int)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((unsigned int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (signed int *)malloc(sizeof(signed int)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((signed int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (long long *)malloc(sizeof(long long)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((long long *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_FLOAT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (float *)malloc(sizeof(float)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((float *)nii_ptr->data)[i] = val;
        }
        break;
    case DT_FLOAT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (double *)malloc(sizeof(double)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((double *)nii_ptr->data)[i] = val;
        }
        break;
    }
        
    (void) sprintf(buffer, "%s%s",output_filename,extension);

    nii_ptr->iname = NULL;
    nii_ptr->iname = malloc(strlen(buffer)+8);
    strcpy(nii_ptr->iname, buffer);

    if ((nii_ptr->nifti_type == 0) || (nii_ptr->nifti_type == 2)) {
        (void) sprintf(buffer, "%s%s",output_filename,".hdr");
    }
    nii_ptr->fname = NULL;
    nii_ptr->fname = malloc(strlen(buffer)+8);
    strcpy(nii_ptr->fname, buffer);

    nifti_image_write(nii_ptr);

    free(nii_ptr->data);
    
    return(1);
}

int
write_nifti_float(const char *output_filename, float image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr)
{
    nifti_image *nii_ptr, nii_rec;
    char *extension, buffer[1024];
    int i;
    double img_range[2], dt_range[2], val;
    
    if ((data_type != DT_UINT8) && 
        (data_type != DT_INT8) && 
        (data_type != DT_UINT16) && 
        (data_type != DT_INT16) && 
        (data_type != DT_UINT32) && 
        (data_type != DT_INT32) && 
        (data_type != DT_INT64) && 
        (data_type != DT_FLOAT32) && 
        (data_type != DT_FLOAT64)) {
        fprintf(stderr,"Datatype %d not supported to write data.\n",data_type);
        return(0);
    }
    
    if (!in_ptr) {
        nii_ptr = &nii_rec;
        init_nifti_header(nii_ptr);
    } else nii_ptr = in_ptr;
    
    extension = nifti_find_file_extension(output_filename);
    
    /* if no valid extension was found use .nii */
    if (!extension) {
        fprintf(stderr,"No valid extension found for output filename %s.\n",output_filename);
        return(0);
    }

    nii_ptr->nifti_type = 1;
    
    if (strcmp(extension,".img") == 0)
        nii_ptr->nifti_type = 2;
    
    if (strcmp(extension,".hdr") == 0) {
        nii_ptr->nifti_type = 2;
        strcpy(extension,".img");
    }
    
    output_filename = nifti_makebasename(output_filename);

    nii_ptr->nx = dim[0];
    nii_ptr->ny = dim[1];
    nii_ptr->nz = dim[2];
    nii_ptr->pixdim[0] = vox[0];
    nii_ptr->pixdim[1] = vox[1];
    nii_ptr->pixdim[2] = vox[2];
    nii_ptr->dx = vox[0];
    nii_ptr->dy = vox[1];
    nii_ptr->dz = vox[2];

    nii_ptr->data = NULL;
    nii_ptr->scl_slope = slope;
    nii_ptr->scl_inter = 0;
    nii_ptr->nvox = dim[0]*dim[1]*dim[2];
    nii_ptr->ndim = 3;

    /* for floating data no rescaling is necessary */
    if ((slope == 0.0) && ((data_type == DT_FLOAT32) || (data_type == DT_FLOAT64)))
        slope = 1.0;
         
    if (slope == 0.0) {
        /* find min and max of image */
        img_range[0] = FLT_MAX;
        img_range[1] = -FLT_MAX;
        for (i = 0; i < nii_ptr->nvox; i++) {
            img_range[0] = MIN((double)image[i], img_range[0]);
            img_range[1] = MAX((double)image[i], img_range[1]);
        }
    }

    switch (data_type) {
    case DT_UINT8:
        dt_range[0] = 0.0; dt_range[1] = 255.0;
        break;
    case DT_INT8:
        dt_range[0] = -128.0; dt_range[1] = 127.0;
        break;
    case DT_INT16:
        dt_range[0] = -32768.0; dt_range[1] = 32767.0;
        break;
    case DT_UINT16:
        dt_range[0] = 0.0; dt_range[1] = 65535.0;
        break;
    case DT_INT32:
        dt_range[0] = -2147483648.0; dt_range[1] = 2147483647.0;
        break;
    case DT_UINT32:
        dt_range[0] = 0.0; dt_range[1] = 4294967295.0;
        break;
    case DT_INT64:
        dt_range[0] = -9223372036854775808.0; dt_range[1] = 9223372036854775807.0;
        break;
    }

    /* rescale data to fit into data range */
    if (slope == 0.0) {
        nii_ptr->scl_slope = ((img_range[1] - img_range[0]) / (dt_range[1] - dt_range[0]));
        nii_ptr->scl_inter = img_range[0] - (dt_range[0] * nii_ptr->scl_slope);
    }

    nii_ptr->datatype = data_type;
    switch (data_type) {
    case DT_UINT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (unsigned char *)malloc(sizeof(unsigned char)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((unsigned char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (signed char *)malloc(sizeof(signed char)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((signed char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (unsigned short *)malloc(sizeof(unsigned short)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((unsigned short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (signed short *)malloc(sizeof(signed short)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((signed short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (unsigned int *)malloc(sizeof(unsigned int)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((unsigned int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (signed int *)malloc(sizeof(signed int)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((signed int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (long long *)malloc(sizeof(long long)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((long long *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_FLOAT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (float *)malloc(sizeof(float)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((float *)nii_ptr->data)[i] = val;
        }
        break;
    case DT_FLOAT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (double *)malloc(sizeof(double)*nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data) {
            fprintf(stderr,"Memory allocation error\n");
            return(0);
        }
        for (i = 0; i < nii_ptr->nvox; i++) {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter)/(float)nii_ptr->scl_slope;
            else val = (double)image[i];
            ((double *)nii_ptr->data)[i] = val;
        }
        break;
    }
        
    (void) sprintf(buffer, "%s%s",output_filename,extension);

    nii_ptr->iname = malloc(strlen(buffer)+8);
    strcpy(nii_ptr->iname, buffer);

    if ((nii_ptr->nifti_type == 0) || (nii_ptr->nifti_type == 2)) {
        (void) sprintf(buffer, "%s%s",output_filename,".hdr");
    }
    nii_ptr->fname = malloc(strlen(buffer)+8);
    strcpy(nii_ptr->fname, buffer);

    nifti_image_write(nii_ptr);

    free(nii_ptr->data);
    
    return(1);
}

nifti_image
*read_nifti_double(const char *input_filename, double *image[], int read_data)
{
    nifti_image *nii_ptr;
    double tmp;
    int i;
    
    nii_ptr = nifti_image_read(input_filename, 1);
    if (!nii_ptr) {
        fprintf(stderr,"read_nifti_double: Error reading %s.\n", input_filename);
        return(NULL);
    }
    
    /* read as double format */
    *image = (double *)malloc(sizeof(double)*nii_ptr->nvox*nii_ptr->nt);
    
    /* check for memory */
    if (!image) {
        fprintf(stderr,"read_nifti_double: Memory allocation error\n");
        return(NULL);
    }
    
    for (i = 0; i < nii_ptr->nvox; i++) {
            switch (nii_ptr->datatype) {
            case DT_INT8:
                tmp = (double) ((signed char *)nii_ptr->data)[i];
                break;
            case DT_UINT8:
                tmp = (double) ((unsigned char *)nii_ptr->data)[i];
                break;
            case DT_INT16:
                tmp = (double) ((signed short *)nii_ptr->data)[i];
                break;
            case DT_UINT16:
                tmp = (double) ((unsigned short *)nii_ptr->data)[i];
                break;
            case DT_INT32:
                tmp = (double) ((signed int *)nii_ptr->data)[i];
                break;
            case DT_UINT32:
                tmp = (double) ((unsigned int *)nii_ptr->data)[i];
                break;
            case DT_INT64:
                tmp = (double) ((long long *)nii_ptr->data)[i];
                break;
            case DT_FLOAT32:
                tmp = (double) ((float *)nii_ptr->data)[i];
                break;
            case DT_FLOAT64:
                tmp = (double) ((double *)nii_ptr->data)[i];
                break;
            default:
                fprintf(stderr,"read_nifti_double: Unknown datatype\n");
                return(NULL);
                break;
            }
            /* check whether scaling is needed */
            if (nii_ptr->scl_slope == 0)
                (*image)[i] = tmp;
            else
                (*image)[i] = (nii_ptr->scl_slope * tmp) + nii_ptr->scl_inter;
        }    
        
        /* ensure that nvox is that of a 3D image */
        if (nii_ptr->nt > 1) {
            nii_ptr->nvox  /= nii_ptr->nt;
            nii_ptr->dim[0] = nii_ptr->ndim = 4;
        }
        
        if (!read_data) free(nii_ptr->data);

        return(nii_ptr);
}

nifti_image
*read_nifti_float(const char *input_filename, float *image[], int read_data)
{
    nifti_image *nii_ptr;
    float tmp;
    int i;
    
    nii_ptr = nifti_image_read(input_filename, 1);
    if (!nii_ptr) {
        fprintf(stderr,"read_nifti_float: Error reading %s.\n", input_filename);
        return(NULL);
    }
    
    /* read as float format */
    *image = (float *)malloc(sizeof(float)*nii_ptr->nvox*nii_ptr->nt);
    
    /* check for memory */
    if (!image) {
        fprintf(stderr,"read_nifti_float: Memory allocation error\n");
        return(NULL);
    }
    
    for (i = 0; i < nii_ptr->nvox; i++) {
            switch (nii_ptr->datatype) {
            case DT_INT8:
                tmp = (float) ((signed char *)nii_ptr->data)[i];
                break;
            case DT_UINT8:
                tmp = (float) ((unsigned char *)nii_ptr->data)[i];
                break;
            case DT_INT16:
                tmp = (float) ((signed short *)nii_ptr->data)[i];
                break;
            case DT_UINT16:
                tmp = (float) ((unsigned short *)nii_ptr->data)[i];
                break;
            case DT_INT32:
                tmp = (float) ((signed int *)nii_ptr->data)[i];
                break;
            case DT_UINT32:
                tmp = (float) ((unsigned int *)nii_ptr->data)[i];
                break;
            case DT_INT64:
                tmp = (float) ((long long *)nii_ptr->data)[i];
                break;
            case DT_FLOAT32:
                tmp = (float) ((float *)nii_ptr->data)[i];
                break;
            case DT_FLOAT64:
                tmp = (float) ((double *)nii_ptr->data)[i];
                break;
            default:
                fprintf(stderr,"read_nifti_float: Unknown datatype\n");
                return(NULL);
                break;
            }
            /* check whether scaling is needed */
            if (nii_ptr->scl_slope == 0)
                (*image)[i] = tmp;
            else
                (*image)[i] = ((float)nii_ptr->scl_slope * tmp) + (float)nii_ptr->scl_inter;
        }    
        
        /* ensure that nvox is that of a 3D image */
        if (nii_ptr->nt > 1) {
            nii_ptr->nvox  /= nii_ptr->nt;
            nii_ptr->dim[0] = nii_ptr->ndim = 4;
        }
        
        if (!read_data) free(nii_ptr->data);
        
        return(nii_ptr);
}

