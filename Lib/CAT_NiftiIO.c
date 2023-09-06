/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Most of the functions are used from nii2mnc.c
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#include "CAT_NiftiIO.h"

/* private function from from libminc2.  This function is private partially
 because it's parameters are somewhat bizarre.  It would be a good idea to
 rework them into a more rational and easily described form. 
*/

extern void restructure_array(int ndims,
                              unsigned char *array,
                              const unsigned long *lengths_perm,
                              int el_size,
                              const int *map,
                              const int *dir);

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
                        nii_ptr->qto_xyz.m[i][j] = 0.0;
                        nii_ptr->qto_ijk.m[i][j] = 0.0;
                        nii_ptr->sto_xyz.m[i][j] = 0.0;
                        nii_ptr->sto_ijk.m[i][j] = 0.0;
                }
        }
}    

static void
find_data_range(int datatype, long nvox, void *data, double range[2])
{
        int i;
        double tmp;

        range[0] = DBL_MAX;
        range[1] = -DBL_MAX;

        for (i = 0; i < nvox; i++) {
                switch (datatype) {
                case DT_INT8:
                        tmp = (float) ((char *)data)[i];
                        break;
                case DT_UINT8:
                        tmp = (float) ((unsigned char *)data)[i];
                        break;
                case DT_INT16:
                        tmp = (float) ((short *)data)[i];
                        break;
                case DT_UINT16:
                        tmp = (float) ((unsigned short *)data)[i];
                        break;
                case DT_INT32:
                        tmp = (float) ((int *)data)[i];
                        break;
                case DT_UINT32:
                        tmp = (float) ((unsigned int *)data)[i];
                        break;
                case DT_FLOAT32:
                        tmp = (float) ((float *)data)[i];
                        break;
                case DT_FLOAT64:
                        tmp = (float) ((double *)data)[i];
                        break;
                default:
                        fprintf(stderr, "Data type %d not handled\n", datatype);
                        break;
                }
                if (tmp < range[0]) {
                        range[0] = tmp;
                }
                if (tmp > range[1]) {
                        range[1] = tmp;
                }
        }
}

/* not working yet, needs some more work */
Status
output_nifti(char *filename, nc_type volume_nc_data_type, 
                 BOOLEAN volume_signed_flag, double volume_voxel_min, 
                 double volume_voxel_max, Volume volume, char *history,
                 minc_output_options *options)
{

    /* NIFTI stuff */
    nifti_image *nii_ptr;
    nifti_image nii_rec;
    int nii_dimids[MAX_NII_DIMS];
    int nii_ndims;
    static int nifti_filetype;
    static int nifti_datatype;
    static int nifti_signed = 1;

    /* MINC stuff */
    int mnc_fd;                 /* MINC file descriptor */
    nc_type mnc_type;           /* MINC data type as read */
    int mnc_ndims;              /* MINC image dimension count */
    int mnc_dimids[MAX_VAR_DIMS]; /* MINC image dimension identifiers */
    int mnc_icv;                /* MINC image conversion variable */
    int mnc_vid;                /* MINC Image variable ID */
    long mnc_start[MAX_VAR_DIMS]; /* MINC data starts */
    long mnc_count[MAX_VAR_DIMS]; /* MINC data counts */
    int mnc_signed;             /* MINC if output voxels are signed */
    double mnc_rrange[2];       /* MINC real range (min, max) */
    double mnc_vrange[2];       /* MINC valid range (min, max) */

    /* Other stuff */
    char att_str[1024];         /* Big string for attribute values */
    int i, j;                      
    int  dims[MAX_VAR_DIMS];   /* dimensions */
    double steps[MAX_VAR_DIMS];
    char *str_ptr;              /* Generic ASCIZ string pointer */
    int r;                      /* Result code. */
    static int qflag = 0;       /* Quiet flag (default is non-quiet) */

    ncopts = 0;                 /* Clear global netCDF error reporting flag */

    /* Default NIfTI file type is "NII", single binary file
     */
    nifti_filetype = FT_UNSPECIFIED;
    nifti_datatype = DT_UNKNOWN;

    /* still experimental */
    fprintf(stderr, "Function output_nifti not yet working.\n");
    return (-1);

    if (!nifti_signed) {
        switch (nifti_datatype) {
        case DT_INT8:
            nifti_datatype = DT_UINT8;
            break;
        case DT_INT16:
            nifti_datatype = DT_UINT16;
            break;
        case DT_INT32:
            nifti_datatype = DT_UINT32;
            break;
        }
    }
    switch (nifti_datatype){
    case DT_INT8:
    case DT_UINT8:
        mnc_type = NC_BYTE;
        break;
    case DT_INT16:
    case DT_UINT16:
        mnc_type = NC_SHORT;
        break;
    case DT_INT32:
    case DT_UINT32:
        mnc_type = NC_INT;
        break;
    case DT_FLOAT32:
        mnc_type = NC_FLOAT;
        break;
    case DT_FLOAT64:
        mnc_type = NC_DOUBLE;
        break;
    }

    mnc_ndims = 3;
    mnc_dimids[0] = 0; mnc_dimids[1] = 1; mnc_dimids[2] = 2;
    nii_ndims = 3;
    get_volume_sizes(volume, dims);
    get_volume_separations(volume, steps);
/*    get_volume_real_range(volume, mnc_rrange);
    get_volume_voxel_range(volume, mnc_vrange);
*/
    /* Initialize the NIfTI structure 
     */
    nii_ptr = &nii_rec;

    init_nifti_header(nii_ptr);

    /* For now we just use the mnc2nii command line as the description
     * field.  Probably we should use something better, perhaps a
     * combination of some other standard MINC fields that might
     * provide more information.
     */
    str_ptr = nii_ptr->descrip;

    nii_ptr->fname = malloc(strlen(filename) + 4);
    nii_ptr->iname = malloc(strlen(filename) + 4);
    strcpy(nii_ptr->fname, filename);
    strcpy(nii_ptr->iname, filename);

    if (!strcmp(filename, ".nii")) {
        fprintf(stderr,"Only single file .nii allowed.");
        return(-1);
    }

    nifti_filetype = FT_NIFTI_SINGLE;

    if (mnc_vrange[1] != mnc_vrange[0] && mnc_rrange[1] != mnc_rrange[0]) {
        nii_ptr->scl_slope = ((mnc_rrange[1] - mnc_rrange[0]) / 
                              (mnc_vrange[1] - mnc_vrange[0]));
        nii_ptr->scl_inter = mnc_rrange[0] - (mnc_vrange[0] * nii_ptr->scl_slope);
    }
    else {
        nii_ptr->scl_slope = 0.0;
    }

    nii_ptr->nvox = 1;          /* Initial value for voxel count */

    /* Find all of the dimensions of the MINC file, in the order they 
     * will be listed in the NIfTI-1/Analyze file.  We use this to build
     * a map for restructuring the data according to the normal rules
     * of NIfTI-1.
     */

    for (i = 0; i < MAX_NII_DIMS; i++) {

        nii_ptr->dim[i] = dims[i];
        nii_ptr->nvox *=  dims[i];
        nii_ptr->pixdim[i] = (float) steps[i];
    }

    nii_ptr->ndim = nii_ndims; /* Total number of dimensions in file */
    nii_ptr->nx = nii_ptr->dim[0];
    nii_ptr->ny = nii_ptr->dim[1];
    nii_ptr->nz = nii_ptr->dim[2];
    nii_ptr->nt = 1;
    nii_ptr->nu = 1;

    nii_ptr->dx = nii_ptr->pixdim[0];
    nii_ptr->dy = nii_ptr->pixdim[1];
    nii_ptr->dz = nii_ptr->pixdim[2];
    nii_ptr->dt = 1;
    nii_ptr->du = 1; /* MINC files don't define a sample size for a vector_dimension */
    nii_ptr->nifti_type = nifti_filetype;

    if (nifti_datatype == DT_UNKNOWN) {
        nii_ptr->datatype = DT_FLOAT32; /* Default */
        mnc_type = NC_FLOAT;
        mnc_signed = 1;
    }
    else {
        nii_ptr->datatype = nifti_datatype;
    }


    /* Load the direction_cosines and start values into the NIfTI-1 
     * sform structure.
     *
     */
    for (i = 0; i < MAX_SPACE_DIMS; i++) {
        int id = ncvarid(mnc_fd, mnc_spatial_names[i]);
        double start;
        double step;
        double dircos[MAX_SPACE_DIMS];
        int tmp;

        if (id < 0) {
            continue;
        }

/*        start = 0.0;
        step = 1.0;
        dircos[DIM_X] = dircos[DIM_Y] = dircos[DIM_Z] = 0.0;
        dircos[i] = 1.0;

        miattget(mnc_fd, id, MIstart, NC_DOUBLE, 1, &start, &tmp);
        miattget(mnc_fd, id, MIstep, NC_DOUBLE, 1, &step, &tmp);
        miattget(mnc_fd, id, MIdirection_cosines, NC_DOUBLE, MAX_SPACE_DIMS, 
                 dircos, &tmp);
        ncdiminq(mnc_fd, ncdimid(mnc_fd, mnc_spatial_names[i]), NULL, 
                 &mnc_dlen);

        if (step < 0) {
            step = -step;
            start = start - step * (mnc_dlen - 1);
        }

        nii_ptr->sto_xyz.m[0][i] = step * dircos[0];
        nii_ptr->sto_xyz.m[1][i] = step * dircos[1];
        nii_ptr->sto_xyz.m[2][i] = step * dircos[2];

        nii_ptr->sto_xyz.m[0][3] += start * dircos[0];
        nii_ptr->sto_xyz.m[1][3] += start * dircos[1];
        nii_ptr->sto_xyz.m[2][3] += start * dircos[2];

        miattgetstr(mnc_fd, id, MIspacetype, sizeof(att_str), att_str);

        if (!strcmp(att_str, MI_NATIVE)) {
            nii_ptr->sform_code = NIFTI_XFORM_SCANNER_ANAT;
        }
        if (!strcmp(att_str, MI_TALAIRACH)) {
            nii_ptr->sform_code = NIFTI_XFORM_TALAIRACH;
        }
*/
    }

    /* So the last row is right... */
    nii_ptr->sto_xyz.m[3][0] = 0.0;
    nii_ptr->sto_xyz.m[3][1] = 0.0;
    nii_ptr->sto_xyz.m[3][2] = 0.0;
    nii_ptr->sto_xyz.m[3][3] = 1.0;

//    nii_ptr->sto_ijk = mat44_inverse(nii_ptr->sto_xyz);

    nifti_datatype_sizes(nii_ptr->datatype, 
                         &nii_ptr->nbyper, &nii_ptr->swapsize);


    if (!qflag) {
        nifti_image_infodump(nii_ptr);
    }

    /* Now load the actual MINC data. */

    nii_ptr->data = malloc(nii_ptr->nbyper * nii_ptr->nvox);
    if (nii_ptr->data == NULL) {
        fprintf(stderr, "Out of memory.\n");
        return (-1);
    }

    if (!qflag) {
        fprintf(stderr, "MINC type %d signed %d\n", mnc_type, mnc_signed);
    }

    mnc_icv = miicv_create();
    miicv_setint(mnc_icv, MI_ICV_TYPE, mnc_type);
    miicv_setstr(mnc_icv, MI_ICV_SIGN, (mnc_signed) ? MI_SIGNED : MI_UNSIGNED);
    miicv_setdbl(mnc_icv, MI_ICV_VALID_MAX, mnc_vrange[1]);
    miicv_setdbl(mnc_icv, MI_ICV_VALID_MIN, mnc_vrange[0]);
    miicv_setint(mnc_icv, MI_ICV_DO_NORM, 1);

    miicv_attach(mnc_icv, mnc_fd, mnc_vid);

    /* Read in the entire hyperslab from the file.
     */
    for (i = 0; i < mnc_ndims; i++) {
        ncdiminq(mnc_fd, mnc_dimids[i], NULL, &mnc_count[i]);
        mnc_start[i] = 0;
    }

    r = miicv_get(mnc_icv, mnc_start, mnc_count, nii_ptr->data);
    if (r < 0) {
        fprintf(stderr, "Read error\n");
        return (-1);
    }

    /* Shut down the MINC stuff now that it has done its work. 
     */
    miicv_detach(mnc_icv);
    miicv_free(mnc_icv);
    miclose(mnc_fd);

    nifti_image_write(nii_ptr);

    return (0);

}

Status
input_nifti(char *filename, int n_dimensions, char *dim_names[],
            nc_type volume_nc_data_type, BOOLEAN volume_signed_flag,
            double volume_voxel_min, double volume_voxel_max,
            BOOLEAN create_volume_flag, Volume *volume,
            minc_input_options *options)
{
        Status               status = OK;

        /* NIFTI stuff */
        nifti_image *nii_ptr;

        /* MINC stuff */
        nc_type mnc_mtype;          /* MINC memory data type */
        int mnc_msign;              /* MINC !0 if signed data */
        static nc_type mnc_vtype;   /* MINC voxel data type */
        static int mnc_vsign;       /* MINC !0 if signed data */
        int mnc_ndims;              /* MINC image dimension count */
        long mnc_start[MAX_VAR_DIMS]; /* MINC data starts */
        long mnc_count[MAX_VAR_DIMS]; /* MINC data counts */
        double mnc_vrange[2];       /* MINC valid min/max */
        double mnc_srange[2];       /* MINC image min/max */
        double mnc_time_step;
        double mnc_time_start;
        int mnc_spatial_axes[MAX_NII_DIMS], dims[MAX_VAR_DIMS];;
        double mnc_starts[MAX_SPACE_DIMS];
        double mnc_steps[MAX_SPACE_DIMS];
        double mnc_dircos[MAX_SPACE_DIMS][MAX_SPACE_DIMS];
        Transform mnc_xform;
        General_transform mnc_linear_xform;
        struct analyze75_hdr ana_hdr;
        static char *mnc_ordered_dim_names[MAX_SPACE_DIMS];

        /* Other stuff */
        int i, j, k, l, d;
        long ind;
        float tmp;
        static int oflag = DIMORDER_XYZ;
        static int flip[MAX_SPACE_DIMS] = {0, 0, 0}; /* not used */
        void *voxels;

        mnc_vtype = NC_NAT;
        
        if(n_dimensions == -1) n_dimensions = 3;

        /* Read in the entire NIfTI file. */
        nii_ptr = nifti_image_read(filename, 1);

        if(nii_ptr == NULL) {
                fprintf(stderr,"input_nifti: Error reading %s.\n", filename);
                return(-1);
        }

        if (nii_ptr->nifti_type == 0) { /* Analyze file!!! */
                FILE *fp;
                int ss, must_swap;

                fp = fopen(filename, "r");
                if (fp != NULL) {
                        fread(&ana_hdr, sizeof (ana_hdr), 1, fp);
                        fclose(fp);

                        must_swap = 0;
                        ss = ana_hdr.dime.dim[0];
                        if (ss != 0) {
                                if (ss < 0 || ss > 7) {
                                        nifti_swap_2bytes(1, &ss);
                                        if (ss < 0 || ss > 7) {
                                                /* We should never get here!! */
                                                fprintf(stderr,
                                                     "Bad dimension count!!\n");
                                        } else {
                                                must_swap = 1;
                                        }
                                }
                        } else {
                                ss = ana_hdr.hk.sizeof_hdr;
                                if (ss != sizeof(ana_hdr)) {
                                        nifti_swap_4bytes(1, &ss);
                                        if (ss != sizeof(ana_hdr)) {
                                                /* We should never get here!! */
                                                fprintf(stderr,
                                                        "Bad header size!!\n");
                                        } else {
                                                must_swap = 1;
                                        }
                                }
                        }

                        if (must_swap) {
                                nifti_swap_4bytes(1, &ana_hdr.hk.sizeof_hdr);
                                nifti_swap_4bytes(1, &ana_hdr.hk.extents);
                                nifti_swap_2bytes(1, &ana_hdr.hk.session_error);

                                nifti_swap_4bytes(1, &ana_hdr.dime.compressed);
                                nifti_swap_4bytes(1, &ana_hdr.dime.verified);
                                nifti_swap_4bytes(1, &ana_hdr.dime.glmax); 
                                nifti_swap_4bytes(1, &ana_hdr.dime.glmin);
                                nifti_swap_2bytes(8, &ana_hdr.dime.dim);
                                nifti_swap_4bytes(8, &ana_hdr.dime.pixdim);
                                nifti_swap_2bytes(1, &ana_hdr.dime.datatype);
                                nifti_swap_2bytes(1, &ana_hdr.dime.bitpix);
                                nifti_swap_4bytes(1, &ana_hdr.dime.vox_offset);
                                nifti_swap_4bytes(1, &ana_hdr.dime.cal_max); 
                                nifti_swap_4bytes(1, &ana_hdr.dime.cal_min);
                        }
                }
        }

        /* Set up the spatial axis correspondence for the call to 
         * convert_transform_to_starts_and_steps()
         */
        switch (oflag) {
        default:
        case DIMORDER_ZYX:
                mnc_ordered_dim_names[DIM_X] = MIxspace;
                mnc_ordered_dim_names[DIM_Y] = MIyspace;
                mnc_ordered_dim_names[DIM_Z] = MIzspace;
                mnc_spatial_axes[DIM_X] = DIM_X;
                mnc_spatial_axes[DIM_Y] = DIM_Y;
                mnc_spatial_axes[DIM_Z] = DIM_Z;
                break;
        case DIMORDER_ZXY:
                mnc_ordered_dim_names[DIM_X] = MIyspace;
                mnc_ordered_dim_names[DIM_Y] = MIxspace;
                mnc_ordered_dim_names[DIM_Z] = MIzspace;
                mnc_spatial_axes[DIM_X] = DIM_Y;
                mnc_spatial_axes[DIM_Y] = DIM_X;
                mnc_spatial_axes[DIM_Z] = DIM_Z;
                break;
        case DIMORDER_XYZ:
                mnc_ordered_dim_names[DIM_X] = MIzspace;
                mnc_ordered_dim_names[DIM_Y] = MIyspace;
                mnc_ordered_dim_names[DIM_Z] = MIxspace;
                mnc_spatial_axes[DIM_X] = DIM_Z;
                mnc_spatial_axes[DIM_Y] = DIM_Y;
                mnc_spatial_axes[DIM_Z] = DIM_X;
                break;
        case DIMORDER_XZY:
                mnc_ordered_dim_names[DIM_X] = MIyspace;
                mnc_ordered_dim_names[DIM_Y] = MIzspace;
                mnc_ordered_dim_names[DIM_Z] = MIxspace;
                mnc_spatial_axes[DIM_X] = DIM_Y;
                mnc_spatial_axes[DIM_Y] = DIM_Z;
                mnc_spatial_axes[DIM_Z] = DIM_X;
                break;
        case DIMORDER_YZX:
                mnc_ordered_dim_names[DIM_X] = MIxspace;
                mnc_ordered_dim_names[DIM_Y] = MIzspace;
                mnc_ordered_dim_names[DIM_Z] = MIyspace;
                mnc_spatial_axes[DIM_X] = DIM_X;
                mnc_spatial_axes[DIM_Y] = DIM_Z;
                mnc_spatial_axes[DIM_Z] = DIM_Y;
                break;
        case DIMORDER_YXZ:
                mnc_ordered_dim_names[DIM_X] = MIzspace;
                mnc_ordered_dim_names[DIM_Y] = MIxspace;
                mnc_ordered_dim_names[DIM_Z] = MIyspace;
                mnc_spatial_axes[DIM_X] = DIM_Z;
                mnc_spatial_axes[DIM_Y] = DIM_X;
                mnc_spatial_axes[DIM_Z] = DIM_Y;
                break;
        }

        switch (nii_ptr->datatype) {
        case DT_INT8:
                mnc_msign = 1;
                mnc_mtype = NC_BYTE;
                mnc_vrange[0] = CHAR_MIN;
                mnc_vrange[1] = CHAR_MAX;
                break;
        case DT_UINT8:
                mnc_msign = 0;
                mnc_mtype = NC_BYTE;
                mnc_vrange[0] = 0;
                mnc_vrange[1] = UCHAR_MAX;
                break;
        case DT_INT16:
                mnc_msign = 1;
                mnc_mtype = NC_SHORT;
                mnc_vrange[0] = SHRT_MIN;
                mnc_vrange[1] = SHRT_MAX;
                break;
        case DT_UINT16:
                mnc_msign = 0;
                mnc_mtype = NC_SHORT;
                mnc_vrange[0] = 0;
                mnc_vrange[1] = USHRT_MAX;
                break;
        case DT_INT32:
                mnc_msign = 1;
                mnc_mtype = NC_INT;
                mnc_vrange[0] = INT_MIN;
                mnc_vrange[1] = INT_MAX;
                break;
        case DT_UINT32:
                mnc_msign = 0;
                mnc_mtype = NC_INT;
                mnc_vrange[0] = 0;
                mnc_vrange[1] = UINT_MAX;
                break;
        case DT_FLOAT32:
                mnc_msign = 1;
                mnc_mtype = NC_FLOAT;
                mnc_vrange[0] = -FLT_MAX;
                mnc_vrange[1] = FLT_MAX;
                break;
        case DT_FLOAT64:
                mnc_msign = 1;
                mnc_mtype = NC_DOUBLE;
                mnc_vrange[0] = -DBL_MAX;
                mnc_vrange[1] = DBL_MAX;
                break;
        default:
                fprintf(stderr, "Data type %d not handled\n",
                                nii_ptr->datatype);
                break;
        }

        if (mnc_vtype == NC_NAT) {
                mnc_vsign = mnc_msign;
                mnc_vtype = mnc_mtype;
        }

        mnc_ndims = 0;

        if (nii_ptr->nt > 1) {

                fprintf(stderr, "ERROR: Input of 4D data is not yet supported!\n");
                return (-1);
                
                mnc_start[mnc_ndims] = 0;
                mnc_count[mnc_ndims] = nii_ptr->nt;
                mnc_ndims++;
        
                switch (nii_ptr->time_units) {
                case NIFTI_UNITS_UNKNOWN:
                case NIFTI_UNITS_SEC:
                        mnc_time_step = nii_ptr->dt;
                        mnc_time_start = nii_ptr->toffset;
                        break;
                case NIFTI_UNITS_MSEC:
                        mnc_time_step = nii_ptr->dt / 1000;
                        mnc_time_start = nii_ptr->toffset / 1000;
                        break;
                case NIFTI_UNITS_USEC:
                        mnc_time_step = nii_ptr->dt / 1000000;
                        mnc_time_step = nii_ptr->dt / 1000000;
                        mnc_time_start = nii_ptr->toffset / 1000000;
                        break;
                default:
                        fprintf(stderr, "Unknown time units value %d\n",
                                        nii_ptr->time_units);
                        break;
                }
                
                /* increase dimensions and add 4th dimension */
                n_dimensions++;
                nii_ptr->dim[mnc_spatial_axes[3]+1] = nii_ptr->nt;
        }

        if (nii_ptr->nz > 1) {
                mnc_start[mnc_ndims] = 0;
                mnc_count[mnc_ndims] = nii_ptr->nz;
                mnc_ndims++;
        }

        if (nii_ptr->ny > 1) {
                mnc_start[mnc_ndims] = 0;
                mnc_count[mnc_ndims] = nii_ptr->ny;
                mnc_ndims++;
        }

        if (nii_ptr->nx > 1) {
                mnc_start[mnc_ndims] = 0;
                mnc_count[mnc_ndims] = nii_ptr->nx;
                mnc_ndims++;
        }

        if (nii_ptr->nu > 1) {
                mnc_start[mnc_ndims] = 0;
                mnc_count[mnc_ndims] = nii_ptr->nu;
                mnc_ndims++;
        }

        /* Calculate the starts, steps, and direction cosines. This only
         * be done properly if the file is NIfTI-1 file.  If it is an Analyze
         * file we have to resort to other methods...
         */

        if (nii_ptr->nifti_type != 0 &&
            (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN ||
             nii_ptr->qform_code != NIFTI_XFORM_UNKNOWN)) {

                make_identity_transform(&mnc_xform);

                if (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN) { /* s-form */
                        for (i = 0; i < 4; i++) {
                                for (j = 0; j < 4; j++) {
                                        Transform_elem(mnc_xform, i, j) =
                                                       nii_ptr->sto_xyz.m[i][j];
                                }
                        }
                } else { /* q-xform */
                        for (i = 0; i < 4; i++) {
                                for (j = 0; j < 4; j++) {
                                        Transform_elem(mnc_xform, i, j) =
                                                       nii_ptr->qto_xyz.m[i][j];
                                }
                        }
                }

                create_linear_transform(&mnc_linear_xform, &mnc_xform);

                convert_transform_to_starts_and_steps(&mnc_linear_xform,
                                                      MAX_SPACE_DIMS, NULL,
                                                      mnc_spatial_axes,
                                                      mnc_starts, mnc_steps,
                                                      mnc_dircos);

        } else {
                /* No official transform was found (possibly this is an Analyze
                 * file).  Just use some reasonable defaults.  */
                mnc_steps[mnc_spatial_axes[DIM_X]] = nii_ptr->dx;
                mnc_steps[mnc_spatial_axes[DIM_Y]] = nii_ptr->dy;
                mnc_steps[mnc_spatial_axes[DIM_Z]] = nii_ptr->dz;
                mnc_starts[mnc_spatial_axes[DIM_X]] = -(nii_ptr->dx *
                                                        nii_ptr->nx) / 2;
                mnc_starts[mnc_spatial_axes[DIM_Y]] = -(nii_ptr->dy *
                                                        nii_ptr->ny) / 2;
                mnc_starts[mnc_spatial_axes[DIM_Z]] = -(nii_ptr->dz *
                                                        nii_ptr->nz) / 2;

                /* Unlike the starts and steps, the direction cosines do NOT
                 * change based upon the data orientation.  */
                for (i = 0; i < MAX_SPACE_DIMS; i++) {
                        for (j = 0; j < MAX_SPACE_DIMS; j++) {
                                mnc_dircos[i][j] = (i == j) ? 1.0 : 0.0;
                        }
                }
        }

        switch (nii_ptr->xyz_units) {
        case NIFTI_UNITS_METER:
                for (i = 0; i < MAX_SPACE_DIMS; i++) {
                        mnc_starts[i] *= 1000;
                        mnc_steps[i] *= 1000;
                }
                break;
        case NIFTI_UNITS_MM:
                break;
        case NIFTI_UNITS_MICRON:
                for (i = 0; i < MAX_SPACE_DIMS; i++) {
                        mnc_starts[i] /= 1000;
                        mnc_steps[i] /= 1000;
                }
                break;
        default:
                break;
        }

        /* Find the valid minimum and maximum of the data, in order to set the
         * global image minimum and image maximum properly.
         */
        find_data_range(nii_ptr->datatype, nii_ptr->nvox, nii_ptr->data,
                        mnc_vrange);

        if (nii_ptr->scl_slope != 0.0) {
                /* Convert slope/offset to min/max */
                mnc_srange[0] = (mnc_vrange[0] * nii_ptr->scl_slope) +
                                nii_ptr->scl_inter;
                mnc_srange[1] = (mnc_vrange[1] * nii_ptr->scl_slope) +
                                nii_ptr->scl_inter;
        } else {
                mnc_srange[0] = mnc_vrange[0];
                mnc_srange[1] = mnc_vrange[1];
        }

        *volume = create_volume(n_dimensions, mnc_ordered_dim_names,
                                mnc_vtype, mnc_msign,
                                mnc_vrange[0], mnc_vrange[1]);

        for (d = 0; d < n_dimensions; d++) {
                (*volume)->spatial_axes[d] = mnc_spatial_axes[d];
                set_volume_direction_cosine(*volume, d, mnc_dircos[d]);
                if (flip[d]) { /* not used */
                        mnc_starts[d] = mnc_starts[d] +
                                        ((mnc_count[d] - 1) * mnc_steps[d]);
                }
        }     

        set_volume_starts(*volume, mnc_starts);
        set_volume_real_range(*volume, mnc_srange[0], mnc_srange[1]);
        set_volume_separations(*volume, mnc_steps);
        
        /* set scaling and offset */
        (*volume)->real_value_scale = nii_ptr->scl_slope;
        (*volume)->real_value_translation = nii_ptr->scl_inter;
        
        if (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN)
                set_voxel_to_world_transform(*volume, &mnc_linear_xform);

        for (d = 0; d < n_dimensions; d++)
                dims[d] = nii_ptr->dim[mnc_spatial_axes[d]+1];
                
        set_volume_sizes(*volume, dims);
        alloc_volume_data(*volume);

        if ((0) && (! (*volume)->is_cached_volume)) { // was not working properly
                GET_VOXEL_PTR(voxels, *volume, 0, 0, 0, 0, 0);
                memcpy(voxels, nii_ptr->data,
                      (size_t) get_volume_total_n_voxels(*volume) *
                      (size_t) get_type_size(get_volume_data_type(*volume)));
        } else {
                ind = 0;
                for (i=0; i<dims[0]; i++) for (j=0; j<dims[1]; j++) for (k=0; k<dims[2]; k++) { 
                        switch (nii_ptr->datatype) {
                        case DT_INT8:
                                tmp = (float) ((char *)nii_ptr->data)[ind];
                                break;
                        case DT_UINT8:
                                tmp = (float) ((unsigned char *)nii_ptr->data)[ind];
                                break;
                        case DT_INT16:
                                tmp = (float) ((short *)nii_ptr->data)[ind];
                                break;
                        case DT_UINT16:
                                tmp = (float) ((unsigned short *)nii_ptr->data)[ind];
                                break;
                        case DT_INT32:
                                tmp = (float) ((long *)nii_ptr->data)[ind];
                                break;
                        case DT_UINT32:
                                tmp = (float) ((unsigned long *)nii_ptr->data)[ind];
                                break;
                        case DT_FLOAT32:
                                tmp = (float) ((float *)nii_ptr->data)[ind];
                                break;
                        case DT_FLOAT64:
                                tmp = (float) ((double *)nii_ptr->data)[ind];
                                break;
                        default:
                                fprintf(stderr, "Data type %d not handled\n",nii_ptr->datatype);
                                break;
                        }
                        ind++;
                        /* check whether scaling is needed */
                        if (nii_ptr->scl_slope == 0)
                                set_volume_real_value( *volume, i, j, k, 0, 0, tmp);
                        else    set_volume_real_value( *volume, i, j, k, 0, 0, (float)(nii_ptr->scl_slope * tmp) + (float)nii_ptr->scl_inter);

                }
        }

        nifti_image_free(nii_ptr);

        return(status);
}

Status
output_volume_all(char *filename, nc_type volume_nc_data_type, 
                 BOOLEAN volume_signed_flag, double volume_voxel_min, 
                 double volume_voxel_max, Volume volume, char *history,
                 minc_output_options *options)
{
        char *str_ptr;

        str_ptr = strrchr(filename, '.');

        if (!strcmp(str_ptr, ".nii") || !strcmp(str_ptr, ".hdr")) {
                if (output_nifti(filename, volume_nc_data_type, volume_signed_flag,
                                volume_voxel_min, volume_voxel_max,
                                volume, history, options) != OK)
                        return(ERROR);
        } else {
                if (output_volume(filename, volume_nc_data_type, volume_signed_flag,
                                volume_voxel_min, volume_voxel_max,
                                volume, history, options) != OK)
                        return(ERROR);
        }
        return(OK);
}

Status
input_volume_all(char *filename, int n_dimensions, char *dim_names[],
                 nc_type volume_nc_data_type, BOOLEAN volume_signed_flag,
                 double volume_voxel_min, double volume_voxel_max,
                 BOOLEAN create_volume_flag, Volume *volume,
                 minc_input_options *options)
{
        char *str_ptr;

        str_ptr = strrchr(filename, '.');

        if (!strcmp(str_ptr, ".nii") || !strcmp(str_ptr, ".hdr")) {
                if (input_nifti(filename, n_dimensions, dim_names,
                                volume_nc_data_type, volume_signed_flag,
                                volume_voxel_min, volume_voxel_max,
                                create_volume_flag, volume, options) != OK)
                        return(ERROR);
        } else {
                if (input_volume(filename, n_dimensions, dim_names,
                                 volume_nc_data_type, volume_signed_flag,
                                 volume_voxel_min, volume_voxel_max,
                                 create_volume_flag, volume, options) != OK)
                        return(ERROR);
        }
        return(OK);
}
