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
                        tmp = (double) ((char *)data)[i];
                        break;
                case DT_UINT8:
                        tmp = (double) ((unsigned char *)data)[i];
                        break;
                case DT_INT16:
                        tmp = (double) ((short *)data)[i];
                        break;
                case DT_UINT16:
                        tmp = (double) ((unsigned short *)data)[i];
                        break;
                case DT_INT32:
                        tmp = (double) ((int *)data)[i];
                        break;
                case DT_UINT32:
                        tmp = (double) ((unsigned int *)data)[i];
                        break;
                case DT_FLOAT32:
                        tmp = (double) ((float *)data)[i];
                        break;
                case DT_FLOAT64:
                        tmp = (double) ((double *)data)[i];
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
        int i, j, d;
        static int oflag = DIMORDER_XYZ;
        static int flip[MAX_SPACE_DIMS] = {0, 0, 0}; /* not used */
        void *voxels;

        mnc_vtype = NC_NAT;

        /* Read in the entire NIfTI file. */
        nii_ptr = nifti_image_read(filename, 1);

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
        
        /* set scaling and offset */
        (*volume)->real_value_scale = nii_ptr->scl_slope;
        (*volume)->real_value_translation = nii_ptr->scl_inter;
        
        if (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN)
                set_voxel_to_world_transform(*volume, &mnc_linear_xform);

        for (d = 0; d < n_dimensions; d++)
                dims[d] = nii_ptr->dim[mnc_spatial_axes[d]+1];
                
        set_volume_sizes(*volume, dims);
        alloc_volume_data(*volume);

        GET_VOXEL_PTR(voxels, *volume, 0, 0, 0, 0, 0);

        memcpy(voxels, nii_ptr->data,
               (size_t) get_volume_total_n_voxels(*volume) *
               (size_t) get_type_size(get_volume_data_type(*volume)));

        nifti_image_free(nii_ptr);

        return(status);
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
