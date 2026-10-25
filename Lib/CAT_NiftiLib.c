/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_NiftiLib.h"

/**
 * \brief Compare image dimensions and voxel sizes of two NIfTI images.
 *
 * Checks whether two NIfTI images have identical spatial dimensions (nx, ny, nz)
 * and voxel spacings (dx, dy, dz). Prints detailed mismatch warnings if differences
 * exceed tolerance (0.001 mm). Useful for verifying spatial compatibility before
 * voxel-wise operations like subtraction or correlation.
 *
 * \param nii_ptr  (in) first NIfTI image
 * \param nii_ptr2 (in) second NIfTI image to compare
 * \return 1 if dimensions match (within tolerance), 0 if mismatched
 */
int equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2)
{

    if ((nii_ptr->nx != nii_ptr2->nx) ||
        (nii_ptr->ny != nii_ptr2->ny) ||
        (nii_ptr->nz != nii_ptr2->nz) ||
        (fabs(nii_ptr->dx - nii_ptr2->dx) > 0.001) ||
        (fabs(nii_ptr->dy - nii_ptr2->dy) > 0.001) ||
        (fabs(nii_ptr->dz - nii_ptr2->dz) > 0.001))
    {
        fprintf(stderr, "Warning: Image %s and image %s differ.\n", nii_ptr->fname, nii_ptr2->fname);
        if (nii_ptr->nx != nii_ptr2->nx)
            fprintf(stderr, "nx: %d vs. %d.\n", nii_ptr->nx, nii_ptr2->nx);
        if (nii_ptr->ny != nii_ptr2->ny)
            fprintf(stderr, "ny: %d vs. %d.\n", nii_ptr->ny, nii_ptr2->ny);
        if (nii_ptr->nz != nii_ptr2->nz)
            fprintf(stderr, "nz: %d vs. %d.\n", nii_ptr->nz, nii_ptr2->nz);
        if (nii_ptr->dx != nii_ptr2->dx)
            fprintf(stderr, "dx: %g vs. %g.\n", nii_ptr->dx, nii_ptr2->dx);
        if (nii_ptr->dy != nii_ptr2->dy)
            fprintf(stderr, "dy: %g vs. %g.\n", nii_ptr->dy, nii_ptr2->dy);
        if (nii_ptr->dz != nii_ptr2->dz)
            fprintf(stderr, "dz: %g vs. %g.\n", nii_ptr->dz, nii_ptr2->dz);
        return (0);
    }
    return (1);
}

/**
 * \brief Initialize NIfTI image header structure to reasonable defaults.
 *
 * Sets all NIfTI header fields to standard values for a 3D volume image.
 * Ensures consistent header state before writing or modified after reading.
 * This function provides explicit control over all header parameters for
 * clarity and maintainability.\n *
 * \param nii_ptr (in/out) NIfTI image structure to initialize
 */
void init_nifti_header(nifti_image *nii_ptr)
{
    int i, j;
    double val;

    /* the writers pass a stack record: fields not set below, such as the
       extension list nifti_image_write() walks, must not be left as garbage */
    memset(nii_ptr, 0, sizeof(*nii_ptr));

    nii_ptr->ndim = 0;

    nii_ptr->nx = nii_ptr->ny = nii_ptr->nz = nii_ptr->nt = nii_ptr->nu =
        nii_ptr->nv = nii_ptr->nw = 0;

    for (i = 0; i < MAX_NII_DIMS; i++)
    {
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

    nii_ptr->xyz_units = NIFTI_UNITS_MM;   /* Default spatial units */
    nii_ptr->time_units = NIFTI_UNITS_SEC; /* Default time units */

    nii_ptr->nifti_type = FT_ANALYZE;
    nii_ptr->intent_code = 0;
    nii_ptr->intent_p1 = 0.0;
    nii_ptr->intent_p2 = 0.0;
    nii_ptr->intent_p3 = 0.0;
    memset(nii_ptr->intent_name, 0, sizeof(nii_ptr->intent_name));

    memset(nii_ptr->descrip, 0, sizeof(nii_ptr->descrip));

    memset(nii_ptr->aux_file, 0, sizeof(nii_ptr->aux_file));

    nii_ptr->fname = NULL;
    nii_ptr->iname = NULL;
    nii_ptr->iname_offset = 0;
    nii_ptr->swapsize = 0;
    nii_ptr->byteorder = 1; /* Use LSB first by default. */
    nii_ptr->data = NULL;

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if (i == j)
                val = 1.0;
            else
                val = 0.0;
            nii_ptr->qto_xyz.m[i][j] = val;
            nii_ptr->qto_ijk.m[i][j] = val;
            nii_ptr->sto_xyz.m[i][j] = val;
            nii_ptr->sto_ijk.m[i][j] = val;
        }
    }
}

/**
 * \brief Make sto_xyz the voxel-to-world matrix FSL and nibabel use.
 *
 * libCAT maps voxels to world coordinates with sto_xyz, but nifti1_io fills it
 * only from an sform: for a file that has just a qform it stays zero, and
 * every registration or surface position computed from it is garbage. This
 * applies the rule of FSL and nibabel: the sform when sform_code > 0,
 * otherwise the qform (which nifti1_io sets from pixdim when qform_code is 0
 * as well). sform_code is left alone, so a header written back still stores
 * only the transforms the file had.
 *
 * \param nii_ptr (in/out) header whose sto_xyz and sto_ijk are set
 */
void set_world_matrix(nifti_image *nii_ptr)
{
    if (nii_ptr->sform_code > 0)
        return;
    nii_ptr->sto_xyz = nii_ptr->qto_xyz;
    nii_ptr->sto_ijk = nii_ptr->qto_ijk;
}

/* C = A * B for 4x4 row-major matrices; C may not alias A or B */
static void mat4_mul(const double A[16], const double B[16], double C[16])
{
    int r, c, k;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
        {
            double s = 0.0;
            for (k = 0; k < 4; k++)
                s += A[4 * r + k] * B[4 * k + c];
            C[4 * r + c] = s;
        }
}

/* inverse of a 4x4 row-major affine matrix (last row 0 0 0 1) */
static void affine_inverse(const double A[16], double inv[16])
{
    double a = A[0], b = A[1], c = A[2], d = A[4], e = A[5], f = A[6];
    double g = A[8], h = A[9], k = A[10];
    double det = a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
    double R[9];
    int r;

    R[0] = (e * k - f * h) / det;
    R[1] = (c * h - b * k) / det;
    R[2] = (b * f - c * e) / det;
    R[3] = (f * g - d * k) / det;
    R[4] = (a * k - c * g) / det;
    R[5] = (c * d - a * f) / det;
    R[6] = (d * h - e * g) / det;
    R[7] = (b * g - a * h) / det;
    R[8] = (a * e - b * d) / det;
    for (r = 0; r < 3; r++)
    {
        inv[4 * r] = R[3 * r];
        inv[4 * r + 1] = R[3 * r + 1];
        inv[4 * r + 2] = R[3 * r + 2];
        inv[4 * r + 3] = -(R[3 * r] * A[3] + R[3 * r + 1] * A[7] +
                           R[3 * r + 2] * A[11]);
    }
    inv[12] = inv[13] = inv[14] = 0.0;
    inv[15] = 1.0;
}

/* sto_xyz of a header as a row-major double matrix */
static void vox2world_matrix(const nifti_image *nii_ptr, double V[16])
{
    int r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            V[4 * r + c] = (double)nii_ptr->sto_xyz.m[r][c];
}

/**
 * \brief Matrix from voxel indices to FSL's scaled voxel coordinates.
 *
 * FSL states transforms in voxel index times voxel size (pixdim), with the x
 * axis reversed for an image whose voxel-to-world matrix has a positive
 * determinant (neurological storage), so that FSL coordinates always run
 * radiologically. A world shift to the right is a negative FSL x either way.
 *
 * \param nii_ptr (in)  image header (sto_xyz, nx, dx, dy, dz)
 * \param M       (out) 4x4 row-major matrix, FSL coordinate = M * voxel index
 */
void vox2fsl_matrix(const nifti_image *nii_ptr, double M[16])
{
    double V[16], det;
    int i;

    vox2world_matrix(nii_ptr, V);
    det = V[0] * (V[5] * V[10] - V[6] * V[9]) - V[1] * (V[4] * V[10] - V[6] * V[8]) +
          V[2] * (V[4] * V[9] - V[5] * V[8]);
    for (i = 0; i < 16; i++)
        M[i] = (i % 5 == 0) ? 1.0 : 0.0;
    M[0] = fabs(nii_ptr->dx);
    M[5] = fabs(nii_ptr->dy);
    M[10] = fabs(nii_ptr->dz);
    if (det > 0.0)
    {
        M[0] = -M[0];
        M[3] = (nii_ptr->nx - 1) * fabs(nii_ptr->dx);
    }
}

/**
 * \brief Convert a world-space transform between two images to FSL's convention.
 *
 * fsl = F_ref * V_ref^-1 * ras * V_src * F_src^-1, with V the voxel-to-world
 * matrix (sto_xyz) and F the vox2fsl_matrix() of an image.
 *
 * \param ras (in)  4x4 row-major world transform, source -> reference
 * \param src (in)  header of the source image
 * \param ref (in)  header of the reference image
 * \param fsl (out) 4x4 row-major FSL transform, source -> reference
 */
void ras_to_fsl_matrix(const double ras[16], const nifti_image *src,
                       const nifti_image *ref, double fsl[16])
{
    double Vs[16], Fs[16], Fs_inv[16], Vr[16], Vr_inv[16], Fr[16];
    double t1[16], t2[16];

    vox2world_matrix(src, Vs);
    vox2fsl_matrix(src, Fs);
    affine_inverse(Fs, Fs_inv);
    vox2world_matrix(ref, Vr);
    affine_inverse(Vr, Vr_inv);
    vox2fsl_matrix(ref, Fr);

    mat4_mul(Vs, Fs_inv, t1);    /* FSL source -> world source */
    mat4_mul(ras, t1, t2);       /* -> world reference */
    mat4_mul(Vr_inv, t2, t1);    /* -> voxel reference */
    mat4_mul(Fr, t1, fsl);       /* -> FSL reference */
}

/**
 * \brief Convert an FSL transform between two images to world space.
 *
 * ras = V_ref * F_ref^-1 * fsl * F_src * V_src^-1, the inverse of
 * ras_to_fsl_matrix().
 *
 * \param fsl (in)  4x4 row-major FSL transform, source -> reference
 * \param src (in)  header of the source image
 * \param ref (in)  header of the reference image
 * \param ras (out) 4x4 row-major world transform, source -> reference
 */
void fsl_to_ras_matrix(const double fsl[16], const nifti_image *src,
                       const nifti_image *ref, double ras[16])
{
    double Vs[16], Vs_inv[16], Fs[16], Vr[16], Fr[16], Fr_inv[16];
    double t1[16], t2[16];

    vox2world_matrix(src, Vs);
    affine_inverse(Vs, Vs_inv);
    vox2fsl_matrix(src, Fs);
    vox2world_matrix(ref, Vr);
    vox2fsl_matrix(ref, Fr);
    affine_inverse(Fr, Fr_inv);

    mat4_mul(Fs, Vs_inv, t1);    /* world source -> FSL source */
    mat4_mul(fsl, t1, t2);       /* -> FSL reference */
    mat4_mul(Fr_inv, t2, t1);    /* -> voxel reference */
    mat4_mul(Vr, t1, ras);       /* -> world reference */
}

/**
 * \brief Write a 4x4 matrix as FSL writes it (flirt -omat, mcflirt -mats).
 *
 * \param filename (in) output file
 * \param m        (in) 4x4 row-major matrix
 * \return 0 on success, -1 if the file cannot be written
 */
int write_fsl_matrix(const char *filename, const double m[16])
{
    FILE *fp = fopen(filename, "w");
    int r;

    if (!fp)
    {
        fprintf(stderr, "write_fsl_matrix: cannot open '%s'\n", filename);
        return -1;
    }
    for (r = 0; r < 4; r++)
        fprintf(fp, "%f %f %f %f \n", m[4 * r], m[4 * r + 1], m[4 * r + 2],
                m[4 * r + 3]);
    return fclose(fp) == 0 ? 0 : -1;
}

/**
 * \brief Make a NIfTI header describe a single 3D volume.
 *
 * The writers always store nx*ny*nz voxels. A template header taken from a 4D
 * file still carries its time dimension, and written as it is the file would
 * claim nt volumes while holding one -- which readers such as FSL take at their
 * word (fslinfo reports dim4 = nt). Sets ndim to 3 and dimensions 4 to 7 to 1.
 *
 * \param nii_ptr (in/out) header of the volume about to be written; nx, ny and nz must be set
 */
static void set_single_volume(nifti_image *nii_ptr)
{
    int i;

    nii_ptr->ndim = 3;
    nii_ptr->nt = nii_ptr->nu = nii_ptr->nv = nii_ptr->nw = 1;
    nii_ptr->dim[0] = 3;
    nii_ptr->dim[1] = nii_ptr->nx;
    nii_ptr->dim[2] = nii_ptr->ny;
    nii_ptr->dim[3] = nii_ptr->nz;
    for (i = 4; i < MAX_NII_DIMS; i++)
        nii_ptr->dim[i] = 1;
}

/**
 * \brief Make a NIfTI header describe nt 3D volumes of a time series.
 *
 * Counterpart of set_single_volume() for the 4D writer. The time step of a 4D
 * template header is kept; a header without one gets 1.
 *
 * \param nii_ptr (in/out) header of the series about to be written; nx, ny and nz must be set
 * \param nt      (in)     number of volumes, at least 2
 */
static void set_volumes(nifti_image *nii_ptr, int nt)
{
    int i;

    nii_ptr->ndim = 4;
    nii_ptr->nt = nt;
    nii_ptr->nu = nii_ptr->nv = nii_ptr->nw = 1;
    nii_ptr->dim[0] = 4;
    nii_ptr->dim[1] = nii_ptr->nx;
    nii_ptr->dim[2] = nii_ptr->ny;
    nii_ptr->dim[3] = nii_ptr->nz;
    nii_ptr->dim[4] = nt;
    for (i = 5; i < MAX_NII_DIMS; i++)
        nii_ptr->dim[i] = 1;
    if (!(nii_ptr->dt > 0.0))
        nii_ptr->dt = 1.0;
    nii_ptr->pixdim[4] = nii_ptr->dt;
}

/**
 * \brief Set the image and header file names of an image about to be written.
 *
 * The names are allocated to fit, so an output path of any length works; a
 * fixed 1024-byte buffer used to overflow on longer paths.
 *
 * \param nim       (in/out) image whose iname and fname are set
 * \param basename  (in)     output file name without its extension
 * \param extension (in)     extension of the image file (".nii", ".img", ...)
 */
static void set_output_filenames(nifti_image *nim, const char *basename,
                                 const char *extension)
{
    const char *hdr_extension = extension;
    size_t len;

    if ((nim->nifti_type == 0) || (nim->nifti_type == 2))
        hdr_extension = ".hdr";

    len = strlen(basename) + strlen(extension) + 1;
    nim->iname = malloc(len);
    (void)snprintf(nim->iname, len, "%s%s", basename, extension);
    len = strlen(basename) + strlen(hdr_extension) + 1;
    nim->fname = malloc(len);
    (void)snprintf(nim->fname, len, "%s%s", basename, hdr_extension);
}

/**
 * \brief Write double-precision volume data to NIfTI file with type conversion.
 *
 * Exports a 3D volume array (stored as linear double array) to NIfTI format.
 * Converts from double precision to specified output datatype (int8, int16, int32,
 * float32, float64, etc.) with optional slope-based scaling. Automatically scales
 * data range to fit target datatype. If no input header provided, initializes
 * default header with given dimensions and voxel sizes. The file always holds a
 * single 3D volume, also when the template header comes from a 4D file.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) linear array of voxel values (size dim[0]*dim[1]*dim[2])
 * \param data_type       (in) target NIfTI datatype (DT_UINT8, DT_INT16, DT_FLOAT32, etc.)
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim            (in) array of 3 image dimensions [nx, ny, nz]
 * \param vox            (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr         (in) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_double(const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr)
{
    nifti_image *nii_ptr, nii_rec;
    char *extension;
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
        (data_type != DT_FLOAT64))
    {
        fprintf(stderr, "Datatype %d not supported to write data.\n", data_type);
        return (0);
    }

    if (!in_ptr)
    {
        nii_ptr = &nii_rec;
        init_nifti_header(nii_ptr);
    }
    else
        nii_ptr = in_ptr;

    extension = nifti_find_file_extension(output_filename);

    /* if no valid extension was found use .nii */
    if (!extension)
    {
        fprintf(stderr, "No valid extension found for output filename %s.\n", output_filename);
        return (0);
    }

    nii_ptr->nifti_type = 1;

    if (strcmp(extension, ".img") == 0)
        nii_ptr->nifti_type = 2;

    if (strcmp(extension, ".hdr") == 0)
    {
        nii_ptr->nifti_type = 2;
        strcpy(extension, ".img");
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
    nii_ptr->nvox = dim[0] * dim[1] * dim[2];
    set_single_volume(nii_ptr);

    /* for floating data no rescaling is necessary */
    if ((slope == 0.0) && ((data_type == DT_FLOAT32) || (data_type == DT_FLOAT64)))
        slope = 1.0;

    if (slope == 0.0)
    {
        /* find min and max of image */
        img_range[0] = FLT_MAX;
        img_range[1] = -FLT_MAX;
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            img_range[0] = MIN(image[i], img_range[0]);
            img_range[1] = MAX(image[i], img_range[1]);
        }
    }

    switch (data_type)
    {
    case DT_UINT8:
        dt_range[0] = 0.0;
        dt_range[1] = 255.0;
        break;
    case DT_INT8:
        dt_range[0] = -128.0;
        dt_range[1] = 127.0;
        break;
    case DT_INT16:
        dt_range[0] = -32768.0;
        dt_range[1] = 32767.0;
        break;
    case DT_UINT16:
        dt_range[0] = 0.0;
        dt_range[1] = 65535.0;
        break;
    case DT_INT32:
        dt_range[0] = -2147483648.0;
        dt_range[1] = 2147483647.0;
        break;
    case DT_UINT32:
        dt_range[0] = 0.0;
        dt_range[1] = 4294967295.0;
        break;
    case DT_INT64:
        dt_range[0] = -9223372036854775808.0;
        dt_range[1] = 9223372036854775807.0;
        break;
    }

    /* rescale data to fit into data range */
    if (slope == 0.0)
    {
        nii_ptr->scl_slope = ((img_range[1] - img_range[0]) /
                              (dt_range[1] - dt_range[0]));
        nii_ptr->scl_inter = img_range[0] - (dt_range[0] * nii_ptr->scl_slope);
    }

    nii_ptr->datatype = data_type;
    switch (data_type)
    {
    case DT_UINT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (unsigned char *)malloc(sizeof(unsigned char) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((unsigned char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (signed char *)malloc(sizeof(signed char) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((signed char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (unsigned short *)malloc(sizeof(unsigned short) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((unsigned short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (signed short *)malloc(sizeof(signed short) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((signed short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (unsigned int *)malloc(sizeof(unsigned int) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((unsigned int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (signed int *)malloc(sizeof(signed int) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((signed int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (long long *)malloc(sizeof(long long) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((long long *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_FLOAT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (float *)malloc(sizeof(float) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((float *)nii_ptr->data)[i] = val;
        }
        break;
    case DT_FLOAT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (double *)malloc(sizeof(double) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((double *)nii_ptr->data)[i] = val;
        }
        break;
    }

    set_output_filenames(nii_ptr, output_filename, extension);

    nifti_image_write(nii_ptr);

    free(nii_ptr->data);
    nii_ptr->data = NULL;

    return (1);
}

/**
 * \brief Write nt float volumes to a NIfTI file with type conversion.
 *
 * Shared body of write_nifti_float() (nt = 1) and write_nifti_float_4d(). The
 * volumes follow one another in image, as the readers return them. Float32
 * output is written straight from image instead of from a converted copy: the
 * values would be identical, and a copy doubles the memory a 4D series needs.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) nt volumes of dim[0]*dim[1]*dim[2] float values
 * \param data_type       (in) target NIfTI datatype
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim             (in) array of 3 image dimensions [nx, ny, nz]
 * \param nt              (in) number of volumes; 1 writes a 3D file
 * \param vox             (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr          (in) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
static int write_nifti_float_volumes(const char *output_filename, float image[],
                                     int data_type, double slope, int dim[], int nt,
                                     double vox[], nifti_image *in_ptr)
{
    nifti_image *nii_ptr, nii_rec;
    char *extension;
    size_t i;
    int borrowed = 0;
    double img_range[2], dt_range[2], val;

    if ((data_type != DT_UINT8) &&
        (data_type != DT_INT8) &&
        (data_type != DT_UINT16) &&
        (data_type != DT_INT16) &&
        (data_type != DT_UINT32) &&
        (data_type != DT_INT32) &&
        (data_type != DT_INT64) &&
        (data_type != DT_FLOAT32) &&
        (data_type != DT_FLOAT64))
    {
        fprintf(stderr, "Datatype %d not supported to write data.\n", data_type);
        return (0);
    }

    if (!in_ptr)
    {
        nii_ptr = &nii_rec;
        init_nifti_header(nii_ptr);
    }
    else
        nii_ptr = in_ptr;

    extension = nifti_find_file_extension(output_filename);

    /* if no valid extension was found use .nii */
    if (!extension)
    {
        fprintf(stderr, "No valid extension found for output filename %s.\n", output_filename);
        return (0);
    }

    nii_ptr->nifti_type = 1;

    if (strcmp(extension, ".img") == 0)
        nii_ptr->nifti_type = 2;

    if (strcmp(extension, ".hdr") == 0)
    {
        nii_ptr->nifti_type = 2;
        strcpy(extension, ".img");
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
    nii_ptr->nvox = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2] * (size_t)nt;
    if (nt > 1)
        set_volumes(nii_ptr, nt);
    else
        set_single_volume(nii_ptr);

    /* for floating data no rescaling is necessary */
    if ((slope == 0.0) && ((data_type == DT_FLOAT32) || (data_type == DT_FLOAT64)))
        slope = 1.0;

    if (slope == 0.0)
    {
        /* find min and max of image */
        img_range[0] = FLT_MAX;
        img_range[1] = -FLT_MAX;
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            img_range[0] = MIN((double)image[i], img_range[0]);
            img_range[1] = MAX((double)image[i], img_range[1]);
        }
    }

    switch (data_type)
    {
    case DT_UINT8:
        dt_range[0] = 0.0;
        dt_range[1] = 255.0;
        break;
    case DT_INT8:
        dt_range[0] = -128.0;
        dt_range[1] = 127.0;
        break;
    case DT_INT16:
        dt_range[0] = -32768.0;
        dt_range[1] = 32767.0;
        break;
    case DT_UINT16:
        dt_range[0] = 0.0;
        dt_range[1] = 65535.0;
        break;
    case DT_INT32:
        dt_range[0] = -2147483648.0;
        dt_range[1] = 2147483647.0;
        break;
    case DT_UINT32:
        dt_range[0] = 0.0;
        dt_range[1] = 4294967295.0;
        break;
    case DT_INT64:
        dt_range[0] = -9223372036854775808.0;
        dt_range[1] = 9223372036854775807.0;
        break;
    }

    /* rescale data to fit into data range */
    if (slope == 0.0)
    {
        nii_ptr->scl_slope = ((img_range[1] - img_range[0]) / (dt_range[1] - dt_range[0]));
        nii_ptr->scl_inter = img_range[0] - (dt_range[0] * nii_ptr->scl_slope);
    }

    nii_ptr->datatype = data_type;
    switch (data_type)
    {
    case DT_UINT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (unsigned char *)malloc(sizeof(unsigned char) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((unsigned char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT8:
        nii_ptr->nbyper = 1;
        nii_ptr->data = (signed char *)malloc(sizeof(signed char) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((signed char *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (unsigned short *)malloc(sizeof(unsigned short) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((unsigned short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT16:
        nii_ptr->nbyper = 2;
        nii_ptr->data = (signed short *)malloc(sizeof(signed short) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((signed short *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_UINT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (unsigned int *)malloc(sizeof(unsigned int) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((unsigned int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT32:
        nii_ptr->nbyper = 4;
        nii_ptr->data = (signed int *)malloc(sizeof(signed int) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((signed int *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_INT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (long long *)malloc(sizeof(long long) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((long long *)nii_ptr->data)[i] = round(val);
        }
        break;
    case DT_FLOAT32:
        /* slope is never 0 here (reset above), so every value is image[i]
           itself: write the caller's buffer instead of a copy of it */
        nii_ptr->nbyper = 4;
        nii_ptr->data = image;
        borrowed = 1;
        break;
    case DT_FLOAT64:
        nii_ptr->nbyper = 8;
        nii_ptr->data = (double *)malloc(sizeof(double) * nii_ptr->nvox);
        /* check for memory */
        if (!nii_ptr->data)
        {
            fprintf(stderr, "Memory allocation error\n");
            return (0);
        }
        for (i = 0; i < nii_ptr->nvox; i++)
        {
            if (slope == 0.0)
                val = (double)(image[i] - (float)nii_ptr->scl_inter) / (float)nii_ptr->scl_slope;
            else
                val = (double)image[i];
            ((double *)nii_ptr->data)[i] = val;
        }
        break;
    }

    set_output_filenames(nii_ptr, output_filename, extension);

    nifti_image_write(nii_ptr);

    if (!borrowed)
        free(nii_ptr->data);
    nii_ptr->data = NULL;

    return (1);
}

/**
 * \brief Write single-precision volume data to NIfTI file with type conversion.
 *
 * Exports a 3D volume array (stored as linear float array) to NIfTI format.
 * Converts from float to specified output datatype. Similar to write_nifti_double
 * but accepts input as single-precision floats. Useful when source data is already
 * float and double conversion overhead should be avoided. The file always holds a
 * single 3D volume, also when the template header comes from a 4D file.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) linear array of float voxel values
 * \param data_type       (in) target NIfTI datatype
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim            (in) array of 3 image dimensions [nx, ny, nz]
 * \param vox            (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr         (in) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_float(const char *output_filename, float image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr)
{
    return write_nifti_float_volumes(output_filename, image, data_type, slope,
                                     dim, 1, vox, in_ptr);
}

/**
 * \brief Write a single-precision 4D time series to a NIfTI file.
 *
 * Like write_nifti_float(), but the file holds nt volumes, stored one after
 * the other in image as read_nifti_float() returns them. The header describes
 * a 4D image (dim[4] = nt); the time step of a 4D template is kept, otherwise
 * it is set to 1. With nt = 1 this is write_nifti_float().
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) nt volumes of dim[0]*dim[1]*dim[2] float values
 * \param data_type       (in) target NIfTI datatype
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim             (in) array of 3 image dimensions [nx, ny, nz]
 * \param nt              (in) number of volumes (at least 1)
 * \param vox             (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr          (in) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_float_4d(const char *output_filename, float image[], int data_type,
                         double slope, int dim[], int nt, double vox[],
                         nifti_image *in_ptr)
{
    if (nt < 1)
    {
        fprintf(stderr, "write_nifti_float_4d: invalid number of volumes %d.\n", nt);
        return (0);
    }
    return write_nifti_float_volumes(output_filename, image, data_type, slope,
                                     dim, nt, vox, in_ptr);
}

/* nifti_image_read() and nifti_image_load() find the image file by
   stripping the extension and trying .nii before .nii.gz, so x.nii.gz with
   an x.nii next to it came back with the header of the one and the voxels
   of the other (upstream nifti1_io does the same).  A single-file NIfTI
   keeps its voxels in the file the header was read from, so read them from
   there; .hdr/.img pairs keep the library's own search. */
static nifti_image *
nifti_image_read_exact(const char *filename)
{
    nifti_image *nim = nifti_image_read(filename, 0);
    znzFile fp;
    size_t ntot;
    char *same;

    if (!nim)
        return NULL;
    if (nim->nifti_type != NIFTI_FTYPE_NIFTI1_1 || nim->fname == NULL ||
        nim->iname_offset < 0)
    {
        if (nifti_image_load(nim) < 0)
        {
            nifti_image_free(nim);
            return NULL;
        }
        return nim;
    }

    /* keep the header consistent with the file the voxels come from */
    same = (char *)malloc(strlen(nim->fname) + 1);
    if (!same)
    {
        nifti_image_free(nim);
        return NULL;
    }
    strcpy(same, nim->fname);
    free(nim->iname);
    nim->iname = same;

    /* the steps of nifti_image_load(), on the file that was asked for */
    ntot = nifti_get_volsize(nim);
    fp = znzopen(nim->fname, "rb", nifti_is_gzfile(nim->fname));
    if (znz_isnull(fp))
    {
        nifti_image_free(nim);
        return NULL;
    }
    nim->data = calloc(1, ntot);
    if (!nim->data || znzseek(fp, (long)nim->iname_offset, SEEK_SET) < 0 ||
        nifti_read_buffer(fp, nim->data, ntot, nim) < ntot)
    {
        znzclose(fp);
        nifti_image_free(nim);
        return NULL;
    }
    znzclose(fp);
    return nim;
}

/**
 * \brief Read NIfTI image file and load data into double-precision array.
 *
 * Parses NIfTI (.nii or .nii.gz) format files and converts voxels to double
 * precision. Handles all NIfTI integer and floating-point datatypes. Automatically
 * allocates output array and applies datatype scaling/slope corrections. Non-finite
 * values of float images are read as 0 (the NIfTI library replaces NaN and Inf).
 * A single-file NIfTI is read from the file named, also when one with the other
 * extension (x.nii next to x.nii.gz) exists.
 *
 * A 4D file is returned as its nt volumes one after the other: volume v starts at
 * image[v * nvox], where nvox of the returned header is the voxel count of ONE
 * volume (nx*ny*nz) and nt the number of volumes.
 *
 * \param input_filename (in)  path to NIfTI file
 * \param image          (out) pointer to allocated double array (size nx*ny*nz*nt)
 * \param read_data      (in)  if non-zero, keep the raw NIfTI data in the returned
 *                             header; if 0, release it (image is returned either way)
 * \return pointer to NIfTI image structure on success; NULL on error
 */
nifti_image *read_nifti_double(const char *input_filename, double *image[], int read_data)
{
    nifti_image *nii_ptr;
    double tmp;
    size_t i;

    nii_ptr = nifti_image_read_exact(input_filename);
    if (!nii_ptr)
    {
        fprintf(stderr, "read_nifti_double: Error reading %s.\n", input_filename);
        return (NULL);
    }

    /* read as double format; nvox already counts all dimensions, also the
       volumes of a 4D file */
    *image = (double *)malloc(sizeof(double) * nii_ptr->nvox);

    /* check for memory */
    if (!*image)
    {
        fprintf(stderr, "read_nifti_double: Memory allocation error\n");
        nifti_image_free(nii_ptr);
        return (NULL);
    }

    for (i = 0; i < nii_ptr->nvox; i++)
    {
        switch (nii_ptr->datatype)
        {
        case DT_INT8:
            tmp = (double)((signed char *)nii_ptr->data)[i];
            break;
        case DT_UINT8:
            tmp = (double)((unsigned char *)nii_ptr->data)[i];
            break;
        case DT_INT16:
            tmp = (double)((signed short *)nii_ptr->data)[i];
            break;
        case DT_UINT16:
            tmp = (double)((unsigned short *)nii_ptr->data)[i];
            break;
        case DT_INT32:
            tmp = (double)((signed int *)nii_ptr->data)[i];
            break;
        case DT_UINT32:
            tmp = (double)((unsigned int *)nii_ptr->data)[i];
            break;
        case DT_INT64:
            tmp = (double)((long long *)nii_ptr->data)[i];
            break;
        case DT_FLOAT32:
            tmp = (double)((float *)nii_ptr->data)[i];
            break;
        case DT_FLOAT64:
            tmp = (double)((double *)nii_ptr->data)[i];
            break;
        default:
            fprintf(stderr, "read_nifti_double: Unknown datatype\n");
            free(*image);
            *image = NULL;
            nifti_image_free(nii_ptr);
            return (NULL);
            break;
        }
        /* check whether scaling is needed */
        if (nii_ptr->scl_slope == 0)
            (*image)[i] = tmp;
        else
            (*image)[i] = (nii_ptr->scl_slope * tmp) + nii_ptr->scl_inter;
    }

    /* sto_xyz is what libCAT maps voxels with; a qform-only file left it zero */
    set_world_matrix(nii_ptr);

    /* ensure that nvox is that of a 3D image */
    if (nii_ptr->nt > 1)
    {
        nii_ptr->nvox /= nii_ptr->nt;
        nii_ptr->dim[0] = nii_ptr->ndim = 4;
    }

    if (!read_data)
    {
        free(nii_ptr->data);
        nii_ptr->data = NULL;
    }

    return (nii_ptr);
}

/**
 * \brief Read NIfTI image file and load data into single-precision array.
 *
 * Parses NIfTI format files and converts voxels to float (single precision).
 * Similar to read_nifti_double but saves memory for large volumes where double
 * precision is not required. Handles all NIfTI datatypes and applies scaling corrections.
 * Like read_nifti_double(), it reads a single-file NIfTI from the file named.
 *
 * A 4D file is returned as its nt volumes one after the other: volume v starts at
 * image[v * nvox], where nvox of the returned header is the voxel count of ONE
 * volume (nx*ny*nz) and nt the number of volumes.
 *
 * \param input_filename (in)  path to NIfTI file
 * \param image          (out) pointer to allocated float array (size nx*ny*nz*nt)
 * \param read_data      (in)  if non-zero, keep the raw NIfTI data in the returned
 *                             header; if 0, release it (image is returned either way)
 * \return pointer to NIfTI image structure on success; NULL on error
 */
nifti_image *read_nifti_float(const char *input_filename, float *image[], int read_data)
{
    nifti_image *nii_ptr;
    float tmp;
    size_t i;

    nii_ptr = nifti_image_read_exact(input_filename);
    if (!nii_ptr)
    {
        fprintf(stderr, "read_nifti_float: Error reading %s.\n", input_filename);
        return (NULL);
    }

    /* read as float format; nvox already counts all dimensions, also the
       volumes of a 4D file -- multiplying by nt again asked a 300-volume
       64x64x36 series for 53 GB */
    *image = (float *)malloc(sizeof(float) * nii_ptr->nvox);

    /* check for memory */
    if (!*image)
    {
        fprintf(stderr, "read_nifti_float: Memory allocation error\n");
        nifti_image_free(nii_ptr);
        return (NULL);
    }

    for (i = 0; i < nii_ptr->nvox; i++)
    {
        switch (nii_ptr->datatype)
        {
        case DT_INT8:
            tmp = (float)((signed char *)nii_ptr->data)[i];
            break;
        case DT_UINT8:
            tmp = (float)((unsigned char *)nii_ptr->data)[i];
            break;
        case DT_INT16:
            tmp = (float)((signed short *)nii_ptr->data)[i];
            break;
        case DT_UINT16:
            tmp = (float)((unsigned short *)nii_ptr->data)[i];
            break;
        case DT_INT32:
            tmp = (float)((signed int *)nii_ptr->data)[i];
            break;
        case DT_UINT32:
            tmp = (float)((unsigned int *)nii_ptr->data)[i];
            break;
        case DT_INT64:
            tmp = (float)((long long *)nii_ptr->data)[i];
            break;
        case DT_FLOAT32:
            tmp = (float)((float *)nii_ptr->data)[i];
            break;
        case DT_FLOAT64:
            tmp = (float)((double *)nii_ptr->data)[i];
            break;
        default:
            fprintf(stderr, "read_nifti_float: Unknown datatype\n");
            free(*image);
            *image = NULL;
            nifti_image_free(nii_ptr);
            return (NULL);
            break;
        }
        /* check whether scaling is needed */
        if (nii_ptr->scl_slope == 0)
            (*image)[i] = tmp;
        else
            (*image)[i] = ((float)nii_ptr->scl_slope * tmp) + (float)nii_ptr->scl_inter;
    }

    /* sto_xyz is what libCAT maps voxels with; a qform-only file left it zero */
    set_world_matrix(nii_ptr);

    /* ensure that nvox is that of a 3D image */
    if (nii_ptr->nt > 1)
    {
        nii_ptr->nvox /= nii_ptr->nt;
        nii_ptr->dim[0] = nii_ptr->ndim = 4;
    }

    if (!read_data)
    {
        free(nii_ptr->data);
        nii_ptr->data = NULL;
    }

    return (nii_ptr);
}
