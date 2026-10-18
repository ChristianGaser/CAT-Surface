#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_NiftiLib.h"

/*
 * 4D files and the template header of the writers.
 *
 * read_nifti_double() returns the volumes of a 4D file one after the other and
 * reports nvox of ONE volume, which is what the tools index with.  It used to
 * allocate nt times what it needed (nvox already counts every dimension), and
 * a header taken from a 4D file made the writers store dim[4] = nt in a file
 * with a single volume -- FSL reads that as nt volumes.
 */

#define NX 3
#define NY 2
#define NZ 2
#define NVOX (NX * NY * NZ)
#define NT 5

/** \brief Write an NT-volume float32 file; voxel i of volume t holds 100 t + i. */
static void
write_4d(const char *file)
{
    int dims[8] = {4, NX, NY, NZ, NT, 1, 1, 1};
    nifti_image *nim = nifti_make_new_nim(dims, DT_FLOAT32, 1);
    int t, i;

    for (t = 0; t < NT; t++)
        for (i = 0; i < NVOX; i++)
            ((float *)nim->data)[t * NVOX + i] = (float)(100 * t + i);

    nifti_set_filenames(nim, file, 0, 1);
    nifti_image_write(nim);
    nifti_image_free(nim);
}

/** \brief Read the dim[] array of the 348-byte header and the file size. */
static int
read_raw_header(const char *file, short dim[8], long *size)
{
    FILE *fp = fopen(file, "rb");

    if (!fp)
        return 0;
    fseek(fp, 40, SEEK_SET);
    if (fread(dim, sizeof(short), 8, fp) != 8)
    {
        fclose(fp);
        return 0;
    }
    fseek(fp, 0, SEEK_END);
    *size = ftell(fp);
    fclose(fp);
    return 1;
}

/** \brief A 4D file comes back as NT volumes of nvox voxels each. */
static void
test_4d_read_layout(void)
{
    const char *file = "test_nifti_4d.nii";
    double *img = NULL;
    nifti_image *nii;
    int t, i, ok = 1;

    write_4d(file);
    nii = read_nifti_double(file, &img, 0);

    MU_ASSERT("the 4D file is read", nii != NULL && img != NULL);
    if (nii && img)
    {
        MU_ASSERT("nt is the number of volumes", nii->nt == NT);
        MU_ASSERT("nvox is that of one volume", nii->nvox == NVOX);
        MU_ASSERT("the raw data are released with read_data = 0", nii->data == NULL);
        for (t = 0; t < NT; t++)
            for (i = 0; i < NVOX; i++)
                if (img[t * NVOX + i] != (double)(100 * t + i))
                    ok = 0;
        MU_ASSERT("volume t starts at t * nvox", ok);

        /* data is NULL, so freeing the header must not free it a second time */
        nifti_image_free(nii);
    }
    free(img);
    remove(file);
}

/** \brief read_nifti_float() returns the same layout as read_nifti_double(). */
static void
test_4d_read_float_layout(void)
{
    const char *file = "test_nifti_4d_float.nii";
    float *img = NULL;
    nifti_image *nii;
    int t, i, ok = 1;

    write_4d(file);
    nii = read_nifti_float(file, &img, 0);

    MU_ASSERT("float: the 4D file is read", nii != NULL && img != NULL);
    if (nii && img)
    {
        MU_ASSERT("float: nt is the number of volumes", nii->nt == NT);
        MU_ASSERT("float: nvox is that of one volume", nii->nvox == NVOX);
        MU_ASSERT("float: the raw data are released with read_data = 0",
                  nii->data == NULL);
        for (t = 0; t < NT; t++)
            for (i = 0; i < NVOX; i++)
                if (img[t * NVOX + i] != (float)(100 * t + i))
                    ok = 0;
        MU_ASSERT("float: volume t starts at t * nvox", ok);
        nifti_image_free(nii);
    }
    free(img);

    /* read_data keeps the raw data, which the header then owns */
    nii = read_nifti_float(file, &img, 1);
    MU_ASSERT("float: read_data = 1 keeps the raw data",
              nii != NULL && nii->data != NULL);
    if (nii)
        nifti_image_free(nii);
    free(img);
    remove(file);
}

/** \brief Both writers store one 3D volume although the template is 4D. */
static void
test_writers_drop_the_time_dimension(void)
{
    const char *src = "test_nifti_src.nii";
    const char *dbl = "test_nifti_dbl.nii";
    const char *flt = "test_nifti_flt.nii";
    int dims[3] = {NX, NY, NZ};
    double vox[3] = {2.0, 2.0, 2.0};
    float fimg[NVOX];
    double *img = NULL;
    nifti_image *nii;
    short dim[8];
    long size;
    int i;

    write_4d(src);
    nii = read_nifti_double(src, &img, 0);
    MU_ASSERT("the 4D template is read", nii != NULL && img != NULL);
    if (!nii || !img)
        return;
    for (i = 0; i < NVOX; i++)
        fimg[i] = (float)img[i];

    MU_ASSERT("write_nifti_double succeeds",
              write_nifti_double(dbl, img, DT_FLOAT32, 1.0, dims, vox, nii));
    MU_ASSERT("double: the header describes a 3D image",
              read_raw_header(dbl, dim, &size) && dim[0] == 3);
    MU_ASSERT("double: dim[4] is 1", dim[4] == 1 && dim[5] == 1);
    MU_ASSERT("double: the file holds one volume",
              size == 352 + (long)sizeof(float) * NVOX);

    /* the template is reusable and freeable after a write */
    MU_ASSERT("write_nifti_float succeeds",
              write_nifti_float(flt, fimg, DT_FLOAT32, 1.0, dims, vox, nii));
    MU_ASSERT("float: the header describes a 3D image",
              read_raw_header(flt, dim, &size) && dim[0] == 3);
    MU_ASSERT("float: dim[4] is 1", dim[4] == 1 && dim[5] == 1);
    MU_ASSERT("float: the file holds one volume",
              size == 352 + (long)sizeof(float) * NVOX);

    nifti_image_free(nii);
    free(img);
    remove(src);
    remove(dbl);
    remove(flt);
}

/** \brief write_nifti_float_4d() stores every volume and reads back unchanged. */
static void
test_4d_writer_round_trip(void)
{
    const char *src = "test_nifti_4d_src.nii";
    const char *out = "test_nifti_4d_out.nii";
    const char *out_int = "test_nifti_4d_int16.nii";
    const char *out_new = "test_nifti_4d_new.nii";
    int dims[3] = {NX, NY, NZ};
    double vox[3] = {2.0, 2.0, 2.0};
    float *img = NULL, *back = NULL;
    nifti_image *nii, *nii2;
    short dim[8];
    long size;
    int t, i, ok = 1, intact = 1, close = 1;

    write_4d(src);
    nii = read_nifti_float(src, &img, 0);
    MU_ASSERT("4D writer: the template is read", nii != NULL && img != NULL);
    if (!nii || !img)
        return;
    nii->dt = 2.5f; /* the repetition time has to survive */

    MU_ASSERT("write_nifti_float_4d succeeds",
              write_nifti_float_4d(out, img, DT_FLOAT32, 0.0, dims, NT, vox, nii));
    MU_ASSERT("4D writer: the header describes NT volumes",
              read_raw_header(out, dim, &size) && dim[0] == 4 && dim[4] == NT);
    MU_ASSERT("4D writer: the file holds every volume",
              size == 352 + (long)sizeof(float) * NVOX * NT);
    for (t = 0; t < NT; t++)
        for (i = 0; i < NVOX; i++)
            if (img[t * NVOX + i] != (float)(100 * t + i))
                intact = 0;
    MU_ASSERT("4D writer: the caller's buffer is left as it was", intact);

    nii2 = read_nifti_float(out, &back, 0);
    MU_ASSERT("4D writer: the output is read", nii2 != NULL && back != NULL);
    if (nii2 && back)
    {
        MU_ASSERT("4D writer: nt survives", nii2->nt == NT);
        MU_ASSERT("4D writer: the time step survives", fabs(nii2->dt - 2.5) < 1e-6);
        for (i = 0; i < NVOX * NT; i++)
            if (back[i] != img[i])
                ok = 0;
        MU_ASSERT("4D writer: float32 values survive exactly", ok);
        nifti_image_free(nii2);
    }
    free(back);

    /* integer output is scaled over all volumes together */
    MU_ASSERT("write_nifti_float_4d int16 succeeds",
              write_nifti_float_4d(out_int, img, DT_INT16, 0.0, dims, NT, vox, nii));
    nii2 = read_nifti_float(out_int, &back, 0);
    MU_ASSERT("4D writer: the int16 output is read", nii2 != NULL && back != NULL);
    if (nii2 && back)
    {
        for (i = 0; i < NVOX * NT; i++)
            if (fabs(back[i] - img[i]) > 0.5 * nii2->scl_slope + 1e-6)
                close = 0;
        MU_ASSERT("4D writer: int16 values within half a step", close);
        nifti_image_free(nii2);
    }
    free(back);

    /* without a template the series gets a time step of 1 */
    MU_ASSERT("write_nifti_float_4d without template succeeds",
              write_nifti_float_4d(out_new, img, DT_FLOAT32, 0.0, dims, NT, vox, NULL));
    nii2 = read_nifti_float(out_new, &back, 0);
    MU_ASSERT("4D writer: a new series has nt volumes and dt 1",
              nii2 != NULL && nii2->nt == NT && fabs(nii2->dt - 1.0) < 1e-6);
    if (nii2)
        nifti_image_free(nii2);
    free(back);

    MU_ASSERT("write_nifti_float_4d rejects nt = 0",
              !write_nifti_float_4d(out_new, img, DT_FLOAT32, 0.0, dims, 0, vox, NULL));

    nifti_image_free(nii);
    free(img);
    remove(src);
    remove(out);
    remove(out_int);
    remove(out_new);
}

/**
 * \brief Write a 3D float file with a qform, and an sform if requested.
 *
 * The qform turns 180 degrees about z with voxel sizes 2, 3 and 4 mm, so its
 * matrix is diag(-2, -3, 4) with offset (10, -20, 30); the sform is diag(2, 3,
 * 4) with offset (1, 2, 3).
 */
static void
write_xform(const char *file, int with_sform)
{
    int dims[8] = {3, NX, NY, NZ, 1, 1, 1, 1};
    nifti_image *nim = nifti_make_new_nim(dims, DT_FLOAT32, 1);
    int i, j;

    nim->dx = nim->pixdim[1] = 2.0f;
    nim->dy = nim->pixdim[2] = 3.0f;
    nim->dz = nim->pixdim[3] = 4.0f;
    nim->qform_code = NIFTI_XFORM_SCANNER_ANAT;
    nim->quatern_b = 0.0f;
    nim->quatern_c = 0.0f;
    nim->quatern_d = 1.0f;
    nim->qoffset_x = 10.0f;
    nim->qoffset_y = -20.0f;
    nim->qoffset_z = 30.0f;
    nim->qfac = 1.0f;
    nim->sform_code = NIFTI_XFORM_UNKNOWN;
    if (with_sform)
    {
        nim->sform_code = NIFTI_XFORM_ALIGNED_ANAT;
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                nim->sto_xyz.m[i][j] = 0.0f;
        nim->sto_xyz.m[0][0] = 2.0f;
        nim->sto_xyz.m[1][1] = 3.0f;
        nim->sto_xyz.m[2][2] = 4.0f;
        nim->sto_xyz.m[0][3] = 1.0f;
        nim->sto_xyz.m[1][3] = 2.0f;
        nim->sto_xyz.m[2][3] = 3.0f;
        nim->sto_xyz.m[3][3] = 1.0f;
    }
    nifti_set_filenames(nim, file, 0, 1);
    nifti_image_write(nim);
    nifti_image_free(nim);
}

/** \brief Without an sform, sto_xyz is the qform; with one, it is the sform. */
static void
test_world_matrix(void)
{
    const char *qonly = "test_nifti_qform.nii";
    const char *both = "test_nifti_sform.nii";
    const char *back = "test_nifti_qform_back.nii";
    int dims[3] = {NX, NY, NZ};
    double vox[3] = {2.0, 3.0, 4.0};
    float *img = NULL;
    nifti_image *nii, *nii2;
    mat44 m;

    write_xform(qonly, 0);
    nii = read_nifti_float(qonly, &img, 0);
    MU_ASSERT("qform only: the file is read", nii != NULL);
    if (nii)
    {
        m = nii->sto_xyz;
        MU_ASSERT("qform only: sto_xyz is the qform",
                  m.m[0][0] == -2.0f && m.m[1][1] == -3.0f && m.m[2][2] == 4.0f &&
                  m.m[0][3] == 10.0f && m.m[1][3] == -20.0f && m.m[2][3] == 30.0f);
        MU_ASSERT("qform only: sto_ijk is its inverse",
                  fabs(nii->sto_ijk.m[0][0] + 0.5) < 1e-6 &&
                  fabs(nii->sto_ijk.m[0][3] - 5.0) < 1e-5);
        MU_ASSERT("qform only: sform_code stays 0", nii->sform_code == 0);

        /* written back, the header still has no sform */
        MU_ASSERT("qform only: written back",
                  write_nifti_float(back, img, DT_FLOAT32, 0.0, dims, vox, nii));
        nii2 = nifti_image_read(back, 0);
        MU_ASSERT("qform only: no sform is invented on writing",
                  nii2 != NULL && nii2->sform_code == 0 && nii2->qform_code > 0);
        if (nii2)
            nifti_image_free(nii2);
        nifti_image_free(nii);
    }
    free(img);

    write_xform(both, 1);
    nii = read_nifti_float(both, &img, 0);
    MU_ASSERT("sform and qform: the file is read", nii != NULL);
    if (nii)
    {
        m = nii->sto_xyz;
        MU_ASSERT("sform and qform: sto_xyz is the sform",
                  m.m[0][0] == 2.0f && m.m[1][1] == 3.0f && m.m[2][2] == 4.0f &&
                  m.m[0][3] == 1.0f && m.m[1][3] == 2.0f && m.m[2][3] == 3.0f);
        nifti_image_free(nii);
    }
    free(img);
    remove(qonly);
    remove(both);
    remove(back);
}

int main(void)
{
    MU_RUN_TEST(test_4d_read_layout);
    MU_RUN_TEST(test_4d_read_float_layout);
    MU_RUN_TEST(test_writers_drop_the_time_dimension);
    MU_RUN_TEST(test_4d_writer_round_trip);
    MU_RUN_TEST(test_world_matrix);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
