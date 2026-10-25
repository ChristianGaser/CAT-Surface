#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_SurfaceIO.h"

/*
 * The surface and annotation readers must turn a malformed file into ERROR.
 * They used to trust the file: a face index past the last vertex was
 * dereferenced by compute_polygon_normals(), and the counts in a header went
 * straight into allocations, so a 44-byte file could ask for gigabytes and a
 * failed allocation was written through as NULL plus a file offset.
 */

/** \brief Write a big-endian int, the byte order of FreeSurfer files. */
static void
put_int(FILE *fp, int v)
{
    unsigned char b[4];

    b[0] = (unsigned char)((unsigned int)v >> 24);
    b[1] = (unsigned char)((unsigned int)v >> 16);
    b[2] = (unsigned char)((unsigned int)v >> 8);
    b[3] = (unsigned char)v;
    fwrite(b, 1, 4, fp);
}

/** \brief Write a big-endian float. */
static void
put_float(FILE *fp, float f)
{
    int v;

    memcpy(&v, &f, sizeof(v));
    put_int(fp, v);
}

/**
 * \brief Write a FreeSurfer triangle surface of a tetrahedron.
 *
 * The header declares n_points and n_faces; n_write_points points and
 * n_write_faces faces are actually written, so a header can promise more
 * than the file holds.  bad_index replaces the last face index if >= 0.
 */
static void
write_fs_surface(const char *file, int n_points, int n_faces,
                 int n_write_points, int n_write_faces, int bad_index)
{
    static const float v[4][3] = {{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {0, 0, 10}};
    static const int f[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    unsigned char magic[3] = {0xff, 0xff, 0xfe};
    FILE *fp = fopen(file, "wb");
    int i, j;

    fwrite(magic, 1, 3, fp);
    fputs("created by test_readers\n\n", fp);
    put_int(fp, n_points);
    put_int(fp, n_faces);
    for (i = 0; i < n_write_points; i++)
        for (j = 0; j < 3; j++)
            put_float(fp, v[i % 4][j]);
    for (i = 0; i < n_write_faces; i++)
        for (j = 0; j < 3; j++)
            put_int(fp, (bad_index >= 0 && i == n_write_faces - 1 && j == 2)
                            ? bad_index : f[i % 4][j]);
    fclose(fp);
}

/** \brief Read a surface; returns the reader's status and the point count. */
static Status
read_surface(char *file, int *n_points)
{
    File_formats format;
    int n_objects = 0;
    object_struct **objects = NULL;
    Status status;

    *n_points = -1;
    status = input_graphics_any_format(file, &format, &n_objects, &objects);
    if (status == OK && n_objects == 1)
        *n_points = get_polygons_ptr(objects[0])->n_points;
    if (n_objects > 0)
        delete_object_list(n_objects, objects);
    return status;
}

static void
test_freesurfer_surface(void)
{
    char file[] = "/tmp/test_readers_lh.surf";
    int n;

    write_fs_surface(file, 4, 4, 4, 4, -1);
    MU_ASSERT("a valid FreeSurfer surface is read",
              read_surface(file, &n) == OK && n == 4);

    write_fs_surface(file, 4, 4, 4, 4, 1000000000);
    MU_ASSERT("a face index past the last vertex is an error",
              read_surface(file, &n) == ERROR);

    write_fs_surface(file, -5, -5, 0, 0, -1);
    MU_ASSERT("negative counts are an error", read_surface(file, &n) == ERROR);

    write_fs_surface(file, 1073741824, 1, 0, 0, -1);
    MU_ASSERT("counts the file does not hold are an error",
              read_surface(file, &n) == ERROR);

    write_fs_surface(file, 1, 715827883, 0, 0, -1);
    MU_ASSERT("a face count past INT_MAX / 3 is an error",
              read_surface(file, &n) == ERROR);

    write_fs_surface(file, 4, 4, 4, 2, -1);
    MU_ASSERT("a truncated surface is an error", read_surface(file, &n) == ERROR);

    remove(file);
}

static void
test_off_surface(void)
{
    char file[] = "/tmp/test_readers.off";
    FILE *fp;
    int n;

    fp = fopen(file, "w");
    fputs("OFF\n4 4 0\n0 0 0\n10 0 0\n0 10 0\n0 0 10\n"
          "3 0 2 1\n3 0 1 3\n3 0 3 2\n3 1 2 3\n", fp);
    fclose(fp);
    MU_ASSERT("a valid OFF surface is read", read_surface(file, &n) == OK && n == 4);

    fp = fopen(file, "w");
    fputs("OFF\n4 4 0\n0 0 0\n10 0 0\n0 10 0\n0 0 10\n"
          "3 0 2 1\n3 0 1 3\n3 0 3 2\n3 1 2 999\n", fp);
    fclose(fp);
    MU_ASSERT("a face index past the last vertex is an error in OFF too",
              read_surface(file, &n) == ERROR);

    remove(file);
}

/** \brief Write an annot: n_vertices pairs, then the given colour-table head. */
static void
write_annot(const char *file, int n_vertices, int n_written, int table_count)
{
    FILE *fp = fopen(file, "wb");
    int i;

    put_int(fp, n_vertices);
    for (i = 0; i < n_written; i++)
    {
        put_int(fp, i);
        put_int(fp, 0x0a0b0c);
    }
    if (table_count != 0)
    {
        put_int(fp, 1);             /* tag: a colour table follows */
        put_int(fp, table_count);   /* old layout: number of entries */
        put_int(fp, 2);
        fwrite("x", 1, 2, fp);      /* table name with its NUL */
    }
    fclose(fp);
}

static void
test_annotation(void)
{
    char file[] = "/tmp/test_readers.annot";
    FILE *fp;
    int n_array, n_labels, *array = NULL;
    ATABLE *atable = NULL;

    fp = fopen(file, "wb");
    put_int(fp, 3);
    put_int(fp, 0);
    put_int(fp, 0x0a0b0c);
    put_int(fp, 1);
    put_int(fp, 0x0a0b0c);
    put_int(fp, 2);
    put_int(fp, 0x0a0b0c);
    put_int(fp, 1);                 /* tag: a colour table follows */
    fclose(fp);
    MU_ASSERT("a colour-table tag without a count is an error",
              read_annotation_table(file, &n_array, &array, &n_labels, &atable) == ERROR &&
              array == NULL);

    write_annot(file, -1, 0, 0);
    MU_ASSERT("a negative vertex count is an error",
              read_annotation_table(file, &n_array, &array, &n_labels, &atable) == ERROR);

    write_annot(file, 100000000, 3, 0);
    MU_ASSERT("a vertex count the file does not hold is an error",
              read_annotation_table(file, &n_array, &array, &n_labels, &atable) == ERROR);

    write_annot(file, 3, 3, 2147483647);
    MU_ASSERT("a colour table the file does not hold is an error",
              read_annotation_table(file, &n_array, &array, &n_labels, &atable) == ERROR);

    n_labels = -7;
    atable = (ATABLE *)1;
    write_annot(file, 3, 3, 0);
    MU_ASSERT("codes without a colour table are read",
              read_annotation_table(file, &n_array, &array, &n_labels, &atable) == OK &&
              n_array == 3 && array[2] == 0x0a0b0c);
    MU_ASSERT("and report an empty colour table", n_labels == 0 && atable == NULL);

    free(array);
    remove(file);
}

int main(void)
{
    MU_RUN_TEST(test_freesurfer_surface);
    MU_RUN_TEST(test_off_surface);
    MU_RUN_TEST(test_annotation);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
