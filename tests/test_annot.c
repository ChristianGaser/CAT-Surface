#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_SurfaceIO.h"

/*
 * Annotation tables must survive a write/read round trip with their names
 * intact.  Before cat-surf 1.0.29 the writer padded each name with a space in
 * place of the terminating NUL and the reader did not terminate what it read
 * into a malloc'd table.  Resampling an annot that CAT itself had written
 * (T1Prep's DK40 and a2009s templates) therefore produced region names like
 * "bankssts M(C\xa0\x17\"C..." -- the heap bytes following the name.
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

/** \brief Read a big-endian int at byte offset pos of buf. */
static int
get_int(const unsigned char *buf, long pos)
{
    return (int)(((unsigned int)buf[pos] << 24) | ((unsigned int)buf[pos + 1] << 16) |
                 ((unsigned int)buf[pos + 2] << 8) | (unsigned int)buf[pos + 3]);
}

/**
 * \brief Write a two-label annot the way write_annotation_table used to.
 *
 * The first name is stored with the given length, which may exceed the
 * reader's 1000-byte buffer, and without a terminating NUL: its last byte is
 * a space, as the old writer left it.
 */
static void
write_legacy_annot(const char *file, int first_len)
{
    FILE *fp = fopen(file, "wb");
    int i;

    put_int(fp, 3);                     /* vertices */
    put_int(fp, 0); put_int(fp, 1 + 2 * 256 + 3 * 65536);
    put_int(fp, 1); put_int(fp, 4 + 5 * 256 + 6 * 65536);
    put_int(fp, 2); put_int(fp, 4 + 5 * 256 + 6 * 65536);
    put_int(fp, 1);                     /* colortable follows */
    put_int(fp, -2);                    /* version 2 */
    put_int(fp, 2);                     /* max structure */
    put_int(fp, 4); fwrite("orig", 1, 4, fp);
    put_int(fp, 2);                     /* entries */

    put_int(fp, 0);
    put_int(fp, first_len);
    for (i = 0; i < first_len - 1; i++)
        fputc(i < 7 ? "unknown"[i] : 'x', fp);
    fputc(' ', fp);
    put_int(fp, 1); put_int(fp, 2); put_int(fp, 3); put_int(fp, 0);

    put_int(fp, 1);
    put_int(fp, 9); fwrite("bankssts ", 1, 9, fp);
    put_int(fp, 4); put_int(fp, 5); put_int(fp, 6); put_int(fp, 0);
    fclose(fp);
}

/** \brief Names written by CAT's own old writer come back terminated. */
static void
test_legacy_names_are_terminated(void)
{
    char file[] = "/tmp/test_annot_legacy.annot";
    int n_array, n_labels, *array = NULL;
    ATABLE *atable = NULL;

    write_legacy_annot(file, 8);        /* "unknown" + space */
    read_annotation_table(file, &n_array, &array, &n_labels, &atable);

    MU_ASSERT("three vertices are read", n_array == 3);
    MU_ASSERT("two labels are read", n_labels == 2);
    MU_ASSERT("the padding space is dropped",
              strcmp(atable[0].name, "unknown") == 0);
    MU_ASSERT("the second name is terminated",
              strcmp(atable[1].name, "bankssts") == 0);
    MU_ASSERT("the colour after the name is intact",
              atable[1].r == 4 && atable[1].g == 5 && atable[1].b == 6);

    free(array);
    free(atable);
    remove(file);
}

/** \brief A name longer than the buffer is cut, and the file stays in step. */
static void
test_oversized_name_is_bounded(void)
{
    char file[] = "/tmp/test_annot_long.annot";
    int n_array, n_labels, *array = NULL;
    ATABLE *atable = NULL;

    write_legacy_annot(file, 1500);
    read_annotation_table(file, &n_array, &array, &n_labels, &atable);

    MU_ASSERT("the long name fills the buffer and no more",
              strlen(atable[0].name) == sizeof(atable[0].name) - 1);
    MU_ASSERT("the rest of the long name is skipped",
              atable[0].r == 1 && atable[0].g == 2 && atable[0].b == 3);
    MU_ASSERT("the entry after it is read in step",
              strcmp(atable[1].name, "bankssts") == 0 && atable[1].b == 6);

    free(array);
    free(atable);
    remove(file);
}

/** \brief The writer produces FreeSurfer's layout and leaves its input alone. */
static void
test_round_trip_is_freesurfer_standard(void)
{
    char file[] = "/tmp/test_annot_roundtrip.annot";
    ATABLE table[2];
    int labels[3], n_array, n_labels, *array = NULL;
    ATABLE *atable = NULL;
    unsigned char buf[512];
    long size, pos;
    FILE *fp;

    memset(table, 0, sizeof(table));
    strcpy(table[0].name, "unknown");
    strcpy(table[1].name, "bankssts");
    table[0].r = 1; table[0].g = 2; table[0].b = 3;
    table[1].r = 4; table[1].g = 5; table[1].b = 6;
    labels[0] = 1 + 2 * 256 + 3 * 65536;
    labels[1] = labels[2] = 4 + 5 * 256 + 6 * 65536;

    write_annotation_table(file, 3, labels, 2, table);
    MU_ASSERT("writing leaves the caller's names unchanged",
              strcmp(table[1].name, "bankssts") == 0);

    fp = fopen(file, "rb");
    size = (long)fread(buf, 1, sizeof(buf), fp);
    fclose(fp);
    /* 3 vertices -> 4 + 3 * 8 bytes, then tag, version, max structure,
     * the 22-byte original file name and the entry count */
    pos = 4 + 3 * 8 + 4 + 4 + 4 + 4 + 22 + 4;
    MU_ASSERT("the file is complete", size > pos + 8);
    MU_ASSERT("the first length counts the NUL", get_int(buf, pos + 4) == 8);
    MU_ASSERT("the first name ends in a NUL",
              memcmp(buf + pos + 8, "unknown\0", 8) == 0);

    read_annotation_table(file, &n_array, &array, &n_labels, &atable);
    MU_ASSERT("the names come back", strcmp(atable[0].name, "unknown") == 0 &&
              strcmp(atable[1].name, "bankssts") == 0);
    MU_ASSERT("the labels come back", array[0] == labels[0] &&
              array[1] == labels[1] && array[2] == labels[2]);

    free(array);
    free(atable);
    remove(file);
}

int main(void)
{
    MU_RUN_TEST(test_legacy_names_are_terminated);
    MU_RUN_TEST(test_oversized_name_is_bounded);
    MU_RUN_TEST(test_round_trip_is_freesurfer_standard);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
