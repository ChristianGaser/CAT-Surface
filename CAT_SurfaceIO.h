/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#define TRIANGLE_FILE_MAGIC_NUMBER  16777214
#define QUAD_FILE_MAGIC_NUMBER      16777215
#define NEW_VERSION_MAGIC_NUMBER    16777215

#define TYPE_REAL 0
#define TYPE_DOUBLE 1
#define TYPE_INTEGER 2
#define TYPE_CHAR 3

Status input_values_any_format(char *, int *, Real **);
Status input_values_integer(char *, int *, int **);
Status output_values_any_format(char *, int, void *, int);
Status input_graphics_any_format(char *, File_formats *, int *,
                                 object_struct  ***); 
Status output_graphics_any_format(char *, File_formats, int, object_struct  **);
Status output_txt(char *, int, Real *);
int    input_oogl(char *, File_formats *, int *, object_struct  ***);
int    output_oogl(char *, File_formats, int, object_struct * []);
int    output_freesurfer(char *, File_formats, int, object_struct * []);
int    input_freesurfer(char *, File_formats *, int *, object_struct  ***);
int    input_freesurfer_curv(char *, int *, Real **);
int    input_dx(char *, File_formats *, int *, object_struct  ***);
int    input_dfs(char *, File_formats *, int *, object_struct  ***);
Real * read_pgm(char *, int *, int *);
int    write_pgm(char *, Real *, int, int);
