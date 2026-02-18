/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#ifndef CAT_SURFACEIO_H
#define CAT_SURFACEIO_H

#include <bicpl.h>
#include <stdio.h>
#include <string.h>
#include <gifti_io.h>
#include <nifti1.h>
#include <nifti1_io.h>
#include <quadric.h>

#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)
#define NEW_VERSION_MAGIC_NUMBER    16777215

#define TYPE_DOUBLE 1
#define TYPE_INTEGER 2
#define TYPE_CHAR 3

/* allow override */
#ifndef BYTE_ORDER

/////////////Linux////////////////////////////
#ifdef linux
#include <endian.h>

#ifndef BYTE_ORDER
#define BYTE_ORDER __BYTE_ORDER
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN __LITTLE_ENDIAN
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN __BIG_ENDIAN
#endif

#endif

/////////////Windows Cygwin////////////////////////////
#ifdef WIN32

#define BIG_ENDIAN  4321
#define LITTLE_ENDIAN 1234
#define BYTE_ORDER  LITTLE_ENDIAN

#endif

////////////MacOS X and BSD ////////////////////////////
#if defined(__APPLE__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <machine/endian.h>
#endif

////////////Solaris 2.5.1//////////////////////
#ifdef sun

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN    4321
#endif

#include <sys/isa_defs.h>
/* only defines one of _LITTLE_ENDIAN or _BIG_ENDIAN */
#ifdef _LITTLE_ENDIAN
#define BYTE_ORDER LITTLE_ENDIAN
#endif

#ifdef _BIG_ENDIAN
#define BYTE_ORDER BIG_ENDIAN
#endif

#endif

/////////////IRIX  ////////////////////////////
#if defined(__sgi) || defined(Mips)
#include <sys/endian.h>
#endif

///////////////////////////////////////////////////
#endif /* BYTE_ORDER */
///////////////////////////////////////////////////

typedef struct
{
  int    r, g, b ;
  int    annotation ;
  char   name[1000] ;
}
ATABLE ;

int polygons_to_tri_arrays(const polygons_struct *poly,
                                      vec3d **V, vec3i **F,
                                      int *nv, int *nf, int *fan_used);
int tri_arrays_to_polygons(polygons_struct *poly,
                                      const vec3d *V, const vec3i *F,
                                      int nv, int nf);
int qem_target(int nf_total, int target);
Status bicpl_to_facevertexdata(polygons_struct *, double **, double **);
/**
 * \brief Public API for input_values_any_format.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of input_values_any_format.
 * \param param (in/out) Parameter of input_values_any_format.
 * \param param (in/out) Parameter of input_values_any_format.
 * \return Return value of input_values_any_format.
 */
Status input_values_any_format(char *, int *, double **);
Status input_values_integer(char *, int *, int **);
/**
 * \brief Public API for output_values_any_format.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of output_values_any_format.
 * \param int (in/out) Parameter of output_values_any_format.
 * \param param (in/out) Parameter of output_values_any_format.
 * \param int (in/out) Parameter of output_values_any_format.
 * \return Return value of output_values_any_format.
 */
Status output_values_any_format(const char *, int, void *, int);
/**
 * \brief Public API for input_graphics_any_format.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of input_graphics_any_format.
 * \param param (in/out) Parameter of input_graphics_any_format.
 * \param param (in/out) Parameter of input_graphics_any_format.
 * \param param (in/out) Parameter of input_graphics_any_format.
 * \return Return value of input_graphics_any_format.
 */
Status input_graphics_any_format(char *, File_formats *, int *, object_struct  ***); 
/**
 * \brief Public API for output_graphics_any_format.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of output_graphics_any_format.
 * \param File_formats (in/out) Parameter of output_graphics_any_format.
 * \param int (in/out) Parameter of output_graphics_any_format.
 * \param param (in/out) Parameter of output_graphics_any_format.
 * \param param (in/out) Parameter of output_graphics_any_format.
 * \return Return value of output_graphics_any_format.
 */
Status output_graphics_any_format(char *, File_formats, int, object_struct **, double *);
int    input_gifti_mesh_and_texture(char *, File_formats *, int *, object_struct ***, int *, double **);
Status output_txt(char *, int, double *);
int    input_oogl(char *, File_formats *, int *, object_struct  ***);
int    output_oogl(char *, File_formats, int, object_struct * []);
int    output_freesurfer(char *, File_formats, int, object_struct * []);
int    output_freesurfer_curv(char *, int, double *);
int    input_freesurfer(char *, File_formats *, int *, object_struct  ***);
/**
 * \brief Public API for input_freesurfer_curv.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of input_freesurfer_curv.
 * \param param (in/out) Parameter of input_freesurfer_curv.
 * \param param (in/out) Parameter of input_freesurfer_curv.
 * \return Return value of input_freesurfer_curv.
 */
int    input_freesurfer_curv(char *, int *, double **);
int    output_gifti(char *, File_formats, int, object_struct * [], double *);
int    output_gifti_curv(char *, int, double *);
int    input_gifti(char *, File_formats *, int *, object_struct  ***, int *, double **);
int    input_gifti_curv(char *, int *, double **);
int    input_dx(char *, File_formats *, int *, object_struct  ***);
int    input_dfs(char *, File_formats *, int *, object_struct  ***);
/**
 * \brief Public API for read_pgm.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of read_pgm.
 * \param param (in/out) Parameter of read_pgm.
 * \param param (in/out) Parameter of read_pgm.
 * \return Return value of read_pgm.
 */
double * read_pgm(char *, int *, int *);
/**
 * \brief Public API for write_pgm.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of write_pgm.
 * \param param (in/out) Parameter of write_pgm.
 * \param int (in/out) Parameter of write_pgm.
 * \param int (in/out) Parameter of write_pgm.
 * \return Return value of write_pgm.
 */
int    write_pgm(char *, double *, int, int);
/**
 * \brief Public API for read_annotation_table.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of read_annotation_table.
 * \param param (in/out) Parameter of read_annotation_table.
 * \param param (in/out) Parameter of read_annotation_table.
 * \param param (in/out) Parameter of read_annotation_table.
 * \param param (in/out) Parameter of read_annotation_table.
 * \return Return value of read_annotation_table.
 */
int    read_annotation_table(char *, int *, int **, int *, ATABLE **);
/**
 * \brief Public API for write_annotation_table.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of write_annotation_table.
 * \param int (in/out) Parameter of write_annotation_table.
 * \param param (in/out) Parameter of write_annotation_table.
 * \param int (in/out) Parameter of write_annotation_table.
 * \param param (in/out) Parameter of write_annotation_table.
 * \return Return value of write_annotation_table.
 */
int    write_annotation_table(char *, int, int *, int, ATABLE *);

#endif
