/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#include <bicpl.h>
#include <stdio.h>
#include <string.h>
#include "gifticlib/gifti_io.h"
#include "niftilib/nifti1.h"
#include "niftilib/nifti1_io.h"

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

#define BIG_ENDIAN	4321
#define LITTLE_ENDIAN	1234
#define BYTE_ORDER	LITTLE_ENDIAN

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

Status bicpl_to_facevertexdata(polygons_struct *, double **, double **);
Status input_values_any_format(char *, int *, double **);
Status input_values_integer(char *, int *, int **);
Status output_values_any_format(char *, int, void *, int);
Status input_graphics_any_format(char *, File_formats *, int *, object_struct  ***); 
Status output_graphics_any_format(char *, File_formats, int, object_struct **, double *);
Status output_txt(char *, int, double *);
int    input_oogl(char *, File_formats *, int *, object_struct  ***);
int    output_oogl(char *, File_formats, int, object_struct * []);
int    output_freesurfer(char *, File_formats, int, object_struct * []);
int    output_freesurfer_curv(char *, int, double *);
int    input_freesurfer(char *, File_formats *, int *, object_struct  ***);
int    input_freesurfer_curv(char *, int *, double **);
int    output_gifti(char *, File_formats, int, object_struct * [], double *);
int    output_gifti_curv(char *, int, double *);
int    input_gifti(char *, File_formats *, int *, object_struct  ***);
int    input_gifti_curv(char *, int *, double **);
int    input_dx(char *, File_formats *, int *, object_struct  ***);
int    input_dfs(char *, File_formats *, int *, object_struct  ***);
double * read_pgm(char *, int *, int *);
int    write_pgm(char *, double *, int, int);
