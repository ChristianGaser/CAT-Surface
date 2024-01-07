#ifndef  DEF_BIC_PROGRAMMING_LIBRARY
#define  DEF_BIC_PROGRAMMING_LIBRARY

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BICAPI
#if defined(_MSC_VER)
/* If we are building on the Microsoft C compiler, we want to
 * explicitly import all public functions from the DLL
 */
#define BICAPI __declspec(dllimport)
#else
#define BICAPI
#endif /* _MSC_VER not defined */
#endif /* BICAPI not defined */

#include  <bicpl/data_structures.h>
#include  <bicpl/geom.h>
#include  <bicpl/trans.h>
#include  <bicpl/numerical.h>
#include  <bicpl/objects.h>
#include  <bicpl/prog_utils.h>
#include  <alloc.h>
//#include  <bicpl/vols.h>
//#include  <bicpl/images.h>

#define  for_less( i, start, end )  for( (i) = (start);  (i) < (end);  ++(i) )

/*! \mainpage The BIC Programming Library

BICPL is the BIC Programming Library, a set of functions for
manipulating volumes, objects such as polygons and lines, some basic
geometry, and some general programming utilities.  It is based on top
of the MNI Volume IO library, which is described in a separate
document.

BICPL is available at the McConnell Brain Imaging Centre (BIC) at the
Montreal Neurological Institute.  It was developed as part of a
medical imaging software testbed by David MacDonald, with source code
and considerable input from Peter Neelin, Louis Collins, and others at
the Centre.


\section compiling Compiling and Linking with bicpl

The library name is bicpl, and will usually reside in a
file named libbicpl.a (static) or libbicpl.so (shared).  The
main header file for the library is <bicpl.h>.

A typical program using bicpl will need to link against several
libraries: 
  cc test_bicpl.c  -o test_bicpl -lbicpl -lvolume_io -lminc -lnetcdf -lm

See the documentation for the BIC Volume IO Library for more
information on the list of libraries on the link line.

Note that you may need to specify where to find the headers and library
file by using -I and -L options, respectively.  For example,
  cc test_bicpl.c  -o test_bicpl -I/usr/local/mni/include -L/usr/local/mni/lib
  -lbicpl -lvolume_io -lminc -lnetcdf -lm


*/

#ifdef __cplusplus
}
#endif

#endif
