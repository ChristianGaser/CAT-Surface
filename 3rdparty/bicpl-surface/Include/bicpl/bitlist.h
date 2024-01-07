#ifndef  DEF_BITLIST
#define  DEF_BITLIST

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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : bitlist.h
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: A header file for routines that manipulate bitlists, which
            : are simply arrays of individually addressable bits.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */


#include  <volume_io.h>

typedef  unsigned int   bitlist_type;

#define  BITS_PER_BITLIST_WORD               32
#define  BITS_PER_BITLIST_WORD_MINUS_1       31
#define  LOG_BITS_PER_BITLIST_WORD            5

typedef  struct
{
    int            n_words;
    bitlist_type   *bits;
}  bitlist_struct;

typedef  struct
{
    int            nx, ny, nz;
    int            n_z_words;
    bitlist_type   ***bits;
}  bitlist_3d_struct;

#endif
