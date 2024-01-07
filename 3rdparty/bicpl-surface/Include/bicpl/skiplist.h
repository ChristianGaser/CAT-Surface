#ifndef  DEF_SKIPLIST
#define  DEF_SKIPLIST

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
@NAME       : skiplist.h
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Maintains a skiplist structure to handle ordered lists, searches.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */


#include  <volume_io.h>

typedef  struct  skip_struct
{
    float                   key;
    void                    *data_ptr;
    struct  skip_struct     *forward[1];
} skip_struct;

#define  LESS_THAN     -1
#define  EQUAL_TO       0
#define  GREATER_THAN   1

typedef  struct
{
    skip_struct    *header;
    int            level;
} skiplist_struct;

#endif
