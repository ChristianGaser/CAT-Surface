#ifndef  DEF_GLOBAL_LOOKUP
#define  DEF_GLOBAL_LOOKUP

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


#include  <volume_io.h>

typedef  enum {
                 BOOLEAN_type,
                 int_type,
                 Real_type,
                 STRING_type,
                 Point_type,
                 Vector_type,
                 Colour_type,
                 Surfprop_type
              } Variable_types;

typedef  struct 
{
    void             *ptr_to_global;
    STRING           variable_name;
    Variable_types   type;
    Smallest_int     set_flag;
} global_struct;

#endif
