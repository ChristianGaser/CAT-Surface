#ifndef  DEF_HASH
#define  DEF_HASH

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
@NAME       : hash.h
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Structures and definitions for the generic hash table routines.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */


#include  <volume_io.h>

typedef  struct  hash_entry_struct
{
    int                          key;
    struct   hash_entry_struct   *next;
    char                         data[1];
} hash_entry_struct;

typedef  struct
{
    int                 data_size;
    int                 size;
    int                 n_entries;
    Real                enlarge_threshold;
    Real                new_density;
    hash_entry_struct   **table;
} hash_table_struct;

typedef  struct
{
    int                current_index;
    hash_entry_struct  *current_entry;
} hash_table_pointer;

/*---- 2 key hash */ 

typedef  struct  hash2_entry_struct
{
    int                           key1;
    int                           key2;
    struct   hash2_entry_struct   *next;
    char                          data[1];
} hash2_entry_struct;

typedef  struct
{
    int                 data_size;
    int                 size;
    int                 n_entries;
    Real                enlarge_threshold;
    Real                new_density;
    hash2_entry_struct   **table;
} hash2_table_struct;

typedef  struct
{
    int                 current_index;
    hash2_entry_struct  *current_entry;
} hash2_table_pointer;

#endif
