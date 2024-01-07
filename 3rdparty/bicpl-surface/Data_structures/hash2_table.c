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

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Data_structures/hash2_table.c,v 1.4 2005/08/17 22:31:12 bert Exp $";
#endif

#define  HASH1_FUNCTION_CONSTANT          0.6180339887498948482
#define  HASH2_FUNCTION_CONSTANT          0.2794536923672642321

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_hash2_table
@INPUT      : 
            : size
            : enlarge_threshold
            : new_density
@OUTPUT     : hash_table
@RETURNS    : Status
@DESCRIPTION: Initializes a hash table to empty, storing the
            : initial size, and parameters controlling automatic growth.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void  initialize_hash2_table(
    hash2_table_struct  *hash_table,
    int                 size,
    int                 data_size,
    Real                enlarge_threshold,
    Real                new_density )
{
    int        i;

    hash_table->size = size;
    hash_table->data_size = data_size;
    hash_table->n_entries = 0;
    hash_table->enlarge_threshold = enlarge_threshold;
    hash_table->new_density = new_density;

    ALLOC( hash_table->table, size );

    for( i = 0;  i < size;  ++i )
        hash_table->table[i] = NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_hash2_table_list
@INPUT      : hash_table
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Deletes the list of pointers in the hash table.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   void  delete_hash2_table_list(
    hash2_table_struct  *hash_table )
{
    if( hash_table->size > 0 )
        FREE( hash_table->table );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_hash2_table
@INPUT      : hash_table
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Deletes the given hash table.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void  delete_hash2_table(
    hash2_table_struct  *hash_table )
{
    int                  i;
    hash2_entry_struct   *entry, *next;

    for( i = 0;  i < hash_table->size;  ++i )
    {
        entry = hash_table->table[i];

        while( entry != NULL )
        {
            next = entry->next;
            FREE( entry );
            entry = next;
        }
    }

    delete_hash2_table_list( hash_table );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : hash2_function
@INPUT      : hash_table
            : key1
            : key2
@OUTPUT     : 
@RETURNS    : an index into the hash table.
@DESCRIPTION: Hashes the two keys to create an index into the table.  Uses
              the multiplicative hash method from Corman book on algorithms.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   int  hash2_function(
    hash2_table_struct   *hash_table,
    int                 key1,
    int                 key2 )
{
    int    index;
    Real   v;

    v = (Real) key1 * HASH1_FUNCTION_CONSTANT +
        (Real) key2 * HASH2_FUNCTION_CONSTANT;

    index = (int) (FRACTION(v) * (Real) hash_table->size);

    return( index );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : lookup
@INPUT      : hash_table
            : key1
            : key2
@OUTPUT     : 
@RETURNS    : a pointer to a pointer to the hash table entry
@DESCRIPTION: Looks up the given key in the hash table, and returns a pointer
            : to a pointer to the entry, suitable for inserting or deleting.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   hash2_entry_struct  **lookup(
    hash2_table_struct  *hash_table,
    int                 key1,
    int                 key2 )
{
    int                  i;
    hash2_entry_struct   **ptr_to_entry;

    i = hash2_function( hash_table, key1, key2 );

    ptr_to_entry = &hash_table->table[i];

    while( *ptr_to_entry != NULL && ((*ptr_to_entry)->key1 != key1 ||
                                     (*ptr_to_entry)->key2 != key2) )
        ptr_to_entry = &(*ptr_to_entry)->next;

    return( ptr_to_entry );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : insert_in_hash2_table
@INPUT      : hash_table
            : key1
            : key2
            : data_ptr
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Inserts the data_ptr in the hash table.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  insert_in_hash2_table(
    hash2_table_struct  *hash_table,
    int                 key1,
    int                 key2,
    void                *data_ptr )
{
    hash2_entry_struct   **ptr_to_entry;
    hash2_entry_struct   *entry;

    ptr_to_entry = lookup( hash_table, key1, key2 );

    entry = *ptr_to_entry;

    if( entry == NULL )
    {
        ALLOC_VAR_SIZED_STRUCT( entry, char, hash_table->data_size );

        entry->key1 = key1;
        entry->key2 = key2;
        entry->next = *ptr_to_entry;
        (void) memcpy( (void *) entry->data, data_ptr,
                       (size_t) hash_table->data_size );

        *ptr_to_entry = entry; 

        ++hash_table->n_entries;

        if( (Real) hash_table->n_entries / (Real) hash_table->size >
            hash_table->enlarge_threshold )
        {
            int   new_size;

            new_size = ROUND( (Real) hash_table->n_entries /
                              hash_table->new_density );
            increase_hash2_table_size( hash_table, new_size );
        }
    }
    else
    {
        print_error( "Insert in hash2 table: entry already present: %d %d\n",
                     key1, key2 );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : lookup_in_hash2_table
@INPUT      : hash_table
            : key1
            : key2
@OUTPUT     : data_ptr
@RETURNS    : TRUE if found
@DESCRIPTION: Lookups the given key in the hash table, copying back the
              data if found.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  lookup_in_hash2_table(
    hash2_table_struct  *hash_table,
    int                 key1,
    int                 key2,
    void                *data_ptr )
{
    BOOLEAN              found;
    hash2_entry_struct   **ptr_to_entry;
    hash2_entry_struct   *entry;

    ptr_to_entry = lookup( hash_table, key1, key2 );

    entry = *ptr_to_entry;

    if( entry == NULL )
    {
        found = FALSE;
    }
    else
    {
        if( data_ptr != NULL )
        {
            (void) memcpy( data_ptr, (void *) entry->data,
                           (size_t) hash_table->data_size );
        }

        found = TRUE;
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_from_hash2_table
@INPUT      : hash_table
            : key1
            : key2
@OUTPUT     : : data_ptr
@RETURNS    : TRUE if it existed in the table
@DESCRIPTION: Finds and removes the given data from the hash table, copying
              the data back.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  remove_from_hash2_table(
    hash2_table_struct  *hash_table,
    int                 key1,
    int                 key2,
    void                *data_ptr )
{
    BOOLEAN              removed;
    hash2_entry_struct   **ptr_to_entry;
    hash2_entry_struct   *entry;

    ptr_to_entry = lookup( hash_table, key1, key2 );

    entry = *ptr_to_entry;

    if( entry == NULL )
    {
        removed = FALSE;
    }
    else
    {
        if( data_ptr != NULL )
        {
            (void) memcpy( data_ptr, (void *) entry->data,
                           (size_t) hash_table->data_size );
        }

        *ptr_to_entry = entry->next;

        FREE( entry );

        removed = TRUE;
        --hash_table->n_entries;
    }

    return( removed );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : move_hash2_entries_to_new_table
@INPUT      : src           hash_table
@OUTPUT     : dest          hash_table
@RETURNS    : 
@DESCRIPTION: Moves the entries of one hash table to the other, without
            : reallocating them.  Used when the hash table size is increased.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   void  move_hash2_entries_to_new_table(
    hash2_table_struct  *dest,
    hash2_table_struct  *src )
{
    int                  i, hash;
    hash2_entry_struct   *entry, *next;

    for( i = 0;  i < src->size;  ++i )
    {
        entry = src->table[i];

        while( entry != NULL )
        {
            next = entry->next;

            hash = hash2_function( dest, entry->key1, entry->key2 );
            entry->next = dest->table[hash];
            dest->table[hash] = entry;
            ++dest->n_entries;

            entry = next;
        }

        src->table[i] = NULL;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : increase_hash2_table_size
@INPUT      : hash_table
            : new_size
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Increases the size of the hash table, by creating a new hash
            : table, moving over the entries, then assigning the new table to
            : the old one.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void  increase_hash2_table_size(
    hash2_table_struct   *hash_table,
    int                  new_size )
{
    hash2_table_struct   new_table;

    initialize_hash2_table( &new_table, new_size, hash_table->data_size,
                            hash_table->enlarge_threshold,
                            hash_table->new_density );

    move_hash2_entries_to_new_table( &new_table, hash_table );

    delete_hash2_table_list( hash_table );

    *hash_table = new_table;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_hash2_pointer
@INPUT      : ptr
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes a pointer to the beginning of the hash table.  This is
            : used to step through the hash table to visit every entry.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_hash2_pointer(
    hash2_table_pointer  *ptr )
{
    ptr->current_index = -1;
    ptr->current_entry = NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_next_hash2_entry
@INPUT      : hash_table
            : ptr
@OUTPUT     : data_ptr
@RETURNS    : TRUE if one is found
@DESCRIPTION: Steps to the next entry in the hash table, updating the ptr,
            : and passing back the data ptr stored in that entry.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_next_hash2_entry(
    hash2_table_struct   *hash_table,
    hash2_table_pointer  *ptr,
    void                 *data_ptr )
{
    BOOLEAN   found;

    if( ptr->current_entry != NULL )
        ptr->current_entry = ptr->current_entry->next;

    if( ptr->current_entry == NULL )
    {
        do
        {
            ++ptr->current_index;
        }
        while( ptr->current_index < hash_table->size &&
               hash_table->table[ptr->current_index] == NULL );

        if( ptr->current_index < hash_table->size )
            ptr->current_entry = hash_table->table[ptr->current_index];
    }

    found = (ptr->current_entry != NULL);

    if( found && data_ptr != NULL )
    {
        (void) memcpy( data_ptr, (void *) ptr->current_entry->data,
                       (size_t) hash_table->data_size );
    }

    return( found );
}
