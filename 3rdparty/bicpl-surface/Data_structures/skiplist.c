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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Data_structures/skiplist.c,v 1.13 2005/08/17 22:31:13 bert Exp $";
#endif

static  int  get_random_level( void );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : skip_list.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Routines for managing an ordered list, using the skiplist
            : technique.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  MAX_SKIP_LEVELS   50
#define  SKIP_P            0.5

typedef  struct
{
    skip_struct   *update[MAX_SKIP_LEVELS];
} update_struct;

#define  ALLOC_SKIP_STRUCT( ptr, n_level )                                    \
         ALLOC_VAR_SIZED_STRUCT( ptr, skip_struct *, n_level )

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_skiplist
@INPUT      : skiplist
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes the skiplist to empty.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void  initialize_skiplist(
    skiplist_struct  *skiplist )
{
    int       i;

    skiplist->level = 1;

    ALLOC_SKIP_STRUCT( skiplist->header, MAX_SKIP_LEVELS );

    for_less( i, 0, MAX_SKIP_LEVELS )
        skiplist->header->forward[i] = NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_data_position
@INPUT      : skiplist
            : key
@OUTPUT     : update
@RETURNS    : TRUE if found
@DESCRIPTION: Searches the skiplist for the given key, and sets the update
            : struct so that it can provide an insert.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  BOOLEAN  find_data_position(
    skiplist_struct    *skiplist,
    float              key,
    update_struct      *update )
{
    int           i;
    skip_struct   *x;
    BOOLEAN       found;

    x = skiplist->header;

    for( i = skiplist->level-1;  i >= 0;  --i )
    {
        while( x->forward[i] != NULL && x->forward[i]->key < key )
            x = x->forward[i];

        update->update[i] = x;
    }

    x = update->update[0]->forward[0];

    found = x != NULL && x->key == key;

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : insert_data_in_skiplist
@INPUT      : skiplist
            : update           - the set of pointers indicating where to insert
            : key              - search key
            : data_ptr         - the data to insert
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Inserts the data ptr in the skiplist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   void  insert_data_in_skiplist(
    skiplist_struct   *skiplist,
    update_struct     *update,
    float             key,
    void              *data_ptr )
{
    int           i, new_level;
    skip_struct   *x;
#ifdef DEBUG
    int       prev_size = -1;
    void      test_skiplist_integrity();
#endif

    new_level = get_random_level();

    if( new_level > skiplist->level )
    {
        for( i = skiplist->level;  i < new_level;  ++i )
            update->update[i] = skiplist->header;

#ifdef DEBUG
        prev_size = skiplist->level;
#endif
        skiplist->level = new_level;
    }

    ALLOC_SKIP_STRUCT( x, new_level );

    x->data_ptr = data_ptr;
    x->key = key;

    for( i = 0;  i < new_level;  ++i )
    {
        x->forward[i] = update->update[i]->forward[i];
        update->update[i]->forward[i] = x;
    }

#ifdef DEBUG
    test_skiplist_integrity( skiplist );
#endif
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_entry_from_skiplist
@INPUT      : skiplist
            : update
@OUTPUT     : 
@RETURNS    : TRUE if it existed
@DESCRIPTION: Deletes the entry in the skip list pointed to by update.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   void  delete_entry_from_skiplist(
    skiplist_struct   *skiplist,
    update_struct     *update )
{
    int           i;
    skip_struct   *x;

    x = update->update[0]->forward[0];

    for( i = 0;  i < skiplist->level;  ++i )
    {
        if( update->update[i]->forward[i] != x )
            break;
        update->update[i]->forward[i] = x->forward[i];
    }

    FREE( x );

    while( skiplist->level > 1 &&
           skiplist->header->forward[skiplist->level-1] == NULL )
    {
        --skiplist->level;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_random_level
@INPUT      : 
@OUTPUT     : 
@RETURNS    : a random level between 1 and MAX_LEVELS
@DESCRIPTION: Determines a random level with exponential probability of higher
            : levels.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_random_level( void )
{
    int    level;

    level = 1;

    while( get_random_0_to_1() < SKIP_P && level < MAX_SKIP_LEVELS-1 )
        ++level;

    return( level );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_skiplist
@INPUT      : skiplist
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Deletes the memory associated with the skiplist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void  delete_skiplist(
    skiplist_struct  *skiplist )
{
    skip_struct   *ptr, *deleting;

    ptr = skiplist->header;

    while( ptr != NULL )
    {
        deleting = ptr;
        ptr = ptr->forward[0];
        FREE( deleting );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : search_skiplist
@INPUT      : skiplist
            : key
@OUTPUT     : data_ptr
@RETURNS    : TRUE if found, FALSE if not
@DESCRIPTION: Searches the skiplist for a data record that matches the
            : key.  This function will find the matching record,
            : and return it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  search_skiplist(
    skiplist_struct          *skiplist,
    float                    key,
    void                     **data_ptr )
{
    BOOLEAN        found;
    update_struct  update_ptrs;

    found = find_data_position( skiplist, key, &update_ptrs );

    if( found )
        *data_ptr = update_ptrs.update[0]->forward[0]->data_ptr;

    return( found );
}

BICAPI  BOOLEAN  search_skiplist_and_return_pointer(
    skiplist_struct          *skiplist,
    float                    key,
    skip_struct              **entry_ptr,
    void                     **data_ptr )
{
    BOOLEAN        found;
    update_struct  update_ptrs;

    found = find_data_position( skiplist, key, &update_ptrs );

    if( found )
    {
        *entry_ptr = update_ptrs.update[0]->forward[0];
        *data_ptr = (*entry_ptr)->data_ptr;
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : insert_in_skiplist
@INPUT      : skiplist
            : key
            : data_ptr
@OUTPUT     : 
@RETURNS    : TRUE if inserted, FALSE if already there
@DESCRIPTION: Inserts an entry into a skiplist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  insert_in_skiplist(
    skiplist_struct          *skiplist,
    float                    key,
    void                     *data_ptr )
{
    BOOLEAN        already_there;
    update_struct  update_ptrs;

    already_there = find_data_position( skiplist, key, &update_ptrs );

    if( !already_there )
        insert_data_in_skiplist( skiplist, &update_ptrs, key, data_ptr );

    return( !already_there );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_from_skiplist
@INPUT      : skiplist
            : key
@OUTPUT     : data_ptr
@RETURNS    : TRUE if key was in list and successfully deleted.
@DESCRIPTION: Deletes the entry for the given data_ptr from the skiplist,
            : passing back the pointer that was stored.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  delete_from_skiplist(
    skiplist_struct  *skiplist,
    float            key,
    void             **data_ptr )
{
    BOOLEAN        in_skiplist;
    update_struct  update_ptrs;

    in_skiplist = find_data_position( skiplist, key, &update_ptrs );

    if( in_skiplist )
    {
        *data_ptr = update_ptrs.update[0]->forward[0]->data_ptr;

        delete_entry_from_skiplist( skiplist, &update_ptrs );
    }

    return( in_skiplist );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_first_skiplist_entry
@INPUT      : skiplist
@OUTPUT     : entry_ptr
            : key
            : data_ptr
@RETURNS    : True if an entry exists in the skiplist.
@DESCRIPTION: Initializes the entry_ptr to point to the first entry,
            : and the key and data_ptr to the data for the first entry.
            : This function is used to traverse the list sequentially.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_first_skiplist_entry(
    skiplist_struct   *skiplist,
    skip_struct       **entry_ptr,
    float             *key,
    void              **data_ptr )
{
    *entry_ptr = skiplist->header->forward[0];

    if( *entry_ptr != NULL )
    {
        *key = (*entry_ptr)->key;
        *data_ptr = (*entry_ptr)->data_ptr;
    }

    return( *entry_ptr != NULL );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_next_skiplist_entry
@INPUT      : entry_ptr
@OUTPUT     : key
            : data_ptr
@RETURNS    : True if another entry is found.
@DESCRIPTION: Advances the entry_ptr to the next entry, and copies that entry's
            : key and data ptr to the function arguments.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_next_skiplist_entry(
    skip_struct       **entry_ptr,
    float             *key,
    void              **data_ptr )
{
    *entry_ptr = (*entry_ptr)->forward[0];

    if( *entry_ptr != NULL )
    {
        *key = (*entry_ptr)->key;
        *data_ptr = (*entry_ptr)->data_ptr;
    }

    return( *entry_ptr != NULL );
}

#ifdef DEBUG

static  BOOLEAN  test_skiplist_integrity(
    skiplist_struct  *skiplist )
{
    int            i;
    update_struct  update;
    BOOLEAN        okay;

    for_less( i, 0, skiplist->level )
        update.update[i] = skiplist->header->forward[i];

    while( update.update[0] != NULL )
    {
        i = 1;
        while( i < skiplist->level && update.update[0] == update.update[i] )
        {
            update.update[i] = update.update[i]->forward[i];
            ++i;
        }
        update.update[0] = update.update[0]->forward[0];
    }

    okay = TRUE;

    for_less( i, 0, skiplist->level )
        if( update.update[i] != NULL )
            okay = FALSE;

    if( !okay )
        handle_internal_error( "Skiplist integrity" );

    return( okay );
}
#endif
