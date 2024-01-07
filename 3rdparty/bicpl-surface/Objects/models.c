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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/models.c,v 1.11 2005/08/17 22:28:27 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_model
@INPUT      : 
@OUTPUT     : model
@RETURNS    : 
@DESCRIPTION: Initializes the model to be an empty list of objects.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_model(
    model_struct  *model )
{
    model->filename = create_string( NULL );
    model->n_objects = 0;
    model->extra_ptr = (void *) NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_model
@INPUT      : model
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the model.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_model(
    model_struct   *model )
{
    int   i;

    delete_string( model->filename );

    for_less( i, 0, model->n_objects )
        delete_object( model->objects[i] );

    if( model->n_objects > 0 )
        FREE( model->objects );

    if( model->extra_ptr != (void *) NULL )
        FREE( model->extra_ptr );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : assign_model_extra_ptr
@INPUT      : model
              ptr
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Assigns the general-purpose void pointer associated with a model.
              This is where application-specific data may be stored.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  assign_model_extra_ptr(
    model_struct   *model,
    void           *ptr )
{
    model->extra_ptr = ptr;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_model_extra_ptr
@INPUT      : model
@OUTPUT     : 
@RETURNS    : ptr
@DESCRIPTION: Returns the general-purpose void pointer.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  *get_model_extra_ptr(
    model_struct   *model )
{
    return( model->extra_ptr );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_object_to_list
@INPUT      : n_objects
              object_list
              new_object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Adds an object to a list of objects.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  add_object_to_list(
    int            *n_objects,
    object_struct  ***object_list,
    object_struct  *new_object )
{
    ADD_ELEMENT_TO_ARRAY( *object_list, *n_objects, new_object, 1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_object_from_list
@INPUT      : n_objects
              object_list
              obj_index
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Removes the obj_index'th object from the list, sliding the
              remaining ones down.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  remove_object_from_list(
    int            *n_objects,
    object_struct  **object_list[],
    int            obj_index )
{
    DELETE_ELEMENT_FROM_ARRAY( *object_list, *n_objects, obj_index, 1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_object_list
@INPUT      : n_objects
              object_list
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the list of objects.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_object_list(
    int              n_objects,
    object_struct    *object_list[] )
{
    int      i;

    if( n_objects > 0 )
    {
        for_less( i, 0, n_objects )
        {
            delete_object( object_list[i] );
        }

        FREE( object_list );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_object_index_in_model
@INPUT      : model
              object
@OUTPUT     : 
@RETURNS    : index
@DESCRIPTION: Finds the position of the object in the model.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  find_object_index_in_model(
    model_struct   *model,
    object_struct  *object )
{
    int            obj_index;

    for_less( obj_index, 0, model->n_objects )
    {
        if( model->objects[obj_index] == object )
            return( obj_index );
    }

    return( -1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_object_to_model
@INPUT      : model
              new_object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Adds an object to a model.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  add_object_to_model(
    model_struct   *model,
    object_struct  *new_object )
{
    add_object_to_list( &model->n_objects, &model->objects, new_object );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : insert_object_in_model
@INPUT      : model
              new_object
              obj_index
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Inserts an object into a model at the specified location.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  insert_object_in_model(
    model_struct   *model,
    object_struct  *new_object,
    int            obj_index )
{
    int    i;

    SET_ARRAY_SIZE( model->objects, model->n_objects, model->n_objects+1, 1 );
    ++model->n_objects;

    for( i = model->n_objects-1;  i > obj_index;  --i )
        model->objects[i] = model->objects[i-1];
    model->objects[obj_index] = new_object;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_ith_object_from_model
@INPUT      : model
              obj_index
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Removes an object from a list, without deleting the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  remove_ith_object_from_model(
    model_struct   *model,
    int            obj_index )
{
    remove_object_from_list( &model->n_objects, &model->objects, obj_index );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_object_from_model
@INPUT      : model
              object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Removes the object from the model.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  remove_object_from_model(
    model_struct   *model,
    object_struct  *object )
{
    remove_ith_object_from_model( model,
                                  find_object_index_in_model(model,object) );
}
