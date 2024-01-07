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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/objects.c,v 1.24 2005/08/17 22:28:27 bert Exp $";
#endif

static  void  advance_object_traverse(
    object_traverse_struct  *object_traverse );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_object
@INPUT      : object_type
@OUTPUT     : 
@RETURNS    : an object
@DESCRIPTION: Creates an object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  object_struct  *create_object(
    Object_types   object_type )
{
    object_struct   *object;

    ALLOC( object, 1 );

    object->object_type = object_type;
    object->visibility = TRUE;

    return( object );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_type
@INPUT      : object
@OUTPUT     : 
@RETURNS    : type
@DESCRIPTION: Returns the type of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Object_types  get_object_type(
    object_struct   *object )
{
    return( object->object_type );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_visibility
@INPUT      : object
@OUTPUT     : 
@RETURNS    : visibility
@DESCRIPTION: Returns the visibility of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_object_visibility(
    object_struct  *object )
{
    return( object->visibility );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_object_visibility
@INPUT      : object
              visibility
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the visibility of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_object_visibility(
    object_struct  *object,
    BOOLEAN        visibility )
{
    object->visibility = visibility;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_lines_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : lines
@DESCRIPTION: Returns a pointer to the lines of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  lines_struct  *get_lines_ptr(
    object_struct  *object )
{
    if( object->object_type != LINES )
    {
        handle_internal_error( "get_lines_ptr" );
        return( (lines_struct *) NULL );
    }
    else
        return( &object->specific.lines );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_marker_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : marker
@DESCRIPTION: Returns a pointer to the marker of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  marker_struct  *get_marker_ptr(
    object_struct  *object )
{
    if( object->object_type != MARKER )
    {
        handle_internal_error( "get_marker_ptr" );
        return( (marker_struct *) NULL );
    }
    else
        return( &object->specific.marker );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_model_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : model
@DESCRIPTION: Returns a pointer to the model of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  model_struct  *get_model_ptr(
    object_struct  *object )
{
    if( object->object_type != MODEL )
    {
        handle_internal_error( "get_model_ptr" );
        return( (model_struct *) NULL );
    }
    else
        return( &object->specific.model );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_pixels_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : pixels
@DESCRIPTION: Returns a pointer to the pixels of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  pixels_struct  *get_pixels_ptr(
    object_struct  *object )
{
    if( object->object_type != PIXELS )
    {
        handle_internal_error( "get_pixels_ptr" );
        return( (pixels_struct *) NULL );
    }
    else
        return( &object->specific.pixels );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : polygons
@DESCRIPTION: Returns a pointer to the polygons of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  polygons_struct  *get_polygons_ptr(
    object_struct  *object )
{
    if( object->object_type != POLYGONS )
    {
        handle_internal_error( "get_polygons_ptr" );
        return( (polygons_struct *) NULL );
    }
    else
        return( &object->specific.polygons );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : quadmesh
@DESCRIPTION: Returns a pointer to the quadmesh of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  quadmesh_struct  *get_quadmesh_ptr(
    object_struct  *object )
{
    if( object->object_type != QUADMESH )
    {
        handle_internal_error( "get_quadmesh_ptr" );
        return( (quadmesh_struct *) NULL );
    }
    else
        return( &object->specific.quadmesh );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_text_ptr
@INPUT      : object
@OUTPUT     : 
@RETURNS    : text
@DESCRIPTION: Returns a pointer to the text of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  text_struct  *get_text_ptr(
    object_struct  *object )
{
    if( object->object_type != TEXT )
    {
        handle_internal_error( "get_text_ptr" );
        return( (text_struct *) NULL );
    }
    else
        return( &object->specific.text );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_lines_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : number of points
@DESCRIPTION: Passes back the points of the object and returns the number
              of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_lines_points(
    object_struct   *object,
    Point           *points[] )
{
    *points = get_lines_ptr(object)->points;

    return( get_lines_ptr(object)->n_points );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_marker_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : number of points
@DESCRIPTION: Passes back the points of the object and returns the number
              of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_marker_points(
    object_struct   *object,
    Point           *points[] )
{
    *points = &get_marker_ptr(object)->position;

    return( 1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_zero_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : number of points
@DESCRIPTION: Passes back the points of objects with no points and returns
              the number of points (0).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  int  get_object_zero_points(
    object_struct   *object,
    Point           *points[] )
{
    return( 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : number of points
@DESCRIPTION: Passes back the points of the object and returns the number
              of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_polygons_points(
    object_struct   *object,
    Point           *points[] )
{
    *points = get_polygons_ptr(object)->points;

    return( get_polygons_ptr(object)->n_points );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : number of points
@DESCRIPTION: Passes back the points of the object and returns the number
              of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_quadmesh_points(
    object_struct   *object,
    Point           *points[] )
{
    *points = get_quadmesh_ptr(object)->points;

    return( get_quadmesh_ptr(object)->m * get_quadmesh_ptr(object)->n );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_text_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : number of points
@DESCRIPTION: Passes back the points of the object and returns the number
              of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_text_points(
    object_struct   *object,
    Point           *points[] )
{
    *points = &get_text_ptr(object)->origin;

    return( 1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_zero_normals
@INPUT      : object
@OUTPUT     : normals
@RETURNS    : number of normals
@DESCRIPTION: Passes back the normals of the object and returns the number
              of normals.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  int  get_object_zero_normals(
    object_struct   *object,
    Vector          *normals[] )
{
    *normals = NULL;
    return( 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_normals
@INPUT      : object
@OUTPUT     : normals
@RETURNS    : number of normals
@DESCRIPTION: Passes back the normals of the object and returns the number
              of normals.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_polygons_normals(
    object_struct   *object,
    Vector          *normals[] )
{
    *normals = get_polygons_ptr(object)->normals;

    return( get_polygons_ptr(object)->n_points );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_normals
@INPUT      : object
@OUTPUT     : normals
@RETURNS    : number of normals
@DESCRIPTION: Passes back the normals of the object and returns the number
              of normals.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_quadmesh_normals(
    object_struct   *object,
    Vector          *normals[] )
{
    *normals = get_quadmesh_ptr(object)->normals;

    return( get_quadmesh_ptr(object)->m * get_quadmesh_ptr(object)->n );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_lines_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour flag
@DESCRIPTION: Passes back the colours of the object and returns the colour
              flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Colour_flags  *get_lines_colours(
    object_struct   *object,
    Colour          *colours[] )
{
    *colours = get_lines_ptr(object)->colours;

    return( &get_lines_ptr(object)->colour_flag );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_marker_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour flag
@DESCRIPTION: Passes back the colours of the object and returns the colour
              flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Colour_flags  *get_marker_colours(
    object_struct    *object,
    Colour           *colours[] )
{
    *colours = &get_marker_ptr(object)->colour;

    return( (Colour_flags *) 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_zero_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour flag
@DESCRIPTION: Passes back the colours of the object and returns the colour
              flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  Colour_flags  *get_object_zero_colours(
    object_struct   *object,
    Colour          *colours[] )
{
    *colours = (Colour *) NULL;

    return( (Colour_flags *) 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour flag
@DESCRIPTION: Passes back the colours of the object and returns the colour
              flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Colour_flags  *get_polygons_colours(
    object_struct   *object,
    Colour          *colours[] )
{
    *colours = get_polygons_ptr(object)->colours;

    return( &get_polygons_ptr(object)->colour_flag );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour flag
@DESCRIPTION: Passes back the colours of the object and returns the colour
              flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Colour_flags  *get_quadmesh_colours(
    object_struct   *object,
    Colour          *colours[] )
{
    *colours = get_quadmesh_ptr(object)->colours;

    return( &get_quadmesh_ptr(object)->colour_flag );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_text_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour flag
@DESCRIPTION: Passes back the colours of the object and returns the colour
              flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Colour_flags  *get_text_colours(
    object_struct   *object,
    Colour          *colours[] )
{
    *colours = &get_text_ptr(object)->colour;

    return( (Colour_flags *) 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_lines_colours
@INPUT      : object
              colours
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Assigns the colours array of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  set_lines_colours(
    object_struct   *object,
    Colour          colours[] )
{
    get_lines_ptr(object)->colours = colours;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_object_no_colours
@INPUT      : object
              colours
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Assigns the colours array of the object, for objects which
              have no colours.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  void  set_object_no_colours(
    object_struct    *object,
    Colour           colours[] )
{
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_polygons_colours
@INPUT      : object
              colours
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Assigns the colours array of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  set_polygons_colours(
    object_struct   *object,
    Colour          colours[] )
{
    get_polygons_ptr(object)->colours = colours;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_quadmesh_colours
@INPUT      : object
              colours
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Assigns the colours array of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  set_quadmesh_colours(
    object_struct   *object,
    Colour          colours[] )
{
    get_quadmesh_ptr(object)->colours = colours;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_no_surfprop
@INPUT      : object
@OUTPUT     : 
@RETURNS    : pointer to surfprop
@DESCRIPTION: Returns a null pointer for objects that do not have surfprops.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  Surfprop  *get_object_no_surfprop(
    object_struct   *object )
{
    return( (Surfprop *) NULL );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_surfprop
@INPUT      : object
@OUTPUT     : 
@RETURNS    : surfprop
@DESCRIPTION: Returns a pointer to the object's surfprop.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Surfprop  *get_polygons_surfprop(
    object_struct   *object )
{
    return( &get_polygons_ptr(object)->surfprop );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_surfprop
@INPUT      : object
@OUTPUT     : 
@RETURNS    : surfprop
@DESCRIPTION: Returns a pointer to the object's surfprop.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Surfprop  *get_quadmesh_surfprop(
    object_struct   *object )
{
    return( &get_quadmesh_ptr(object)->surfprop );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_lines_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_lines_name(
    object_struct   *object )
{
    char   buffer[EXTREMELY_LARGE_STRING_SIZE];

    (void) sprintf( buffer, "Lines (%d:%d)",
                    get_lines_ptr(object)->n_items,
                    get_lines_ptr(object)->n_points );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_marker_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_marker_name(
    object_struct   *object )
{
    STRING   label;
    char     buffer[EXTREMELY_LARGE_STRING_SIZE];

    if( string_length(get_marker_ptr(object)->label) == 0 )
        label = "Marker";
    else
        label = get_marker_ptr(object)->label;

    (void) sprintf( buffer, "%s:%d:%d: (%g %g %g)",
                    label,
                    get_marker_ptr(object)->patient_id,
                    get_marker_ptr(object)->structure_id,
                    Point_x(get_marker_ptr(object)->position),
                    Point_y(get_marker_ptr(object)->position),
                    Point_z(get_marker_ptr(object)->position) );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_model_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_model_name(
    object_struct   *object )
{
    char   buffer[EXTREMELY_LARGE_STRING_SIZE];

    (void) sprintf( buffer, "Model (%s) ->", get_model_ptr(object)->filename );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_pixels_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_pixels_name(
    object_struct   *object )
{
    char   buffer[EXTREMELY_LARGE_STRING_SIZE];

    (void) sprintf( buffer, "Pixels (%d by %d)",
                    get_pixels_ptr(object)->x_size,
                    get_pixels_ptr(object)->y_size );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_polygons_name(
    object_struct   *object )
{
    char   buffer[EXTREMELY_LARGE_STRING_SIZE];

    (void) sprintf( buffer, "Polygons (%d:%d)",
                    get_polygons_ptr(object)->n_items,
                    get_polygons_ptr(object)->n_points );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_quadmesh_name(
    object_struct   *object )
{
    char   buffer[EXTREMELY_LARGE_STRING_SIZE];

    (void) sprintf( buffer, "Quadmesh [%d][%d]",
                    get_quadmesh_ptr(object)->m, get_quadmesh_ptr(object)->n );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_text_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets a name description for the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_text_name(
    object_struct   *object )
{
    char   buffer[EXTREMELY_LARGE_STRING_SIZE];

    (void) sprintf( buffer, "Text (%s)", get_text_ptr(object)->string );

    return( create_string( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_lines_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  delete_lines_object(
    object_struct   *object )
{
    delete_lines( get_lines_ptr(object) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_marker_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  void  delete_marker_object(
    object_struct   *object )
{
    delete_marker( get_marker_ptr(object) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_model_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  delete_model_object(
    object_struct   *object )
{
    delete_model( get_model_ptr(object) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_pixels_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  delete_pixels_object(
    object_struct   *object )
{
    delete_pixels( get_pixels_ptr(object) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_polygons_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  delete_polygons_object(
    object_struct   *object )
{
    delete_polygons( get_polygons_ptr(object) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_quadmesh_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  void  delete_quadmesh_object(
    object_struct   *object )
{
    delete_quadmesh( get_quadmesh_ptr(object) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_text_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  void  delete_text_object(
    object_struct   *object )
{
    delete_text( get_text_ptr(object) );
}

/* -------------------- function lookup table ------------------- */

typedef  struct
{
    int           (*get_points_function)( object_struct *, Point *[] );
    int           (*get_normals_function)( object_struct *, Vector *[] );
    Colour_flags  *(*get_colours_function)( object_struct *, Colour *[] );
    void          (*set_colours_function)( object_struct *, Colour [] );
    Surfprop      *(*get_surfprop_function)( object_struct * );
    STRING        (*get_name_function)( object_struct * );
    void          (*delete_function)( object_struct * );
}
object_functions_list;

static  object_functions_list   object_functions[N_OBJECT_TYPES] =
{
    {
        get_lines_points,
        get_object_zero_normals,
        get_lines_colours,
        set_lines_colours,
        get_object_no_surfprop,
        get_lines_name,
        delete_lines_object
    },                              /* LINES */
    {
        get_marker_points,
        get_object_zero_normals,
        get_marker_colours,
        set_object_no_colours,
        get_object_no_surfprop,
        get_marker_name,
        delete_marker_object
    },                              /* MARKER */
    {
        get_object_zero_points,
        get_object_zero_normals,
        get_object_zero_colours,
        set_object_no_colours,
        get_object_no_surfprop,
        get_model_name,
        delete_model_object
    },                              /* MODEL */
    {
        get_object_zero_points,
        get_object_zero_normals,
        get_object_zero_colours,
        set_object_no_colours,
        get_object_no_surfprop,
        get_pixels_name,
        delete_pixels_object
    },                              /* PIXELS */
    {
        get_polygons_points,
        get_polygons_normals,
        get_polygons_colours,
        set_polygons_colours,
        get_polygons_surfprop,
        get_polygons_name,
        delete_polygons_object
    },                              /* POLYGONS */
    {
        get_quadmesh_points,
        get_quadmesh_normals,
        get_quadmesh_colours,
        set_quadmesh_colours,
        get_quadmesh_surfprop,
        get_quadmesh_name,
        delete_quadmesh_object
    },                              /* QUADMESH */
    {
        get_text_points,
        get_object_zero_normals,
        get_text_colours,
        set_object_no_colours,
        get_object_no_surfprop,
        get_text_name,
        delete_text_object
    }                               /* TEXT */
};

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_object
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_object(
    object_struct  *object )
{
    object_functions[object->object_type].delete_function( object );
    FREE( object );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_points
@INPUT      : object
@OUTPUT     : points
@RETURNS    : n_points
@DESCRIPTION: Passes back a pointer to the points of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_object_points(
    object_struct   *object,
    Point           *points[] )
{
    return( object_functions[object->object_type].get_points_function
                   ( object, points ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_normals
@INPUT      : object
@OUTPUT     : normals
@RETURNS    : n_normals
@DESCRIPTION: Passes back a pointer to the normals of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_object_normals(
    object_struct   *object,
    Vector          *normals[] )
{
    return( object_functions[object->object_type].get_normals_function
                   ( object, normals ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_colours
@INPUT      : object
@OUTPUT     : colours
@RETURNS    : colour_flag
@DESCRIPTION: Passes back a pointer to the colours array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Colour_flags  *get_object_colours(
    object_struct   *object,
    Colour          *colours[] )
{
    return( object_functions[object->object_type].get_colours_function
                   ( object, colours ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_object_colours
@INPUT      : object
              colours
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the colours array of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_object_colours(
    object_struct   *object,
    Colour          colours[] )
{
    object_functions[object->object_type].set_colours_function(object,colours);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_surfprop
@INPUT      : object
@OUTPUT     : 
@RETURNS    : surfprop
@DESCRIPTION: Returns a pointer to the surfprop of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Surfprop  *get_object_surfprop(
    object_struct   *object )
{
    return( object_functions[object->object_type].get_surfprop_function
                   ( object ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_name
@INPUT      : object
@OUTPUT     : name
@RETURNS    : 
@DESCRIPTION: Gets the name description of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  STRING  get_object_name(
    object_struct  *object )
{
    return( object_functions[object->object_type].get_name_function( object ) );
}

static  void   pop_object_stack(
    object_traverse_struct  *object_traverse )
{
    if( object_traverse->top_of_stack > 0 )
        --object_traverse->top_of_stack;
    else
        handle_internal_error( "pop_object_stack" );
}

static  void  push_object_stack(
    object_traverse_struct  *object_traverse,
    object_stack_struct     *entry )
{
    if( object_traverse->top_of_stack+1 >= object_traverse->n_stack_alloced )
    {
        SET_ARRAY_SIZE( object_traverse->alloced_stack,
                        object_traverse->top_of_stack,
                        object_traverse->top_of_stack+1,
                        DEFAULT_CHUNK_SIZE );

        object_traverse->stack = object_traverse->alloced_stack;
    }

    object_traverse->stack[object_traverse->top_of_stack] = *entry;
    ++object_traverse->top_of_stack;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_object_traverse
@INPUT      : object_traverse
              n_objects
              object_list
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes a structure for traversing an object, which could
              potentially have a hierarchy, if models are present.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_object_traverse(
    object_traverse_struct  *object_traverse,
    BOOLEAN                 visible_ones_only,
    int                     n_objects,
    object_struct           *object_list[] )
{
    object_stack_struct   push_entry;

    object_traverse->visible_ones_only = visible_ones_only;
    object_traverse->n_stack_alloced = MAX_OBJECT_TRAVERSE;
    object_traverse->top_of_stack = 0;
    object_traverse->stack = object_traverse->static_stack;

    if( n_objects > 0 )
    {
        push_entry.index = 0;
        push_entry.n_objects = n_objects;
        push_entry.object_list = object_list;

        while( push_entry.index < n_objects &&
               visible_ones_only &&
               !get_object_visibility( object_list[push_entry.index] ) )
            ++push_entry.index;

        if( push_entry.index < n_objects )
            push_object_stack( object_traverse, &push_entry );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_next_object_traverse
@INPUT      : object_traverse
@OUTPUT     : object
@RETURNS    : TRUE if object found
@DESCRIPTION: Gets the next object from the traversal.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_next_object_traverse(
    object_traverse_struct  *object_traverse,
    object_struct           **object )
{
    BOOLEAN               object_found;
    object_stack_struct   *top_entry;

    if( object_traverse->top_of_stack > 0 )
    {
        top_entry = &object_traverse->stack[object_traverse->top_of_stack-1];
        *object = top_entry->object_list[top_entry->index];

        advance_object_traverse( object_traverse );
        object_found = TRUE;
    }
    else
    {
        terminate_object_traverse( object_traverse );
        object_found = FALSE;
    }

    return( object_found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : advance_object_traverse
@INPUT      : object_traverse
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Advances the object_traverse struct to point to the next object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  advance_object_traverse(
    object_traverse_struct  *object_traverse )
{
    object_stack_struct   *top_entry, push_entry;
    object_struct         *object;
    model_struct          *model;

    top_entry = &object_traverse->stack[object_traverse->top_of_stack-1];
    object = top_entry->object_list[top_entry->index];
    ++top_entry->index;

    while( top_entry->index < top_entry->n_objects &&
           object_traverse->visible_ones_only &&
           !get_object_visibility( top_entry->object_list[top_entry->index]) )
    {
        ++top_entry->index;
    }

    if( object->object_type == MODEL )
    {
        model = get_model_ptr( object );
        push_entry.index = 0;
        push_entry.n_objects = model->n_objects;
        push_entry.object_list = model->objects;

        while( push_entry.index < push_entry.n_objects &&
               object_traverse->visible_ones_only &&
               !get_object_visibility(push_entry.object_list[push_entry.index]))
            ++push_entry.index;

        if( push_entry.index < push_entry.n_objects )
            push_object_stack( object_traverse, &push_entry );
    }

    while( object_traverse->top_of_stack > 0 &&
           object_traverse->stack[object_traverse->top_of_stack-1].index >=
           object_traverse->stack[object_traverse->top_of_stack-1].n_objects )
    {
        pop_object_stack( object_traverse );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : terminate_object_traverse
@INPUT      : object_traverse
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the object traversal structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  terminate_object_traverse(
    object_traverse_struct   *object_traverse )
{
    if( object_traverse->n_stack_alloced > MAX_OBJECT_TRAVERSE )
        FREE( object_traverse->alloced_stack );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_range_of_object
@INPUT      : object
              visible_ones_only
@OUTPUT     : min_corner
              max_corner
@RETURNS    : TRUE if range found.
@DESCRIPTION: Returns the range of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_range_of_object(
    object_struct   *object,
    BOOLEAN         visible_ones_only,
    Point           *min_corner,
    Point           *max_corner )
{
    BOOLEAN                 found_flag;
    Point                   min_obj, max_obj, *points;
    int                     n_points;
    object_struct           *current_object;
    object_traverse_struct  object_traverse;

    found_flag = FALSE;

    initialize_object_traverse( &object_traverse, visible_ones_only, 1,&object);

    while( get_next_object_traverse( &object_traverse, &current_object ) )
    {
        n_points = get_object_points( current_object, &points );

        if( n_points > 0 )
        {
            get_range_points( n_points, points, &min_obj, &max_obj );

            if( !found_flag )
            {
                found_flag = TRUE;
                *min_corner = min_obj;
                *max_corner = max_obj;
            }
            else
            {
                expand_min_and_max_points( min_corner, max_corner,
                                           &min_obj, &max_obj );
            }
        }
    }

    return( found_flag );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : reverse_object_normals
@INPUT      : object
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Reverses the normals of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  reverse_object_normals(
    object_struct   *object )
{
    int                     n_normals;
    Vector                  *normals;
    object_struct           *current_object;
    object_traverse_struct  object_traverse;

    initialize_object_traverse( &object_traverse, FALSE, 1, &object );

    while( get_next_object_traverse( &object_traverse, &current_object ) )
    {
        n_normals = get_object_normals( current_object, &normals );
        reverse_vectors( n_normals, normals );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_object_colour
@INPUT      : object
@OUTPUT     : colour
@RETURNS    : TRUE if object has single colour
@DESCRIPTION: Returns the global colour of the object, if any. 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_object_colour(
    object_struct  *object,
    Colour         *colour )
{
    BOOLEAN       has_single_colour;
    Colour_flags  *colour_flag;
    Colour        *colours;

    colour_flag = get_object_colours( object, &colours );

    if( (colour_flag == (Colour_flags *) NULL ||
         *colour_flag == ONE_COLOUR) && colours != (Colour *) NULL )
    {
        has_single_colour = TRUE;
        *colour = colours[0];
    }
    else
    {
        has_single_colour = FALSE;
    }

    return( has_single_colour );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_object_colour
@INPUT      : object
              col
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the colour of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_object_colour(
    object_struct  *object,
    Colour         col )
{
    Colour_flags  *colour_flag;
    Colour        *colours;

    colour_flag = get_object_colours( object, &colours );

    if( colour_flag != (Colour_flags *) NULL && *colour_flag != ONE_COLOUR )
    {
        *colour_flag = ONE_COLOUR;

        if( colours != (Colour *) NULL )
            FREE( colours );

        ALLOC( colours, 1 );
        set_object_colours( object, colours );
    }

    if( colours != (Colour *) NULL )
        colours[0] = col;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_object_surfprop
@INPUT      : object
              spr
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the surfprop of the object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_object_surfprop(
    object_struct  *object,
    Surfprop       *spr )
{
    Surfprop   *object_spr;

    object_spr = get_object_surfprop( object );

    if( object_spr != (Surfprop *) NULL )
        *object_spr = *spr;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : free_colours
@INPUT      : colour_flag
              colours
              n_points
              n_items
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Frees the colours.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  free_colours(
    Colour_flags   colour_flag,
    Colour         colours[],
    int            n_points,
    int            n_items )
{
    BOOLEAN  should_free;

    should_free = FALSE;

    switch( colour_flag )
    {
    case  ONE_COLOUR:
        should_free = TRUE;
        break;

    case  PER_ITEM_COLOURS:
        should_free = (n_items > 0);
        break;

    case  PER_VERTEX_COLOURS:
        should_free = (n_points > 0);
        break;
    }

    if( should_free )
        FREE( colours );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : objects_are_same_topology
@INPUT      : n_points1
              n_items1
              end_indices1
              indices1
              n_points2
              n_items2
              end_indices2
              indices2
@OUTPUT     : 
@RETURNS    : TRUE if same topology
@DESCRIPTION: Tests if two objects are the same topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  objects_are_same_topology(
    int    n_points1,
    int    n_items1,
    int    *end_indices1,
    int    *indices1,
    int    n_points2,
    int    n_items2,
    int    *end_indices2,
    int    *indices2 )
{
    int   p;

    if( n_points1 != n_points2 || n_items1 != n_items2 )
        return( FALSE );

    for_less( p, 0, n_items1 )
        if( end_indices1[p] != end_indices2[p] )
            return( FALSE );

    for_less( p, 0, end_indices1[n_items1-1] )
        if( indices1[p] != indices2[p] )
            return( FALSE );

    return( TRUE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_n_colours
@INPUT      : colour_flag
              n_points
              n_items
@OUTPUT     : 
@RETURNS    : number of colours
@DESCRIPTION: Returns the number of colours corresponding to the colour_flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_n_colours(
    Colour_flags  colour_flag,
    int           n_points,
    int           n_items )
{
    int   n_colours;

    switch( colour_flag )
    {
    case ONE_COLOUR:
        n_colours = 1;
        break;

    case PER_ITEM_COLOURS:
        n_colours = n_items;
        break;

    case PER_VERTEX_COLOURS:
        n_colours = n_points;
        break;

    default:
        print_error( "invalid colour flag\n" );
        n_colours = 0;
        break;
    }

    return( n_colours );
}
