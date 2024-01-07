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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/graphics_io.c,v 1.9 2005/08/17 22:28:26 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_graphics_file
@INPUT      : filename
@OUTPUT     : format
              n_objects
              object_list
@RETURNS    : OK or ERROR
@DESCRIPTION: Inputs a file of graphics objects.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status   input_graphics_file(
    STRING         filename,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list )
{
    Status         status;
    FILE           *file;
    BOOLEAN        eof;
    object_struct  *object;
    STRING         current_directory;

    status = open_file_with_default_suffix( filename, "obj", READ_FILE,
                                            BINARY_FORMAT, &file );

    *n_objects = 0;

    if( status == OK )
    {
        current_directory = extract_directory( filename );

        do
        {
            status = input_object( current_directory, file, format,
                                   &object, &eof );

            if( status == OK && !eof )
                add_object_to_list( n_objects, object_list, object );

        } while( status == OK && !eof );

        delete_string( current_directory );
    }

    if( status == OK )
        status = close_file( file );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_graphics_file
@INPUT      : filename
              format
              n_objects
              object_list
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Writes a file of graphics objects.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status   output_graphics_file(
    STRING         filename,
    File_formats   format,
    int            n_objects,
    object_struct  *object_list[] )
{
    Status         status;
    int            i;
    FILE           *file;

    status = open_file_with_default_suffix( filename, "obj", WRITE_FILE,
                                            BINARY_FORMAT, &file );

    if( status == OK )
    {
        for_less( i, 0, n_objects )
        {
            if( status == OK )
                status = output_object( file, format, object_list[i] );
        }
    }

    if( status == OK )
        status = close_file( file );

    return( status );
}

