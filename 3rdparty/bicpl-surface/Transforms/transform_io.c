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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Transforms/transform_io.c,v 1.9 2005/08/17 22:26:48 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_transform_file
@INPUT      : filename
              comments
              transform
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Writes the specified transform to the file, in MNI format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  write_transform_file(
    STRING       filename,
    STRING       comments,
    Transform    *transform )
{
    Status              status;
    FILE                *file;
    General_transform   gen_transform;

    create_linear_transform( &gen_transform, transform );

    status = open_file_with_default_suffix( filename, "xfm",
                                 WRITE_FILE, ASCII_FORMAT, &file );

    if( status == OK )
    {
        status = output_transform( file, filename, NULL, comments,
                                   &gen_transform );

        (void) close_file( file );
    }

    delete_general_transform( &gen_transform );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : read_transform_file
@INPUT      : filename
@OUTPUT     : transform
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads a transform from a file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  read_transform_file(
    STRING       filename,
    Transform    *transform )
{
    Status             status;
    FILE               *file;
    Transform          *lin_trans;
    General_transform  gen_transform;

    status = open_file_with_default_suffix( filename, "xfm",
                         READ_FILE, ASCII_FORMAT, &file );

    if( status == OK )
    {
        if( input_transform( file, filename, &gen_transform ) != OK )
            status = ERROR;

        (void) close_file( file );
    }

    if( status == OK && get_transform_type(&gen_transform) != LINEAR )
    {
        print( "File %s does not contain a linear transform.\n", filename );
        status = ERROR;
    }

    if( status == OK )
    {
        lin_trans = get_linear_transform_ptr( &gen_transform );
        *transform = *lin_trans;
    }

    delete_general_transform( &gen_transform );

    return( status );
}
