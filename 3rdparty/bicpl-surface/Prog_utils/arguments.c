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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Prog_utils/arguments.c,v 1.14 2005/08/17 22:27:56 bert Exp $";
#endif

typedef struct
{
    int   argc;
    char  **argv;
    int   current_arg;
    int   current_offset_within_arg;
} arguments_struct;

static  arguments_struct   arguments;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : arguments_remaining
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE if some arguments remain to be parsed
@DESCRIPTION: Checks if all arguments have been parsed.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  BOOLEAN  arguments_remaining( void )
{
    return( arguments.current_arg < arguments.argc );
}

BICAPI  int  get_n_arguments_remaining( void )
{
    return( arguments.argc - arguments.current_arg );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_current_argument_string
@INPUT      : 
@OUTPUT     : 
@RETURNS    : a pointer to the argument string.
@DESCRIPTION: Returns a point to an argument string.  This may be one
            : of the argv[i] parameters, or may point within one of these,
            : for instance, if the -n part of -nfilename has been parsed.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  get_current_argument_string( void )
{
    STRING   arg;

    arg = arguments.argv[arguments.current_arg];
    return( &arg[arguments.current_offset_within_arg] );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : advance_argument
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Advances the current arg pointer to the next argument.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  advance_argument( void )
{
    ++arguments.current_arg;
    arguments.current_offset_within_arg = 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_argument_processing
@INPUT      : argc
            : argv
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Records the arguments from the main program, and initializes
            : parsing to the first argument.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_argument_processing(
    int     argc,
    char    *argv[] )
{
    arguments.argc = argc;
    arguments.argv = argv;
    arguments.current_arg = 1;
    arguments.current_offset_within_arg = 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_int_argument
@INPUT      : default_value
@OUTPUT     : value
@RETURNS    : TRUE if next argument exists and is an integer
@DESCRIPTION: Checks if arguments remain, and if the next one is an integer,
            : stores it in "value", and advances the current pointer to the
            : next argument.  Otherwise, value is assigned the default_value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_int_argument(
    int   default_value,
    int   *value )
{
    BOOLEAN   found;

    if( arguments_remaining() &&
        sscanf( get_current_argument_string(), "%d", value ) == 1 )
    {
        advance_argument();
        found = TRUE;
    }
    else
    {
        *value = default_value;
        found = FALSE;
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_real_argument
@INPUT      : default_value
@OUTPUT     : value
@RETURNS    : TRUE if next argument exists and is a real
@DESCRIPTION: Checks if arguments remain, and if the next one is a real,
            : stores it in "value", and advances the current pointer to the
            : next argument.  Otherwise, value is assigned the default_value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_real_argument(
    Real   default_value,
    Real   *value )
{
    BOOLEAN   found;

    if( arguments_remaining() &&
        ((real_is_double() &&
	  sscanf( get_current_argument_string(), "%lf", (double *) value ) == 1) ||
         (!real_is_double() &&
	  sscanf( get_current_argument_string(), "%f", (float *) value ) == 1))
	)
    {
        advance_argument();
        found = TRUE;
    }
    else
    {
        *value = default_value;
        found = FALSE;
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_string_argument
@INPUT      : default_value
@OUTPUT     : value
@RETURNS    : TRUE if next argument exists
@DESCRIPTION: Checks if arguments remain, stores the pointer to the next one
            : in "value", and advances the current pointer to the
            : next argument.  Otherwise, value is assigned the default_value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_string_argument(
    STRING   default_value,
    STRING   *value )
{
    BOOLEAN   found;

    if( arguments_remaining() )
    {
        *value = get_current_argument_string();
        advance_argument();
        found = TRUE;
    }
    else
    {
        *value = default_value;
        found = FALSE;
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_prefix_argument
@INPUT      : prefix          - string to look for, e.g., "-n"
@OUTPUT     : 
@RETURNS    : TRUE if next argument begins with the prefix
@DESCRIPTION: Checks if arguments remain, and if the next one starts with
            : the prefix, advances the current pointer to just after the
            : prefix, which may be within the same argument.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_prefix_argument(
    STRING  prefix )
{
    char      *next_str;
    BOOLEAN   found;

    if( arguments_remaining() &&
        string_length(get_current_argument_string()) >= string_length(prefix) &&
        strncmp( get_current_argument_string(), prefix,
                 (size_t) string_length(prefix) ) == 0 )
    {
        arguments.current_offset_within_arg += string_length(prefix);

        next_str = get_current_argument_string();
        while( *next_str == ' ' || *next_str == '\t' )
            ++next_str;

        if( *next_str == END_OF_STRING )
            advance_argument();

        found = TRUE;
    }
    else
    {
        found = FALSE;
    }

    return( found );
}
