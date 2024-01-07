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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Prog_utils/globals.c,v 1.16 2005/08/17 22:27:56 bert Exp $";
#endif

static    Status  input_global_variable( int, global_struct [],
                                         FILE *, BOOLEAN * );
static    STRING    extract_string( STRING );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_globals_file
@INPUT      : n_globals_lookup
              globals_lookup
              filename
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads the global variable values from filename, using the
              globals lookup to find the relevant globals.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  input_globals_file(
    int             n_globals_lookup,
    global_struct   globals_lookup[],
    STRING          filename )
{
    Status  status;
    BOOLEAN eof;
    FILE    *file;

    status = open_file( filename, READ_FILE, ASCII_FORMAT, &file );

    if( status == OK )
    {
        eof = FALSE;

        while( !eof )
            status = input_global_variable( n_globals_lookup,
                                            globals_lookup, file, &eof );
    }

    if( status == OK )
        status = close_file( file );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_global_variable
@INPUT      : n_globals_lookup
              globals_lookup
              file
@OUTPUT     : eof  - TRUE if end of file
@RETURNS    : 
@DESCRIPTION: Inputs a global variable value from the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Status  input_global_variable(
    int             n_globals_lookup,
    global_struct   globals_lookup[],
    FILE            *file,
    BOOLEAN         *eof )
{
    Status   set_status, status;
    STRING   variable_name;
    STRING   value;

    *eof = FALSE;

    value = NULL;

    status = input_string( file, &variable_name, (char) '=' );

    if( status == OK )
    {
        status = input_string( file, &value, (char) ';' );

        if( status == OK )
            set_status = set_global_variable( n_globals_lookup, globals_lookup,
                                              variable_name, value );
        else
            *eof = TRUE;
    }
    else
    {
        status = OK;
        set_status = OK;
        *eof = TRUE;
    }

    if( set_status != OK || status != OK )
    {
        print_error( "Error inputting global.\n" );
        print_error( "Variable name is %s\n", variable_name );
    }

    delete_string( variable_name );
    delete_string( value );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : lookup_global
@INPUT      : n_globals
              global_lookup
              variable_name
@OUTPUT     : ptr    - pointer to global variable
              type
@RETURNS    : OK or ERROR
@DESCRIPTION: Given a name string, looks up the corresponding global variable.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Status  lookup_global(
    int              n_globals,
    global_struct    global_lookup[],
    STRING           variable_name,
    void             **ptr,
    Variable_types   *type,
    Smallest_int     **set_flag )
{
    Status  status;
    STRING  stripped;
    char    *global_name;
    int     i, s, len;

    status = ERROR;

    stripped = strip_outer_blanks( variable_name );

    for_less( i, 0, n_globals )
    {
        global_name = global_lookup[i].variable_name;
        len = string_length( global_name );
        s = 0;

        while( s < len && global_name[s] == ' ' )
            ++s;

        if( equal_strings( &global_name[s], stripped ) )
        {
            *ptr = global_lookup[i].ptr_to_global;
            *type = global_lookup[i].type;
            if( set_flag != NULL )
                *set_flag = &global_lookup[i].set_flag;
            status = OK;
            break;
        }
    }

    delete_string( stripped );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_global_variable
@INPUT      : n_globals_lookup
              globals_lookup
              variable_name
@OUTPUT     : value
@RETURNS    : 
@DESCRIPTION: Looks up the given variable name and places its value in string
              form in the 'value' parameter.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  get_global_variable(
    int              n_globals_lookup,
    global_struct    globals_lookup[],
    STRING           variable_name,
    STRING           *value )
{
    Status             status;
    void               *ptr;
    Surfprop           *surfprop;
    Variable_types     type;
    char               buffer[EXTREMELY_LARGE_STRING_SIZE];
    STRING             *string_ptr;
    STRING             tmp_str;

    status = lookup_global( n_globals_lookup, globals_lookup, variable_name,
                            &ptr, &type, NULL );

    if( status == OK )
    {
        switch( type )
        {
        case BOOLEAN_type:
            if( * ((BOOLEAN *) ptr) )
                (void) strcpy( buffer, "True" );
            else
                (void) strcpy( buffer, "False" );
            break;

        case int_type:
            (void) sprintf( buffer, "%d", * (int *) ptr );
            break;

        case Real_type:
            (void) sprintf( buffer, "%g", * (Real *) ptr );
            break;

        case STRING_type:
            string_ptr = (STRING *) ptr;
            (void) strcpy( buffer, *string_ptr );
            break;

        case Point_type:
            (void) sprintf( buffer, "%g %g %g",
                            Point_x(* (Point *) ptr),
                            Point_y(* (Point *) ptr),
                            Point_z(* (Point *) ptr) );
            break;

        case Vector_type:
            (void) sprintf( buffer, "%g %g %g",
                            Vector_x(* (Vector *) ptr),
                            Vector_y(* (Vector *) ptr),
                            Vector_z(* (Vector *) ptr) );
            break;

        case Colour_type:
            tmp_str = convert_colour_to_string( * (Colour *) ptr );
            (void) strcpy( buffer, tmp_str );
            delete_string( tmp_str );
            break;

        case Surfprop_type:
            surfprop = (Surfprop *) ptr;
            (void) sprintf( buffer, "%g %g %g %g %g",
                            Surfprop_a(*surfprop),
                            Surfprop_d(*surfprop),
                            Surfprop_s(*surfprop),
                            Surfprop_se(*surfprop),
                            Surfprop_t(*surfprop) );
            break;

        default:
            handle_internal_error( "get_global_variable\n" );
            buffer[0] = END_OF_STRING;
            break;
        }
    }

    *value = create_string( buffer );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_global_variable
@INPUT      : n_globals_lookup
              globals_lookup
              variable_name
              value_to_set
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Looks up the given variable name, converts the value_to_set
              string to the correct type, and sets the global variable.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  set_global_variable(
    int              n_globals_lookup,
    global_struct    globals_lookup[],
    STRING           variable_name,
    STRING           value_to_set )
{
    Status             status;
    STRING             value, *string_ptr;
    void               *ptr;
    Variable_types     type;
    int                tmp_int;
    Real               tmp_real;
    Point              tmp_point;
    Vector             tmp_vector;
    Colour             tmp_colour;
    Surfprop           tmp_surfprop;
    Smallest_int       *set_flag;

    value = strip_outer_blanks( value_to_set );

    status = lookup_global( n_globals_lookup, globals_lookup, variable_name,
                            &ptr, &type, &set_flag );

    if( status == OK )
    {
        switch( type )
        {
        case BOOLEAN_type:
            if( (value[0] == 't' || value[0] == 'T') )
            {
                * (BOOLEAN *) ptr = TRUE;
            }
            else if( (value[0] == 'f' || value[0] == 'F') )
            {
                * (BOOLEAN *) ptr = FALSE;
            }
            else
            {
                status = ERROR;
            }
            break;

        case int_type:
            if( sscanf( value, "%d", &tmp_int ) == 1 )
                * (int *) ptr = tmp_int;
            else
                status = ERROR;
            break;

        case Real_type:
            if( real_is_double() )
            {
                tmp_real = 0.0;
                if( sscanf( value, "%lf", (double *) &tmp_real ) == 1 )
                    * (Real *) ptr = tmp_real;
                else
                    status = ERROR;
            }
            else
            {
                if( sscanf( value, "%f", (float *) &tmp_real ) == 1 )
                    * (Real *) ptr = tmp_real;
                else
                    status = ERROR;
            }
            break;

        case STRING_type:
            string_ptr = (STRING *) ptr;
            if( *set_flag )
                delete_string( *string_ptr );
            *string_ptr = extract_string( value );
            *set_flag = TRUE;
            break;

        case Point_type:
            if( sscanf( value, "%f %f %f", &Point_x(tmp_point),
                        &Point_y(tmp_point), &Point_z(tmp_point) ) == 3 )
                * (Point *) ptr = tmp_point;
            else
                status = ERROR;
            break;

        case Vector_type:
            if( sscanf( value, "%f %f %f", &Vector_x(tmp_vector),
                        &Vector_y(tmp_vector), &Vector_z(tmp_vector) ) == 3 )
                * (Vector *) ptr = tmp_vector;
            else
                status = ERROR;
            break;

        case Colour_type:
            tmp_colour = convert_string_to_colour( value );

            if( status == OK )
            {
                * (Colour *) ptr = tmp_colour;
            }
            break;

        case Surfprop_type:
            if( sscanf( value, "%f %f %f %f %f",
                            &Surfprop_a(tmp_surfprop),
                            &Surfprop_d(tmp_surfprop),
                            &Surfprop_s(tmp_surfprop),
                            &Surfprop_se(tmp_surfprop),
                            &Surfprop_t(tmp_surfprop) ) == 5 )
            {
                * (Surfprop *) ptr = tmp_surfprop;
            }
            else
                status = ERROR;
            break;

        default:
            handle_internal_error( "set_global_variable\n" );
            break;
        }
    }

    delete_string( value );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_or_get_global_variable
@INPUT      : n_globals_lookup
              globals_lookup
              input_str
@OUTPUT     : variable_name
              value_string
@RETURNS    : 
@DESCRIPTION: input_str is of the form "variable_name" or 
              "variable_name = value".  In both cases the variable name and
              value of the variable are passed back.  In the second case,
              the variable is first assigned the value.  This routine is
              usually used by programs to allow the user to type in input
              strings to query or modify the values of global variables.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  set_or_get_global_variable(
    int              n_globals_lookup,
    global_struct    globals_lookup[],
    STRING           input_str,
    STRING           *variable_name,
    STRING           *value_string )
{
    Status  status;
    STRING  tmp_var_name, value_to_set;
    int     equal_index;

    status = OK;

    tmp_var_name = create_string( input_str );

    equal_index = find_character( tmp_var_name, (char) '=' );

    if( equal_index >= 0 )
    {
        value_to_set = create_string( &tmp_var_name[equal_index+1] );
        tmp_var_name[equal_index] = END_OF_STRING;
        status = set_global_variable( n_globals_lookup, globals_lookup,
                                      tmp_var_name, value_to_set );
        delete_string( value_to_set );
    }

    *variable_name = strip_outer_blanks( tmp_var_name );

    if( status == OK )
    {
        status = get_global_variable( n_globals_lookup, globals_lookup,
                                      *variable_name, value_string );
    }

    delete_string( tmp_var_name );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : extract_string
@INPUT      : str
@OUTPUT     : extracted
@RETURNS    : 
@DESCRIPTION: Extracts a string from str.  Skips leading and trailing
              white space.  String can be a quoted string or not.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  STRING  extract_string(
    STRING   str )
{
    int     len_str, start, end, i;
    STRING  extracted;

    start = 0;

    len_str = string_length( str );

    while( start < len_str && (str[start] == ' ' || str[start] == '\t') )
    {
        ++start;
    }

    end = len_str - 1;

    while( end >= 0 && (str[end] == ' ' || str[end] == '\t') )
    {
        --end;
    }

    if( start <= end && str[start] == '"' && str[end] == '"' )
    {
        ++start;
        --end;
    }

    extracted = create_string( NULL );

    if( start <= end )
    {
        for_inclusive( i, start, end )
            concat_char_to_string( &extracted, str[i] );
    }

    return( extracted );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_global_variables
@INPUT      : n_globals_lookup
              globals_lookup
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_global_variables(
    int             n_globals_lookup,
    global_struct   globals_lookup[] )
{
    int     i;
    STRING  *ptr;

    for_less( i, 0, n_globals_lookup )
    {
        if( globals_lookup[i].type == STRING_type &&
            globals_lookup[i].set_flag )
        {
            ptr = (STRING *) globals_lookup[i].ptr_to_global;
            delete_string( *ptr );
        }
    }
}
