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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/object_io.c,v 1.30 2005/08/17 22:28:27 bert Exp $";
#endif


/*! \defgroup grp_file File Format

A BIC OBJ file stores one or more geometrical objects with
associated data like colour, and patient identifiers.  The
file name typically ends with ".obj".

The list of objects is stored in the file one after
the other.  There are no separator bytes between objects.

The first charcter signals the type of the next object.  If the object
is stored in BINARY_FORMAT, the type byte is a lower-case character.
If the object is in ASCII_FORMAT, the type is the corresponding
upper-case character.  The following table maps the type byte to its
object type.

- L = LINES
- M = MARKER   
- F = MODEL    
- X = PIXELS   
- P = POLYGONS 
- Q = QUADMESH 
- T = TEXT     

  @{
*/

/*! \defgroup grp_file_lines Lines format

A LINES object file follows the following sequence of data items.

  - 'L'
  - (float) thickness
  - (int) n_points
  - (newline)
  - points 
    - three floats per line
  - (newline)
  - (int) n_items (#lines)
  - (newline)
  - colours
    - (int) colour_flag 0=ONE_COLOUR, 1=PER_ITEM, 2=PER_VERTEX
    - the specified number of colours, one per line
  - (newline)
  - end_indices
  - (newline)
  - indices

*/

/*! \defgroup grp_file_marker Marker format

A MARKER object consists of the following data items.

  - 'M'
  - (int) type    0=box, 1=sphere
  - (real) size
  - (colour) colour
  - (point) position
  - (int) structure_id
  - (int) patient_id
  - (quoted_string) label
  - (newline)

*/

/*! \defgroup grp_file_polygons Polygons format

A POLYGONS object consists of the following data items.

  - 'P'
  - five floats representing 'surfprop'
  - n_points if normal format,
    or -n_items (#polygons) if compressed format
  - newline
  - table of point coordinate, one point per line
  - newline (blank line)
  - if not compressed
    - table of point normals, one per line 
    - newline (blank line)
    - n_items (#polygons)
    - newline
  - colour flag
  - for each allocated colour (1, n_points, or n_items)
    - colour and newline
  - newline (blank line)
  - if not compressed
    - end_indices (two formats, see Objects/o.c)
    - newline (blank line)
    - indices
    - newline (blank line)


In compressed format, the polygons form either a tetrahedron (4
polygons), two tetrahedrons glued together (6 polygons), two pyramids
glued together (8 polygons) or a regular subdivision of one of these.
Surface normals for each point are computed, after reading file.

*/


/*! \defgroup grp_file_quadmesh Quadmesh format

A QUADMESH object consists of the following data items.

  - 'Q'
  - surface properties
    - five floats: a, d, s, se, t
  - (integer) m
  - (integer) n
  - (T or F) m_closed
  - (T or F) n_closed
  - newline
  - colours
    - (integer) colour_flag 0=ONE_COLOUR, 1=PER_ITEM, 2=PER_VERTEX
    - the specified number of colours, one colour per line
    - each colour is four floats: r g b a
  - newline
  - points
    - one point (3 floats) per line
  - newline
  - normals
    - one vector (3 floats) per line
  - newline

*/

/*! \defgroup grp_file_text Text format

A TEXT object consists of the following data items.

  - 'T'
  - font type (0 = fixed, 1 = sized)
  - text size (real)
  - colour (4 floats: rgba)
  - point
  - quoted string
  - newline

*/

/*! 
  @} 
*/


static  Status  io_vectors(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    Vector          *vectors[] );
static  Status  io_line_thickness(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    float           *line_thickness );
static  Status  io_end_indices(
    FILE           *file,
    IO_types       io_flag,
    File_formats   format,
    int            n_items,
    int            *end_indices[],
    int            min_size );
static  Status  io_points(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    Point           *points[] );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_lines
@INPUT      : file
              io_flag
              format
              lines
@OUTPUT     : (lines)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a lines structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_lines(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    lines_struct        *lines )
{
    Status   status;

    status = OK;

    if( io_flag == READ_FILE )
    {
        initialize_lines( lines, WHITE );
        FREE( lines->colours );
    }

    if( io_flag == READ_FILE ||
        (lines->n_points > 0 && lines->n_items > 0) )
    {
        status = io_object_type( file, io_flag, format, LINES );

        if( status == OK )
            status = io_line_thickness( file, io_flag, format,
                                        &lines->line_thickness );

        if( status == OK )
            status = io_int( file, io_flag, format, &lines->n_points );

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
        {
            status = io_points( file, io_flag, format,
                                lines->n_points, &lines->points );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
            status = io_int( file, io_flag, format, &lines->n_items );

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
        {
            status = io_colours( file, io_flag, format, &lines->colour_flag,
                                 lines->n_items, lines->n_points,
                                 &lines->colours );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
        {
            status = io_end_indices( file, io_flag, format,
                                     lines->n_items, &lines->end_indices, 1 );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
        {
            status = io_ints( file, io_flag, format,
                              NUMBER_INDICES(*lines),
                              &lines->indices );
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_marker
@INPUT      : file
              io_flag
              format
              marker
@OUTPUT     : (marker)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a marker structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_marker(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    marker_struct       *marker )
{
    Status   status;

    status = io_object_type( file, io_flag, format, MARKER );

    if( status == OK )
        status = io_int( file, io_flag, format, (int *) &marker->type );

    if( status == OK )
        status = io_real( file, io_flag, format, &marker->size );

    if( status == OK )
        status = io_colour( file, io_flag, format, &marker->colour );

    if( status == OK )
        status = io_point( file, io_flag, format, &marker->position );

    if( status == OK )
        status = io_int( file, io_flag, format, &marker->structure_id );

    if( status == OK )
        status = io_int( file, io_flag, format, &marker->patient_id );

    if( status == OK )
        status = io_quoted_string( file, io_flag, format, &marker->label );

    if( status == OK )
        status = io_newline( file, io_flag, format );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_model
@INPUT      : file
              io_flag
              format
              model
@OUTPUT     : (model)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a model structure in ascii or binary format.
              Rather than writing all the objects in the model, it simply
              writes the name of the file containing the objects.  A
              separate call (e.g. output_graphics_file) must be made to
              input or output the objects in the model.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_model(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    model_struct    *model )
{
    Status   status;

    status = io_object_type( file, io_flag, format, MODEL );

    if( status == OK )
        status = io_quoted_string( file, io_flag, format, &model->filename );

    if( status == OK )
        status = io_newline( file, io_flag, format );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_pixels
@INPUT      : file
              io_flag
              format
              pixels
@OUTPUT     : (pixels)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a pixels structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_pixels(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    pixels_struct       *pixels )
{
    Status       status;
    int          n_pixels;

    status = OK;

    if( io_flag == READ_FILE || (pixels->x_size > 0 && pixels->y_size > 0) )
    {
        status = io_object_type( file, io_flag, format, PIXELS );

        if( status == OK )
            status = io_int( file, io_flag, format, (int *)&pixels->pixel_type);

        if( status == OK )
            status = io_int( file, io_flag, format, &pixels->x_size );

        if( status == OK )
            status = io_int( file, io_flag, format, &pixels->y_size );

        if( status == OK && io_flag == READ_FILE )
        {
            pixels->x_position = 0;
            pixels->y_position = 0;
            pixels->x_zoom = 1.0;
            pixels->y_zoom = 1.0;
        }

        if( status == OK )
        {
            n_pixels = pixels->x_size * pixels->y_size;

            switch( pixels->pixel_type )
            {
            case COLOUR_INDEX_8BIT_PIXEL:
                status = io_unsigned_chars( file, io_flag, format, n_pixels,
                                  &pixels->data.pixels_8bit_colour_index );
                break;

            case COLOUR_INDEX_16BIT_PIXEL:
                if( io_flag == READ_FILE )
                    ALLOC( pixels->data.pixels_16bit_colour_index, n_pixels );

                status = io_binary_data( file, io_flag,
                          (void *) (&pixels->data.pixels_16bit_colour_index),
                          sizeof(pixels->data.pixels_16bit_colour_index[0]),
                          n_pixels );
                break;

            case RGB_PIXEL:
                status = io_pixel_colours( file, io_flag, format, n_pixels,
                                           &pixels->data.pixels_rgb );
                break;

            default:
                print_error( "Error, unrecognized pixel type %d.\n",
                             pixels->pixel_type );
                status = ERROR;
            }
        }
    }

    return( status );
}

static  BOOLEAN  use_compressed = FALSE;
static  BOOLEAN  compressed_initialized = FALSE;

static  BOOLEAN   use_compressed_polygons( void )
{
    if( !compressed_initialized )
    {
        compressed_initialized = TRUE;
        use_compressed = getenv( "USE_COMPRESSED_POLYGONS" ) != NULL;
    }

    return( use_compressed );
}

BICAPI  void  set_use_compressed_polygons_flag(
    BOOLEAN  value )
{
    use_compressed = value;
    compressed_initialized = TRUE;
}

BICAPI  BOOLEAN  get_use_compressed_polygons_flag( void )
{
    return( use_compressed_polygons() );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_polygons
@INPUT      : file
              io_flag
              format
              polygons
@OUTPUT     : (polygons)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a polygons structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Mar. 23, 1997   D. MacDonald - added compressed polygons.
---------------------------------------------------------------------------- */

BICAPI  Status  io_polygons(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    polygons_struct     *polygons )
{
    int      n_items;
    Status   status;
    Surfprop save_surfprop;
    Point    centre;
    BOOLEAN  compressed_format;

    status = OK;

    if( io_flag == READ_FILE )
    {
        initialize_polygons( polygons, WHITE, NULL );
        FREE( polygons->colours );
    }

    if( io_flag == READ_FILE ||
        (polygons->n_points > 0 && polygons->n_items > 0) )
    {
        status = io_object_type( file, io_flag, format, POLYGONS );

        if( status == OK )
            status = io_surfprop( file, io_flag, format, &polygons->surfprop );

        compressed_format = FALSE;

        if( status == OK )
        {
            if( io_flag == WRITE_FILE && use_compressed_polygons() &&
                is_this_tetrahedral_topology(polygons) )
            {
                n_items = -polygons->n_items;
                status = io_int( file, io_flag, format, &n_items );
                compressed_format = TRUE;
            }
            else
                status = io_int( file, io_flag, format, &polygons->n_points );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( io_flag == READ_FILE && polygons->n_points < 0 )
        {
            n_items = -polygons->n_points;
            compressed_format = TRUE;
            fill_Point( centre, 0.0, 0.0, 0.0 );
            save_surfprop = polygons->surfprop;
            create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0, n_items,
                                       polygons );
            polygons->surfprop = save_surfprop;
            FREE( polygons->points );
        }

        if( status == OK )
        {
            status = io_points( file, io_flag, format,
                                polygons->n_points, &polygons->points );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( !compressed_format )
        {
            if( status == OK )
            {
                status = io_vectors( file, io_flag, format,
                                     polygons->n_points, &polygons->normals );
            }

            if( status == OK )
                status = io_newline( file, io_flag, format );

            if( status == OK )
                status = io_int( file, io_flag, format, &polygons->n_items );

            if( status == OK )
                status = io_newline( file, io_flag, format );
        }

        if( status == OK )
        {
            status = io_colours( file, io_flag, format, &polygons->colour_flag,
                                 polygons->n_items, polygons->n_points,
                                 &polygons->colours );

            if( status == OK )
                status = io_newline( file, io_flag, format );
        }

        if( !compressed_format )
        {
            if( status == OK )
            {
                status = io_end_indices( file, io_flag, format,
                                     polygons->n_items, &polygons->end_indices,
                                     3 );
            }

            if( status == OK )
                status = io_newline( file, io_flag, format );

            if( status == OK )
            {
                status = io_ints( file, io_flag, format,
                                  NUMBER_INDICES(*polygons),
                                  &polygons->indices );
            }

            if( status == OK )
                status = io_newline( file, io_flag, format );
        }

        if( io_flag == READ_FILE && compressed_format )
            compute_polygon_normals( polygons );

        if( io_flag == READ_FILE )
            polygons->line_thickness = 1.0f;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_quadmesh
@INPUT      : file
              io_flag
              format
              quadmesh
@OUTPUT     : (quadmesh)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a quadmesh structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Read or write quadmesh structure.
 *
 * \param io_flag one of READ_FILE or WRITE_FILE
 * \param format one of ASCII_FORMAT or BINARY_FORMAT
 */
BICAPI  Status  io_quadmesh(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    quadmesh_struct     *quadmesh )
{
    Status   status;

    status = OK;

    if( io_flag == READ_FILE )
    {
        initialize_quadmesh( quadmesh, WHITE, NULL, 0, 0 );
        FREE( quadmesh->colours );
    }

    if( io_flag == READ_FILE ||
        (quadmesh->m > 1 && quadmesh->n > 1) )
    {
        status = io_object_type( file, io_flag, format, QUADMESH );

        if( status == OK )
            status = io_surfprop( file, io_flag, format, &quadmesh->surfprop );

        if( status == OK )
            status = io_int( file, io_flag, format, &quadmesh->m );

        if( status == OK )
            status = io_int( file, io_flag, format, &quadmesh->n );

        if( status == OK )
            status = io_boolean( file, io_flag, format, &quadmesh->m_closed );

        if( status == OK )
            status = io_boolean( file, io_flag, format, &quadmesh->n_closed );

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
            status = io_colours( file, io_flag, format, &quadmesh->colour_flag,
                                 (quadmesh->m-1) * (quadmesh->n-1),
                                 quadmesh->m * quadmesh->n,
                                 &quadmesh->colours );

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
        {
            status = io_points( file, io_flag, format,
                                quadmesh->m * quadmesh->n, &quadmesh->points );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );

        if( status == OK )
        {
            status = io_vectors( file, io_flag, format,
                                quadmesh->m * quadmesh->n, &quadmesh->normals );
        }

        if( status == OK )
            status = io_newline( file, io_flag, format );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_text
@INPUT      : file
              io_flag
              format
              text
@OUTPUT     : (text)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a text structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_text(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    text_struct     *text )
{
    Status   status;

    status = io_object_type( file, io_flag, format, TEXT );

    if( status == OK )
        status = io_int( file, io_flag, format, (int *) &text->font );

    if( status == OK )
        status = io_real( file, io_flag, format, &text->size );

    if( status == OK )
        status = io_colour( file, io_flag, format, &text->colour );

    if( status == OK )
        status = io_point( file, io_flag, format, &text->origin );

    if( status == OK )
        status = io_quoted_string( file, io_flag, format, &text->string );

    if( status == OK )
        status = io_newline( file, io_flag, format );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_point
@INPUT      : file
              io_flag
              format
              point
@OUTPUT     : (point)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a point structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_point(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Point           *point )
{
    Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        status = io_float( file, io_flag, format, &Point_x(*point) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Point_y(*point) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Point_z(*point) );
    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) point,
                                 sizeof(*point), 1 );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_vector
@INPUT      : file
              io_flag
              format
              vector
@OUTPUT     : (vector)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a vector structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_vector(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Vector          *v )
{
    Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        status = io_float( file, io_flag, format, &Vector_x(*v) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Vector_y(*v));

        if( status == OK )
            status = io_float( file, io_flag, format, &Vector_z(*v));
    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) v, sizeof(*v), 1 );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_colour
@INPUT      : file
              io_flag
              format
              colour
@OUTPUT     : (colour)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a colour in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_colour(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Colour          *colour )
{
    float          r, g, b, a;
    unsigned char  comps[4];
    Status         status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == WRITE_FILE )
        {
            r = (float) get_Colour_r_0_1( *colour );
            g = (float) get_Colour_g_0_1( *colour );
            b = (float) get_Colour_b_0_1( *colour );
            a = (float) get_Colour_a_0_1( *colour );
        }

        status = io_float( file, io_flag, ASCII_FORMAT, &r );
        if( status == OK )
            status = io_float( file, io_flag, ASCII_FORMAT, &g );
        if( status == OK )
            status = io_float( file, io_flag, ASCII_FORMAT, &b );
        if( status == OK )
            status = io_float( file, io_flag, ASCII_FORMAT, &a );

        if( io_flag == READ_FILE )
            *colour = make_rgba_Colour_0_1( (Real) r, (Real) g,
                                            (Real) b, (Real) a );
    }
    else
    {
        if( io_flag == WRITE_FILE )
        {
            comps[3] = (unsigned char) get_Colour_r( *colour );
            comps[2] = (unsigned char) get_Colour_g( *colour );
            comps[1] = (unsigned char) get_Colour_b( *colour );
            comps[0] = (unsigned char) get_Colour_a( *colour );
        }

        status = io_binary_data( file, io_flag, (void *) comps,
                                 sizeof(comps[0]), 4 );

        if( io_flag == READ_FILE )
        {
            *colour = make_rgba_Colour( (int) comps[3], (int) comps[2],
                                        (int) comps[1], (int) comps[0]);
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_colours
@INPUT      : file
              io_flag
              format
              colour_flag
              n_items
              n_points
              colours
@OUTPUT     : (colours)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a colour_flag and list of colours, where the
              number of colours depends on the colour flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_colours(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Colour_flags    *colour_flag,
    int             n_items,
    int             n_points,
    Colour          **colours )
{
    int      i, n_colours;
    Status   status;

    status = io_int( file, io_flag, format, (int *) colour_flag );

    if( status == OK )
    {
        switch( *colour_flag )
        {
        case ONE_COLOUR:          n_colours = 1;   break;
        case PER_ITEM_COLOURS:    n_colours = n_items;   break;
        case PER_VERTEX_COLOURS:  n_colours = n_points;   break;
        default:
            print_error( "Error inputting colour flag.\n" );
            status = ERROR;
            break;
        }
    }

    if( status == OK && io_flag == READ_FILE && n_colours > 0 )
        ALLOC( *colours, n_colours );

    if( status == OK )
    {
        for_less( i, 0, n_colours )
        {
            status = io_colour( file, io_flag, format, &(*colours)[i] );

            if( status == OK )
                status = io_newline( file, io_flag, format );
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_surfprop
@INPUT      : file
              io_flag
              format
              surfprop
@OUTPUT     : (surfprop)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a surfprop structure in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_surfprop(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Surfprop        *surfprop )
{
    Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        status = io_float( file, io_flag, format,
                          &Surfprop_a(*surfprop) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Surfprop_d(*surfprop) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Surfprop_s(*surfprop) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Surfprop_se(*surfprop) );

        if( status == OK )
            status = io_float( file, io_flag, format, &Surfprop_t(*surfprop) );

    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) surfprop,
                                 sizeof(*surfprop), 1 );
    }

    if( io_flag == READ_FILE && Surfprop_t(*surfprop) == 0.0f )
        Surfprop_t(*surfprop) = 1.0f;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_points
@INPUT      : file
              io_flag
              format
              n
              points
@OUTPUT     : (points)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a list of points in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Status  io_points(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    Point           *points[] )
{
    Status   status;
    int      i;

    status = OK;

    if( io_flag == READ_FILE )
    {
        ALLOC( *points, n );
    }

    if( status == OK )
    {
        if( format == ASCII_FORMAT )
        {
            for_less( i, 0, n )
            {
                status = io_point( file, io_flag, format, &(*points)[i]);

                if( status == OK )
                    status = io_newline( file, io_flag, format );

                if( status == ERROR )
                    break;
            }
        }
        else
        {
            status = io_binary_data( file, io_flag, (void *) (*points),
                                     sizeof((*points)[0]), n );
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_vectors
@INPUT      : file
              io_flag
              format
              n
              vectors
@OUTPUT     : (vectors)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a list of vectors in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Status  io_vectors(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    Vector          *vectors[] )
{
    Status   status;
    int      i;

    status = OK;

    if( io_flag == READ_FILE )
    {
        ALLOC( *vectors, n );
    }

    if( format == ASCII_FORMAT )
    {
        for_less( i, 0, n )
        {
            status = io_vector( file, io_flag, format, &(*vectors)[i] );

            if( status == OK )
                status = io_newline( file, io_flag, format );

            if( status == ERROR )
                break;
        }
    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) (*vectors),
                                 sizeof((*vectors)[0]), n );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_object_type
@INPUT      : file
              io_flag
              format
              type
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Either writes an object type or if reading, does not read
              anything.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_object_type(
    FILE           *file,
    IO_types       io_flag,
    File_formats   format,
    Object_types   type )
{
    int      ch;
    Status   status;

    status = OK;

    if( io_flag == WRITE_FILE )
    {
        switch( type )
        {
        case  LINES:       ch = 'l';    break;
        case  MARKER:      ch = 'm';    break;
        case  MODEL:       ch = 'f';    break;
        case  PIXELS:      ch = 'x';    break;
        case  POLYGONS:    ch = 'p';    break;
        case  QUADMESH:    ch = 'q';    break;
        case  TEXT:        ch = 't';    break;
	default: 
            print_error( "io_object_type: object type %d not handled.\n", type );
            return ERROR;
        }

        if( format == ASCII_FORMAT )
        {
            ch += 'A' - 'a';
        }

        if( fputc( ch, file ) == EOF )
        {
            print_error( "Error outputting char.\n" );
            status = ERROR;
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_object_type
@INPUT      : file
@OUTPUT     : type
              format
              eof
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads the object type, passing it back, as well as whether
              it is in ascii or binary format, and an eof flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  input_object_type(
    FILE           *file,
    Object_types   *type,
    File_formats   *format,
    BOOLEAN        *eof )
{
    char     ch;
    Status   status;

    status = OK;
    *eof = FALSE;

    status = input_nonwhite_character( file, &ch );

    if( status == OK )
    {
        if( ch >= 'A' && ch <= 'Z' )
        {
            *format = ASCII_FORMAT;
            ch -= 'A' - 'a';
        }
        else
        {
            *format = BINARY_FORMAT;
        }

        switch( ch )
        {
        case  'l':   *type = LINES;        break;
        case  'm':   *type = MARKER;       break;
        case  'f':   *type = MODEL;        break;
        case  'x':   *type = PIXELS;       break;
        case  'p':   *type = POLYGONS;     break;
        case  'q':   *type = QUADMESH;     break;
        case  't':   *type = TEXT;         break;
/*
        case  'v':   *type = VOLUME;       break;
*/

        default:
            print_error( "Unrecognized object type in file.\n" );
            status = ERROR;
        }
    }
    else
    {
        *eof = TRUE;
        status = OK;
    }

    return( status );
}

#define  PIXELS_PER_LINE    1
#define  BUFFER_SIZE        256

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_pixel_colours
@INPUT      : file
              io_flag
              format
              n
              pixel_colours
@OUTPUT     : (pixel_colours)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a list of pixel colours in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_pixel_colours(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    Colour          *pixel_colours[] )
{
    Status         status;
    unsigned char  buffer[3*BUFFER_SIZE];
    Colour         pixel_colour;
    int            i, block, n_blocks, n_to_do, start_pixel_index;

    status = OK;

    if( io_flag == READ_FILE )
        ALLOC( *pixel_colours, n );

    if( status == OK )
    {
        if( format == ASCII_FORMAT )
        {
            for_less( i, 0, n )
            {
                status = io_pixel_colour( file, io_flag, format,
                                          &(*pixel_colours)[i] );

                if( status == OK )
                {
                    if( i == n - 1 || (i+1) % PIXELS_PER_LINE == 0 )
                    {
                        status = io_newline( file, io_flag, format );
                    }
                }

                if( status != OK )
                    break;
            }
        }
        else
        {
            n_blocks = n / BUFFER_SIZE;
            if( n % BUFFER_SIZE != 0 )
                ++n_blocks;

            for_less( block, 0, n_blocks )
            {
                start_pixel_index = block * BUFFER_SIZE;
                n_to_do = BUFFER_SIZE;
                if( start_pixel_index + n_to_do > n )
                    n_to_do = n - start_pixel_index;

                if( io_flag == WRITE_FILE )
                {
                    for_less( i, 0, n_to_do )
                    {
                        pixel_colour = (*pixel_colours)[start_pixel_index+i];
                        buffer[3*i+0] = (unsigned char)
                                        get_Colour_r(pixel_colour);
                        buffer[3*i+1] = (unsigned char)
                                        get_Colour_g(pixel_colour);
                        buffer[3*i+2] = (unsigned char)
                                        get_Colour_b(pixel_colour);
                    }
                }

                status = io_binary_data( file, io_flag, (void *) buffer,
                                     sizeof(buffer[0]), 3 * n_to_do );

                if( io_flag == READ_FILE )
                {
                    for_less( i, 0, n_to_do )
                    {
                        pixel_colour = make_Colour( (int) buffer[3*i+0],
                                                    (int) buffer[3*i+1],
                                                    (int) buffer[3*i+2] );
                        (*pixel_colours)[start_pixel_index+i] = pixel_colour;
                    }
                }
            }
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_pixel_colour
@INPUT      : file
              io_flag
              format
              pixel_colour
@OUTPUT     : (pixel_colour)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a pixel colour in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_pixel_colour(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Colour          *pixel_colour )
{
    int      r, g, b;
    Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == WRITE_FILE )
        {
            r = get_Colour_r( *pixel_colour );
            g = get_Colour_g( *pixel_colour );
            b = get_Colour_b( *pixel_colour );
        }

        status = io_int( file, io_flag, format, &r );

        if( status == OK )
            status = io_int( file, io_flag, format, &g );

        if( status == OK )
            status = io_int( file, io_flag, format, &b );

        if( io_flag == READ_FILE && status == OK )
            *pixel_colour = make_Colour( r, g, b );
    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) pixel_colour,
                                 sizeof(*pixel_colour), 1 );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_line_thickness
@INPUT      : file
              io_flag
              format
              line_thickness
@OUTPUT     : (line_thickness)
@RETURNS    : OK or ERROR
@DESCRIPTION: Reads or writes a line thickness in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Status  io_line_thickness(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    float           *line_thickness )
{
    Status   status;

    status = io_float( file, io_flag, format, line_thickness );

    if( status == OK && io_flag == READ_FILE && format == BINARY_FORMAT &&
        (*line_thickness <= 0.001f || *line_thickness > 20.0f) )
        *line_thickness = 1.0f;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_object
@INPUT      : directory
              file
@OUTPUT     : format
              object
              eof
@RETURNS    : OK or ERROR
@DESCRIPTION: Inputs an object from a file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  input_object(
    STRING         directory,
    FILE           *file,
    File_formats   *format,
    object_struct  **object,
    BOOLEAN        *eof )
{
    Status         status;
    File_formats   sub_format;
    Object_types   type;
    STRING         abs_filename;

    status = input_object_type( file, &type, format, eof );

    if( status == OK && !(*eof) )
    {
        *object = create_object( type );

        switch( type )
        {
        case LINES:
            status = io_lines( file, READ_FILE, *format,
                               get_lines_ptr(*object) );
            break;

        case MARKER:
            status = io_marker( file, READ_FILE, *format,
                                get_marker_ptr(*object) );
            break;

        case MODEL:
            status = io_model( file, READ_FILE, *format,
                               get_model_ptr(*object) );

            abs_filename = get_absolute_filename(
                              get_model_ptr(*object)->filename,
                              directory );

            if( status == OK )
            {
                status = input_graphics_file( abs_filename, &sub_format,
                                          &get_model_ptr(*object)->n_objects,
                                          &get_model_ptr(*object)->objects );
            }

            delete_string( abs_filename );

            break;

        case PIXELS:
            status = io_pixels( file, READ_FILE, *format,
                                get_pixels_ptr(*object) );
            break;

        case POLYGONS:
            status = io_polygons( file, READ_FILE, *format,
                                   get_polygons_ptr( *object ) );

            break;

        case QUADMESH:
            status = io_quadmesh( file, READ_FILE, *format,
                                  get_quadmesh_ptr( *object ) );
            break;

        case TEXT:
            status = io_text( file, READ_FILE, *format,
                              get_text_ptr( *object ) );
            break;

        default:
            print_error( "Unrecognized object type %d\n", type );
            status = ERROR;
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_object
@INPUT      : file
              format
              object
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Writes an object to a file in ascii or binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  output_object(
    FILE           *file,
    File_formats   format,
    object_struct  *object )
{
    Status         status;

    switch( object->object_type )
    {
    case LINES:
        status = io_lines( file, WRITE_FILE, format, get_lines_ptr(object) );
        break;

    case MARKER:
        status = io_marker( file, WRITE_FILE, format, get_marker_ptr(object) );
        break;

    case MODEL:
        status = io_model( file, WRITE_FILE, format, get_model_ptr(object) );
        break;

    case PIXELS:
        status = io_pixels( file, WRITE_FILE, format, get_pixels_ptr(object) );
        break;

    case POLYGONS:
        status = io_polygons( file, WRITE_FILE, format,
                              get_polygons_ptr(object) );
        break;

    case QUADMESH:
        status = io_quadmesh( file, WRITE_FILE, format,
                              get_quadmesh_ptr(object) );
        break;

    case TEXT:
        status = io_text( file, WRITE_FILE, format, get_text_ptr(object) );
        break;

    default:
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_end_indices
@INPUT      : file
              io_flag
              format
              n_items
              end_indices
              min_size
@OUTPUT     : 
@RETURNS    : ERROR or OK
@DESCRIPTION: Converts between old file format of cumulative sizes of 
              individual polygons to new format of outputting each size.
              It tries to detect on input whether it is new format or old
              format by assuming it is cumulative sizes and checking for
              any individual size less than min_size.  This is not foolproof
              but should work for most files. 
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Jul. 31, 1996    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Status  io_end_indices(
    FILE           *file,
    IO_types       io_flag,
    File_formats   format,
    int            n_items,
    int            *end_indices[],
    int            min_size )
{
    int      *sizes, item;
    Status   status;

    if( io_flag == WRITE_FILE )
    {
        if( getenv( "NEW_OBJ_FORMAT" ) != NULL )
        {
            ALLOC( sizes, n_items );
            sizes[0] = (*end_indices)[0];
            for_less( item, 1, n_items )
                sizes[item] = (*end_indices)[item] - (*end_indices)[item-1];

            status = io_ints( file, io_flag, format, n_items, &sizes );
 
            FREE( sizes );
        }
        else
            status = io_ints( file, io_flag, format, n_items, end_indices );
    }
    else
    {
        status = io_ints( file, io_flag, format, n_items, end_indices );

        if( status != OK )
            return( status );

        for_less( item, 1, n_items )
        {
            if( (*end_indices)[item] - (*end_indices)[item-1] < min_size )
                break;
        }

        if( item < n_items )
        {
            for_less( item, 1, n_items )
                (*end_indices)[item] += (*end_indices)[item-1];
        }
    }

    return( status );
}
