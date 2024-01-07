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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/pixels.c,v 1.27 2005/08/17 22:28:27 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_pixels
@INPUT      : x_position
              y_position
              x_size
              y_size
              x_zoom
              y_zoom
              pixel_type
@OUTPUT     : pixels
@RETURNS    : 
@DESCRIPTION: Initializes the pixels struct and allocates memory as needed.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_pixels(
    pixels_struct  *pixels,
    int            x_position,
    int            y_position,
    int            x_size,
    int            y_size,
    Real           x_zoom,
    Real           y_zoom,
    Pixel_types    pixel_type )
{
    int            n_alloced;

    n_alloced = 0;

    pixels->pixel_type = pixel_type;
    pixels->x_size = -1;
    pixels->y_size = -1;
    pixels->x_position = x_position;
    pixels->y_position = y_position;
    pixels->x_zoom = x_zoom;
    pixels->y_zoom = y_zoom;

    modify_pixels_size( &n_alloced, pixels, x_size, y_size, pixel_type );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_pixels
@INPUT      : pixels
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the pixels.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_pixels( pixels_struct *pixels )
{
    if( pixels->x_size > 0 && pixels->y_size > 0 )
    {
        switch( pixels->pixel_type )
        {
        case COLOUR_INDEX_8BIT_PIXEL:
            FREE( pixels->data.pixels_8bit_colour_index );
            break;

        case COLOUR_INDEX_16BIT_PIXEL:
            FREE( pixels->data.pixels_16bit_colour_index );
            break;

        case RGB_PIXEL:
            FREE( pixels->data.pixels_rgb );
            break;
        }
    }

    pixels->x_size = 0;
    pixels->y_size = 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : modify_pixels_size
@INPUT      : n_pixels_alloced
              pixels
              x_size
              y_size
              pixel_type
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Changes the size of the pixels struct, reallocating as necessary.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  modify_pixels_size(
    int            *n_pixels_alloced,
    pixels_struct  *pixels,
    int            x_size,
    int            y_size,
    Pixel_types    pixel_type )
{
    int                new_n_pixels;
    Colour             *tmp_int_ptr;

    if( pixel_type != pixels->pixel_type )
    {
        if( *n_pixels_alloced > 0 )
            delete_pixels( pixels );
        *n_pixels_alloced = 0;
        pixels->pixel_type = pixel_type;
    }

    if( x_size > 0 && y_size > 0 )
        new_n_pixels = x_size * y_size;
    else
        new_n_pixels = 0;

    pixels->x_size = x_size;
    pixels->y_size = y_size;

    /*--- if the new size is equal to or smaller than the current size */

    if( new_n_pixels > 0 && new_n_pixels <= *n_pixels_alloced )
        return;

    switch( pixel_type )
    {
    case COLOUR_INDEX_8BIT_PIXEL:    /* make sure it's word aligned */
        tmp_int_ptr = (Colour *) (void *)
                          pixels->data.pixels_8bit_colour_index;

        SET_ARRAY_SIZE( tmp_int_ptr, (*n_pixels_alloced + 3) / 4,
                        (new_n_pixels + 3) / 4, DEFAULT_CHUNK_SIZE );
        pixels->data.pixels_8bit_colour_index = (unsigned char *)
                         ((void *) tmp_int_ptr );
        break;

    case COLOUR_INDEX_16BIT_PIXEL:    /* make sure it's word aligned */
        tmp_int_ptr = (Colour *) (void *)
                          pixels->data.pixels_16bit_colour_index;
        SET_ARRAY_SIZE( tmp_int_ptr, (*n_pixels_alloced + 1) / 2,
                        (new_n_pixels + 1) / 2, DEFAULT_CHUNK_SIZE );
        pixels->data.pixels_16bit_colour_index = (unsigned short *)
                         ((void *) tmp_int_ptr );
        break;

    case RGB_PIXEL:
        SET_ARRAY_SIZE( pixels->data.pixels_rgb, *n_pixels_alloced,
                        new_n_pixels, DEFAULT_CHUNK_SIZE );
        break;
    }

    *n_pixels_alloced = new_n_pixels;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_pixels24_to_pixels8
@INPUT      : pixels_rgb
@OUTPUT     : pixels_8
@RETURNS    : 
@DESCRIPTION: Converts the 24 bit pixels to an 8 bit (3-3-2 bit) lookup.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_pixels24_to_pixels8(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8 )
{
    int   x, y;

    initialize_pixels( pixels_8, pixels_rgb->x_position, pixels_rgb->y_position,
                       pixels_rgb->x_size, pixels_rgb->y_size,
                       pixels_rgb->x_zoom, pixels_rgb->y_zoom,
                       COLOUR_INDEX_8BIT_PIXEL );

    for_less( x, 0, pixels_rgb->x_size )
    {
        for_less( y, 0, pixels_rgb->y_size )
        {
            PIXEL_COLOUR_INDEX_8(*pixels_8,x,y) =
                   (unsigned char) convert_rgb_pixel_to_8bit_lookup(
                                PIXEL_RGB_COLOUR(*pixels_rgb,x,y) );
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_pixels24_to_index8
@INPUT      : pixels_rgb
              n_colours
              colour_table
@OUTPUT     : pixels_8
@RETURNS    : 
@DESCRIPTION: Converts the rgb pixels to the 8-bit colour index pixels,
              given a colour table.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_pixels24_to_index8(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8,
    int             n_colours,
    Colour          colour_table[] )
{
    Colour  col;
    int     x, y;

    initialize_pixels( pixels_8, pixels_rgb->x_position, pixels_rgb->y_position,
                       pixels_rgb->x_size, pixels_rgb->y_size,
                       pixels_rgb->x_zoom, pixels_rgb->y_zoom,
                       COLOUR_INDEX_8BIT_PIXEL );

    for_less( x, 0, pixels_rgb->x_size )
    {
        for_less( y, 0, pixels_rgb->y_size )
        {
            col = PIXEL_RGB_COLOUR(*pixels_rgb,x,y);

            PIXEL_COLOUR_INDEX_8(*pixels_8,x,y) =
                   (unsigned char) find_closest_colour( get_Colour_r(col),
                                        get_Colour_g(col),
                                        get_Colour_b(col),
                                        n_colours, colour_table );
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_index8_to_pixels24
@INPUT      : pixels_8
              colour_table
@OUTPUT     : pixels_rgb
@RETURNS    : 
@DESCRIPTION: Creates an rgb pixel map from 8 bit indices and a colour table.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_index8_to_pixels24(
    pixels_struct   *pixels_8,
    Colour          colour_table[],
    pixels_struct   *pixels_rgb )
{
    int     x, y;

    initialize_pixels( pixels_rgb, pixels_8->x_position, pixels_8->y_position,
                       pixels_8->x_size, pixels_8->y_size,
                       pixels_8->x_zoom, pixels_8->y_zoom, RGB_PIXEL );

    for_less( x, 0, pixels_rgb->x_size )
    {
        for_less( y, 0, pixels_rgb->y_size )
        {
            PIXEL_RGB_COLOUR(*pixels_rgb,x,y) = colour_table[
                            (int) PIXEL_COLOUR_INDEX_8(*pixels_8,x,y)];
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_pixels24_to_gray_scale
@INPUT      : pixels_rgb
@OUTPUT     : pixels_8
@RETURNS    : 
@DESCRIPTION: Converts 24 bit pixels to 8 bit gray scale.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_pixels24_to_gray_scale(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8 )
{
    int   x, y;

    initialize_pixels( pixels_8, 0, 0, pixels_rgb->x_size, pixels_rgb->y_size,
                       1.0, 1.0,
                       COLOUR_INDEX_8BIT_PIXEL );

    for_less( x, 0, pixels_rgb->x_size )
    {
        for_less( y, 0, pixels_rgb->y_size )
        {
            PIXEL_COLOUR_INDEX_8(*pixels_8,x,y) = (unsigned char)
                   get_Colour_luminance( PIXEL_RGB_COLOUR(*pixels_rgb,x,y) );
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_pixels24_to_dithered
@INPUT      : pixels_rgb
              n_colours
              colour_table
@OUTPUT     : pixels_8
@RETURNS    : 
@DESCRIPTION: Performs Floyd-Steinberg error diffusion to convert a 24 bit
              image to an 8 bit image with the given colour map.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_pixels24_to_dithered(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8,
    int             n_colours,
    Colour          colour_table[] )
{
    int      x, y, dir, error_row, c, comp[3], used[3];
    int      first_x, last_x, ind, one, three, five, seven, error;
    Colour   col;
    short    ***errors;

    initialize_pixels( pixels_8, 0, 0, pixels_rgb->x_size, pixels_rgb->y_size,
                       1.0, 1.0,
                       COLOUR_INDEX_8BIT_PIXEL );

    ALLOC3D( errors, 2, 3, pixels_rgb->x_size );

    for_less( x, 0, pixels_rgb->x_size )
    {
        errors[0][0][x] = 0;
        errors[0][1][x] = 0;
        errors[0][2][x] = 0;
        errors[1][0][x] = 0;
        errors[1][1][x] = 0;
        errors[1][2][x] = 0;
    }

    error_row = 0;
    dir = 1;

    for_less( y, 0, pixels_rgb->y_size )
    {
        if( dir > 0 )
        {
            first_x = 0;
            last_x = pixels_rgb->x_size;
        }
        else
        {
            first_x = pixels_rgb->x_size-1;
            last_x = -1;
        }
        for( x = first_x;  x != last_x;  x += dir )
        {
            col = PIXEL_RGB_COLOUR(*pixels_rgb,x,y);
            comp[0] = get_Colour_r( col ) + errors[error_row][0][x];
            comp[1] = get_Colour_g( col ) + errors[error_row][1][x];
            comp[2] = get_Colour_b( col ) + errors[error_row][2][x];

            ind = find_closest_colour( comp[0], comp[1], comp[2],
                                       n_colours, colour_table );

            used[0] = get_Colour_r( colour_table[ind] );
            used[1] = get_Colour_g( colour_table[ind] );
            used[2] = get_Colour_b( colour_table[ind] );

            for_less( c, 0, 3 )
            {
                error = comp[c] - used[c];
                seven = ROUND( (Real) error * 7.0/ 16.0 );
                five = ROUND( (Real) error * 5.0 / 16.0 );
                three = ROUND( (Real) error * 3.0 / 16.0 );
                one = error - seven - five - three;

                if( x + dir >= 0 && x + dir < pixels_rgb->x_size )
                {
                    errors[error_row][c][x+dir] += seven;
                    errors[1-error_row][c][x+dir] += one;
                }

                if( x - dir >= 0 && x - dir < pixels_rgb->x_size )
                    errors[1-error_row][c][x-dir] += three;

                errors[1-error_row][c][x] += five;
            }

            PIXEL_COLOUR_INDEX_8(*pixels_8,x,y) = (unsigned char) ind;
        }

        error_row = 1 - error_row;
        dir = - dir;

        for_less( x, 0, pixels_rgb->x_size )
        {
            errors[1-error_row][0][x] = 0;
            errors[1-error_row][1][x] = 0;
            errors[1-error_row][2][x] = 0;
        }
    }

    FREE3D( errors );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : resample_pixels
@INPUT      : pixels
              transform
              n_samples
              background_colour
@OUTPUT     : new_pixels
@RETURNS    : 
@DESCRIPTION: Resamples the pixels in the transformed space.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  resample_pixels(
    pixels_struct   *pixels,
    Transform_2d    *transform,
    int             n_samples,
    pixels_struct   *new_pixels,
    Colour          background_colour )
{
    int           x, y, i, j;
    Real          xs, ys, r_sum, g_sum, b_sum;
    Real          x_trans, y_trans, weight;
    Colour        colour;
    Transform_2d  inverse;

    if( pixels->pixel_type != RGB_PIXEL ||
        new_pixels->pixel_type != RGB_PIXEL )
    {
        print_error( "Can only resample 24 bit rgb pixels.\n" );
        return;
    }

    weight = 1.0 / (Real) n_samples / (Real) n_samples;
    get_inverse_transform_2d( transform, &inverse );

    for_less( x, 0, new_pixels->x_size )
    {
        for_less( y, 0, new_pixels->y_size )
        {
            if( n_samples == 1 )
            {
                transform_point_2d( &inverse, (Real) x, (Real) y,
                                    &x_trans, &y_trans );
                if( x_trans >= -0.5 && x_trans < (Real) pixels->x_size - 0.5 &&
                    y_trans >= -0.5 && y_trans < (Real) pixels->y_size - 0.5 )
                {
                    colour = PIXEL_RGB_COLOUR(*pixels,ROUND(x_trans),
                                                      ROUND(y_trans));
                }
                else
                {
                    colour = background_colour;
                }
            }
            else
            {
                r_sum = 0.0;
                g_sum = 0.0;
                b_sum = 0.0;
                for_less( i, 0, n_samples )
                {
                    xs = (Real) x - 0.5 + ((Real) i + 0.5) / (Real) n_samples;

                    for_less( j, 0, n_samples )
                    {
                        ys = (Real) y - 0.5 + ((Real) j + 0.5) /(Real)n_samples;
                        transform_point_2d( &inverse, xs, ys,
                                            &x_trans, &y_trans );
                        if( x_trans >= -0.5 &&
                            x_trans < (Real) pixels->x_size - 0.5 &&
                            y_trans >= -0.5 &&
                            y_trans < (Real) pixels->y_size - 0.5 )
                        {
                            colour = PIXEL_RGB_COLOUR(*pixels,ROUND(x_trans),
                                                              ROUND(y_trans));
                        }
                        else
                        {
                            colour = background_colour;
                        }
                        r_sum += (Real) get_Colour_r(colour);
                        g_sum += (Real) get_Colour_g(colour);
                        b_sum += (Real) get_Colour_b(colour);
                    }
                }

                colour = make_Colour( (int) (r_sum * weight),
                                      (int) (g_sum * weight),
                                      (int) (b_sum * weight) );
            }

            PIXEL_RGB_COLOUR( *new_pixels, x, y ) = colour;
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : copy_pixel_region
@INPUT      : pixels
              x_min
              x_max
              y_min
              y_max
@OUTPUT     : new_pixels
@RETURNS    : 
@DESCRIPTION: Copies the subset of the pixels, creating new_pixels.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  copy_pixel_region(
    pixels_struct   *pixels,
    int             x_min,
    int             x_max,
    int             y_min,
    int             y_max,
    pixels_struct   *new_pixels )
{
    int           y, y_dest, x_size, y_size;

    x_min -= pixels->x_position;
    y_min -= pixels->y_position;
    x_max -= pixels->x_position;
    y_max -= pixels->y_position;

    if( x_min < 0 )
        x_min = 0;
    if( x_max >= pixels->x_size )
        x_max = pixels->x_size-1;

    if( y_min < 0 )
        y_min = 0;
    if( y_max >= pixels->y_size )
        y_max = pixels->y_size-1;

    x_size = MAX( 0, x_max - x_min + 1 );
    y_size = MAX( 0, y_max - y_min + 1 );

    initialize_pixels( new_pixels, pixels->x_position + x_min,
                       pixels->y_position + y_min,
                       x_size, y_size, pixels->x_zoom,
                       pixels->y_zoom, pixels->pixel_type );

    if( x_size == 0 || y_size == 0 )
        return;

    for_inclusive( y, y_min, y_max )
    {
        y_dest = y - y_min;

        switch( pixels->pixel_type )
        {
        case COLOUR_INDEX_8BIT_PIXEL:
            (void) memcpy( (void *) &PIXEL_COLOUR_INDEX_8(*new_pixels,0,y_dest),
                           (void *) &PIXEL_COLOUR_INDEX_8(*pixels,x_min,y),
                           (size_t) x_size *
                               sizeof(PIXEL_COLOUR_INDEX_8(*new_pixels,0,0)) );
            break;

        case COLOUR_INDEX_16BIT_PIXEL:
            (void) memcpy( (void *)&PIXEL_COLOUR_INDEX_16(*new_pixels,0,y_dest),
                           (void *) &PIXEL_COLOUR_INDEX_16(*pixels,x_min,y),
                           (size_t) x_size *
                               sizeof(PIXEL_COLOUR_INDEX_16(*new_pixels,0,0)) );
            break;

        case RGB_PIXEL:
            (void) memcpy( (void *) &PIXEL_RGB_COLOUR(*new_pixels,0,y_dest),
                           (void *) &PIXEL_RGB_COLOUR(*pixels,x_min,y),
                           (size_t) x_size *
                               sizeof(PIXEL_RGB_COLOUR(*new_pixels,0,0)) );
            break;
        }
    }
}
