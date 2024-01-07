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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/rgb_lookup.c,v 1.13 2005/08/17 22:28:27 bert Exp $";
#endif

#define  NUM_PIXELS8   256

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_8bit_rgb_pixel_lookup
@INPUT      : 
@OUTPUT     : 
@RETURNS    : colour table
@DESCRIPTION: Returns the 8 bit (3-3-2 bit) RGB colour table.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Colour     *get_8bit_rgb_pixel_lookup( void )
{
    static  BOOLEAN       initialized = FALSE;
    static  Colour        lookup_table[NUM_PIXELS8];
    int                   r, g, b, ind, r255, g255, b255;

    if( !initialized )
    {
        initialized = TRUE;

        for_less( r, 0, 8 )
        {
            for_less( g, 0, 8 )
            {
                for_less( b, 0, 4 )
                {
                    ind = ( r | (g << 3) | (b << 6) );

                    r255 = (int) (((Real) r + 0.5) * 256.0 / 8.0);
                    g255 = (int) (((Real) g + 0.5) * 256.0 / 8.0);
                    b255 = (int) (((Real) b + 0.5) * 256.0 / 4.0);

                    lookup_table[ind] = make_Colour( r255, g255, b255 );
                }
            }
        }
    }

    return( lookup_table );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_rgb_pixel_to_8bit_lookup
@INPUT      : colour
@OUTPUT     : 
@RETURNS    : index
@DESCRIPTION: Converts an RGB pixel to an 8 bit index.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  convert_rgb_pixel_to_8bit_lookup(
    Colour    colour )
{
    int   ind, r, g, b, r_ind, g_ind, b_ind;

    r = get_Colour_r( colour );
    g = get_Colour_g( colour );
    b = get_Colour_b( colour );

    r_ind = 8 * r / 256;
    g_ind = 8 * g / 256;
    b_ind = 4 * b / 256;
    ind = r_ind | (g_ind << 3) | (b_ind << 6);

    return( ind );
}
