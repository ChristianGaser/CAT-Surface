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

#include  <internal_volume_io.h>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/minc/volume_io/Geometry/colour.c,v 1.3.2.1 2004/10/04 20:18:40 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_Colour
@INPUT      : r
              g
              b
@OUTPUT     : 
@RETURNS    : Colour
@DESCRIPTION: Packs the three components, which are in the range 0 to 255,
              into a Colour type, unsigned long.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Colour  make_Colour(
    int   r,
    int   g,
    int   b )
{
    return( make_rgba_Colour( r, g, b, 255 ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_Colour_r_0_1
@INPUT      : colour
@OUTPUT     : 
@RETURNS    : red component
@DESCRIPTION: Returns the red component of the colour in the range of 0.0 to 1.0
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Real  get_Colour_r_0_1(
    Colour   colour )
{
    return( (Real) get_Colour_r(colour) / 255.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_Colour_g_0_1
@INPUT      : colour
@OUTPUT     : 
@RETURNS    : green component
@DESCRIPTION: Returns the green component of the colour in the range of
              0.0 to 1.0
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Real  get_Colour_g_0_1(
    Colour   colour )
{
    return( (Real) get_Colour_g(colour) / 255.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_Colour_b_0_1
@INPUT      : colour
@OUTPUT     : 
@RETURNS    : blue component
@DESCRIPTION: Returns the blue component of the colour in the range of
              0.0 to 1.0
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Real  get_Colour_b_0_1(
    Colour   colour )
{
    return( (Real) get_Colour_b(colour) / 255.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_Colour_a_0_1
@INPUT      : colour
@OUTPUT     : 
@RETURNS    : alpha component
@DESCRIPTION: Returns the alpha (opacity) component of the colour in the
              range of 0.0 to 1.0
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Real  get_Colour_a_0_1(
    Colour   colour )
{
    return( (Real) get_Colour_a(colour) / 255.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_Colour_0_1
@INPUT      : r
              g
              b
@OUTPUT     : 
@RETURNS    : Colour
@DESCRIPTION: Takes the three components, each in the range of 0 to 1,
              and packs them into a colour.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Colour  make_Colour_0_1(
    Real   r,
    Real   g,
    Real   b )
{
    return( make_Colour( (int) (r * 255.0 + 0.5),
                         (int) (g * 255.0 + 0.5),
                         (int) (b * 255.0 + 0.5) ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rgba_Colour_0_1
@INPUT      : r
              g
              b
              a       - alpha (opacity)
@OUTPUT     : 
@RETURNS    : Colour
@DESCRIPTION: Takes the four components, each in the range of 0 to 1,
              and packs them into a colour.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 10, 1995        D. MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Colour  make_rgba_Colour_0_1(
    Real   r,
    Real   g,
    Real   b,
    Real   a )
{
    return( make_rgba_Colour( (int) (r * 255.0 + 0.5),
                              (int) (g * 255.0 + 0.5),
                              (int) (b * 255.0 + 0.5),
                              (int) (a * 255.0 + 0.5) ) );
}
