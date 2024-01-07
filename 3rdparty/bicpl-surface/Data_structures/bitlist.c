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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Data_structures/bitlist.c,v 1.13 2005/08/17 22:31:12 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_bitlist
@INPUT      : n_bits
@OUTPUT     : bitlist
@RETURNS    : Status
@DESCRIPTION: Creates a bitlist containing the required number of bits.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_bitlist(
    int             n_bits,
    bitlist_struct  *bitlist )
{
    bitlist->n_words = (n_bits + BITS_PER_BITLIST_WORD - 1) /
                        (BITS_PER_BITLIST_WORD);

    if( bitlist->n_words > 0 )
        ALLOC( bitlist->bits, bitlist->n_words );

    zero_bitlist( bitlist );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : zero_bitlist
@INPUT      : bitlist
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets all bits to zero (FALSE).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  zero_bitlist(
    bitlist_struct  *bitlist )
{
    int   i;

    for_less( i, 0, bitlist->n_words )
        bitlist->bits[i] = (bitlist_type) 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : fill_bitlist
@INPUT      : bitlist
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets all bits to 1 (TRUE).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  fill_bitlist(
    bitlist_struct  *bitlist )
{
    int   i;

    for_less( i, 0, bitlist->n_words )
        bitlist->bits[i] = (bitlist_type) -1;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_bitlist_bit
@INPUT      : bitlist
            : bit_index
            : value
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the bit_index'th bit to the given value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void     set_bitlist_bit(
    bitlist_struct  *bitlist,
    int             bit_index,
    BOOLEAN         value )
{
    int           word_index, bit_sub_index;
    bitlist_type  bit;

    word_index = bit_index >> LOG_BITS_PER_BITLIST_WORD;
    bit_sub_index = bit_index & BITS_PER_BITLIST_WORD_MINUS_1;

    bit = (bitlist_type) 1 << (bitlist_type) bit_sub_index;

    if( value )                                    /* setting bit to 1 */
        bitlist->bits[word_index] |= bit;
    else if( bitlist->bits[word_index] & bit )     /* setting bit to 0 */
            bitlist->bits[word_index] ^= bit;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_bitlist_bit
@INPUT      : bitlist
            : bit_index
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns the bit_index'th bit.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN     get_bitlist_bit(
    bitlist_struct  *bitlist,
    int             bit_index )
{
    int           word_index, bit_sub_index;
    bitlist_type  word_bits;
    BOOLEAN       value;

    word_index = bit_index >> LOG_BITS_PER_BITLIST_WORD;
    bit_sub_index = bit_index & BITS_PER_BITLIST_WORD_MINUS_1;

    word_bits = bitlist->bits[word_index];

    value = ((word_bits >> (bitlist_type) bit_sub_index) & 1) != 0;

    return( value );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_bitlist
@INPUT      : bitlist
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Deletes the given bitlist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_bitlist(
    bitlist_struct  *bitlist )
{
    if( bitlist->n_words > 0 )
    {
        FREE( bitlist->bits );
    }
}

/* ############################### bitlist 3d ############################### */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_bitlist_3d
@INPUT      : nx
            : ny
            : nz
@OUTPUT     : bitlist
@RETURNS    : Status
@DESCRIPTION: Creates a 3D bitlist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_bitlist_3d(
    int                nx,
    int                ny,
    int                nz,
    bitlist_3d_struct  *bitlist )
{
    bitlist->nx = nx;
    bitlist->ny = ny;
    bitlist->nz = nz;
    bitlist->n_z_words = (nz + BITS_PER_BITLIST_WORD - 1) /
                         BITS_PER_BITLIST_WORD;

    if( nx > 0 && ny > 0 && nz > 0 )
        ALLOC3D( bitlist->bits, nx, ny, bitlist->n_z_words );

    zero_bitlist_3d( bitlist );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : zero_bitlist_3d
@INPUT      : bitlist
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets all bits to 0 (FALSE).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  zero_bitlist_3d(
    bitlist_3d_struct  *bitlist )
{
    int   x, y, z;

    for_less( x, 0, bitlist->nx )
    {
        for_less( y, 0, bitlist->ny )
        {
            for_less( z, 0, bitlist->n_z_words )
            {
                bitlist->bits[x][y][z] = (bitlist_type) 0;
            }
        }
    }
}

BICAPI  void  get_bitlist_3d_sizes(
    bitlist_3d_struct  *bitlist,
    int                *nx,
    int                *ny,
    int                *nz )
{
    *nx = bitlist->nx;
    *ny = bitlist->ny;
    *nz = bitlist->nz;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : fill_bitlist_3d
@INPUT      : bitlist
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets all bits to 1 (TRUE).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  fill_bitlist_3d(
    bitlist_3d_struct  *bitlist )
{
    int   x, y, z;

    for_less( x, 0, bitlist->nx )
    {
        for_less( y, 0, bitlist->ny )
        {
            for_less( z, 0, bitlist->n_z_words )
            {
                bitlist->bits[x][y][z] = (bitlist_type) -1;
            }
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_bitlist_bit_3d
@INPUT      : bitlist
            : x
            : y
            : z
            : value
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the bit indicated by x y z to the given value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void     set_bitlist_bit_3d(
    bitlist_3d_struct  *bitlist,
    int                x,
    int                y,
    int                z,
    BOOLEAN            value )
{
    int           word_index, bit_sub_index;
    bitlist_type  bit;

    word_index = z >> LOG_BITS_PER_BITLIST_WORD;
    bit_sub_index = z & BITS_PER_BITLIST_WORD_MINUS_1;

    bit = (bitlist_type) 1 << (bitlist_type) bit_sub_index;

    if( value )
    {
        bitlist->bits[x][y][word_index] |= bit;
    }
    else
    {
        if( bitlist->bits[x][y][word_index] & bit )
        {
            bitlist->bits[x][y][word_index] ^= bit;
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_bitlist_bit_3d
@INPUT      : bitlist
            : x
            : y
            : z
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns the bit indicated by x, y, z.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN     get_bitlist_bit_3d(
    bitlist_3d_struct  *bitlist,
    int                x,
    int                y,
    int                z )
{
    int           word_index, bit_sub_index;
    bitlist_type  word;

    bit_sub_index = z & BITS_PER_BITLIST_WORD_MINUS_1;

    word_index = z >> LOG_BITS_PER_BITLIST_WORD;

    word = bitlist->bits[x][y][word_index];

    if( word == 0 )
        return( FALSE );
    else
        return( ((word >> (bitlist_type) bit_sub_index) & 1) != 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_bitlist_3d
@INPUT      : bitlist
@OUTPUT     : 
@RETURNS    : Status
@DESCRIPTION: Deletes the 3D bitlist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_bitlist_3d(
    bitlist_3d_struct  *bitlist )
{;
    if( bitlist->nx > 0 && bitlist->ny > 0 && bitlist->nz > 0 )
    {
        FREE3D( bitlist->bits );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_bitlist_3d
@INPUT      : file
            : io_type         READ_FILE or WRITE_FILE 
            : bitlist
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Inputs or outputs the 3D bitlist.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status  io_bitlist_3d(
    FILE               *file,
    IO_types           io_type,
    bitlist_3d_struct  *bitlist )
{
    Status  status;
    int     x, y;

    status = OK;

    for_less( x, 0, bitlist->nx )
    {
        for_less( y, 0, bitlist->ny )
        {
            status = io_binary_data( file, io_type,
                                     (void *) (bitlist->bits[x][y]),
                                     sizeof( bitlist->bits[x][y] ),
                                     bitlist->n_z_words );
       
            if( status != OK )
                break;
        }
        if( status != OK )
            break;
    }

    return( status );
}
