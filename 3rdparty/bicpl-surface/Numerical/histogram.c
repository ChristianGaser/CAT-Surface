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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Numerical/histogram.c,v 1.12 2005/08/17 22:28:59 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_histogram
@INPUT      : delta
              offset
@OUTPUT     : histogram
@RETURNS    : 
@DESCRIPTION: Initializes a histogram of boxes whose boundaries are integer
              multiples of delta from offset.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_histogram(
    histogram_struct  *histogram,
    Real              delta,
    Real              offset )
{
    histogram->delta = delta;
    histogram->offset = offset;

    histogram->min_index = 0;
    histogram->max_index = -1;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_histogram
@INPUT      : histogram
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the histogram.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_histogram(
    histogram_struct  *histogram )
{
    if( histogram->max_index >= histogram->min_index )
    {
        FREE( histogram->counts );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_histogram_index
@INPUT      : histogram
              value
@OUTPUT     : 
@RETURNS    : box index
@DESCRIPTION: Converts a real value to a histogram box index, which can be
              positive, zero, or negative.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_histogram_index(
    histogram_struct  *histogram,
    Real              value )
{
    int    ind;

    ind = FLOOR( (value - histogram->offset) / histogram->delta );

    return( ind );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_real_index_to_value
@INPUT      : histogram
              ind
@OUTPUT     : 
@RETURNS    : value
@DESCRIPTION: Converts a histogram box index to a real value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  convert_real_index_to_value(
    histogram_struct  *histogram,
    Real              ind )
{
    return( ind * histogram->delta + histogram->offset );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_index_to_value
@INPUT      : histogram
              ind
@OUTPUT     : 
@RETURNS    : value
@DESCRIPTION: Converts a histogram box index to a real value, which is the
              left edge of the ind'th box.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  convert_index_to_value(
    histogram_struct  *histogram,
    int               ind )
{
    return( convert_real_index_to_value( histogram, (Real) ind ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_to_histogram
@INPUT      : histogram
              value
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Adds a value to the histogram, by finding the box it belongs in,
              expanding the list of histogram boxes, if needed, and incrementing
              the box count for that box.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  add_to_histogram(
    histogram_struct  *histogram,
    Real              value )
{
    int   ind, i, prev_size, new_size;

    ind = get_histogram_index( histogram, value );

    if( histogram->min_index > histogram->max_index )   /* first box */
    {
        ALLOC( histogram->counts, 1 );
        histogram->counts[0] = 0;
        histogram->min_index = ind;
        histogram->max_index = ind;
    }
    else if( ind < histogram->min_index )   /* need to expand below */
    {
        prev_size = histogram->max_index - histogram->min_index + 1;
        new_size = histogram->max_index - ind + 1;
        SET_ARRAY_SIZE( histogram->counts, prev_size, new_size, 1 );

        for_down( i, histogram->max_index, histogram->min_index )
        {
            histogram->counts[i-ind] =
                         histogram->counts[i-histogram->min_index];
        }

        for_less( i, ind, histogram->min_index )
            histogram->counts[i-ind] = 0;

        histogram->min_index = ind;
    }
    else if( ind > histogram->max_index )    /* need to expand above */
    {
        prev_size = histogram->max_index - histogram->min_index + 1;
        new_size = ind - histogram->min_index + 1;
        SET_ARRAY_SIZE( histogram->counts, prev_size, new_size, 1 );

        for_inclusive( i, histogram->max_index + 1, ind )
            histogram->counts[i-histogram->min_index] = 0;

        histogram->max_index = ind;
    }

    ++histogram->counts[ind-histogram->min_index];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_histogram_range
@INPUT      : histogram
@OUTPUT     : min_value
              max_value
@RETURNS    : 
@DESCRIPTION: Passes back the range of the histogram, which is the left
              edge of the first box, and the right edge of the last box.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_histogram_range(
    histogram_struct  *histogram,
    Real              *min_value,
    Real              *max_value )
{
    *min_value = convert_index_to_value( histogram, histogram->min_index );
    *max_value = convert_index_to_value( histogram, histogram->max_index + 1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_histogram_max_count
@INPUT      : histogram
@OUTPUT     : 
@RETURNS    : maximum count
@DESCRIPTION: Finds the maximum number of counts in the histogram boxes.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_histogram_max_count(
    histogram_struct  *histogram )
{
    int   ind, max_count;

    max_count = 0;
    for_inclusive( ind, histogram->min_index, histogram->max_index )
    {
         if( ind == 0 ||
             histogram->counts[ind-histogram->min_index] > max_count )
             max_count = histogram->counts[ind-histogram->min_index];
    }

    return( max_count );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : box_filter_histogram
@INPUT      : n
              counts
              half_width
@OUTPUT     : new_counts
@RETURNS    : 
@DESCRIPTION: Box filters the histogram boxes into another set of boxes
              which corresponds to the average number of counts in the
              range i - half_width to i + half_width.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  box_filter_histogram(
    int          n,
    Real         counts[],
    Real         new_counts[],
    int          half_width )
{
    int    i, window_width, start_index, end_index;
    Real   current_value;

    start_index = - half_width;
    end_index   = half_width;

    current_value = 0.0;
    for_inclusive( i, 0, MIN( end_index, n-1 ) )
        current_value += counts[i];

    for_less( i, 0, n )
    {
        window_width = MIN( end_index, n-1 ) - MAX( start_index, 0 ) + 1;
        new_counts[i] = (Real) current_value / (Real) window_width;
        if( start_index >= 0 )
            current_value -= counts[start_index];
        ++start_index;
        ++end_index;
        if( end_index < n )
            current_value += counts[end_index];
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_histogram_counts
@INPUT      : histogram
              filter_width
@OUTPUT     : counts
              scale_factor
              trans_factor
@RETURNS    : n
@DESCRIPTION: Gets the histogram counts, after filtering as necessary.
              and a scale and trans factor for converting counts indices
              to real values.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_histogram_counts(
    histogram_struct  *histogram,
    Real              *counts[],
    Real              filter_width,
    Real              *scale_factor,
    Real              *trans_factor )
{
    int    i, n, width;
    Real   *tmp_counts;

    n = histogram->max_index - histogram->min_index + 1;

    if( n <= 0 )
        return( 0 );

    ALLOC( tmp_counts, n );
    ALLOC( *counts, n );

    for_less( i, 0, n )
        tmp_counts[i] = (Real) histogram->counts[i];

    width = ROUND( filter_width / histogram->delta / 2.0 );

    box_filter_histogram( n, tmp_counts, *counts, width );

    FREE( tmp_counts );

    *scale_factor = histogram->delta;
    *trans_factor = convert_real_index_to_value( histogram,
                                          (Real) histogram->min_index + 0.5 );

    return( n );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : resample_histogram
@INPUT      : histogram
              x_size
              y_size
@OUTPUT     : x_scale
              y_scale
              height
@RETURNS    : 
@DESCRIPTION: Resamples the histogram to fit in a grid of x_size boxes by
              max y_size counts.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  resample_histogram(
    histogram_struct  *histogram,
    int               x_size,
    int               y_size,
    Real              *x_scale,
    Real              *x_trans,
    Real              height[] )
{
    int   ind, x, max_count, left_ind, right_ind;
    Real  weight, sum_count, min_value, max_value;
    Real  left_side, right_side, left, right;

    get_histogram_range( histogram, &min_value, &max_value );

    if( max_value <= min_value )
    {
        for_less( x, 0, x_size )
            height[x] = 0.0;
        *x_scale = 1.0;
        *x_trans = 0.0;
        return;
    }

    max_count = get_histogram_max_count( histogram );

    *x_scale = 1.0 / (Real) x_size * (max_value - min_value);
    *x_trans = min_value + 0.5 * (*x_scale);

    for_less( x, 0, x_size )
    {
        left = min_value + (Real) x / (Real) x_size * (max_value - min_value);
        right = min_value + (Real) (x+1) / (Real) x_size *
                (max_value - min_value);

        left_ind = get_histogram_index( histogram, left );
        right_ind = get_histogram_index( histogram, right );

        sum_count = 0.0;

        for_inclusive( ind, left_ind, right_ind )
        {
            left_side = convert_index_to_value( histogram, ind );
            right_side = convert_index_to_value( histogram, ind+1 );

            weight = MIN( right_side, right ) - MAX( left_side, left );

            sum_count += weight *
                         (Real) histogram->counts[ind-histogram->min_index];
        }

        height[x] = sum_count / (right - left) * (Real) y_size /
                    (Real) max_count;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : display_histogram
@INPUT      : histogram
              x_size
              y_size
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Displays the histogram as an ascii plot.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  display_histogram(
    histogram_struct  *histogram,
    int               x_size,
    int               y_size )
{
    int   x, y, max_count;
    Real  min_value, max_value, *n_chars, x_scale, x_trans;

    ALLOC( n_chars, x_size );

    resample_histogram( histogram, x_size, y_size, &x_scale, &x_trans,
                        n_chars );

    for( y = y_size-1;  y >= 0;  --y )
    {
        for_less( x, 0, x_size )
        {
            if( y < ROUND(n_chars[x]) )
                 print( "#" );
            else
                 print( " " );
        }

        print( "\n" );
    }

    get_histogram_range( histogram, &min_value, &max_value );
    max_count = get_histogram_max_count( histogram );

    print( "%g to %g with max count = %d\n", min_value, max_value, max_count );

    FREE( n_chars );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_histogram_line
@INPUT      : histogram
              x_size
              y_size
              filter_width
@OUTPUT     : lines
@RETURNS    : 
@DESCRIPTION: Creates a lines structure corresponding to the histogram.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_histogram_line(
    histogram_struct  *histogram,
    int               x_size,
    int               y_size,
    Real              filter_width,
    lines_struct      *lines )
{
    int     x, width;
    Real    *height, *smooth_height, x_scale, x_trans;
    Point   p;

    ALLOC( height, x_size );

    resample_histogram( histogram, x_size, y_size, &x_scale, &x_trans, height );

    width = ROUND( filter_width / x_scale / 2.0 );

    ALLOC( smooth_height, x_size );

    box_filter_histogram( x_size, height, smooth_height, width );

    initialize_lines( lines, WHITE );

    for_less( x, 0, x_size )
    {
        fill_Point( p, (Real) x * x_scale + x_trans, smooth_height[x], 0.0 );
        add_point_to_line( lines, &p );
    }

    FREE( height );
    FREE( smooth_height );
}
