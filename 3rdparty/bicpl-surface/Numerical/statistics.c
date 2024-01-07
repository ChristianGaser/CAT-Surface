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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Numerical/statistics.c,v 1.15 2005/08/17 22:28:59 bert Exp $";
#endif

#define  DEFAULT_N_MEDIAN_BOXES    100000

#define  MAX_SAMPLES_RECORDED      100000

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_statistics
@INPUT      : n
              samples
@OUTPUT     : min_value
              max_value
              mean_value
              std_dev
              median     - may be null pointer if median not desired
@RETURNS    : 
@DESCRIPTION: Computes the most basic statistics on a set of n samples.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  compute_statistics(
    int      n,
    Real     samples[],
    Real     *min_value,
    Real     *max_value,
    Real     *mean_value,
    Real     *std_dev,
    Real     *median )
{
    int                   i;
    Real                  x, min_median, max_median, median_error;
    statistics_struct     stats;
    BOOLEAN               done;

    if( median != NULL )
    {
        min_median = 0.0;
        max_median = 0.0;

        for_less( i, 0, n )
        {
            x = samples[i];
            if( i == 0 )
            {
                min_median = x;
                max_median = x;
            }
            else if( x < min_median )
                min_median = x;
            else if( x > max_median )
                max_median = x;
        }
    }
    else
    {
        min_median = 0.0;
        max_median = -1.0;
    }

    initialize_statistics( &stats, min_median, max_median );

    done = FALSE;

    while( !done )
    {
        for_less( i, 0, n )
            add_sample_to_statistics( &stats, samples[i] );

        get_statistics( &stats, NULL, mean_value, median, &median_error,
                        min_value, max_value, std_dev );

        if( median != NULL && median_error > 0.0 )
            restart_statistics_with_narrower_median_range( &stats );
        else
            done = TRUE;
    }

    terminate_statistics( &stats );
}

BICAPI  void  initialize_statistics(
    statistics_struct  *stats,
    Real               median_lower_bound,
    Real               median_upper_bound )
{
    int    i;

    stats->n_samples = 0;
    stats->sum_x = 0.0;
    stats->sum_xx = 0.0;

    stats->min_median_range = median_lower_bound;
    stats->max_median_range = median_upper_bound;

    if( median_lower_bound < median_upper_bound )
    {
        stats->n_median_boxes = DEFAULT_N_MEDIAN_BOXES;
        ALLOC( stats->median_box_counts, stats->n_median_boxes );
        ALLOC( stats->median_box_values, stats->n_median_boxes );

        for_less( i, 0, stats->n_median_boxes )
            stats->median_box_counts[i] = 0;

        stats->n_below_median_range = 0;
        stats->n_above_median_range = 0;
    }
}

BICAPI  void  add_sample_to_statistics(
    statistics_struct  *stats,
    Real               sample )
{
    int         median_box;

    if( stats->n_samples == 0 )
    {
        stats->min_value = sample;
        stats->max_value = sample;
    }
    else if( sample < stats->min_value )
        stats->min_value = sample;
    else if( sample > stats->max_value )
        stats->max_value = sample;

    if( stats->n_samples < MAX_SAMPLES_RECORDED )
    {
        SET_ARRAY_SIZE( stats->samples, stats->n_samples,
                        stats->n_samples+1, DEFAULT_CHUNK_SIZE );
        stats->samples[stats->n_samples] = sample;
    }
    else if( stats->n_samples == MAX_SAMPLES_RECORDED )
        FREE( stats->samples );

    stats->n_samples += 1;

    stats->sum_x += sample;
    stats->sum_xx += sample * sample;

    if( stats->min_median_range < stats->max_median_range )
    {
        if( sample < stats->min_median_range )
            ++stats->n_below_median_range;
        else if( sample >= stats->max_median_range )
            ++stats->n_above_median_range;
        else
        {
            median_box = (int) ((Real) stats->n_median_boxes *
                           (sample - stats->min_median_range)/
                           (stats->max_median_range - stats->min_median_range));

            ++stats->median_box_counts[median_box];
            if( stats->median_box_counts[median_box] == 1 )
                stats->median_box_values[median_box] = sample;
        }
    }
}

static  void  get_median(
    statistics_struct  *stats,
    Real               *min_range,
    Real               *max_range )
{
    int   box_index, median_index;

    median_index = stats->n_samples / 2;

    if( stats->n_samples <= MAX_SAMPLES_RECORDED )
    {
        int  i, j, best_index;
        Real tmp;

        for_less( i, 0, stats->n_samples-1 )
        {
            best_index = i;
            for_less( j, i+1, stats->n_samples )
            {
                if( stats->samples[j] < stats->samples[best_index] )
                    best_index = j;
            }

            tmp = stats->samples[best_index];
            stats->samples[best_index] = stats->samples[i];
            stats->samples[i] = tmp;
        }

        *min_range = stats->samples[median_index];
        *max_range = stats->samples[median_index];

        return;
    }

    if( stats->min_median_range >= stats->max_median_range )
    {
        *min_range = -1.0e30;
        *max_range = 1.0e30;
        return;
    }

    if( median_index < stats->n_below_median_range )
    {
        *min_range = stats->min_value;
        *max_range = stats->min_median_range;
        return;
    }

    median_index -= stats->n_below_median_range;

    box_index = 0;
    while( box_index < stats->n_median_boxes &&
           median_index >= stats->median_box_counts[box_index] )
    {
        median_index -= stats->median_box_counts[box_index];
        ++box_index;
    }

    if( box_index == stats->n_median_boxes )
    {
        *min_range = stats->max_median_range;
        *max_range = stats->max_value;
    }
    else
    {
        if( stats->median_box_counts[box_index] == 1 )
        {
            *min_range = stats->median_box_values[box_index];
            *max_range = stats->median_box_values[box_index];
        }
        else
        {
            *min_range = stats->min_median_range +
                         (stats->max_median_range - stats->min_median_range) *
                          (Real) box_index / (Real) stats->n_median_boxes;
            *max_range = stats->min_median_range +
                         (stats->max_median_range - stats->min_median_range) *
                          (Real) (box_index+1) / (Real) stats->n_median_boxes;
        }
    }
}

BICAPI  void  restart_statistics_with_narrower_median_range(
    statistics_struct  *stats )
{
    Real   min_median_range, max_median_range;

    get_median( stats, &min_median_range, &max_median_range );

    if( min_median_range == max_median_range )
    {
        min_median_range = stats->min_median_range;
        max_median_range = stats->max_median_range;
        print_error( "Median range already narrow enough.\n" );
    }

    terminate_statistics( stats );

    initialize_statistics( stats, min_median_range, max_median_range );
}

BICAPI  void  get_statistics(
    statistics_struct  *stats,
    int                *n_samples,
    Real               *mean,
    Real               *median,
    Real               *median_error,
    Real               *min_value,
    Real               *max_value,
    Real               *std_deviation )
{
    int      n;
    Real     sum_x, sum_xx, min_median_range, max_median_range, variance;

    if( n_samples != NULL )
        *n_samples = stats->n_samples;

    if( stats->n_samples <= 0 )
    {
        if( median_error != NULL )
            *median_error = 0.0;
        return;
    }

    if( median != NULL )
    {
        get_median( stats, &min_median_range, &max_median_range );

        if( min_median_range == max_median_range )
        {
            *median = min_median_range;
            if( median_error != NULL )
                *median_error = 0.0;
        }
        else
        {
            *median = (min_median_range + max_median_range) / 2.0;
            if( median_error != NULL )
                *median_error = (max_median_range - min_median_range) / 2.0;
        }
    }

    if( min_value != NULL )
        *min_value = stats->min_value;

    if( max_value != NULL )
        *max_value = stats->max_value;

    n = stats->n_samples;
    sum_x = stats->sum_x;
    sum_xx = stats->sum_xx;

    if( mean != NULL )
        *mean = sum_x / (Real) n;

    if( n == 1 )
        variance = 0.0;
    else
        variance = (sum_xx - sum_x * sum_x / (Real) n) / (Real) (n - 1);

    if( std_deviation != NULL )
    {
        if( variance <= 0.0 )
            *std_deviation = 0.0;
        else
            *std_deviation = sqrt( variance );
    }
}

BICAPI  void  terminate_statistics(
    statistics_struct  *stats )
{
    if( stats->n_samples > 0 && stats->n_samples <= MAX_SAMPLES_RECORDED )
        FREE( stats->samples );

    if( stats->min_median_range < stats->max_median_range )
    {
        FREE( stats->median_box_counts );
        FREE( stats->median_box_values );
    }
}

BICAPI  void  compute_mean_and_variance(
    int   n,
    Real  samples[],
    Real  *mean,
    Real  *variance )
{
    int    i;
    Real   sum_x, sum_xx;

    sum_x = 0.0;
    sum_xx = 0.0;

    for_less( i, 0, n )
    {
        sum_x += samples[i];
        sum_xx += samples[i] * samples[i];
    }

    *mean = sum_x / (Real) n;

    if( n == 1 )
        *variance = 0.0;
    else
        *variance = (sum_xx - sum_x * sum_x / (Real) n) / (Real) (n - 1);
}

BICAPI  Real  compute_two_means_t_statistic(
    int    n1,
    Real   samples1[],
    int    n2,
    Real   samples2[] )
{
    Real   mean1, mean2, variance1, variance2, std_dev, std_err;
    Real   t;

    compute_mean_and_variance( n1, samples1, &mean1, &variance1 );
    compute_mean_and_variance( n2, samples2, &mean2, &variance2 );

    std_dev = sqrt( ((Real) n1 * variance1 + (Real) n2 * variance2) /
                    (Real) (n1 + n2 - 2) );

    std_err = std_dev * sqrt( 1.0 / (Real) n1 + 1.0 / (Real) n2 );

    if( std_err == 0.0 )
        t = 0.0;
    else
        t = (mean1 - mean2) / std_err;

    return( t );
}
