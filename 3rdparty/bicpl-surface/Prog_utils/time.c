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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Prog_utils/time.c,v 1.8 2005/08/17 22:27:56 bert Exp $";
#endif

static  Real  cpu_start_time;
static  Real  real_start_time;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : start_timing
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Starts timing a subset of a process.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  start_timing( void )
{
    cpu_start_time = current_cpu_seconds();
    real_start_time = current_realtime_seconds();
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : end_timing
@INPUT      : title
              n     - number of iterations of job that were performed.
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Computes the cpu and real time since the start time,
              prints this out divided by the number of iterations performed.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  end_timing(
    STRING   title,
    int      n )
{
    Real  per_unit_real_time, per_unit_cpu_time;
    Real  cpu_time, real_time;

    if( n < 1 )
        n = 1;

    cpu_time = current_cpu_seconds() - cpu_start_time;
    real_time = current_realtime_seconds() - real_start_time;

    per_unit_cpu_time = cpu_time / (Real) n;
    per_unit_real_time = real_time / (Real) n;

    if( title != NULL )
        print( "%s: ", title );
    print( "per cycle: " );
    print_time( "Cpu:  %g %s\t", per_unit_cpu_time );
    print_time( "Real:  %g %s", per_unit_real_time );

    if( cpu_time > 0.0 && real_time > 0.0 )
        print( "  c/s: %g  %g\n", (Real) n / cpu_time, (Real) n / real_time );
}
