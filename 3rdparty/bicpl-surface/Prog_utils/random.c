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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "bicpl_internal.h"

#if HAVE_SYS_TIME_H
#include  <sys/time.h>
#endif

#include <stdlib.h>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Prog_utils/random.c,v 1.15 2005/08/17 22:27:56 bert Exp $";
#endif

static  BOOLEAN  initialized = FALSE;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_random_seed
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the random seed for the generator
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :      1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_random_seed( int seed )
{
#if HAVE_SRANDOM
    (void) srandom( (unsigned int) seed );
#else
    (void) srand( (unsigned int) seed );
#endif /* HAVE_SRANDOM */

    initialized = TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_initialized
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Ensures that the random seed has been set.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :      1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  check_initialized( void )
{
#if HAVE_GETTIMEOFDAY
    struct   timeval   t;
#endif
    int                seed;

    if( !initialized )
    {
#if HAVE_GETTIMEOFDAY
        (void) gettimeofday( &t, (struct timezone *) 0 );

        seed = (int) t.tv_usec;
#else
        seed = time(NULL);
#endif

        set_random_seed( seed );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_random
@INPUT      : 
@OUTPUT     : 
@RETURNS    : int
@DESCRIPTION: Gets a random number from 0 to 2^30
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :      1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_random( void )
{
    check_initialized();

#if HAVE_RANDOM
    return( (int) random() );
#else
    return( (int) rand() );
#endif /* HAVE_RANDOM */
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_random_int
@INPUT      : n
@OUTPUT     : 
@RETURNS    : int
@DESCRIPTION: Returns a random number from 0 to n-1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :      1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_random_int( int n )
{
    return( get_random() % n );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_random_0_to_1
@INPUT      : 
@OUTPUT     : 
@RETURNS    : random number
@DESCRIPTION: Returns a random number greater than or equal to 0 and and less
            : than 1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :      1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_random_0_to_1( void )
{
    return( (Real) get_random() / (Real) RAND_MAX );
}
