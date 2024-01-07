/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              Peter Neelin, Louis Collins,
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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : safe_compute_xfm.c
@DESCRIPTION: Routine to calculate a General_transform from a pair of tag
              point lists. This routine is safe in that it returns
              an identity matrix if any error occurs.
@METHOD     : 
@GLOBALS    : 
@CREATED    : April 21, 1994 (Peter Neelin)
@MODIFIED   : $Log: safe_compute_xfm.c,v $
@MODIFIED   : Revision 1.13  2005/08/17 22:26:48  bert
@MODIFIED   : Replace public/private with BICAPI/static
@MODIFIED   :
@MODIFIED   : Revision 1.12  2005/06/03 18:08:37  bert
@MODIFIED   : Include config.h and use HAVE_FORK rather than NO_FORK, also omit unistd.h
@MODIFIED   :
@MODIFIED   : Revision 1.11  2000/02/06 15:30:51  stever
@MODIFIED   : rearranged header file structure; add gcc -Wall fixes
@MODIFIED   :
@MODIFIED   : Revision 1.10  2000/02/05 21:27:21  stever
@MODIFIED   : change include lines <foo.h> --> <bicpl/foo.h>
@MODIFIED   :
@MODIFIED   : Revision 1.9  1995/10/19 15:48:59  david
@MODIFIED   : check_in_all
@MODIFIED   :
 * Revision 1.8  1995/09/29  19:06:37  david
 * *** empty log message ***
 *
 * Revision 1.7  1995/07/31  13:46:03  david
 * check_in_all
 *
 * Revision 1.6  1995/07/10  18:02:53  david
 * check_in_all
 *
 * Revision 1.5  1995/07/10  14:36:11  david
 * check_in_all
 *
 * Revision 1.4  1995/04/28  18:30:19  david
 * check_in_all
 *
 * Revision 1.3  1995/03/07  18:54:51  david
 * check_in_all
 *
 * Revision 1.2  94/11/25  14:23:27  david
 * check_in_all
 * 
 * Revision 1.1  94/11/04  14:45:55  david
 * Initial revision
 * 
 * Revision 1.1  94/04/22  08:07:20  neelin
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#if HAVE_FORK
#include <sys/wait.h>
#endif

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Transforms/safe_compute_xfm.c,v 1.13 2005/08/17 22:26:48 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : safe_compute_transform_from_tags
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
              This routine is safe in that it returns an identity matrix 
              if any error occurs.
@METHOD     : 
@GLOBALS    : 
@CALLS      :
@CREATED    : April 21, 1994 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void safe_compute_transform_from_tags(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform )
{
#if !HAVE_FORK
        compute_transform_from_tags( npoints, tag_list1, tag_list2, trans_type,
                                     transform );
#else
    int                 fildes[2];
    FILE                *fpin, *fpout;
    Status              status;
    int                 statptr;
    General_transform   computed_transform;

    /* Create a pipe */

    if (pipe(fildes)) {
        create_linear_transform(transform, NULL);
        return;
    }

    /* Fork */
    if (fork()) {          /* Parent */
        (void) close(fildes[1]);
        fpin = fdopen(fildes[0], "r");
        status = input_transform(fpin, NULL, transform);
        (void) fclose(fpin);
        do {
            (void) wait(&statptr);
        } while (WIFSTOPPED(statptr));
        if (WEXITSTATUS(statptr) || status != OK) {
           if( status == OK )
               delete_general_transform( transform );
           create_linear_transform(transform, NULL);
           return;
        }
    }

    else {                 /* Child */
        (void) close(fildes[0]);
        fpout = fdopen(fildes[1], "w");
        compute_transform_from_tags(npoints, tag_list1, tag_list2, trans_type,
                                    &computed_transform);
        status = output_transform(fpout, NULL, NULL, NULL, &computed_transform);
        delete_general_transform( &computed_transform );
        (void) fclose(fpout);
        if (status != OK) {
            exit(EXIT_FAILURE);
        }
        else {
           exit(EXIT_SUCCESS);
        }
    }

    return;
#endif
}
