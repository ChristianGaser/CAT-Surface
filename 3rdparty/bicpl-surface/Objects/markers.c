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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/markers.c,v 1.5 2005/08/17 22:28:27 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_marker
@INPUT      : type
              colour
@OUTPUT     : marker
@RETURNS    : 
@DESCRIPTION: Initializes the marker structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep 25, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   initialize_marker(
    marker_struct     *marker,
    Marker_types      type,
    Colour            colour )
{
    marker->type = type;
    marker->colour = colour;
    marker->label = NULL;
    marker->size = 1.0;
    marker->structure_id = 0;
    marker->patient_id = 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_marker
@INPUT      : marker
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the marker object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 25, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   delete_marker(
    marker_struct     *marker )
{
    delete_string( marker->label );
}
