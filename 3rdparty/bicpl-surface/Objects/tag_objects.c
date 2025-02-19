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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/tag_objects.c,v 1.8 2005/08/17 22:28:27 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_tag_objects_file
@INPUT      : filename
              marker_colour
              default_size
              default_type
@OUTPUT     : n_objects
              object_list
@RETURNS    : OK or ERROR
@DESCRIPTION: Inputs a tag file into a list of objects (markers).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Status   input_tag_objects_file(
    STRING         filename,
    Colour         marker_colour,
    Real           default_size,
    Marker_types   default_type,
    int            *n_objects,
    object_struct  **object_list[] )
{
    Status             status;
    object_struct      *object;
    marker_struct      *marker;
    int                i, n_volumes, n_tag_points, *structure_ids, *patient_ids;
    STRING             *labels;
    Real             *weights;
    Real             **tags1, **tags2;

    *n_objects = 0;

    status = input_tag_file( filename, &n_volumes, &n_tag_points,
                             &tags1, &tags2, &weights,
                             &structure_ids, &patient_ids, &labels );

    if( status == OK )
    {
        for_less( i, 0, n_tag_points )
        {
            object = create_object( MARKER );
            marker = get_marker_ptr( object );
            fill_Point( marker->position, tags1[i][X], tags1[i][Y],tags1[i][Z]);
            marker->label = create_string( labels[i] );

            if( structure_ids[i] >= 0 )
                marker->structure_id = structure_ids[i];
            else
                marker->structure_id = -1;

            if( patient_ids[i] >= 0 )
                marker->patient_id = patient_ids[i];
            else
                marker->patient_id = -1;

            if( weights[i] > 0.0 )
                marker->size = weights[i];
            else
                marker->size = default_size;

            marker->colour = marker_colour;
            marker->type = default_type;

            add_object_to_list( n_objects, object_list, object );
        }

        free_tag_points( n_volumes, n_tag_points, tags1, tags2, weights,
                         structure_ids, patient_ids, labels );
    }

    return( status );
}
