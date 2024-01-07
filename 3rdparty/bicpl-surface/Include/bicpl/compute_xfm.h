#ifndef  _COMPUTE_XFM_H
#define  _COMPUTE_XFM_H

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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_xfm.h
@DESCRIPTION: Header file for compute_transform_from_tags
@GLOBALS    :
@CALLS      :
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : $Log: compute_xfm.h,v $
@MODIFIED   : Revision 1.2  2002/05/15 13:17:55  stever
@MODIFIED   : Remove rcsid strings.
@MODIFIED   :
@MODIFIED   : Revision 1.1  2000/02/06 15:30:35  stever
@MODIFIED   : rearranged header file structure; add gcc -Wall fixes
@MODIFIED   :
@MODIFIED   : Revision 1.3  1995/07/31 13:44:40  david
@MODIFIED   : check_in_all
@MODIFIED   :
 * Revision 1.2  1995/07/12  18:07:27  david
 * *** empty log message ***
 *
 * Revision 1.1  94/11/04  14:09:05  david
 * Initial revision
 * 
 * Revision 1.1  94/11/04  13:11:36  david
 * Initial revision
 * 
 * Revision 1.2  93/09/08  15:47:44  neelin
 * Added 9 and 10 parameter linear fitting.
 * 
 * Revision 1.1  93/09/01  15:25:37  neelin
 * Initial revision
 * 
 * Revision 1.1  93/09/01  13:22:24  neelin
 * Initial revision
 * 
---------------------------------------------------------------------------- */


/* Transformation types */
typedef enum {
   TRANS_LSQ6, TRANS_LSQ7, TRANS_LSQ9, TRANS_LSQ10, TRANS_LSQ12, TRANS_TPS
} Trans_type;

/* Minimum number of points needed */
#define MIN_POINTS_LINEAR 4
#define MIN_POINTS_TPS    5

#endif
