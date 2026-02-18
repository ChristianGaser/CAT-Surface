/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_PATCH_H_
#define _CAT_PATCH_H_

struct patchinfo {
        int num;
        int pts[3];

        struct patchinfo *next;
};

object_struct ** extract_patch_around_polygon(polygons_struct *, int, int);
object_struct ** extract_patch_around_point(polygons_struct *, int, int);
object_struct ** extract_patch_polys(polygons_struct *, int *, int);
/**
 * \brief Public API for extract_patch_points.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of extract_patch_points.
 * \param param (in/out) Parameter of extract_patch_points.
 * \param int (in/out) Parameter of extract_patch_points.
 * \return Return value of extract_patch_points.
 */
object_struct ** extract_patch_points(polygons_struct *, int *, int);

#endif
