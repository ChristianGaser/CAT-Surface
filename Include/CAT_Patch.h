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
object_struct ** extract_patch_points(polygons_struct *, int *, int);

#endif