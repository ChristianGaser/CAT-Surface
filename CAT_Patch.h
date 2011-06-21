/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

struct patchinfo {
        int num;
        int pts[3];

        struct patchinfo *next;
};


object_struct ** extract_patch_around_polygon(polygons_struct *, int, int);
object_struct ** extract_patch_around_point(polygons_struct *, int, int);
object_struct ** extract_patch_polys(polygons_struct *, int *, int);
object_struct ** extract_patch_points(polygons_struct *, int *, int);
