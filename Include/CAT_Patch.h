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

#include <bicpl.h>

struct patchinfo {
        int num;
        int pts[3];

        struct patchinfo *next;
};

/**
 * \brief Extract a polygon patch around a seed triangle.
 *
 * \param polygons (in)  source mesh
 * \param poly     (in)  seed polygon index
 * \param level    (in)  neighborhood depth
 * \return Allocated object list containing the patch mesh
 */
object_struct ** extract_patch_around_polygon(polygons_struct *polygons,
                                              int poly, int level);
/**
 * \brief Extract a polygon patch around a seed vertex.
 *
 * \param polygons (in)  source mesh
 * \param point    (in)  seed vertex index
 * \param level    (in)  neighborhood depth
 * \return Allocated object list containing the patch mesh
 */
object_struct ** extract_patch_around_point(polygons_struct *polygons,
                                            int point, int level);
/**
 * \brief Extract a patch from a polygon selection mask.
 *
 * \param polygons (in)  source mesh
 * \param polys    (in)  polygon selection mask
 * \param num      (in)  selection value (0 to include all non-zero)
 * \return Allocated object list containing the patch mesh
 */
object_struct ** extract_patch_polys(polygons_struct *polygons, int *polys,
                                     int num);
/**
 * \brief Extract a patch from a point selection mask.
 *
 * \param polygons (in)  source mesh
 * \param points   (in)  point selection mask
 * \param num      (in)  selection value (0 to include all non-zero)
 * \return Allocated object list containing the patch mesh
 */
object_struct ** extract_patch_points(polygons_struct *polygons, int *points,
                                      int num);

#endif
