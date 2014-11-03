/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Gyrification.c 312 2014-05-07 14:28:54Z gaser $
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Smooth.h"
#include "CAT_Gyrification.h"
#include "CAT_Resample.h"
#include "CAT_SPH.h"
#include "CAT_ConvexHull.h"


double
gyrification_index_sph(polygons_struct *surface, polygons_struct *sphere,
                      char *file, int n_triangles, polygons_struct *reparam)
{
        polygons_struct *polygons, *convex, *convex_sphere, *convex_resampled, *surface_resampled;
        object_struct **object;
        int i, n_points;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz, *spectral_power;
        double *convex_areas, *areas, convex_area, area, *local_gi;

        n_points = 81920;
        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        
        /* get convex hull */
        object = surface_get_convex_hull(surface, sphere);
        convex = get_polygons_ptr(*object);

output_graphics_any_format( "convex0.obj", ASCII_FORMAT, 1, object, NULL);
        /* get sphere of convex hull */
        convex_sphere = get_polygons_ptr(create_object(POLYGONS));
        copy_polygons(convex, convex_sphere);
        surf_to_sphere(convex_sphere, 6);

        convex_areas = (double *) malloc(sizeof(double) * surface->n_points);
        convex_area = get_area_of_points_normalized_to_sphere(convex, convex_sphere, convex_areas);

        areas = (double *) malloc(sizeof(double) * surface->n_points);
        area = get_area_of_points_normalized_to_sphere(surface, sphere, areas);

        /* estimate ratio */
        for (i = 0; i < surface->n_points; i++)
                areas[i] /= convex_areas[i];

        output_values_any_format("gi.txt", surface->n_points, areas,
                                 TYPE_DOUBLE);

        free(convex_areas);
        free(areas);
        free(object);

        return area/convex_area;
}


