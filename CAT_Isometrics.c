/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* algorithms to make a map more isometric: smoothing, distortion correction,
 * and stretch optimization */

#include "CAT_Isometrics.h"

struct metricdata *
getmetricdata(polygons_struct *polygons)
{
        struct metricdata *brain;
        struct pointdata **ptr;
        int *n_neighbours, **neighbours;
        int p, n, n2;
        Point pts[3];
        Vector dir;

        brain = (struct metricdata *) malloc(sizeof(struct metricdata));

        compute_polygon_normals(polygons);
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                         &neighbours, NULL, NULL);

        ptr = (struct pointdata **)
              malloc(sizeof(struct pointdata *) * polygons->n_points);

        for (p = 0; p < polygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                ptr[p] = (struct pointdata *) malloc(sizeof(struct pointdata));

                ptr[p]->areas = (double *) malloc(sizeof(double) *
                                                        n_neighbours[p]);

                ptr[p]->lengths = (double *) malloc(sizeof(double) *
                                                          n_neighbours[p]);

                /* allocate this but don't calculating yet */
                ptr[p]->norm = (Vector *) malloc(sizeof(Vector) *
                                                       n_neighbours[p]);


                /* Get 2 consecutive neighbors of this node */
                pts[0] = polygons->points[p];
                for (n = 0; n < n_neighbours[p]; n++) {
                        n2 = (n + 1) % n_neighbours[p];

                        /* area of the triangle */
                        pts[1] = polygons->points[ neighbours[p][n] ];
                        pts[2] = polygons->points[ neighbours[p][n2] ];
                        ptr[p]->areas[n] = get_polygon_surface_area(3, pts);

                        /* lengths to neighbouring nodes */
                        SUB_POINTS(dir, pts[1], polygons->points[p]);
                        ptr[p]->lengths[n] = MAGNITUDE(dir);
                }
        }
        brain->polygons = polygons;
        brain->n_neigh = n_neighbours;
        brain->neigh = neighbours;
        brain->ptdata = ptr;

        return brain;
}

int
smooth(struct metricdata *brain, polygons_struct *map, int iters,
       int selectflag, int quietflag)
{
        int                i, n, it, p;
        Point              pts[3], newcenter;
        double             a, areas[128], centers[128], totalArea;
        double             ad, new_ad, weight, radius, xyz[3];
        int                n1, n2, count;
        progress_struct    progress;

        if (brain->polygons->n_points != map->n_points
            || brain->polygons->n_items != map->n_items) {
                printf("Input mesh and conformal map mesh do not match.\n");
                return(-1);
        }
    
        /* map the spherical map to a sphere of the same surface area */
        radius = sqrt(get_polygons_surface_area(brain->polygons) / (4.0 * PI));
        for (p = 0; p < map->n_points; p++)
                set_vector_length(&map->points[p], radius);

        if (quietflag == 0)
                initialize_progress_report(&progress, FALSE, iters, "Smooth");
        
        for (it = 1; it <= iters; it++) {
                count = 0;
                for (p = 0; p < map->n_points; p++) {
                        if (brain->n_neigh[p] <= 1)
                                continue; /* skip this point */

                        ad = 0.0;
                        totalArea = 0.0;

                        /* lengths to neighbouring nodes */
                        pts[0] = map->points[p];
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                n2 = (n + 1) % brain->n_neigh[p];

                                /* area of the triangle */
                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                areas[n] = get_polygon_surface_area(3, pts);

                                totalArea += areas[n];

                                a = brain->ptdata[p]->areas[n];
                                if (a > 0)
                                        ad += log10(areas[n] / a);

                                /* Save center of this tile */
                                centers[n*3    ] = (Point_x(pts[0]) +
                                                    Point_x(pts[1]) +
                                                    Point_x(pts[2])) / 3.0;
                                centers[n*3 + 1] = (Point_y(pts[0]) +
                                                    Point_y(pts[1]) +
                                                    Point_y(pts[2])) / 3.0;
                                centers[n*3 + 2] = (Point_z(pts[0]) +
                                                    Point_z(pts[1]) +
                                                    Point_z(pts[2])) / 3.0;
                        }
                        ad = fabs(ad / brain->n_neigh[p]);

                        /* Algorithm #1 - Area Smoothing */
                        for (i = 0; i < 3; i++)
                                xyz[i] = 0.0;
                        for (n = 0; n <  brain->n_neigh[p]; n++) {
                                if (totalArea > 0) {
                                        weight = areas[n] / totalArea;
                                        for (i = 0; i < 3; i++)
                                                xyz[i] += weight *
                                                          centers[n*3 + i];
                                }
                        }

                        if (selectflag == SELECT_OFF) {
                                count++;
                                fill_Point(map->points[p],
                                           xyz[0], xyz[1], xyz[2]);
                                continue;
                        }

                        /* see if the new center reduces area distortion! */
                        fill_Point(newcenter, xyz[0], xyz[1], xyz[2]);
                        new_ad = 0;

                        pts[0] = newcenter;
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                n2 = (n + 1) % brain->n_neigh[p];

                                /* area of the triangle */
                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                areas[n] = get_polygon_surface_area(3, pts);

                                a = brain->ptdata[p]->areas[n];
                                if (a > 0)
                                        new_ad += log10(areas[n] / a);
                        }
                        new_ad = fabs(new_ad / brain->n_neigh[p]);

                        if (new_ad < ad) {
                                count++;
                                fill_Point(map->points[p],
                                           xyz[0], xyz[1], xyz[2]);
                        }

                }

                for (p = 0; p < map->n_points; p++)
                        set_vector_length(&map->points[p], radius);

                if (quietflag == 0)
                        update_progress_report(&progress, it);

                if (count == 0)
                        break; /* done! */
        }
        if (quietflag == 0)
                terminate_progress_report(&progress);

        return it;
}

int
distortcorrect(struct metricdata *brain, polygons_struct *map, int iters,
       int selectflag, int quietflag)
{
        int                i, n, it, p;
        Point              pts[3], newcenter;
        double             a, areas[128], centers[128];
        double             ad, new_ad, total_ad, weight, radius, xyz[3];
        int                n1, n2, count;
        progress_struct    progress;

        if (brain->polygons->n_points != map->n_points
            || brain->polygons->n_items != map->n_items) {
                printf("Input mesh and conformal map mesh do not match.\n");
                return(-1);
        }
    
        /* map the spherical map to a sphere of the same surface area */
        radius = sqrt(get_polygons_surface_area(brain->polygons) / (4.0 * PI));
        for (p = 0; p < map->n_points; p++)
                set_vector_length(&map->points[p], radius);

        if (quietflag == 0)
                initialize_progress_report(&progress, FALSE, iters,
                                           "DistortCorrect");
        
        for (it = 1; it <= iters; it++) {
                count = 0;
                for (p = 0; p < map->n_points; p++) {
                        if (brain->n_neigh[p] <= 1)
                                continue; /* skip this point */

                        ad = 0.0;
                        total_ad = 0.0;

                        /* Get 2 consecutive neighbors of this node */
                        pts[0] = map->points[p];
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                n2 = (n + 1) % brain->n_neigh[p];

                                /* area of the triangle */
                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                areas[n] = get_polygon_surface_area(3, pts);

                                a = brain->ptdata[p]->areas[n];
                                if (a > 0) {
                                        ad += log10(areas[n] / a);
                                        total_ad += areas[n] / a;
                                }

                                /* Save center of this tile */
                                centers[n*3    ] = (Point_x(pts[0]) +
                                                    Point_x(pts[1]) +
                                                    Point_x(pts[2])) / 3.0;
                                centers[n*3 + 1] = (Point_y(pts[0]) +
                                                    Point_y(pts[1]) +
                                                    Point_y(pts[2])) / 3.0;
                                centers[n*3 + 2] = (Point_z(pts[0]) +
                                                    Point_z(pts[1]) +
                                                    Point_z(pts[2])) / 3.0;
                        }
                        ad = fabs(ad / brain->n_neigh[p]);

                        /* Algorithm #2 - Area Distortion */
                        for (i = 0; i < 3; i++)
                                xyz[i] = 0.0;
                        for (n = 0; n <  brain->n_neigh[p]; n++) {
                                a = brain->ptdata[p]->areas[n];
                                if (a > 0) {
                                        weight = (areas[n]/a) / total_ad;
                                        for (i = 0; i < 3; i++)
                                                xyz[i] += weight *
                                                          centers[n*3 + i];
                                }
                        }

                        if (selectflag == SELECT_OFF) {
                                count++;
                                fill_Point(map->points[p],
                                           xyz[0], xyz[1], xyz[2]);
                                continue;
                        }

                        /* see if the new center reduces area distortion! */
                        fill_Point(newcenter, xyz[0], xyz[1], xyz[2]);
                        new_ad = 0;

                        pts[0] = newcenter;
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                n2 = (n + 1) % brain->n_neigh[p];

                                /* area of the triangle */
                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                areas[n] = get_polygon_surface_area(3, pts);

                                a = brain->ptdata[p]->areas[n];
                                if (a > 0)
                                        new_ad += log10(areas[n] / a);
                        }
                        new_ad = fabs(new_ad / brain->n_neigh[p]);

                        if (new_ad < ad) {
                                count++;
                                fill_Point(map->points[p],
                                           xyz[0], xyz[1], xyz[2]);
                        }

                }

                for (p = 0; p < map->n_points; p++)
                        set_vector_length(&map->points[p], radius);

                if (quietflag == 0)
                        update_progress_report(&progress, it);

                if (count == 0)
                        break; /* done! */
        }
        if (quietflag == 0)
                terminate_progress_report(&progress);

        return it;
}

int
stretch(struct metricdata *brain, polygons_struct *map, int iters,
        int selectflag, int quietflag, int largeonly)
{
        int                i, n, it, p;
        Point              pts[3], npt, dir, newcenter;
        Vector             norm;
        double             dist, radius;
        double             ad, new_ad, a, areas[128];
        int                n2, count, ok;
        progress_struct    progress;

        if (brain->polygons->n_points != map->n_points
            || brain->polygons->n_items != map->n_items) {
                printf("Input mesh and conformal map mesh do not match.\n");
                return(-1);
        }

        /* map the spherical map to a sphere of twice the surface area */
        radius = sqrt(get_polygons_surface_area(brain->polygons) / (2.0 * PI));
        for (p = 0; p < map->n_points; p++)
                set_vector_length(&map->points[p], radius);

        /* get the original normals */
        for (p = 0; p < map->n_points; p++) {
                if (brain->n_neigh[p] <= 1)
                        continue; /* skip this point */

                /* Get 2 consecutive neighbors of this node */
                pts[0] = map->points[p];
                for (n = 0; n < brain->n_neigh[p]; n++) {
                        n2 = (n + 1) % brain->n_neigh[p];

                        /* normal of the triangle */
                        pts[1] = map->points[brain->neigh[p][n]];
                        pts[2] = map->points[brain->neigh[p][n2]];
                        find_polygon_normal(3, pts,
                                            &brain->ptdata[p]->norm[n]);
                }
        }

        if (quietflag == 0)
                initialize_progress_report(&progress, FALSE, iters, "Stretch");

        for (it = 1; it <= iters; it++) {
                count = 0;

                for (p = 0; p < map->n_points; p++) {
                        if (brain->n_neigh[p] <= 1)
                                continue; /* skip this point */

                        ad = 0.0;

                        pts[0] = map->points[p];
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                n2 = (n + 1) % brain->n_neigh[p];

                                /* area of the triangle */
                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                areas[n] = get_polygon_surface_area(3, pts);

                                a = brain->ptdata[p]->areas[n] * 2;
                                if (a > 0)
                                        ad += log10(areas[n] / a);
                        }

                        if (largeonly == 1 && ad + 1 < 0)
                                continue; /* skip this point */

                        ad = fabs(ad / brain->n_neigh[p]);

                        /* Algorithm - Stretch Optimization */
                        fill_Point(newcenter, 0.0, 0.0, 0.0);

                        for (n = 0; n <  brain->n_neigh[p]; n++) {
                                npt = map->points[ brain->neigh[p][n] ];

                                SUB_POINTS(dir, npt, map->points[p]);

                                dist = MAGNITUDE(dir) -
                                       brain->ptdata[p]->lengths[n];
                                dist /= MAGNITUDE(dir);
                                SCALE_VECTOR(dir, dir, dist);
                                ADD_POINT_VECTOR(newcenter, newcenter, dir);
                        }

                        Point_x(newcenter) /= brain->n_neigh[p];
                        Point_y(newcenter) /= brain->n_neigh[p];
                        Point_z(newcenter) /= brain->n_neigh[p];
                        ADD_POINT_VECTOR(newcenter, map->points[p], newcenter);
                        set_vector_length(&newcenter, radius);

                        /* check for flips */
                        ok = 1;
                        pts[0] = newcenter;
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                n2 = (n + 1) % brain->n_neigh[p];

                                /* normal of the triangle */
                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                find_polygon_normal(3, pts, &norm);
                                if (acos(DOT_VECTORS(norm,
                                           brain->ptdata[p]->norm[n])) > 0.01) {
                                        ok = 0;
                                        break;
                                }
                        }
                        if (ok == 0)
                                continue; /* flipped, so skip it */

                        if (selectflag == SELECT_ON) {
                                /* see if new center reduces area distortion */
                                new_ad = 0; 
                                pts[0] = newcenter;
                                for (n = 0; n < brain->n_neigh[p]; n++) {
                                         n2 = (n + 1) % brain->n_neigh[p];

                                        /* area of the triangle */
                                        pts[1]=map->points[brain->neigh[p][n]];
                                        pts[2]=map->points[brain->neigh[p][n2]];
                                        areas[n] =
                                               get_polygon_surface_area(3, pts);

                                        a = brain->ptdata[p]->areas[n] * 2;
                                        if (a > 0)
                                                new_ad += log10(areas[n]/a);
                                }
                                new_ad = fabs(new_ad / brain->n_neigh[p]);
                                if (new_ad >= ad) {
                                        continue; /* skip this point */
                                }
                        }

                        count++;
                        fill_Point(map->points[p], Point_x(newcenter),
                                   Point_y(newcenter), Point_z(newcenter));
                        pts[0] = map->points[p];
                        for (n = 0; n < brain->n_neigh[p]; n++) {
                                 n2 = (n + 1) % brain->n_neigh[p];

                                pts[1] = map->points[brain->neigh[p][n]];
                                pts[2] = map->points[brain->neigh[p][n2]];
                                find_polygon_normal(3, pts,
                                                    &brain->ptdata[p]->norm[n]);
                        }
                }
                if (quietflag == 0)
                        update_progress_report(&progress, it);

                if (count == 0)
                        break; /* done! */
        }

        if (quietflag == 0)
                terminate_progress_report(&progress);

        radius = sqrt(get_polygons_surface_area(brain->polygons) / (4.0 * PI));
        for (p = 0; p < map->n_points; p++)
                set_vector_length(&map->points[p], radius);

        return it;
}
