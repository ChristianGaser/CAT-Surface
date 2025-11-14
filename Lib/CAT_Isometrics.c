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
#include "CAT_SafeAlloc.h"

struct metricdata *
getmetricdata(polygons_struct *polygons)
{
    struct metricdata *brain;
    struct pointdata **ptr;
    int *n_neighbours, **neighbours;
    int p, n, n2;
    Point pts[3];
    Vector dir;

    brain = SAFE_MALLOC(struct metricdata, 1);

    compute_polygon_normals(polygons);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                     &neighbours, NULL, NULL);

    ptr = SAFE_MALLOC(struct pointdata*, polygons->n_points);

    for (p = 0; p < polygons->n_points; p++) {
        if (n_neighbours[p] <= 1)
            continue; /* skip this point */

    ptr[p] = SAFE_MALLOC(struct pointdata, 1);

    ptr[p]->areas = SAFE_MALLOC(double, n_neighbours[p]);

    ptr[p]->lengths = SAFE_MALLOC(double, n_neighbours[p]);

        /* allocate this but don't calculating yet */
    ptr[p]->norm = SAFE_MALLOC(Vector, n_neighbours[p]);


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

    free(n_neighbours);
    for (p = 0; p < polygons->n_points; p++) {
        free(neighbours[p]);
    }
    free(neighbours);

    return brain;
}

double
areadistortion(struct metricdata *brain, polygons_struct *map)
{
    double ratio, area, area2, value;
    double distortion = 0;
    int poly, size;
    Point pts[3], pts2[3];

    ratio = get_polygons_surface_area(brain->polygons) /
        get_polygons_surface_area(map);

    for (poly = 0; poly < brain->polygons->n_items; poly++) {
        size = get_polygon_points(brain->polygons, poly, pts);
        area = get_polygon_surface_area(size, pts);

        size = get_polygon_points(map, poly, pts2);
        area2 = get_polygon_surface_area(size, pts2);

        size = GET_OBJECT_SIZE(*brain->polygons, poly);

        value = log10(ratio * area2 / area);

        if (value < PINF && value > NINF)
            distortion += fabs(value);
    }

    return (distortion / brain->polygons->n_items);
}

int
smooth(struct metricdata *brain, polygons_struct *map, int maxiters,
     int selectflag, double tolerance)
{
    int          i, n, it, p;
    Point        pts[3], *oldpts, *newpts, newcenter;
    double         a, areas[128], centers[128], totalArea;
    double         ad, new_ad, weight, radius, xyz[3];
    double         metric, newmetric;
    int          n1, n2, count, stepsize;

    if (brain->polygons->n_points != map->n_points
      || brain->polygons->n_items != map->n_items) {
        printf("Input mesh and spherical map mesh do not match.\n");
        return(EXIT_FAILURE);
    }
  
    /* map the spherical map to a sphere of the same surface area */
    radius = sqrt(get_polygons_surface_area(brain->polygons) / (4.0 * PI));
    for (p = 0; p < map->n_points; p++)
        set_vector_length(&map->points[p], radius);

    newpts = SAFE_MALLOC(Point, map->n_points);
    for (p = 0; p < map->n_points; p++)
        newpts[p] = map->points[p];

    metric = areadistortion(brain, map);
    printf("  0: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", metric);

    stepsize = 100;
    for (it = 1; it <= maxiters; it++) {
        count = 0;
        for (p = 0; p < map->n_points; p++) {
            if (brain->n_neigh[p] <= 1)
                continue; /* skip this point */

            ad = 0.0;
            totalArea = 0.0;

            /* Get 2 consecutive neighbors of this node */
            pts[0] = newpts[p];
            for (n = 0; n < brain->n_neigh[p]; n++) {
                n2 = (n + 1) % brain->n_neigh[p];

                /* area of the triangle */
                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
                areas[n] = get_polygon_surface_area(3, pts);

                totalArea += areas[n];

                a = brain->ptdata[p]->areas[n];
                if (a > 0) {
                    ad += log10(areas[n] / a);
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
            if (totalArea == 0) continue;

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
                fill_Point(newpts[p], xyz[0], xyz[1], xyz[2]);
                continue;
            }

            /* see if the new center reduces area distortion! */
            fill_Point(newcenter, xyz[0], xyz[1], xyz[2]);
            new_ad = 0;

            pts[0] = newcenter;
            for (n = 0; n < brain->n_neigh[p]; n++) {
                n2 = (n + 1) % brain->n_neigh[p];

                /* area of the triangle */
                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
                areas[n] = get_polygon_surface_area(3, pts);

                a = brain->ptdata[p]->areas[n];
                if (a > 0)
                    new_ad += log10(areas[n] / a);
            }
            new_ad = fabs(new_ad / brain->n_neigh[p]);

            if (new_ad < ad) {
                count++;
                fill_Point(newpts[p], xyz[0], xyz[1], xyz[2]);
            }
        }

        for (p = 0; p < map->n_points; p++)
            set_vector_length(&newpts[p], radius);

        if (it % stepsize == 0) {
            oldpts = map->points;
            map->points = newpts;
            newmetric = areadistortion(brain, map);

            if (selectflag == SELECT_ON) {
                if (newmetric >= metric) {
                    /* don't update, decrement step */
                    map->points = oldpts;
                    for (p = 0; p < map->n_points; p++)
                        newpts[p] = map->points[p];
                    it -= stepsize;
                    stepsize = round(stepsize / 2);
                    if (stepsize == 0) break;
                } else if (metric - newmetric < tolerance) {
                    metric = newmetric;
                    newpts = oldpts;
                    printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                    break;
                } else {
                    metric = newmetric;
                    newpts = oldpts;
                    printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                    if (count == 0) break;
                    for (p = 0; p < map->n_points; p++)
                        newpts[p] = map->points[p];
                }
            } else {
                metric = newmetric;
                newpts = oldpts;
                printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                for (p = 0; p < map->n_points; p++)
                    newpts[p] = map->points[p];
            }
        }
    }
    if (it > maxiters) it = maxiters;
    if (it < 0) it = 0;

    printf("Smooth %d iters, Area distortion %f\n", it, areadistortion(brain, map));

    free(newpts);
    return it;
}

int
distortcorrect(struct metricdata *brain, polygons_struct *map, int maxiters,
         int selectflag, double tolerance)
{
    int          i, n, it, p;
    Point        pts[3], *oldpts, *newpts, newcenter;
    double         a, areas[128], centers[128];
    double         ad, new_ad, total_ad, weight, radius, xyz[3];
    double         metric, newmetric;
    int          n1, n2, count, stepsize;
    progress_struct    progress;

    if (brain->polygons->n_points != map->n_points
      || brain->polygons->n_items != map->n_items) {
        printf("Input mesh and spherical map mesh do not match.\n");
        return(-1);
    }
  
    /* map the spherical map to a sphere of the same surface area */
    radius = sqrt(get_polygons_surface_area(brain->polygons) / (4.0 * PI));
    for (p = 0; p < map->n_points; p++)
        set_vector_length(&map->points[p], radius);

    newpts = SAFE_MALLOC(Point, map->n_points);
    for (p = 0; p < map->n_points; p++)
        newpts[p] = map->points[p];

    metric = areadistortion(brain, map);
    printf("  0: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", metric);

    stepsize = 100;
    for (it = 1; it <= maxiters; it++) {
        count = 0;
        for (p = 0; p < map->n_points; p++) {
            if (brain->n_neigh[p] <= 1)
                continue; /* skip this point */

            ad = 0.0;
            total_ad = 0.0;

            /* Get 2 consecutive neighbors of this node */
            pts[0] = newpts[p];
            for (n = 0; n < brain->n_neigh[p]; n++) {
                n2 = (n + 1) % brain->n_neigh[p];

                /* area of the triangle */
                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
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
            if (total_ad == 0) continue;

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
                fill_Point(newpts[p], xyz[0], xyz[1], xyz[2]);
                continue;
            }

            /* see if the new center reduces area distortion! */
            fill_Point(newcenter, xyz[0], xyz[1], xyz[2]);
            new_ad = 0;

            pts[0] = newcenter;
            for (n = 0; n < brain->n_neigh[p]; n++) {
                n2 = (n + 1) % brain->n_neigh[p];

                /* area of the triangle */
                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
                areas[n] = get_polygon_surface_area(3, pts);

                a = brain->ptdata[p]->areas[n];
                if (a > 0)
                    new_ad += log10(areas[n] / a);
            }
            new_ad = fabs(new_ad / brain->n_neigh[p]);

            if (new_ad < ad) {
                count++;
                fill_Point(newpts[p], xyz[0], xyz[1], xyz[2]);
            }
        }

        for (p = 0; p < map->n_points; p++)
            set_vector_length(&newpts[p], radius);

        if (it % stepsize == 0) {
            oldpts = map->points;
            map->points = newpts;
            newmetric = areadistortion(brain, map);

            if (newmetric >= metric) {
                /* don't update, decrement step */
                map->points = oldpts;
                for (p = 0; p < map->n_points; p++)
                    newpts[p] = map->points[p];
                it -= stepsize;
                stepsize = round(stepsize / 2);
                if (stepsize == 0) break;
            } else if (metric - newmetric < tolerance) {
                metric = newmetric;
                newpts = oldpts;
                printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                break;
            } else {
                metric = newmetric;
                newpts = oldpts;
                printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                if (count == 0) break;
                for (p = 0; p < map->n_points; p++)
                    newpts[p] = map->points[p];
            }
        }
    }
    if (it > maxiters) it = maxiters;
    if (it < 0) it = 0;

    printf("Distortion correction %d iters, Area distortion %f\n", it, metric);

    free(newpts);
    return it;
}

int
stretch(struct metricdata *brain, polygons_struct *map, int maxiters,
    int selectflag, int largeonly, double tolerance)
{
    int          i, n, it, p;
    Point        pts[3], *oldpts, *newpts, npt, dir, newcenter;
    Vector         norm;
    double         dist, radius;
    double         metric, newmetric;
    double         ad, new_ad, a, areas[128];
    int          n2, count, ok, stepsize;

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

    newpts = SAFE_MALLOC(Point, map->n_points);
    for (p = 0; p < map->n_points; p++)
        newpts[p] = map->points[p];

    metric = areadistortion(brain, map);
    printf("  0: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", metric);

    stepsize = 10;
    for (it = 1; it <= maxiters; it++) {
        count = 0;

        for (p = 0; p < map->n_points; p++) {
            if (brain->n_neigh[p] <= 1)
                continue; /* skip this point */

            ad = 0.0;

            pts[0] = newpts[p];
            for (n = 0; n < brain->n_neigh[p]; n++) {
                n2 = (n + 1) % brain->n_neigh[p];

                /* area of the triangle */
                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
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
                npt = newpts[ brain->neigh[p][n] ];

                SUB_POINTS(dir, npt, newpts[p]);

                dist = MAGNITUDE(dir) -
                     brain->ptdata[p]->lengths[n];
                dist /= MAGNITUDE(dir);
                SCALE_VECTOR(dir, dir, dist);
                ADD_POINT_VECTOR(newcenter, newcenter, dir);
            }

            Point_x(newcenter) /= brain->n_neigh[p];
            Point_y(newcenter) /= brain->n_neigh[p];
            Point_z(newcenter) /= brain->n_neigh[p];
            ADD_POINT_VECTOR(newcenter, newpts[p], newcenter);
            set_vector_length(&newcenter, radius);

            if (isnan(Point_x(newcenter)) ||
              isnan(Point_y(newcenter)) ||
              isnan(Point_z(newcenter)))
                continue;

            /* check for flips */
            ok = 1;
            pts[0] = newcenter;
            for (n = 0; n < brain->n_neigh[p]; n++) {
                n2 = (n + 1) % brain->n_neigh[p];

                /* normal of the triangle */
                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
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
                    pts[1] = newpts[brain->neigh[p][n]];
                    pts[2] = newpts[brain->neigh[p][n2]];
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
            fill_Point(newpts[p], Point_x(newcenter),
                   Point_y(newcenter), Point_z(newcenter));
            pts[0] = newpts[p];
            for (n = 0; n < brain->n_neigh[p]; n++) {
                 n2 = (n + 1) % brain->n_neigh[p];

                pts[1] = newpts[brain->neigh[p][n]];
                pts[2] = newpts[brain->neigh[p][n2]];
                find_polygon_normal(3, pts,
                          &brain->ptdata[p]->norm[n]);
            }
        }

        if (it % stepsize == 0) {
            oldpts = map->points;
            map->points = newpts;
            newmetric = areadistortion(brain, map);

            if (newmetric >= metric) {
                /* don't update, decrement step */
                map->points = oldpts;
                for (p = 0; p < map->n_points; p++)
                    newpts[p] = map->points[p];
                it -= stepsize;
                stepsize = round(stepsize / 2);
                if (stepsize == 0) break;
            } else if (metric - newmetric < tolerance) {
                metric = newmetric;
                newpts = oldpts;
                printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                break;
            } else {
                metric = newmetric;
                newpts = oldpts;
                printf("%d: %f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it, metric);
                if (count == 0) break;
                for (p = 0; p < map->n_points; p++)
                    newpts[p] = map->points[p];
            }
        }
    }
    if (it > maxiters) it = maxiters;
    if (it < 0) it = 0;

    radius = sqrt(get_polygons_surface_area(brain->polygons) / (4.0 * PI));
    for (p = 0; p < map->n_points; p++)
        set_vector_length(&map->points[p], radius);

    printf("Stretch %d iters, Are distortion %f\n", it, areadistortion(brain, map));

    free(newpts);
    return it;
}
