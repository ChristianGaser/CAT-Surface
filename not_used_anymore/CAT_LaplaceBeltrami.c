/*****************************************************************************/
/*
 * Laplace-Beltrami Conformal Mapping implemented with OpenNL linear solver
 *
 *  OpenNL: Numerical Library
 *  Copyright (C) 2004 Bruno Levy
 *
 *  Laplace-Beltrami: from itkConformalFlatteningFilter.txx
 *  Copyright (C) 2006 John Melonakos
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ISA Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */
/****************************************************************************/

#include <bicpl.h>
#include <float.h>
#include "CAT_Surf.h"

#include "nl/nl.h"

struct ptmap {
        int *cells;
        int size;
};

int xpt, ypt, zpt;

/**
 * project_to_sphere(): calculate the projection from the complex plane
 * to the sphere.  npts = number of points, f = scaling factor,
 * zRf = translation of real value, zIf = translation of imaginary value.
 */
void
project_to_sphere(Point *points, int npts, double f,
                double *zR, double zRf, double *zI, double zIf) {
        double zRit, zIit, r2;
        int p;

        for (p = 0; p < npts; p++) {
                zRit = f * (zR[p] + zRf);
                zIit = f * (zI[p] + zIf);

                r2 = zRit*zRit + zIit*zIit;
                Point_x(points[p]) = 2*zRit / (1 + r2);
                Point_y(points[p]) = 2*zIit / (1 + r2);
                Point_z(points[p]) = (2*r2 / (1 + r2)) - 1;
        }
}

/* calc_polygon_areas(): calculate the polygon areas, store in areas */
void
calc_polygon_areas(polygons_struct *polygons, Point *pts, double *areas)
{
        Point tpts[3];
        int t;

        for (t = 0; t < polygons->n_items; t++) {
                tpts[0] = pts[polygons->indices[
                               POINT_INDEX(polygons->end_indices,t,0)]];
                tpts[1] = pts[polygons->indices[
                               POINT_INDEX(polygons->end_indices,t,1)]];
                tpts[2] = pts[polygons->indices[
                               POINT_INDEX(polygons->end_indices,t,2)]];
                areas[t] = get_polygon_surface_area(3, tpts);
        }
}

/* calc_area_distortion(): calculate the area distortion.  calculates
 * log difference of distortion between area1 and area2.
 */
double
calc_area_distortion(polygons_struct *polygons, double *area1, double *area2)
{
        double area_distortion = 0.0;
        int i;

        for (i = 0; i < polygons->n_items; i++) {
                if (area2[i] > 0 && area1[i] > 0)
                        area_distortion += fabs(log10(area1[i] / area2[i]));
        }
        return area_distortion;
}

/**
 * moebius_transform(): optimizes the solution for area distortion using
 * a simplified moebius transformation: f(x) = ax + b
 */
void
moebius_transform(polygons_struct *polygons, double *zR, double *zI, double *factor,
                  double *zRf, double *zIf)
{
        Point *pts, tpts[MAX_POINTS_PER_POLYGON];
        double *mareas, *areas;
        int size, f, i;
        double areadist = FLT_MAX, ad_lo, ad_hi, new_ad;
        double fstep;

        *factor = 0.0; *zRf = 0.0; *zIf = 0.0;

        mareas = (double *) malloc(sizeof(double) * polygons->n_items);
        areas = (double *) malloc(sizeof(double) * polygons->n_items);
        pts = (Point *) malloc(sizeof(Point) * polygons->n_points);

        /* first get all of the areas of the original triangles */
        for (i = 0; i < polygons->n_items; i++) {
                size = get_polygon_points(polygons, i, tpts);
                mareas[i] = get_polygon_surface_area(size, tpts);
        }

        /* get an estimate of the optimal factor value */
        for (f = 400; f <= 10000; ) { // get in the ballpark first
                project_to_sphere(pts, polygons->n_points, f, zR, 0.0, zI, 0.0);
                calc_polygon_areas(polygons, pts, areas);
                new_ad = calc_area_distortion(polygons, areas, mareas);
                if (new_ad < areadist) {
                        *factor = f;
                        areadist = new_ad;
                }
                f += 200;
        }

        fstep = 100;
        for (i = 0; i < 5; i++) { /* optimize factor 5x more */
                project_to_sphere(pts, polygons->n_points, *factor - fstep,
                                  zR, 0.0, zI, 0.0);
                calc_polygon_areas(polygons, pts, areas);
                ad_lo = calc_area_distortion(polygons, areas, mareas);

                project_to_sphere(pts, polygons->n_points, *factor + fstep,
                                  zR, 0.0, zI, 0.0);
                calc_polygon_areas(polygons, pts, areas);
                ad_hi = calc_area_distortion(polygons, areas, mareas);

                if (ad_lo < areadist && ad_lo <= ad_hi && fstep < *factor) {
                        *factor -= fstep;
                        areadist = ad_lo;
                } else if (ad_hi < areadist) {
                        *factor += fstep;
                        areadist = ad_hi;
                }
                fstep /= 2;
        }

        /* shift */
        fstep = 0.001;
        for (i = 0; i < 50; i++) { /* optimize 50x */
                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf - fstep, zI, *zIf);
                calc_polygon_areas(polygons, pts, areas);
                ad_lo = calc_area_distortion(polygons, areas, mareas);

                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf + fstep, zI, *zIf);
                calc_polygon_areas(polygons, pts, areas);
                ad_hi = calc_area_distortion(polygons, areas, mareas);

                if (ad_lo < areadist && ad_lo <= ad_hi) {
                        *zRf -= fstep;
                        areadist = ad_lo;
                } else if (ad_hi < areadist) {
                        *zRf += fstep;
                        areadist = ad_hi;
                }

                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf, zI, *zIf - fstep);
                calc_polygon_areas(polygons,pts, areas);
                ad_lo = calc_area_distortion(polygons, areas, mareas);

                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf, zI, *zIf + fstep);
                calc_polygon_areas(polygons,pts, areas);
                ad_hi = calc_area_distortion(polygons, areas, mareas);

                if (ad_lo < areadist && ad_lo <= ad_hi) {
                        *zIf -= fstep;
                        areadist = ad_lo;
                } else if (ad_hi < areadist) {
                        *zIf += fstep;
                        areadist = ad_hi;
                }
                fstep /= 2;
        }
        printf("factor = %f, zRf = %f, zIf = %f\n", *factor, *zRf, *zIf);
        free(mareas);
        free(areas);
        free(pts);
}

/**
 * solver_to_mesh(): gets the solution from the solver, optimizes
 * the solution, & stereographically projects solution into spherical
 * coordinates via stereographic projection
 */
void
solver_to_mesh(polygons_struct *polygons)
{
        Point points[polygons->n_points];
        double factor = 0.0, zRf = 0.0, zIf = 0.0;
        double *zR, *zI;
        double xmin, xmax, ymin, ymax;
        double phi, theta, xx, yy, zz, x, y, z, ct, st, cp, sp;
        int i, p;
        double dist, sphereRadius, xyz[3];

        zR = (double *) malloc(sizeof(double) * polygons->n_points);
        zI = (double *) malloc(sizeof(double) * polygons->n_points);

        for (p = 0; p < polygons->n_points; p++) {
                zR[p] = nlGetVariable(p*2);
                zI[p] = nlGetVariable(p*2 + 1);
        }

        /* optimize the solution */
        moebius_transform(polygons, zR, zI, &factor, &zRf, &zIf);

        project_to_sphere(points, polygons->n_points, factor, zR, zRf, zI, zIf);

        /* map it to a sphere */
        sphereRadius = sqrt(get_polygons_surface_area(polygons) / (4.0 * PI));

        for (p = 0; p < polygons->n_points; p++) {
                xyz[0] = Point_x(points[p]);
                xyz[1] = Point_y(points[p]);
                xyz[2] = Point_z(points[p]);

                dist = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
                if (dist != 0.0) {
                        xyz[0] /= dist;
                        xyz[1] /= dist;
                        xyz[2] /= dist;
                }

                /* Push coordinate onto the sphere */
                Point_x(points[p]) = (sphereRadius * xyz[0]);
                Point_y(points[p]) = (sphereRadius * xyz[1]);
                Point_z(points[p]) = (sphereRadius * xyz[2]);
        }

        /* rough alignment */
        xx = 1; yy = 1; zz = 1;
        if (Point_x(points[xpt]) < 0)
                xx = -1;
        if (Point_y(points[ypt]) < 0)
                yy = -1;
        if (Point_z(points[zpt]) < 0)
                zz = -1;
        for (p = 0; p < polygons->n_points; p++) {
                Point_x(points[p]) *= xx;
                Point_y(points[p]) *= yy;
                Point_z(points[p]) *= zz;
        }

        /* align the maximum z-point (phi) */
        phi = -acos(Point_z(points[zpt]) / sphereRadius);
        cp = cos(phi); sp = sin(phi);
        for (p = 0; p < polygons->n_points; p++) {
                xx = Point_x(points[p]);
                yy = cp*Point_y(points[p]) - sp*Point_z(points[p]);
                zz = sp*Point_y(points[p]) + cp*Point_z(points[p]);

                Point_x(points[p]) = xx;
                Point_y(points[p]) = yy;
                Point_z(points[p]) = zz;
        }

        /* align the maximum x-point (theta & phi) */
        theta = PI + atan(Point_y(points[xpt]) / Point_x(points[xpt]));
        phi = PI/2 - acos(Point_z(points[xpt]) / sphereRadius);
        ct = cos(theta); st = sin(theta);
        cp = cos(phi); sp = sin(phi);
        for (p = 0; p < polygons->n_points; p++) {
                xx = ct*Point_x(points[p]) - st*Point_y(points[p]);
                yy = st*Point_x(points[p]) + ct*Point_y(points[p]);
                zz = Point_z(points[p]);

                x = xx;
                y = cp*yy - sp*zz;
                z = sp*yy + cp*zz;

                Point_x(polygons->points[p]) = x;
                Point_y(polygons->points[p]) = y;
                Point_z(polygons->points[p]) = z;
        }

}

/**
 * finds pointP at the top of the brain, as well as other coordinates
 * for rotating the sphere later on.
 */
int
findPointP(polygons_struct *polygons)
{
        double max_x, max_y, max_z;
        int i;

        xpt = 0; ypt = 0; zpt = 0;
        max_x = Point_coord(polygons->points[0], 0);
        max_y = Point_coord(polygons->points[0], 1);
        max_z = Point_coord(polygons->points[0], 2);

        for (i = 1; i < polygons->n_points; i++) {
                if (Point_x(polygons->points[i]) > max_x) {
                        xpt = i;
                        max_x = Point_x(polygons->points[i]);
                }
                if (Point_y(polygons->points[i]) > max_y) {
                        ypt = i;
                        max_y = Point_y(polygons->points[i]);
                }
                if (Point_z(polygons->points[i]) > max_z) {
                        zpt = i;
                        max_z = Point_z(polygons->points[i]);
                }
        }

        return(zpt);
}


/**
 * generates the matrices for computing Dx=b
 */
void
mesh_to_solver(polygons_struct *polygons)
{
        int npts = polygons->n_points;
        int nfacets = polygons->n_items;
        int i, j, f, n, p, q, r, s, pidx, size, ptP;
        int  **cellinfo;
        int  *n_neighbours, **neighbours;
        struct ptmap *ptinfo;

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        /* find point closest to center-of-mass */
        ptP = findPointP(polygons);

        // 2. store the relationship from point to facet, i.e. for each
        // point, which facets contain it?  For each point in the mesh,
        // generate a vector, storing the id number of triangles containing
        // this point.  The vectors of all the points form a vector:
        // ptinfo

        // 3. store the relationship from cell to point, i.e. for each cell,
        // which points does it contains? store in vector: cellinfo

        ptinfo = (struct ptmap *) malloc(sizeof(struct ptmap) * npts);
        cellinfo = (int **) malloc(sizeof(int *) * nfacets);

        for (p = 0; p < npts; p++)
                ptinfo[p].size = 0;
        for (f = 0; f < nfacets; f++)
                cellinfo[f] = (int *) malloc(sizeof(int) * 3);

        for (f = 0; f < nfacets; f++) {
            size = GET_OBJECT_SIZE(*polygons, f);
            for (p = 0; p < size; p++) {
                pidx = polygons->indices[POINT_INDEX(polygons->end_indices,
                                                     f, p)];

                if (ptinfo[pidx].size == 0) /* initialize */
                        ptinfo[pidx].cells = (int *) malloc(sizeof(int) *
                                                         n_neighbours[pidx]);
                ptinfo[pidx].cells[ptinfo[pidx].size] = f;
                ptinfo[pidx].size++;
                cellinfo[f][p] = pidx;
            }
        }

        // compute b = bR + i*bI separately
        double *bR, *bI;
        int facet = ptinfo[ptP].cells[0];

        bR = (double *) malloc(sizeof(double) * npts);
        bI = (double *) malloc(sizeof(double) * npts);

        for (i = 0; i < npts; i++) {
                bR[i] = bI[i] = 0;
        }

        double A[3], B[3], C[3];
        for (i = 0; i < 3; i++) {
              A[i] = polygons->points[ cellinfo[facet][0] ].coords[i];
              B[i] = polygons->points[ cellinfo[facet][1] ].coords[i];
              C[i] = polygons->points[ cellinfo[facet][2] ].coords[i];
        }
        double ABnorm, CA_BAip; // the inner product of vector C-A and B-A;
        ABnorm = (A[0] - B[0]) * (A[0] - B[0]) +
                 (A[1] - B[1]) * (A[1] - B[1]) +
                 (A[2] - B[2]) * (A[2] - B[2]);
        
        CA_BAip = (C[0] - A[0]) * (B[0] - A[0]) +
                  (C[1] - A[1]) * (B[1] - A[1]) +
                  (C[2] - A[2]) * (B[2] - A[2]);

        double theta = CA_BAip / ABnorm;
        // Here ABnorm is actually the square of AB's norm, which is what we
        // want. So don't bother square the square root.
        
        ABnorm = sqrt(ABnorm); // This is real norm of vector AB.
        
        double E[3];
        for (i = 0; i < 3; i++)
                E[i] = A[i] + theta*(B[i] - A[i]);
        
        double CEnorm;
        CEnorm = (C[0] - E[0]) * (C[0] - E[0]) +
                 (C[1] - E[1]) * (C[1] - E[1]) +
                 (C[2] - E[2]) * (C[2] - E[2]);
        CEnorm = sqrt(CEnorm); // This is real norm of vector CE.
        
        bR[cellinfo[facet][0]] = -1 / ABnorm;
        bR[cellinfo[facet][1]] = 1 / ABnorm;
        
        bI[cellinfo[facet][0]] = (1 - theta) / CEnorm;
        bI[cellinfo[facet][1]] = theta / CEnorm;
        bI[cellinfo[facet][2]] = -1 / CEnorm;

        // 1. Iterate point P from 0 to the last point in the mesh.
        // 2. For each P, find its neighbors
        // 3. For each of P's neighbors, Q, calculate R, S
        // 4. Write the value in matrix.

        progress_struct progress;
        double *Dr;

        Dr = (double *) malloc(sizeof(double) * npts * 2);

        initialize_progress_report(&progress, FALSE, polygons->n_points,
                                   "Mesh2Solver");

        for (p = 0; p < polygons->n_points; p++) {
                for (i = 0; i < npts * 2; i++)
                        Dr[i] = 0;
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                for (n = 0; n < n_neighbours[p]; n++) {
                        q = neighbours[p][n];

                        /* find the two triangles containing both P and Q */
                        int facets[2] = {-1, -1};
                        for (i = 0; i < ptinfo[p].size; i++) {
                                for (j = 0; j < ptinfo[q].size; j++) {
                                        if (ptinfo[p].cells[i] ==
                                            ptinfo[q].cells[j]) {
                                                if (facets[0] < 0)
                                                        facets[0] =
                                                             ptinfo[p].cells[i];
                                                else
                                                        facets[1] =
                                                             ptinfo[p].cells[i];
                                        }
                                }
                        }
        
                        /* get R and S */
                        size = GET_OBJECT_SIZE(*polygons, facets[0]);
                        for (i = 0; i < size; i++) {
                                pidx = polygons->indices[
                                              POINT_INDEX(polygons->end_indices,
                                                          facets[0], i)];
                                if (pidx != p && pidx != q)
                                        r = pidx;
                        }
                        size = GET_OBJECT_SIZE(*polygons, facets[1]);
                        for (i = 0; i < size; i++) {
                                pidx = polygons->indices[
                                             POINT_INDEX(polygons->end_indices,
                                                         facets[1], i)];
                                if (pidx != p && pidx != q)
                                        s = pidx;
                        }

                        double P[3], Q[3], R[3], S[3];

                        for (i = 0; i < 3; i++) {
                                P[i] = polygons->points[p].coords[i];
                                Q[i] = polygons->points[q].coords[i];
                                R[i] = polygons->points[r].coords[i];
                                S[i] = polygons->points[s].coords[i];
                        }
                
                        double SP[3], SQ[3], RP[3], RQ[3];
                        double SPnorm = 0, SQnorm = 0, RPnorm = 0;
                        double RQnorm = 0, SPSQinnerProd = 0, RPRQinnerProd = 0;

                        for (i = 0; i < 3; ++i) {
                                SP[i] = P[i] - S[i];
                                SPnorm += SP[i] * SP[i];

                                SQ[i] = Q[i] - S[i];
                                SQnorm += SQ[i] * SQ[i];
                                SPSQinnerProd += SP[i] * SQ[i];

                                RP[i] = P[i] - R[i];
                                RPnorm += RP[i] * RP[i];

                                RQ[i] = Q[i] - R[i];
                                RQnorm += RQ[i] * RQ[i];
                                RPRQinnerProd += RP[i] * RQ[i];
                        }
                        SPnorm = sqrt(SPnorm);
                        SQnorm = sqrt(SQnorm);
                        RPnorm = sqrt(RPnorm);
                        RQnorm = sqrt(RQnorm);
                
                        double cosS = SPSQinnerProd / (SPnorm * SQnorm);
                        double cosR = RPRQinnerProd / (RPnorm * RQnorm);
                        double ctgS = cosS / sqrt(1 - cosS*cosS);
                        double ctgR = cosR / sqrt(1 - cosR*cosR);

                        Dr[q] = -0.5 * (ctgS + ctgR);
                        Dr[p] += 0.5 * (ctgS + ctgR);
                }

                // write the real values
                nlRowParameterd(NL_RIGHT_HAND_SIDE, bR[p]);
                nlBegin(NL_ROW);

                for (i = 0; i < npts; i++) {
                        if (Dr[i] != 0) { // store it
                                nlCoefficient(i*2, Dr[i]);
                        }
                }
                nlEnd(NL_ROW);
            
                // write the imaginary values
                nlRowParameterd(NL_RIGHT_HAND_SIDE, -bI[p]);
                nlBegin(NL_ROW);
                for (i = 0; i < npts; i++) {
                        if (Dr[i] != 0) {
                                nlCoefficient(i*2 + 1, Dr[i]);
                        }
                }
                nlEnd(NL_ROW);
                update_progress_report(&progress, p);
        }
        terminate_progress_report(&progress);
}

/* Apply the Laplace-Beltrami operator */
void find_conformal_map(polygons_struct *polygons) {
        int i;
        double time;

        int eulerNum = euler_characteristic(polygons);
        if (eulerNum != 2) {
                printf("WARNING: Euler characteristics is %d, not 2!  Not genus 0 surface...\n", eulerNum);
        }

        nlNewContext();
        if (nlInitExtension("SUPERLU")) {
                nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
                printf("Using SuperLU, cool !\n");
        } else {
                nlSolverParameteri(NL_SOLVER, NL_CG);
                nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_NONE);
                printf("Using conjugate gradient\n");
        }
        nlSolverParameteri(NL_NB_VARIABLES, 2*polygons->n_points);
        nlSolverParameteri(NL_LEAST_SQUARES, NL_FALSE);
        nlSolverParameteri(NL_MAX_ITERATIONS, 100*polygons->n_points);
        /*nlSolverParameteri(NL_SYMMETRIC, NL_FALSE); */
        nlSolverParameterd(NL_THRESHOLD, 1e-10);

        nlBegin(NL_SYSTEM);
        for (i = 0; i < 2*polygons->n_points; i++)
                nlSetVariable(i, 0);

        nlBegin(NL_MATRIX);
        mesh_to_solver(polygons); /* calculate the matrices and send it to solver */
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);

        printf("Solving ...\n");
        nlSolve();
        solver_to_mesh(polygons);

        nlGetDoublev(NL_ELAPSED_TIME, &time);
        printf("Solver time: %f\n", time);
        nlDeleteContext(nlGetCurrent());

        compute_polygon_normals(polygons);
}
