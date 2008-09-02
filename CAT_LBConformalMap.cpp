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

#undef MOEBIUS

extern "C"
{   
        #include "nl/nl.h"
        #include <volume_io/internal_volume_io.h>
        #include <bicpl.h>
        #include <CAT_Curvature.h>
        #include <CAT_Blur2d.h>
        #include <CAT_Surf.h>
        #include <CAT_SurfaceIO.h>
}   

#include <ParseArgv.h>

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cfloat>

#include <math.h>
#include <assert.h>

struct ptmap {
        int *cells;
        int size;
};

double mapScale = -1.0;
int ptP = -1;
polygons_struct *polygons;

static ArgvInfo argTable[] = {
  {"-pointP", ARGV_INT, (char *) 0, (char *) &ptP,
       "point index of the locked point P [default: optimized using COM]."},
  {"-mapscale", ARGV_FLOAT, (char *) 1, (char *) &mapScale,
       "scaling factor [default: optimized]. Set mapScale to turn off optimization, 0 to estimate empirically."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


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
calc_polygon_areas(Point *pts, double *areas)
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
calc_area_distortion(double *area1, double *area2)
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
moebius_transform(double *zR, double *zI, double *factor,
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
        for (f = 0; f < polygons->n_items; f++) {
                int size = get_polygon_points(polygons, f, tpts);
                mareas[f] = get_polygon_surface_area(size, tpts);
        }

        /* get an estimate of the optimal factor value */
        for (f = 400; f <= 10000; ) { // get in the ballpark first
                project_to_sphere(pts, polygons->n_points, f, zR, 0.0, zI, 0.0);
                calc_polygon_areas(pts, areas);
                new_ad = calc_area_distortion(areas, mareas);
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
                calc_polygon_areas(pts, areas);
                ad_lo = calc_area_distortion(areas, mareas);

                project_to_sphere(pts, polygons->n_points, *factor + fstep,
                                  zR, 0.0, zI, 0.0);
                calc_polygon_areas(pts, areas);
                ad_hi = calc_area_distortion(areas, mareas);

                if (ad_lo < areadist && ad_lo <= ad_hi && fstep < *factor) {
                        *factor -= fstep;
                        areadist = ad_lo;
                } else if (ad_hi < areadist) {
                        *factor += fstep;
                        areadist = ad_hi;
                }
                fstep /= 2;
        }
        printf("f = %f, ad: %f\n", *factor, areadist);

        /* shift */
        fstep = 0.001;
        for (i = 0; i < 50; i++) { /* optimize 50x */
                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf - fstep, zI, *zIf);
                calc_polygon_areas(pts, areas);
                ad_lo = calc_area_distortion(areas, mareas);

                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf + fstep, zI, *zIf);
                calc_polygon_areas(pts, areas);
                ad_hi = calc_area_distortion(areas, mareas);

                if (ad_lo < areadist && ad_lo <= ad_hi) {
                        *zRf -= fstep;
                        areadist = ad_lo;
                } else if (ad_hi < areadist) {
                        *zRf += fstep;
                        areadist = ad_hi;
                }

                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf, zI, *zIf - fstep);
                calc_polygon_areas(pts, areas);
                ad_lo = calc_area_distortion(areas, mareas);

                project_to_sphere(pts, polygons->n_points, *factor,
                                  zR, *zRf, zI, *zIf + fstep);
                calc_polygon_areas(pts, areas);
                ad_hi = calc_area_distortion(areas, mareas);

                if (ad_lo < areadist && ad_lo <= ad_hi) {
                        *zIf -= fstep;
                        areadist = ad_lo;
                } else if (ad_hi < areadist) {
                        *zIf += fstep;
                        areadist = ad_hi;
                }
                fstep /= 2;
        }
        printf("zRf = %f, zIf = %f, ad: %f\n", *zRf, *zIf, areadist);
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
solver_to_mesh()
{
        Point points[polygons->n_points];
        double factor = 0.0, zRf = 0.0, zIf = 0.0;
        double *zR, *zI;
        double xmin, xmax, ymin, ymax;
        int i, p;
        double dist, sphereRadius, xyz[3];

        zR = (double *) malloc(sizeof(double) * polygons->n_points);
        zI = (double *) malloc(sizeof(double) * polygons->n_points);

        for (p = 0; p < polygons->n_points; p++) {
                zR[p] = nlGetVariable(p*2);
                zI[p] = nlGetVariable(p*2 + 1);
        }

        /* empirical best guess for mapScale */
        /* mapScale = round(pow(polygons->n_points, 0.5719)); */

        if (mapScale >= 0) {
                printf("MapScale = %f\n", mapScale);
                xmin =  zR[0];
                xmax =  zR[0];
                ymin =  zI[0];
                ymax =  zI[0];

                /* get approximate scale */
                for (p = 1; p < polygons->n_points;  p++) {
                        xmin = (xmin < zR[p]) ? xmin : zR[p];
                        xmax = (xmax > zR[p]) ? xmax : zR[p];
                        ymin = (ymin < zI[p]) ? ymin : zI[p];
                        ymax = (ymax > zI[p]) ? ymax : zI[p];
                }

                xmax = ( fabs(xmin) > fabs(xmax) ) ? fabs(xmin) : fabs(xmax);
                ymax = ( fabs(ymin) > fabs(ymax) ) ? fabs(ymin) : fabs(ymax);

                /* the factor is used to re-scale the points in the plane */
                factor = mapScale / ( (xmax > ymax) ? xmax : ymax);
        } else { /* optimize the solution */
                moebius_transform(zR, zI, &factor, &zRf, &zIf);
        }

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
                xyz[0] = (sphereRadius * xyz[0]);
                xyz[1] = (sphereRadius * xyz[1]);
                xyz[2] = (sphereRadius * xyz[2]);
                for (i = 0; i < 3; i++)
                        Point_coord(polygons->points[p], i) = xyz[i];
        }
}

/**
 * finds pointP as close to the COM as possible
 */
int
findPointP()
{
        Point p;
        Real dist, closest_dist = 0.0;
        double c[3] = {0.0, 0.0, 0.0};
        int i, j, closest_pt;

        for (i = 0; i < polygons->n_points; i++) {
                for (j = 0; j < 3; j++)
                        c[j] += Point_coord(polygons->points[i], j);
        }

        Point_x(p) = c[0] / polygons->n_points;
        Point_y(p) = c[1] / polygons->n_points;
        Point_z(p) = c[2] / polygons->n_points;

        for (i = 0; i < polygons->n_points; i++) {
                dist = sq_distance_between_points(&p, &polygons->points[i]);
                if (i == 0 || dist < closest_dist) {
                        closest_pt = i;
                        closest_dist = dist;
                }
        }
        return(closest_pt);
}


/**
 * generates the matrices for computing Dx=b
 */
void
mesh_to_solver()
{
        int npts = polygons->n_points;
        int nfacets = polygons->n_items;
        int i, f, n, p, pidx, size;
        int  **cellinfo;
        int  *n_neighbours, **neighbours;
        struct ptmap *ptinfo;

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        // 1. store the points coordinates: pointXYZ
        std::vector< std::vector<double> > pointXYZ(npts, std::vector<
                                                          double>(3, 0) );

        for (p = 0; p < npts; p++) {
              pointXYZ[p][0] = polygons->points[p].coords[0];
              pointXYZ[p][1] = polygons->points[p].coords[1];
              pointXYZ[p][2] = polygons->points[p].coords[2];            
        }
          
        /* find point closest to center-of-mass */
        if (ptP < 0)
                ptP = findPointP();

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
        std::vector <std::vector<int> > pointCell(npts);

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
                pointCell[pidx].push_back(f);
                cellinfo[f][p] = pidx;
            }
        }

        std::cerr << "  Checking the existence of boundary...";
        for (p = 0; p < npts; p++) {
            if (pointCell[p].size() < 3) {
                    // If one node has <= 2 neighbors, it's on the boundary.
                    // This is the sufficient condition, i.e., it may still be a
                    // boundary point even having more than 3 neighbors.
                    // So, what's the equivalent expression of being a boundary point?
                    std::cerr<<"There is boundary in mesh! exiting..."<<std::endl;
                    exit(-1);
            }
        }
        std::cerr<<"No boundary found!"<<std::endl;
        
        // compute b = bR + i*bI separately
        std::vector<double> bR(npts), bI(npts);
        int facet = pointCell[ptP][0];

        std::vector<double> A( pointXYZ[ cellinfo[facet][0] ] ),
                            B( pointXYZ[ cellinfo[facet][1] ] ),
                            C( pointXYZ[ cellinfo[facet][2] ] );
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
        
        std::vector<double> E(3);
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
        
        std::vector< std::vector<int> >::iterator itP, itPEnd = pointCell.end();
        int idP = 0;
        unsigned long numOfEdges = 0;

        progress_struct progress;
        initialize_progress_report(&progress, FALSE, npts, "Mesh2Solver");
        for (itP = pointCell.begin(); itP != itPEnd; ++itP, ++idP) {
            update_progress_report(&progress, idP);
            std::vector<int> neighborOfP;
            std::vector<double> Dr(npts*2); // matrix row values
            // for each point P, traverse all cells containing it.
            std::vector<int>::iterator itCellEnd = (*itP).end();
            
            for (f = 0; f < ptinfo[idP].size; f++) {
                facet = ptinfo[idP].cells[f];
                // for each cell containing P, store the point with larger point Id.
                // only three points, don't use for-loop to save time.
                if (cellinfo[facet][0] != idP)
                    neighborOfP.push_back(cellinfo[facet][0]);
                if (cellinfo[facet][1] != idP)
                    neighborOfP.push_back(cellinfo[facet][1]);
                if (cellinfo[facet][2] != idP)
                    neighborOfP.push_back(cellinfo[facet][2]);
            } // for itCell. Ok, now all neighbors of P is stored in neighborOfP

            sort(neighborOfP.begin(), neighborOfP.end());
            std::vector<int>::iterator it;
            it = unique(neighborOfP.begin(), neighborOfP.end());
            neighborOfP.erase(it, neighborOfP.end());
            

            // //-----------------------------------------------
            // // print out the neighbors
            // std::vector<int>::iterator itNeighbor = neighborOfP.begin();
            // std::vector<int>::iterator itNeighborEnd = neighborOfP.end();
            // std::cerr<<"The neighbors of "<<idP<<" are: ";
            // for (; itNeighbor != itNeighborEnd; ++itNeighbor) {
            //     std::cerr<<*itNeighbor<<" , ";
            // }
            // std::cerr<<std::endl;
            // // ----------------------------------------------------
            
            
            // next, from P to each neighbor...
            // note: itP and itQ point at different type of vectors...
            // *itP is a vector containing a list of cell Ids, all of which contains point P
            // idP is the point Id of P
            // *itQ is the point Id of Q (so idP and *itQ are same type)
            std::vector<int>::iterator itQ, itQEnd = neighborOfP.end();
            for (itQ = neighborOfP.begin(); itQ != itQEnd; ++itQ) {
                if (*itQ > idP)
                    numOfEdges++;

                // first check whether PQ is a boundary edge:
                std::vector<int> cellsContainingP(*itP), cellsContainingQ(pointCell[*itQ]);
                std::vector<int> cells(cellsContainingP.size() + cellsContainingQ.size());
                std::vector<int>::iterator itv, endIter;
                
                sort(cellsContainingP.begin(), cellsContainingP.end());
                sort(cellsContainingQ.begin(), cellsContainingQ.end());
                
                endIter = set_intersection(cellsContainingP.begin(), cellsContainingP.end(),
                                           cellsContainingQ.begin(), cellsContainingQ.end(),
                                           cells.begin());
                cells.erase(endIter, cells.end());
                if (cells.size() != 2) continue;
                // If P and Q are not shared by two triangles, i.e. 1: are not
                // connected by and edge, or, 2: are on the surface boundary
                // thus only shared by one triangle. then skip.  However, in
                // this paper the surface is closed thus there is not boundary.
                
                // If passed test above, then P and Q are two valid points.
                // i.e. PQ is a valid edge.  i.e. cells now contain two int's,
                // which are the Id of the triangles containing P and Q
                
                
                // //------------------------------------------------------------
                // //print out valid edge
                // std::cerr<<idP<<" and "<<*itQ<<" are two valid points"<<std::endl;
                // std::cerr<<(endIter == cells.end())<<std::endl;
                // //-----------------------------------------------------------
                
                
                // Next we extract R and S from cells
                int itS, itR; // the Id of point S and R;

                for (int it = 0; it < 3; ++it) {
                    if (cellinfo[cells[0]][it] != idP && cellinfo[cells[0]][it] != *itQ)
                        itS = cellinfo[cells[0]][it];
                    if (cellinfo[cells[1]][it] != idP && cellinfo[cells[1]][it] != *itQ) 
                        itR = cellinfo[cells[1]][it];
                }

                std::vector<double> P(pointXYZ[idP]),
                        Q(pointXYZ[*itQ]),
                        R(pointXYZ[itR]),
                        S(pointXYZ[itS]);
                
                std::vector<double> SP(3), SQ(3), RP(3), RQ(3);
                double SPnorm = 0, SQnorm = 0, RPnorm = 0, RQnorm = 0, SPSQinnerProd = 0, RPRQinnerProd = 0;
                for (int it = 0; it < 3; ++it) {
                    SP[it] = P[it] - S[it];
                    SPnorm += SP[it] * SP[it];

                    SQ[it] = Q[it] - S[it];
                    SQnorm += SQ[it] * SQ[it];
                    SPSQinnerProd += SP[it] * SQ[it];

                    RP[it] = P[it] - R[it];
                    RPnorm += RP[it] * RP[it];

                    RQ[it] = Q[it] - R[it];
                    RQnorm += RQ[it] * RQ[it];
                    RPRQinnerProd += RP[it] * RQ[it];
                }
                SPnorm = sqrt(SPnorm);
                SQnorm = sqrt(SQnorm);
                RPnorm = sqrt(RPnorm);
                RQnorm = sqrt(RQnorm);
                
                double cosS = SPSQinnerProd / (SPnorm * SQnorm);
                double cosR = RPRQinnerProd / (RPnorm * RQnorm);
                double ctgS = cosS / sqrt(1 - cosS*cosS);
                double ctgR = cosR / sqrt(1 - cosR*cosR);
                
                Dr[*itQ] = -0.5 * (ctgS + ctgR);
                Dr[idP] += 0.5 * (ctgS + ctgR);
                // add to the diagonal element of this line.
                
                //D(*itQ, idP) = -0.5*(ctgS + ctgR); // symmetric
                //D(*itQ, *itQ) += 0.5*(ctgS + ctgR);
            }
            
            // write the real values
            nlRowParameterd(NL_RIGHT_HAND_SIDE, bR[idP]);
            nlBegin(NL_ROW);

            for (int i = 0; i < npts; i++) {
                if (Dr[i] != 0) { // store it
                    nlCoefficient(i*2, Dr[i]);
                }
            }
            nlEnd(NL_ROW);
            
            // write the imaginary values
            nlRowParameterd(NL_RIGHT_HAND_SIDE, -bI[idP]);
            nlBegin(NL_ROW);
            for (int i = 0; i < npts; i++) {
                if (Dr[i] != 0) {
                    nlCoefficient(i*2 + 1, Dr[i]);
                }
            }
            nlEnd(NL_ROW);
        }

        terminate_progress_report(&progress);

        /* calc Euler number to test whether mesh is genus 0 (Euler num = 2) */
        int eulerNum = npts - numOfEdges + nfacets;

        fprintf(stderr, "  Calculating Euler characteristics......\n");
        fprintf(stderr, "    Euler Characteristics = %d\n", eulerNum);
        fprintf(stderr, "    genus = %d\n", (2.0 - eulerNum)/2);
}

/* Apply the Laplace-Beltrami operator */
void find_conformal_map() {
        int i;
        double time;

        nlNewContext();
        if (nlInitExtension("SUPERLU")) {
                nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
                fprintf(stderr, "Using SuperLU, cool !\n");
        } else {
                nlSolverParameteri(NL_SOLVER, NL_CG);
                nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_NONE);
                fprintf(stderr, "Using conjugate gradient\n");
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
        mesh_to_solver(); /* calculate the matrices and send it to solver */
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);

        fprintf(stderr, "Solving ...\n");
        nlSolve();
        solver_to_mesh();

        nlGetDoublev(NL_ELAPSED_TIME, &time);
        fprintf(stderr, "Solver time: %f\n", time);
        nlDeleteContext(nlGetCurrent());
}

int
main(int argc, char** argv)
{
        char *input_file, *output_file;
        object_struct **objects;
        int n_objects;
        Real *curvatures;
        File_formats format;
    
        /* Get arguments */
        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
                fprintf(stderr, "\n");
                fprintf(stderr, "Usage: %s [options] infile.obj outfile.obj\n",
                        argv[0]);
                fprintf(stderr, "\nGenerate the Laplace-Beltrami conformal ");
                fprintf(stderr, "map of the mesh specified in\ninfile.obj.  ");
                fprintf(stderr, "Results are saved in outfile.obj.\n\n");
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);
        get_string_argument(NULL, &input_file);
        get_string_argument(NULL, &output_file);

        if (input_graphics_any_format(input_file, &format,
                                      &n_objects, &objects) != OK) {
                printf("Error reading input file\n");
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        int eulerNum = euler_characteristic(polygons);
        if (eulerNum != 2) {
                fprintf(stderr, "WARNING: Euler characteristics is %d, not 2!  Not genus 0 surface...\n", eulerNum);
        }

        find_conformal_map();

        compute_polygon_normals(polygons);
        output_graphics_any_format(output_file, format, 1, objects);
}
