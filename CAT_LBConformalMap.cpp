
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

double ms = -1.0;

static ArgvInfo argTable[] = {
  {"-mapscale", ARGV_FLOAT, (char *) 1, (char *) &ms,
       "scaling factor [default: optimized]. Set mapScale to turn off optimization, 0 to estimate empirically."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

class LSCM {
public:

    LSCM(polygons_struct &p, double ms) : polygons(&p), mapScale(ms) { }

    // Apply the Laplace-Beltrami operator
    void apply() {
        int nb_vertices = polygons->n_points;
        nlNewContext();
        if(nlInitExtension("SUPERLU")) {
            nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
            std::cerr << "Using SuperLU, cool !" << std::endl;
        } else {
            nlSolverParameteri(NL_SOLVER, NL_CG) ;
            nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_NONE ); //NL_PRECOND_JACOBI) ;
            std::cerr << "Using Jacobi pre-conditioned conjugate gradient" 
                      << std::endl;
        }
        nlSolverParameteri(NL_NB_VARIABLES, 2*nb_vertices);
        nlSolverParameteri(NL_LEAST_SQUARES, NL_FALSE);
        nlSolverParameteri(NL_MAX_ITERATIONS, 100*nb_vertices);
        //nlSolverParameteri(NL_SYMMETRIC, NL_FALSE);
        nlSolverParameterd(NL_THRESHOLD, 1e-10);
        nlBegin(NL_SYSTEM);
        for(unsigned int i=0; i<nb_vertices*2; i++) {
            nlSetVariable(i, 0);
        }
        nlBegin(NL_MATRIX);
        mesh_to_solver(); // calculate the matrices and send it to the solver
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        std::cerr << "Solving ..." << std::endl;
        nlSolve();
        solver_to_mesh();
        double time;
        nlGetDoublev(NL_ELAPSED_TIME, &time);
        std::cerr << "Solver time: " << time << std::endl;
        nlDeleteContext(nlGetCurrent());
    }

protected:
    /**
     * calc_projection(): calculate the projection from the complex plane
     * to the sphere.  np = number of points, f = scaling factor,
     * zRf = translation of real value, zIf = translation of imaginary value.
     */
    void calc_projection(Point *points, int np, double f,
                         std::vector<double> zR, double zRf,
                         std::vector<double> zI, double zIf) {
        for (int it = 0; it < np; ++it) {
            double zRit = f * (zR[it] + zRf);
            double zIit = f * (zI[it] + zIf);

            double r2 = zRit*zRit + zIit*zIit;
            Point_x(points[it]) = 2*zRit/(1+r2);
            Point_y(points[it]) = 2*zIit/(1+r2);
            Point_z(points[it]) = 2*r2/(1+r2) - 1;
        }
    }


    /* calc_polygon_areas(): calculate the polygon areas, store in areas */
    void calc_polygon_areas(Point *points, double *areas) {
        Point tpoints[3];

        for (int poly = 0; poly < polygons->n_items; poly++) {
            tpoints[0] = points[polygons->indices[
                                POINT_INDEX(polygons->end_indices,poly,0)]];
            tpoints[1] = points[polygons->indices[
                                POINT_INDEX(polygons->end_indices,poly,1)]];
            tpoints[2] = points[polygons->indices[
                                POINT_INDEX(polygons->end_indices,poly,2)]];
            areas[poly] = get_polygon_surface_area(3, tpoints);
        }
    }

    /* calc_area_distortion(): calculate the area distortion.  calculates
     * log difference of distortion between area1 and area2.
     */
    double calc_area_distortion(double *area1, double *area2) {
        double area_distortion = 0.0;
        for (int i = 0; i < polygons->n_items; i++) {
            if (area2[i] > 0 && area1[i] > 0) { // get sure that areas are > 0
                area_distortion += fabs(log10(area1[i]/area2[i]));
            }
        }
        return area_distortion;
    }

    /**
     * solver_to_mesh(): gets the solution from the solver, optimizes
     * the solution, & stereographically projects solution into spherical
     * coordinates via stereographic projection
     */
    void solver_to_mesh() {
        Point points[polygons->n_points];
        double areas[polygons->n_items];
        double factor = 0.0;
        double zRf = 0.0; double zIf = 0.0;
        std::vector<double> zR(polygons->n_points), zI(polygons->n_points);

        for (int i = 0; i < polygons->n_points; i++) {
            zR[i] = nlGetVariable(i*2);
            zI[i] = nlGetVariable(i*2 + 1);
        }

        if (mapScale >= 0) {
            std::cout << "MapScale = " << mapScale << std::endl;
            double xmin =  zR[0]; double ymin =  zI[0];
            double xmax =  zR[0]; double ymax =  zI[0];

            // get approximate scale
            for (int i = 1; i < polygons->n_points;  i++) {
                xmin = (xmin<zR[i]) ? xmin : zR[i];
                xmax = (xmax>zR[i]) ? xmax : zR[i];
                ymin = (ymin<zI[i]) ? ymin : zI[i];
                ymax = (ymax>zI[i]) ? ymax : zI[i];
            }

            xmax = ( fabs(xmin)>fabs(xmax) ) ? fabs(xmin) : fabs(xmax);
            ymax = ( fabs(ymin)>fabs(ymax) ) ? fabs(ymin) : fabs(ymax);

            // the factor is used to re-scale the points in the plane.
            factor = mapScale/( (xmax>ymax) ? xmax : ymax);
        } else { // optimize the solution
            // first get all of the areas of the original triangles
            double mareas[polygons->n_items];
            for (int i = 0; i < polygons->n_items; i++) {
                Point tpoints[3];
                int size = get_polygon_points(polygons, i, tpoints);
                mareas[i] = get_polygon_surface_area(size, tpoints);
            }

            double area_distortion = FLT_MAX;

            // get an estimate of the optimal factor value
            for (int f = 400; f <= 10000; ) { // get in the ballpark first
                calc_projection(points, polygons->n_points, f, zR, 0.0, zI, 0.0);
                calc_polygon_areas(points, areas);
                double new_ad = calc_area_distortion(areas, mareas);
                if (new_ad < area_distortion) {
                    factor = f;
                    area_distortion = new_ad;
                }
                f += 200;
            }

            double fstep = 100;
            for (int i = 0; i < 5; i++) { // optimize factor 5x more
                calc_projection(points, polygons->n_points, factor - fstep, zR, 0.0, zI, 0.0);
                calc_polygon_areas(points, areas);
                double ad_lo = calc_area_distortion(areas, mareas);
                calc_projection(points, polygons->n_points, factor + fstep, zR, 0.0, zI, 0.0);
                calc_polygon_areas(points, areas);
                double ad_hi = calc_area_distortion(areas, mareas);
                if (ad_lo < area_distortion && ad_lo <= ad_hi && fstep < factor) {
                    factor -= fstep;
                    area_distortion = ad_lo;
                } else if (ad_hi < area_distortion) {
                    factor += fstep;
                    area_distortion = ad_hi;
                }
                fstep /= 2;
            }
            std::cout<<"f = "<<factor<<", ad: "<<area_distortion<<std::endl;


            // shift
            fstep = 0.001;
            for (int i = 0; i < 50; i++) { // optimize 5x
                calc_projection(points, polygons->n_points, factor, zR, zRf - fstep, zI, zIf);
                calc_polygon_areas(points, areas);
                double ad_lo = calc_area_distortion(areas, mareas);
                calc_projection(points, polygons->n_points, factor, zR, zRf + fstep, zI, zIf);
                calc_polygon_areas(points, areas);
                double ad_hi = calc_area_distortion(areas, mareas);
                if (ad_lo < area_distortion && ad_lo <= ad_hi) {
                    zRf -= fstep;
                    area_distortion = ad_lo;
                } else if (ad_hi < area_distortion) {
                    zRf += fstep;
                    area_distortion = ad_hi;
                }

                calc_projection(points, polygons->n_points, factor, zR, zRf, zI, zIf - fstep);
                calc_polygon_areas(points, areas);
                ad_lo = calc_area_distortion(areas, mareas);
                calc_projection(points, polygons->n_points, factor, zR, zRf, zI, zIf + fstep);
                calc_polygon_areas(points, areas);
                ad_hi = calc_area_distortion(areas, mareas);
                if (ad_lo < area_distortion && ad_lo <= ad_hi) {
                    zIf -= fstep;
                    area_distortion = ad_lo;
                } else if (ad_hi < area_distortion) {
                    zIf += fstep;
                    area_distortion = ad_hi;
                }
                fstep /= 2;
            }
            std::cout << "zRf = " << zRf << ", zIf =" << zIf << ", ad: " << area_distortion << std::endl;

        }

        calc_projection(points, polygons->n_points, factor, zR, zRf, zI, zIf);

        // map it to a sphere
        double original_area = get_polygons_surface_area(polygons);
        double sphereRadius = sqrt(original_area / (4.0 * PI));
        double xyz[3] = { 0.0, 0.0, 0.0 };        

        for (int i = 0; i < polygons->n_points; i++) {
            xyz[0] = Point_x(points[i]);
            xyz[1] = Point_y(points[i]);
            xyz[2] = Point_z(points[i]);

            double dist = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
            if (dist != 0.0) {
                xyz[0] /= dist;
                xyz[1] /= dist;
                xyz[2] /= dist;
            }

            // Push coordinate onto the sphere
            xyz[0] = (sphereRadius * xyz[0]);
            xyz[1] = (sphereRadius * xyz[1]);
            xyz[2] = (sphereRadius * xyz[2]);
            for (int j = 0; j < 3; j++) {
                Point_coord(polygons->points[i],j) = xyz[j];
            }
        }

    }

    /**
     * finds pointP as close to the COM as possible
     */
    int findPointP() {
        Point p;
        Real dist, closest_dist = 0.0;
        double c[3] = {0.0, 0.0, 0.0};
        int closest_pt;

        for (int i = 0; i < polygons->n_points; i++) {
            for (int j = 0; j< 3; j++)
                c[j] += Point_coord(polygons->points[i], j);
        }

        Point_x(p) = c[0]/polygons->n_points;
        Point_y(p) = c[1]/polygons->n_points;
        Point_z(p) = c[2]/polygons->n_points;

        for (int i = 0; i < polygons->n_points; i++) {
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
    void mesh_to_solver() {
        int numOfPoints = polygons->n_points;
        int numOfFacets = polygons->n_items;
        
        
        // 1. store the points coordinates: pointXYZ
        std::vector< std::vector<double> > pointXYZ(numOfPoints, std::vector<
                                                          double>(3, 0) );

        for (int it = 0; it < numOfPoints; ++it) {
              double x = polygons->points[it].coords[0];
              double y = polygons->points[it].coords[1];
              double z = polygons->points[it].coords[2];            

              pointXYZ[it][0] = x;
              pointXYZ[it][1] = y;
              pointXYZ[it][2] = z;
        }
          
        // find point closest to center-of-mass
        int pointP = findPointP();

        // 2. store the relationship from point to facet, i.e. for each
        // point, which facets contain it?  For each point in the mesh,
        // generate a vector, storing the id number of triangles containing
        // this point.  The vectors of all the points form a vector:
        // pointCell

        // 3. store the relationship from cell to point, i.e. for each cell,
        // which points does it contains? store in vector: cellPoint

        std::vector <std::vector<int> > pointCell(numOfPoints);
        std::vector <std::vector<int> > cellPoint(numOfFacets, std::vector<int>(3,0));

        int size;
        for (unsigned int itCell = 0; itCell < numOfFacets; itCell++) {
            size = GET_OBJECT_SIZE(*polygons, itCell);
            for (unsigned int itPntInCell = 0; itPntInCell < size; itPntInCell++) {
                int p = polygons->indices[POINT_INDEX(polygons->end_indices, itCell, itPntInCell)];

                pointCell[p].push_back(itCell);
                cellPoint[itCell][itPntInCell] = p;
            }
        }

        // //--------------------------------------------------------------
        // // print out the result for debuging
        // std::cout<<std::endl;
        // std::cout<<std::endl;
        // std::cout<<std::endl;
        // std::cout<<std::endl;
        // 
        // for (int it = 0; it < numOfPoints; ++it) {
        //     std::cout<<"point# "<<it<<" :"                 <<"       X: "<<pointXYZ[it][0]
        //     <<"       Y: "<<pointXYZ[it][1]
        //     <<"       Z: "<<pointXYZ[it][2]<<std::endl;
        // }   
        // 
        // for (int it = 0; it < numOfPoints; ++it) {
        //     std::cout<<"point# "<<it<<" is contained by    "<<pointCell[it].size()<<"   cells:"<<std::endl;
        //     for (std::vector<int>::const_iterator vi = pointCell[it].begin();
        //          vi != pointCell[it].end();
        //          ++vi) {
        //         std::cout<<*vi<<"     ";
        //     }
        //     std::cout<<std::endl;
        // }
        // 
        // for (int it = 0; it < numOfFacets; ++it) {
        //     std::cout<<"cell# "<<it<<" has points: "
        //     <<cellPoint[it][0]<<"  "
        //     <<cellPoint[it][1]<<"  "
        //     <<cellPoint[it][2]<<std::endl;
        // }
        //---------------------------------------------------------------


                 
        std::cerr << "  Checking the existence of boundary...";
        for (int itPointCell = 0; itPointCell < numOfPoints; itPointCell++) {
            if (pointCell[itPointCell].size() < 3) {
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
        std::vector<double> bR(numOfPoints), bI(numOfPoints);
        int facet = pointCell[pointP][0];
        //     for (std::vector<int>::const_iterator vi = pointCell[it].begin();
        std::vector<double> A( pointXYZ[ cellPoint[ facet ][ 0 ] ] ),
        B( pointXYZ[ cellPoint[ facet][ 1 ] ] ),
        C( pointXYZ[ cellPoint[ facet ][ 2 ] ] );
        double ABnorm, CA_BAip; // the inner product of vector C-A and B-A;
        ABnorm = (A[0] - B[0]) * (A[0] - B[0])
        + (A[1] - B[1]) * (A[1] - B[1])
        + (A[2] - B[2]) * (A[2] - B[2]);
        
        CA_BAip = (C[0] - A[0]) * (B[0] - A[0])
        + (C[1] - A[1]) * (B[1] - A[1])
        + (C[2] - A[2]) * (B[2] - A[2]);
        
        double theta = CA_BAip / ABnorm;
        // Here ABnorm is actually the square of AB's norm, which is what we
        // want. So don't bother square the square root.
        
        ABnorm = sqrt(ABnorm); // This is real norm of vector AB.
        
        std::vector<double> E(3);
        for (int it = 0; it < 3; ++it)
            E[it] = A[it] + theta*(B[it] - A[it]);
        
        double CEnorm;
        CEnorm = (C[0] - E[0]) * (C[0] - E[0])
        + (C[1] - E[1]) * (C[1] - E[1])
        + (C[2] - E[2]) * (C[2] - E[2]);
        CEnorm = sqrt(CEnorm); // This is real norm of vector CE.
        
        bR[cellPoint[ facet ][0]] = -1 / ABnorm;
        bR[cellPoint[ facet ][1]] = 1 / ABnorm;
        
        bI[cellPoint[ facet ][0]] = (1-theta)/ CEnorm;
        bI[cellPoint[ facet ][1]] = theta/ CEnorm;
        bI[cellPoint[ facet ][2]] = -1 / CEnorm;
        
        // 1. Iterate point P from 0 to the last point in the mesh.
        // 2. For each P, find its neighbors
        // 3. For each of P's neighbors, Q, calculate R, S
        // 4. Write the value in matrix.
        
        std::vector< std::vector<int> >::iterator itP, itPEnd = pointCell.end();
        int idP = 0;
        unsigned long numOfEdges = 0;

        progress_struct progress;
        initialize_progress_report(&progress, FALSE, numOfPoints, "Mesh2Solver");
        for (itP = pointCell.begin(); itP != itPEnd; ++itP, ++idP) {
            update_progress_report(&progress, idP);	
            std::vector<int> neighborOfP;
            std::vector<double> Dr(numOfPoints*2); // matrix row values
            // for each point P, traverse all cells containing it.
            std::vector<int>::iterator itCellEnd = (*itP).end();
            
            for (std::vector<int>::iterator itCell = (*itP).begin(); itCell != itCellEnd; ++itCell) {
                // for each cell containing P, store the point with larger point Id.
                // only three points, don't use for-loop to save time.
                if ( cellPoint[*itCell][0] != idP )
                    neighborOfP.push_back(cellPoint[*itCell][0]);
                if ( cellPoint[*itCell][1] != idP )
                    neighborOfP.push_back(cellPoint[*itCell][1]);
                if ( cellPoint[*itCell][2] != idP )
                    neighborOfP.push_back(cellPoint[*itCell][2]);
            } // for itCell. Ok, now all neighbors of P is stored in neighborOfP;
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
                if (*itQ > idP) {
                    numOfEdges++;
                }
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
                    if (cellPoint[cells[0]][it] != idP && cellPoint[cells[0]][it] != *itQ)
                        itS = cellPoint[cells[0]][it];
                    if (cellPoint[cells[1]][it] != idP && cellPoint[cells[1]][it] != *itQ) 
                        itR = cellPoint[cells[1]][it];
                }
                
                std::vector<double> P(pointXYZ[idP]),
                        Q(pointXYZ[*itQ]),
                        R(pointXYZ[itR]),
                        S(pointXYZ[itS]);
                
                std::vector<double> SP(3), SQ(3), RP(3), RQ(3);
                double SPnorm = 0, SQnorm = 0, RPnorm = 0, RQnorm = 0, SPSQinnerProd = 0, RPRQinnerProd = 0;
                for (int it = 0; it<3; ++it) {
                    SP[it] = P[it] - S[it]; SPnorm += SP[it]*SP[it];
                    SQ[it] = Q[it] - S[it]; SQnorm += SQ[it]*SQ[it]; SPSQinnerProd += SP[it]*SQ[it];
                    RP[it] = P[it] - R[it]; RPnorm += RP[it]*RP[it];
                    RQ[it] = Q[it] - R[it]; RQnorm += RQ[it]*RQ[it]; RPRQinnerProd += RP[it]*RQ[it];
                } //it
                SPnorm = sqrt(SPnorm);
                SQnorm = sqrt(SQnorm);
                RPnorm = sqrt(RPnorm);
                RQnorm = sqrt(RQnorm);
                
                double cosS = SPSQinnerProd / (SPnorm * SQnorm);
                double cosR = RPRQinnerProd / (RPnorm * RQnorm);
                double ctgS = cosS/sqrt(1-cosS*cosS), ctgR = cosR/sqrt(1-cosR*cosR);
                
                Dr[*itQ] = -0.5*(ctgS + ctgR);
                Dr[idP] += 0.5*(ctgS + ctgR);
                // add to the diagonal element of this line.
                
                //D(*itQ, idP) = -0.5*(ctgS + ctgR); // symmetric
                //D(*itQ, *itQ) += 0.5*(ctgS + ctgR);
                // add to the diagonal element of this line.
            } // itQ
            
            // write the real values
            nlRowParameterd(NL_RIGHT_HAND_SIDE,-bR[idP]);
            nlBegin(NL_ROW);
            for (int i = 0; i < numOfPoints; i++) {
                if (Dr[i] != 0) { // store it
                    nlCoefficient(i*2, Dr[i]);
                }
            }
            nlEnd(NL_ROW);
            
            // write the imaginary values
            nlRowParameterd(NL_RIGHT_HAND_SIDE,-bI[idP]);
            nlBegin(NL_ROW);
            for (int i = 0; i < numOfPoints; i++) {
                if (Dr[i] != 0) {
                    nlCoefficient(i*2 + 1, Dr[i]);
                }
            }
            nlEnd(NL_ROW);
        } // itP

        terminate_progress_report(&progress);
        ////////////////////////////////////////////////////////
        // calculate Euler Number to test whether the mesh is genus 0. i.e. Euler Num is 2;
        //    std::cout<<"Total number of edges: "<<numOfEdges<<std::endl;
        int eulerNum = numOfPoints - numOfEdges + numOfFacets;
        
        std::cerr << "  Calculating Euler characteristics......" << std::endl;
        std::cerr << "    Euler Characteristics = " << eulerNum << std::endl;
        std::cerr << "    genus = "<< (2.0 - eulerNum)/2 << std::endl;        
        
    }


    polygons_struct *polygons;
    double mapScale;
};

int main(int argc, char** argv) {
    STRING input_file, output_file;
    object_struct **objects;
    polygons_struct *polygons;
    int n_objects;
    Real *curvatures;
    File_formats format;
    
    /* Get arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
      (void) fprintf(stderr, 
      "\nUsage: %s [options] infile.obj outfile.obj\n\n\
         Generate the Laplace-Beltrami conformal map of the mesh specified in\n\
         infile.obj.  Results are saved in outfile.obj.\n\n",
                     argv[0]);
      (void) fprintf(stderr, 
        "       %s -help\n\n", argv[0]);
        return(1);
    }

    initialize_argument_processing(argc, argv);
    get_string_argument(NULL, &input_file);
    get_string_argument(NULL, &output_file);

    if (input_graphics_any_format(input_file, &format, &n_objects, &objects) != OK) {
        print("Error reading input file\n");
        return(1);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        print("Input file must contain one polygon object.\n");
        return(1);
    }

    polygons = get_polygons_ptr(objects[0]);
    
    double mapScale;
    // guess for mapScale, no optimization
    if (ms == 0) mapScale = round(pow(polygons->n_points,0.5719));
    else mapScale = ms;

    int eulerNum = euler_characteristic(polygons);
    if (eulerNum != 2) {
        print_error("WARNING: Euler characteristics is %d, not 2! Not genus 0 surface...\n", eulerNum);
    }

    LSCM lscm(*polygons, mapScale);
    lscm.apply();

    compute_polygon_normals(polygons);
    output_graphics_any_format(output_file, format, 1, objects);
}


