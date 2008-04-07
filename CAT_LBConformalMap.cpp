
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
        #include  <bicpl.h>
        #include <CAT_Curvature.h>
        #include <CAT_Surf.h>
        #include <CAT_SurfaceIO.h>
}   

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
#include <algorithm>

#include <math.h>
#include <assert.h>


class LSCM {
public:

    LSCM(polygons_struct &p, int ms, STRING ppl) : polygons(&p), mapScale(ms) {
        if (ppl[0] == '-') {
            minFlag = TRUE;
            pointPloc = ppl[1];
        } else if (ppl[0] == '+') {
            minFlag = FALSE;
            pointPloc = ppl[1];
        } else {
            pointPloc = ppl[0];
        }
    }


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
     * copies x + iy to mesh, transforms it into spherical coordinates via stereographic projection
     */
    void solver_to_mesh() {
        const unsigned int numberOfPoints = polygons->n_points;
        std::vector<double> zR(numberOfPoints), zI(numberOfPoints);
        for(unsigned int i=0; i<numberOfPoints; i++) {
            zR[i] = nlGetVariable(i*2);
            zI[i] = nlGetVariable(i*2 + 1);
        }
        
        double xmin =  zR[0];
        double ymin =  zI[0];
        double xmax =  zR[0];
        double ymax =  zI[0];
        
        for (int it = 0; it < numberOfPoints;  ++it) {
            xmin = (xmin<zR[it])?xmin:zR[it];
            xmax = (xmax>zR[it])?xmax:zR[it];
            ymin = (ymin<zI[it])?ymin:zI[it];
            ymax = (ymax>zI[it])?ymax:zI[it];
        } // for it

        double temp1 = ( fabs(xmin)>fabs(xmax) )?fabs(xmin):fabs(xmax);
        double temp2 = ( fabs(ymin)>fabs(ymax) )?fabs(ymin):fabs(ymax);
        //    std::cout<<std::max( temp1, temp2 )<<std::endl;
        double factor = mapScale/( ( temp1>temp2 )?temp1:temp2 );
        
        // the factor is used to re-scale the points in the plane.
        
        
        std::vector<double> x(numberOfPoints), y(numberOfPoints), z(numberOfPoints);
        std::vector<double>::iterator
        itX = x.begin(),
        itY = y.begin(),
        itZ = z.begin(),
        itXend = x.end();
        for (int it = 0; itX != itXend; ++itX, ++itY, ++itZ, ++it) {
            double r2 = factor*zR[it]*factor*zR[it] + factor*zI[it]*factor*zI[it];
            *itX = 2*factor*zR[it]/(1+r2);
            *itY = 2*factor*zI[it]/(1+r2);
            *itZ = 2*r2/(1+r2) - 1;
            polygons->points[it].coords[0] = *itX;
            polygons->points[it].coords[1] = *itY;
            polygons->points[it].coords[2] = *itZ;
        } // for it
        
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
        
        double minval = 1e30;
        double maxval = -1e30;
        int pointP = 0;
        
        for (int it = 0; it < numOfPoints; ++it) {
            double x = polygons->points[it].coords[0];
            double y = polygons->points[it].coords[1];
            double z = polygons->points[it].coords[2];
            
            switch (pointPloc) {
                case 'x':
                    if (minFlag && x < minval) {
                        pointP = it;
                        minval = x;
                    } else if (x > maxval) {
                        pointP = it;
                        maxval = x;
                    }
                case 'y':
                    if (minFlag && y < minval) {
                        pointP = it;
                        minval = y;
                    } else if (y > maxval) {
                        pointP = it;
                        maxval = y;
                    }
                case 'z':
                    if (minFlag && z < minval) {
                        pointP = it;
                        minval = z;
                    } else if (z > maxval) {
                        pointP = it;
                        maxval = z;
                    }
            }
            
            pointXYZ[it][0] = x;
            pointXYZ[it][1] = y;
            pointXYZ[it][2] = z;
        } // for it
        
        
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
                    // If one node has two or less neighbors, it's on the boundary.
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
        std::vector<double> A( pointXYZ[ cellPoint[ pointP ][ 0 ] ] ),
        B( pointXYZ[ cellPoint[ pointP ][ 1 ] ] ),
        C( pointXYZ[ cellPoint[ pointP ][ 2 ] ] );
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
        
        bR[cellPoint[ pointP ][0]] = -1 / ABnorm;
        bR[cellPoint[ pointP ][1]] = 1 / ABnorm;
        
        bI[cellPoint[ pointP ][0]] = (1-theta)/ CEnorm;
        bI[cellPoint[ pointP ][1]] = theta/ CEnorm;
        bI[cellPoint[ pointP ][2]] = -1 / CEnorm;
        
        // 1. Iterate point P from 0 to the last point in the mesh.
        // 2. For each P, find its neighbors, each neighbor must have at least two triangles containing P and itself ---not the boundary.
        // 3. For each of P's neighbors, Q, calculate R, S
        // 4. Write the value in matrix.
        
        std::vector< std::vector<int> >::iterator itP, itPEnd = pointCell.end();
        int idP = 0;
        unsigned long numOfEdges = 0;
	    progress_struct      progress;
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
                
                Dr[*itQ]= -0.5*(ctgS + ctgR);
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
    char pointPloc;
    bool minFlag;
    double mapScale;
};


void usage(STRING executable) {
    STRING usage_str = "\n\
Usage: %s infile.obj outfile.obj [mapScale] [pointPlocation]\n\n\
     Generate the Laplace-Beltrami conformal map of the mesh specified in\n\
     infile.obj.  Results are saved in outfile.obj.\n\n\
     mapScale = scaling factor [default: 100]\n\
     pointPlocation = +x,-x,+y,-y,+z,-z,0 [default: 0]\n";

     print_error(usage_str, executable);
}

int main(int argc, char** argv) {
    STRING ifname, ofname;
    FILE *file;
    object_struct **objects;
    polygons_struct *polygons;
    int n_objects;
    Real *curvatures;
    File_formats format;
    STRING pointPloc;
    int mapScale;
    
    initialize_argument_processing(argc, argv);
    if (!get_string_argument(NULL, &ifname) ||
                !get_string_argument(NULL, &ofname) ) {
        usage(argv[0]);
        return(1);
    }

    get_int_argument(100, &mapScale);
    get_string_argument("0", &pointPloc);

    if (input_graphics_any_format(ifname, &format, &n_objects, &objects) != OK) {
        print("Error reading input file\n");
        return(1);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        print("Input file must contain one polygon object.\n");
        return(1);
    }

    polygons = get_polygons_ptr(objects[0]);
    
    int eulerNum = euler_characteristic(polygons);
    if (eulerNum != 2) {
        print_error("    Euler characteristics is %d, not 2! Not genus 0 surface, exiting...\n", eulerNum);
        exit(1);
    }

    LSCM lscm(*polygons, mapScale, pointPloc);
    lscm.apply();

    compute_polygon_normals(polygons);
    output_graphics_any_format(ofname, format, 1, objects);
}

