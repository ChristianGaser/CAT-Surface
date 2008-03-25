/*************************************************************************************/
/*
 * OpenNL example:
 * Least Squares Conformal Maps implemented with OpenNL
 *
 *  OpenNL: Numerical Library
 *  Copyright (C) 2004 Bruno Levy
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
/*************************************************************************************/


#include "nl.h"

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>

#include <math.h>
#include <assert.h>

/*************************************************************************************/
/* Basic functions */

template <class T> inline T nl_min(T x, T y) { return x < y ? x : y ; }
template <class T> inline T nl_max(T x, T y) { return x > y ? x : y ; }

/*************************************************************************************/
/* Basic geometric types */

class Vector2 {
public:
    Vector2(double x_in, double y_in) : x(x_in), y(y_in) { }
    Vector2() : x(0), y(0) { }
    double x ;
    double y ;
} ;

class Vector3 {
public:
    Vector3(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) { }
    Vector3() : x(0), y(0), z(0) { }
    double length() const { 
        return sqrt(x*x + y*y + z*z) ;
    }
    void normalize() {
        double l = length() ;
        x /= l ; y /= l ; z /= l ;
    }
    double x ;
    double y ;
    double z ;
} ;

Vector3 inverseStereoProjection(double x, double y, double m_fRadius)
{
    double fRadiusSqr = m_fRadius * m_fRadius;
	double fRSqr = x*x + y*y;
	double fMult = 1.0/(fRSqr + fRadiusSqr);
	double fX = 2.0*fMult*fRadiusSqr*x;
    double fY = 2.0*fMult*fRadiusSqr*y;
    double fZ = fMult*m_fRadius*(fRSqr-fRadiusSqr);

	return Vector3(fX/m_fRadius, fY/m_fRadius, fZ/m_fRadius);
}

Vector3 inverseSphericalProjection(double theta, double phi, double r)
{
	// theta [0..pi]
	// phi [0..2pi]
	double x = r*sin(theta)*cos(phi);
	double y = r*sin(theta)*sin(phi);
	double z = r*cos(theta);

	return Vector3(x,y,z);
}

// I/O

std::ostream& operator<<(std::ostream& out, const Vector2& v) {
    return out << v.x << " " << v.y ;
}

std::ostream& operator<<(std::ostream& out, const Vector3& v) {
    return out << v.x << " " << v.y << " " << v.z ;
}

std::istream& operator>>(std::istream& in, Vector2& v) {
    return in >> v.x >> v.y ;
}

std::istream& operator>>(std::istream& in, Vector3& v) {
    return in >> v.x >> v.y >> v.z ;
}

/* dot product */
double operator*(const Vector3& v1, const Vector3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z ;
}

/* cross product */
Vector3 operator^(const Vector3& v1, const Vector3& v2) {
    return Vector3(
        v1.y*v2.z - v2.y*v1.z,
        v1.z*v2.x - v2.z*v1.x,
        v1.x*v2.y - v2.x*v1.y
    ) ;
}

Vector3 operator+(const Vector3& v1, const Vector3& v2) {
    return Vector3(
        v1.x + v2.x,
        v1.y + v2.y,
        v1.z + v2.z
    ) ;
}

Vector3 operator-(const Vector3& v1, const Vector3& v2) {
    return Vector3(
        v1.x - v2.x,
        v1.y - v2.y,
        v1.z - v2.z
    ) ;
}

Vector2 operator+(const Vector2& v1, const Vector2& v2) {
    return Vector2(
        v1.x + v2.x,
        v1.y + v2.y
    ) ;
}

Vector2 operator-(const Vector2& v1, const Vector2& v2) {
    return Vector2(
        v1.x - v2.x,
        v1.y - v2.y
    ) ;
}

/***********************************************************************************/
/* Mesh class */

// Note1: this is a minimum mesh class, it does not have facet adjacency information
// (we do not need it for LSCM). This is just like an Inventor indexed face set.
//

class Vertex {
public:
    Vertex() : locked(false), id(-1) { }
    Vertex(
        const Vector3& p, const Vector2& t
    ) : point(p), tex_coord(t), locked(false), id(-1) { }
    Vector3 point ;
    Vector2 tex_coord ;
    bool    locked ;
    int     id ;
} ;

class Facet : public std::vector<int> {
public:
} ;

class IndexedMesh {
public:
    IndexedMesh() : in_facet(false) { }

    Vertex* add_vertex() {
        vertex.push_back(Vertex()) ;
        vertex.rbegin()->id = vertex.size() - 1 ;
        return &*(vertex.rbegin()) ;
    }

    Vertex* add_vertex(const Vector3& p, const Vector2& t) {
        vertex.push_back(Vertex(p,t)) ;
        vertex.rbegin()->id = vertex.size() - 1 ;
        return &*(vertex.rbegin()) ;
    }

    void begin_facet() {
        assert(!in_facet) ;
        facet.push_back(Facet()) ;
        in_facet = true ;
    }

    void end_facet() {
        assert(in_facet) ;
        in_facet = false ;
    }

    void add_vertex_to_facet(unsigned int i) {
        assert(in_facet) ;
        assert(i < vertex.size()) ;
        facet.rbegin()->push_back(i) ;
    }

    void clear() {
        vertex.clear() ;
        facet.clear() ;
    }

    void load_bicobj(const std::string& file_name) {
        std::ifstream input(file_name.c_str()) ;
        char            t, line[1024];
        std::string     s;
        unsigned int    i, n_items, n_points, c;
        double          dummy;
        Vector3         p;
        
        clear() ;
        input >> t >> dummy >> dummy >> dummy >> dummy >> dummy >> n_points;
        if(! t == 'P') {
            std::cerr << "Wrong BIC obj format" << std::endl ;
            return;
        }
        // blank line
        input.getline(line, 1024) ;
        // read points
        for (i=0; i<n_points; ++i) {
            input.getline(line, 1024) ;
            std::istrstream line_input(line) ;                
            line_input >> p;
            add_vertex(p,Vector2(0,0));
        }
        // blank line
        input.getline(line, 1024) ;
        // skip normals
        for (i=0; i<n_points; ++i) {
            input.getline(line, 1024) ;
        }
        // blank line
        input.getline(line, 1024) ;
        // # of polygons
        input >> n_items;
        // number of colors
        input.getline(line, 1024) ;
        input >> c;
        if(! c == 0) {
            std::cerr << "Only single color allowed." << std::endl ;
            return;
        }
        input.getline(line, 1024) ;
        for (i=0; i<n_items; ++i) {
            int index ;
            input >> index ;
        }
        input.getline(line, 1024) ;
        
       int index ;
       for (i=0; i<n_items; ++i) {
            begin_facet() ;
            for (int j=0; j<3; ++j) {
                input >> index ;
                if (index > 0) {
                    add_vertex_to_facet(index) ;
                } else {
                    add_vertex_to_facet(0) ;
                }
            }
            end_facet() ;
        }

        std::cerr << "Loaded " << vertex.size() << " vertices and " 
                  << facet.size() << " facets" << std::endl ;
    }

    void save_bicobj(const std::string& file_name) {
        unsigned int i,j ;
        std::ofstream out(file_name.c_str()) ;
 
        out << "P 0.3 0.3 0.4 10 1 " << vertex.size() << std::endl;
        // vertices
        for(i=0; i<vertex.size(); i++) {
            out << inverseStereoProjection(vertex[i].tex_coord.x,vertex[i].tex_coord.y,1.0) << std::endl;
        }
        out << std::endl;
        // normals
        for(i=0; i<vertex.size(); i++) {
            out << "1.0 1.0 1.0" << std::endl ;
        }
        out << std::endl;
        out << facet.size() << std::endl;
        out << "0 1 1 1 1" << std::endl;
        out << std::endl;
        int k = 1;
        for(i=1; i<(facet.size()+1); i++) {
            out << i*3 << " ";
            if (k==8) {
                out << std::endl;
                k = 0;
            }
            k++;
        }

        out << std::endl;
        for(i=0; i<facet.size(); i++) {
            const Facet& F = facet[i] ;
            for(j=0; j<F.size(); j++) {
                out << (F[j] + 0) << " " ;
            }
            out << std::endl ;
        }

    }
    std::vector<Vertex> vertex ;
    std::vector<Facet>  facet ;
    bool in_facet ;
    
} ;


class LSCM {
public:

    LSCM(IndexedMesh& m) : mesh_(&m) {
    }


    // Outline of the algorithm:

    // 1) Find an initial solution by projecting on a plane
    // 2) Lock two vertices of the mesh
    // 3) Copy the initial u,v coordinates to OpenNL
    // 3) Construct the LSCM equation with OpenNL
    // 4) Solve the equation with OpenNL
    // 5) Copy OpenNL solution to the u,v coordinates


    void apply() {
        int nb_vertices = mesh_->vertex.size() ;
        project() ;
        nlNewContext() ;
        if(nlInitExtension("SUPERLU")) {
            nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT) ;
            std::cerr << "Using SuperLU, cool !" << std::endl ;
        } else {
            nlSolverParameteri(NL_SOLVER, NL_CG) ;
            nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_JACOBI) ;
            std::cerr << "Using Jacobi pre-conditioned conjugate gradient" 
                      << std::endl ;
        }
        nlSolverParameteri(NL_NB_VARIABLES, 2*nb_vertices) ;
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE) ;
        nlSolverParameteri(NL_MAX_ITERATIONS, 5*nb_vertices) ;
        nlSolverParameterd(NL_THRESHOLD, 1e-10) ;
        nlBegin(NL_SYSTEM) ;
        mesh_to_solver() ;
        nlBegin(NL_MATRIX) ;
        setup_lscm() ;
        nlEnd(NL_MATRIX) ;
        nlEnd(NL_SYSTEM) ;
        std::cerr << "Solving ..." << std::endl ;
        nlSolve() ;
        solver_to_mesh() ;
        double time ;
        nlGetDoublev(NL_ELAPSED_TIME, &time) ;
        std::cerr << "Solver time: " << time << std::endl ;
        nlDeleteContext(nlGetCurrent()) ;
    }

protected:

    void setup_lscm() {
        for(unsigned int i=0; i<mesh_->facet.size(); i++) {
            const Facet& F = mesh_->facet[i] ;
            setup_lscm(F) ;
        }
    }

    // Note: no-need to triangulate the facet,
    // we can do that "virtually", by creating triangles
    // radiating around vertex 0 of the facet.
    // (however, this may be invalid for concave facets)
    void setup_lscm(const Facet& F) {
        for(unsigned int i=1; i < F.size() - 1; i++) {
            setup_conformal_map_relations(
                mesh_->vertex[F[0]], mesh_->vertex[F[i]], mesh_->vertex[F[i+1]]
            ) ;
        }
    }

    // Computes the coordinates of the vertices of a triangle
    // in a local 2D orthonormal basis of the triangle's plane.
    static void project_triangle(
        const Vector3& p0, 
        const Vector3& p1, 
        const Vector3& p2,
        Vector2& z0,
        Vector2& z1,
        Vector2& z2
    ) {
        Vector3 X = p1 - p0 ;
        X.normalize() ;
        Vector3 Z = X ^ (p2 - p0) ;
        Z.normalize() ;
        Vector3 Y = Z ^ X ;
        const Vector3& O = p0 ;
        
        double x0 = 0 ;
        double y0 = 0 ;
        double x1 = (p1 - O).length() ;
        double y1 = 0 ;
        double x2 = (p2 - O) * X ;
        double y2 = (p2 - O) * Y ;        
            
        z0 = Vector2(x0,y0) ;
        z1 = Vector2(x1,y1) ;
        z2 = Vector2(x2,y2) ;        
    }

    // LSCM equation, geometric form :
    // (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
    // Where Uk = uk + i.vk is the complex number 
    //                       corresponding to (u,v) coords
    //       Zk = xk + i.yk is the complex number 
    //                       corresponding to local (x,y) coords
    // cool: no divide with this expression,
    //  makes it more numerically stable in
    //  the presence of degenerate triangles.
    
    void setup_conformal_map_relations(
        const Vertex& v0, const Vertex& v1, const Vertex& v2
    ) {
        
        int id0 = v0.id ;
        int id1 = v1.id ;
        int id2 = v2.id ;
            
        const Vector3& p0 = v0.point ; 
        const Vector3& p1 = v1.point ;
        const Vector3& p2 = v2.point ;
            
        Vector2 z0,z1,z2 ;
        project_triangle(p0,p1,p2,z0,z1,z2) ;
        Vector2 z01 = z1 - z0 ;
        Vector2 z02 = z2 - z0 ;
        double a = z01.x ;
        double b = z01.y ;
        double c = z02.x ;
        double d = z02.y ;
        assert(b == 0.0) ;

        // Note  : 2*id + 0 --> u
        //         2*id + 1 --> v
        int u0_id = 2*id0     ;
        int v0_id = 2*id0 + 1 ;
        int u1_id = 2*id1     ;
        int v1_id = 2*id1 + 1 ;
        int u2_id = 2*id2     ;
        int v2_id = 2*id2 + 1 ;
        
        // Note : b = 0

        // Real part
        nlBegin(NL_ROW) ;
        nlCoefficient(u0_id, -a+c)  ;
        nlCoefficient(v0_id,  b-d)  ;
        nlCoefficient(u1_id,   -c)  ;
        nlCoefficient(v1_id,    d)  ;
        nlCoefficient(u2_id,    a) ;
        nlEnd(NL_ROW) ;

        // Imaginary part
        nlBegin(NL_ROW) ;
        nlCoefficient(u0_id, -b+d) ;
        nlCoefficient(v0_id, -a+c) ;
        nlCoefficient(u1_id,   -d) ;
        nlCoefficient(v1_id,   -c) ;
        nlCoefficient(v2_id,    a) ;
        nlEnd(NL_ROW) ;
    }

    /**
     * copies u,v coordinates from OpenNL solver to the mesh.
     */
    void solver_to_mesh() {
        for(unsigned int i=0; i<mesh_->vertex.size(); i++) {
            Vertex& it = mesh_->vertex[i] ;
            double u = nlGetVariable(2 * it.id    ) ;
            double v = nlGetVariable(2 * it.id + 1) ;
            it.tex_coord = Vector2(u,v) ;
        }
    }
        
    /**
     * copies u,v coordinates from the mesh to OpenNL solver.
     */
    void mesh_to_solver() {
        for(unsigned int i=0; i<mesh_->vertex.size(); i++) {
            Vertex& it = mesh_->vertex[i] ;
            double u = it.tex_coord.x ;
            double v = it.tex_coord.y ;
            nlSetVariable(2 * it.id    , u) ;
            nlSetVariable(2 * it.id + 1, v) ;
            if(it.locked) {
                nlLockVariable(2 * it.id    ) ;
                nlLockVariable(2 * it.id + 1) ;
            } 
        }
    }


    // Chooses an initial solution, and locks two vertices
    void project() {
        // Get bbox
        unsigned int i ;

        double xmin =  1e30 ;
        double ymin =  1e30 ;
        double zmin =  1e30 ;
        double xmax = -1e30 ;
        double ymax = -1e30 ;
        double zmax = -1e30 ;

        for(i=0; i<mesh_->vertex.size(); i++) {
            const Vertex& v = mesh_->vertex[i] ;
            xmin = nl_min(v.point.x, xmin) ;
            ymin = nl_min(v.point.y, xmin) ;
            zmin = nl_min(v.point.z, xmin) ;

            xmax = nl_max(v.point.x, xmin) ;
            ymax = nl_max(v.point.y, xmin) ;
            zmax = nl_max(v.point.z, xmin) ;
        }

        double dx = xmax - xmin ;
        double dy = ymax - ymin ;
        double dz = zmax - zmin ;
        
        Vector3 V1,V2 ;

        // Find shortest bbox axis
        if(dx < dy && dx < dz) {
            if(dy > dz) {
                V1 = Vector3(0,1,0) ;
                V2 = Vector3(0,0,1) ;
            } else {
                V2 = Vector3(0,1,0) ;
                V1 = Vector3(0,0,1) ;
            }
        } else if(dy < dx && dy < dz) {
            if(dx > dz) {
                V1 = Vector3(1,0,0) ;
                V2 = Vector3(0,0,1) ;
            } else {
                V2 = Vector3(1,0,0) ;
                V1 = Vector3(0,0,1) ;
            }
        } else if(dz < dx && dz < dy) {
            if(dx > dy) {
                V1 = Vector3(1,0,0) ;
                V2 = Vector3(0,1,0) ;
            } else {
                V2 = Vector3(1,0,0) ;
                V1 = Vector3(0,1,0) ;
            }
        }

        // Project onto shortest bbox axis,
        // and lock extrema vertices

        Vertex* vxmin = NULL ;
        double  umin = 1e30 ;
        Vertex* vxmax = NULL ;
        double  umax = -1e30 ;

        for(i=0; i<mesh_->vertex.size(); i++) {
            Vertex& V = mesh_->vertex[i] ;
            double u = V.point * V1 ;
            double v = V.point * V2 ;
            V.tex_coord = Vector2(u,v) ;
            if(u < umin) {
                vxmin = &V ;
                umin = u ;
            } 
            if(u > umax) {
                vxmax = &V ;
                umax = u ;
            } 
        }

        vxmin->locked = true ;
        vxmax->locked = true ;
    }

    IndexedMesh* mesh_ ;
} ;

int main(int argc, char** argv) {
    if(argc != 3) {
        std::cerr << "usage: lscm_bicobj infile.obj outfile.obj" << std::endl ;
        return -1 ;
    }

    IndexedMesh mesh ;
    std::cerr << "Loading ..." << std::endl ;
    mesh.load_bicobj(argv[1]) ;

    LSCM lscm(mesh) ;
    lscm.apply() ;

    std::cerr << "Saving ..." << std::endl ;
    mesh.save_bicobj(argv[2]) ;
}
