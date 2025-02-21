/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "nifti/nifti1_io.h"

#define KAPPA   10

double 
getPixel(double *vol, int *dims, double x0, double y0, double z0)
/* trilinear intensity interpolation */
{
    double value;
    
    int x = (int)x0;
    int y = (int)y0;
    int z = (int)z0;
    
    if (x < 1 || x >= dims[0] - 1)
        return 0.0;
    if (y < 1 || y >= dims[1] - 1)
        return 0.0;
    if (z < 1 || z >= dims[2] - 1)
        return 0.0;

    double dx = x0 - x;
    double dy = y0 - y;
    double dz = z0 - z;
    double xm = 1 - dx;
    double ym = 1 - dy;
    double zm = 1 - dz;
    
    value =  (xm * ym * zm) * (vol[((z+0)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0+1]);
    value += (dx * ym * zm) * (vol[((z+0)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+1+1]);
    value += (xm * dy * zm) * (vol[((z+0)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+0+1]);
    value += (dx * dy * zm) * (vol[((z+0)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+1+1]);
    value += (xm * ym * dz) * (vol[((z+1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0+1]);
    value += (dx * ym * dz) * (vol[((z+1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+1+1]);
    value += (xm * dy * dz) * (vol[((z+1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+0+1]);
    value += (dx * dy * dz) * (vol[((z+1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+1+1]);

    return value;
}


void 
getForce(double *vol, double grad[3], int *dims, double *vx, double x0, double y0, double z0)
{
    
    int x = (int)x0;
    int y = (int)y0;
    int z = (int)z0;
    
    /* set gradient to zero if dimensions exceed in the next step */
    grad[0] = grad[1] = grad[2] = 0.0;
    
    if (x < 1 || x >= dims[0] - 1)
        return;
    if (y < 1 || y >= dims[1] - 1)
        return;
    if (z < 1 || z >= dims[2] - 1)
        return;

    double dx = x0 - x;
    double dy = y0 - y;
    double dz = z0 - z;
    double xm = 1 - dx;
    double ym = 1 - dy;
    double zm = 1 - dz;
    
    grad[0] =    (xm * ym * zm) * (vol[((z+0)*dims[1]*dims[0])   + ((y+0)*dims[0]) + x+0+1] - vol[((z+0)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0-1]) / (2*vx[0]);
    grad[1] =    (xm * ym * zm) * (vol[((z+0)*dims[1]*dims[0])   + ((y+0+1)*dims[0]) + x+0] - vol[((z+0)*dims[1]*dims[0]) + ((y+0-1)*dims[0]) + x+0]) / (2*vx[1]);
    grad[2] =    (xm * ym * zm) * (vol[((z+0+1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0] - vol[((z+0-1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0])   / (2*vx[2]);

    grad[0] += (dx * ym * zm) * (vol[((z+0)*dims[1]*dims[0])     + ((y+0)*dims[0]) + x+1+1] - vol[((z+0)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+1-1]) / (2*vx[0]);
    grad[1] += (dx * ym * zm) * (vol[((z+0)*dims[1]*dims[0])     + ((y+0+1)*dims[0]) + x+1] - vol[((z+0)*dims[1]*dims[0]) + ((y+0-1)*dims[0]) + x+1]) / (2*vx[1]);
    grad[2] += (dx * ym * zm) * (vol[((z+0+1)*dims[1]*dims[0])   + ((y+0)*dims[0]) + x+1] - vol[((z+0-1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+1])   / (2*vx[2]);

    grad[0] += (xm * dy * zm) * (vol[((z+0)*dims[1]*dims[0])     + ((y+1)*dims[0]) + x+0+1] - vol[((z+0)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+0-1]) / (2*vx[0]);
    grad[1] += (xm * dy * zm) * (vol[((z+0)*dims[1]*dims[0])     + ((y+1+1)*dims[0]) + x+0] - vol[((z+0)*dims[1]*dims[0]) + ((y+1-1)*dims[0]) + x+0]) / (2*vx[1]);
    grad[2] += (xm * dy * zm) * (vol[((z+0+1)*dims[1]*dims[0])   + ((y+1)*dims[0]) + x+0] - vol[((z+0-1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+0])   / (2*vx[2]);

    grad[0] += (dx * dy * zm) * (vol[((z+0)*dims[1]*dims[0])     + ((y+1)*dims[0]) + x+1+1] - vol[((z+0)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+1-1]) / (2*vx[0]);
    grad[1] += (dx * dy * zm) * (vol[((z+0)*dims[1]*dims[0])     + ((y+1+1)*dims[0]) + x+1] - vol[((z+0)*dims[1]*dims[0]) + ((y+1-1)*dims[0]) + x+1]) / (2*vx[1]);
    grad[2] += (dx * dy * zm) * (vol[((z+0+1)*dims[1]*dims[0])   + ((y+1)*dims[0]) + x+1] - vol[((z+0-1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+1])   / (2*vx[2]);

    grad[0] += (xm * ym * dz) * (vol[((z+1)*dims[1]*dims[0])     + ((y+0)*dims[0]) + x+0+1] - vol[((z+1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0-1]) / (2*vx[0]);
    grad[1] += (xm * ym * dz) * (vol[((z+1)*dims[1]*dims[0])     + ((y+0+1)*dims[0]) + x+0] - vol[((z+1)*dims[1]*dims[0]) + ((y+0-1)*dims[0]) + x+0]) / (2*vx[1]);
    grad[2] += (xm * ym * dz) * (vol[((z+1+1)*dims[1]*dims[0])   + ((y+0)*dims[0]) + x+0] - vol[((z+1-1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+0])   / (2*vx[2]);

    grad[0] += (dx * ym * dz) * (vol[((z+1)*dims[1]*dims[0])     + ((y+0)*dims[0]) + x+1+1] - vol[((z+1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+1-1]) / (2*vx[0]);
    grad[1] += (dx * ym * dz) * (vol[((z+1)*dims[1]*dims[0])     + ((y+0+1)*dims[0]) + x+1] - vol[((z+1)*dims[1]*dims[0]) + ((y+0-1)*dims[0]) + x+1]) / (2*vx[1]);
    grad[2] += (dx * ym * dz) * (vol[((z+1+1)*dims[1]*dims[0])   + ((y+0)*dims[0]) + x+1] - vol[((z+1-1)*dims[1]*dims[0]) + ((y+0)*dims[0]) + x+1])   / (2*vx[2]);

    grad[0] += (xm * dy * dz) * (vol[((z+0)*dims[1]*dims[0])     + ((y+1)*dims[0]) + x+0+1] - vol[((z+1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+0-1]) / (2*vx[0]);
    grad[1] += (xm * dy * dz) * (vol[((z+0)*dims[1]*dims[0])     + ((y+1+1)*dims[0]) + x+0] - vol[((z+1)*dims[1]*dims[0]) + ((y+1-1)*dims[0]) + x+0]) / (2*vx[1]);
    grad[2] += (xm * dy * dz) * (vol[((z+0+1)*dims[1]*dims[0])   + ((y+1)*dims[0]) + x+0] - vol[((z+1-1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+0])   / (2*vx[2]);

    grad[0] += (dx * dy * dz) * (vol[((z+1)*dims[1]*dims[0])     + ((y+1)*dims[0]) + x+1+1] - vol[((z+1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+1-1]) / (2*vx[0]);
    grad[1] += (dx * dy * dz) * (vol[((z+1)*dims[1]*dims[0])     + ((y+1+1)*dims[0]) + x+1] - vol[((z+1)*dims[1]*dims[0]) + ((y+1-1)*dims[0]) + x+1]) / (2*vx[1]);
    grad[2] += (dx * dy * dz) * (vol[((z+1+1)*dims[1]*dims[0])   + ((y+1)*dims[0]) + x+1] - vol[((z+1-1)*dims[1]*dims[0]) + ((y+1)*dims[0]) + x+1])   / (2*vx[2]);

}

void
usage(char *executable)
{
        fprintf(stderr, "%s volume_file surface_file output_surface_file\n", executable);
}

int
main(int argc, char *argv[])
{
    char              *volume_file = NULL;
    char              *input_file = NULL, *output_surface_file = NULL;
    File_formats      file_format;
    int               i, j, v, dims[3], n_objects, it, pidx;
    int               *n_neighbours, **neighbours;
    object_struct     **object_list;
    polygons_struct   *polygons;
    nifti_image       *nii_ptr;
    double            *input, vx[3], lim;
    Point             points[MAX_POINTS_PER_POLYGON];
    Vector            c, n, p;

    initialize_argument_processing(argc, argv);

    /* read first image to get image parameters */
    nii_ptr = read_nifti_double(argv[1], &input, 0);
    if(nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", volume_file);
        return(EXIT_FAILURE);
    }

    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;

    vx[0] = nii_ptr->dx;
    vx[1] = nii_ptr->dy;
    vx[2] = nii_ptr->dz;

    if (input_graphics_any_format(input_file, &file_format,
              &n_objects, &object_list) == ERROR ||
              n_objects != 1 || object_list[0]->object_type != POLYGONS) {
        fprintf(stderr, "File must contain 1 polygons struct.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);    
    compute_polygon_normals(polygons);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                    &neighbours, NULL, NULL);
    
    it = 20;
    lim = 0.5;
    
    double s = 0.0; int nv = 0;
    for (i = 0; i < it; i++) {
        s = 0.0; nv = 0;
        for (v = 0; v < polygons->n_points; v++) {
            fill_Vector(c, 0.0, 0.0, 0.0);
            for (j = 0; j < n_neighbours[v]; j++) {
                pidx = neighbours[v][j];
                //double x = Point_x(pidx);
            }

/*
Vector p = polygons->points[v]);
Vector n = polygons->normals[v]);
Point_x(polygons_out->points[p]) += extents[p]*length*thickness_values[p]*Point_x(polygons->normals[p]);
x = Point_x(polygons->points[p]) + weight*curvatures[p]*Point_x(polygons->normals[p]);
                        point = polygons->points[p];
                        Point_x(point) += weight*curvatures[p]*Point_x(polygons->normals[p]);
        for (i = 0; i < polygons->n_points; i++) {
                fill_Vector(xyz, Point_x(polygons->points[i]),
                                 Point_y(polygons->points[i]),
                                 Point_z(polygons->points[i]));
                radius[i] = MAGNITUDE(xyz);
        }

        for (i = 0; i < polygons->n_points; i++) {
                dist = 0.0;
                for (j = 0; j < 3; j++)
                        dist += Point_coord(polygons->points[i], j) *
                                Point_coord(polygons->points[i], j);
                radius = MAX(radius, dist);        
        }

*/
/*
      for (nvertex * v = (nvertex *) VGraphFirstNode (vertices); v;
         v = (nvertex *) VGraphNextNode (vertices)) {
    
      // compute internal force
      vec3 c = v->center (vertices);
      vec3 p = v->getPoint ();
      vec3 n = v->getNormal ();
      
      vec3 pv = p / dim;
    
      // compute external forces
      double di = getPixel(src, pv) - lim;
      double f3 = tanh (di / ll.);
      vec3 f = getForce (grad, pv);
      double f2 = f.dot (n);
double f2 = DOT_VECTORS(f, n);

      getForce(input, grad, dims, vx, double x0, double y0, double z0)

      vec3 d = w1 * (c - p) + ((w2 * f2 - w3) * f3) * n;
    
      // move vertex
      p += d;
      v->setPoint (p);
      s += di*di; nv++;
*/
        }  
        compute_polygon_normals(polygons);
        fprintf(stdout, "mesh::deform: iter %03d %8.2lf", i, sqrt(s/(nv-1)));
    }

    if (output_graphics_any_format(output_surface_file, ASCII_FORMAT,
                     n_objects, object_list, NULL) == ERROR)
        exit(EXIT_FAILURE);

    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);

    return(EXIT_SUCCESS);
}
