/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/geometry.h>

#include "CAT_SurfMatrix.h"
#include "CAT_MarchingCubes.h"

void
resel_fun(int *curr, int *prev, /* current and previous planes */
          int m, int n, /* image dimensions */
          int *P,   /* # points */
          int E[3], /* # edges  */
          int F[4], /* # faces  */
          int *C)   /* # cubes  */
{
        int p=0, ex=0, ey=0, ez=0, fxy=0, fxz=0, fyz=0, c=0;
        int i, j, o;

        for (i = 1; i < m; i++) {
                for (j = 1; j < n; j++) {
                        o = i + j*m;

                        if (curr[o]) {
                                p++;

                                /* A simple way of doing it... 
                                if (curr[o-1]) ex++;
                                if (curr[o-m]) ey++;
                                if (prev[o  ]) ez++;
                                if (curr[o-1] && curr[o-m] && curr[o-m-1])fxy++;
                                if (curr[o-1] && prev[o  ] && prev[o  -1])fxz++;
                                if (curr[o-m] && prev[o  ] && prev[o  -m])fyz++;
                                if (curr[o-1] && curr[o-m] && curr[o-m-1] &&
                                    prev[o] && prev[o-1] && prev[o-m] &&
                                    prev[o-m-1]) c++;
                                */

                                /* It could be optimized much further with
                                 * "if else" statements but it was beginning
                                 * to hurt my head */
                                if (curr[o-1]) {
                                        ex++;
                                        if (curr[o-m] && curr[o-m-1]) {
                                                fxy++;
                                                if (prev[o] && prev[o-1] &&
                                                    prev[o-m] &&
                                                    prev[o-m-1]) c++;
                                        }
                                }
                                if (curr[o-m]) {
                                        ey++;
                                        if (prev[o  ] && prev[o  -m]) fyz++;
                                }
                                if (prev[o  ]) {
                                        ez++;
                                        if (curr[o-1] && prev[o  -1]) fxz++;
                                }
                        }
                }
        }
        *P   += p;
        E[0] += ex;
        E[1] += ey;
        E[2] += ez;
        F[0] += fxy;
        F[1] += fxz;
        F[2] += fyz;
        *C   += c;
}

int
volume_matrix_euler(struct voxmatrix *vmat, int val)
{
        int v, x, y, z, E[3], F[3], P, C;
        double R[4], r[3];
        int *curr, *prev, *tmpp; /* planes */

        curr = (int *) calloc(vmat->vsizes[0]*vmat->vsizes[1],sizeof(int));
        prev = (int *) calloc(vmat->vsizes[0]*vmat->vsizes[1],sizeof(int));

        P = C = E[0] = E[1] = E[2] = F[0] = F[1] = F[2] = 0;
        for (z = 0; z < vmat->vsizes[2]; z++) {
                int i1, j1;
                for (x = 0; x < vmat->vsizes[0]; x++) {
                        for (y = 0; y < vmat->vsizes[1]; y++) {
                                v = xyz2mat(vmat, x, y, z);
                                if (vmat->voxels[v] == val)
                                        curr[x + y*vmat->vsizes[0]] = 1;
                                else
                                        curr[x + y*vmat->vsizes[0]] = 0;
                        }
                }

                /* count edges, faces etc */
                resel_fun(curr, prev, vmat->vsizes[0], vmat->vsizes[1],
                          &P, E, F, &C);

                /* make current plane previous */
                tmpp = prev; prev = curr; curr = tmpp;
        }

        free(curr);
        free(prev);

        return(2*(P - (E[0]+E[1]+E[2])+(F[0]+F[1]+F[2])-C));
}

int
volume_euler(Volume volume, int righthemi)
{
        struct voxmatrix *vmat;
        int euler;

        vmat = make_matrix_from_volume(volume, righthemi);
        euler = volume_matrix_euler(vmat, INSIDE);

        delete_voxmatrix(vmat);
        return(euler);
}

/* create a surface from a volume, then calculate the euler */
int
volume_euler_surface(char *filename, Volume volume, int righthemi)
{
        double min, max;
        int euler;
        object_struct **objects;
        polygons_struct *surface;

        if (righthemi) {
                min = 126; max = 128;
        } else {
                min = 254; max = 256;
        }

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = marching_cubes(filename, volume, min, max);

        surface = get_polygons_ptr(*objects);
        output_graphics_any_format("test.obj", ASCII_FORMAT,
                                           1, objects);

        euler = surface_euler(surface);

        delete_object_list(1, objects);
        return(euler);
}

int
surface_euler(polygons_struct *surface)
{
        int n_edges, *n_neighbours, **neighbours;

        create_polygon_point_neighbours(surface, FALSE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        n_edges = count_edges(surface, n_neighbours, neighbours);

        return(surface->n_items + surface->n_points - n_edges);
}
