/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#define BINTREE_FACTOR 0.5

#define HOLE 1
#define HANDLE 2

Vector defect_direction(polygons_struct *, int *, int);
int defect_euler(polygons_struct *, int *, int *, int, int *, int **);

int isedge(polygons_struct *, int *, int *, int *, int **, int);
int split_defects(polygons_struct *, int *, int *, int *, int **, int);
int minimize_defect(polygons_struct *, int *, int *, int *, int **, int);

int find_topological_defects(polygons_struct *, polygons_struct *, int *,
                int *, int **);
int find_artifacts(polygons_struct *, polygons_struct *,
                int *, int *, int **, double );
void expand_defects(polygons_struct *, int *, int *, int, int, int *, int **);
void update_defects(polygons_struct *, int *, int *);
void update_polydefects(polygons_struct *, int *, int *);
Point get_defect_center(polygons_struct *, int *, int);

void get_defect_size(polygons_struct *, int *, int, double *);

void bisect_defects(polygons_struct *, polygons_struct *,int *, int, int *, int *, int);

void remap_defect(polygons_struct *, int *, int *, polygons_struct *, int *,
                  int *);

void inflate_surface_with_topology_defects(polygons_struct *);
