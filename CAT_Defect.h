/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>

#define BINTREE_FACTOR 0.5

#define HOLE 1
#define HANDLE 2
#define LARGE_DEFECT 3
#define VENTRICLE 4


Vector defect_direction(polygons_struct *, int *, int);
int defect_euler(polygons_struct *, int *, int *, int, int *, int **);

int isedge(polygons_struct *, int *, int *, int *, int **, int);
int split_defects(polygons_struct *, int *, int *, int *, int **, int);
int minimize_defect(polygons_struct *, int *, int *, int *, int **, int);

int find_topological_defects(polygons_struct *, polygons_struct *, int *,
                             int *, int **);

void expand_defects(polygons_struct *, int *, int *, int, int, int *, int **);
void update_defects(polygons_struct *, int *, int *);
void update_polydefects(polygons_struct *, int *, int *);
Point get_defect_center(polygons_struct *, int *, int);

double get_holes_handles(polygons_struct *, polygons_struct *, int *, int,
                         int *, Volume, int *, int **);

void bisect_defects(polygons_struct *, int *, int, int *, int *);

void remap_defect(polygons_struct *, int *, int *, polygons_struct *, int *,
                  int *);
