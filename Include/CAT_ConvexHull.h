/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

static int dbg = 0;
static int dbg2 = 0;

#include  <bicpl.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Refine.h"
#include "CAT_Resample.h"

#define  TOLERANCE_2D   1.0e-3
#define  TOLERANCE_DISTANCE   1.0e-6

#define  POINT_USED_IN_CONVEX_HULL  1
#define  POINT_DISCARDED            2

static int KeyFactor = 100000;

object_struct **  surface_get_convex_hull(polygons_struct *, polygons_struct * );
private  int  get_points_of_region(polygons_struct  *, Point ** );
private  void  get_convex_hull(int, Point *, polygons_struct * );
private  int  get_convex_hull_2d(int, Real *, Real *, int *, int, int );
private int get_surface_point_normals( polygons_struct *, int *, Point *[],
                              Vector *[], int *[], int **[] );
private int get_surface_neighbours( polygons_struct *, int *[],
                                    int ** [] );
void find_conformal_map(polygons_struct *polygons);
