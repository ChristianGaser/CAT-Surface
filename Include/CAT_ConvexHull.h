/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_ConvexHull.h 267 2012-06-21 14:57:59Z gaser $
 *
 */

typedef  struct
{
    int  poly;
    int  edge;
} queue_entry;

typedef struct
{
    Smallest_int  ref_count;
} edge_struct;

typedef  QUEUE_STRUCT( queue_entry )  queue_struct;

#define  ENLARGE_THRESHOLD         0.25
#define  NEW_DENSITY               0.125
#define  KEY_FACTOR                1000000

#define  TOLERANCE_2D   1.0e-3
#define  TOLERANCE_DISTANCE   1.0e-6

#define  POINT_USED_IN_CONVEX_HULL  1
#define  POINT_DISCARDED            2

static int KeyFactor = 100000;

int
calculate_convex_hull(polygons_struct *);


double
gyrification_index_sph(polygons_struct *surface, polygons_struct *sphere,
                      char *file, int n_triangles, polygons_struct *reparam,
                      int smoothflag);

