#ifndef  DEF_MARCHING_CUBES_DATA
#define  DEF_MARCHING_CUBES_DATA

#include  <bicpl.h>
#include  <bicpl/marching_cubes.h>

typedef struct
{

    int                  n_triangles;
    voxel_point_type     triangle_points[MAX_POINTS_PER_VOXEL];

}  marching_cube_data_struct;

#endif
