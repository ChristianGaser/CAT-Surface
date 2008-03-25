/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

public void get_all_polygon_point_neighbours(
    polygons_struct  *polygons,
    int              *n_point_neighbours_ptr[],
    int              **point_neighbours_ptr[] );

public  void  heatkernel_blur_points(
    int               n_polygon_points,
    Point             polygon_points[],
    Real              values[],
    int               n_neighbours,
    int               *neighbours,
    int               point_index,
    Real              sigma,
    Point             *smooth_point,
    Real              *value);

