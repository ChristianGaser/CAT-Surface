/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

void  get_polygon_vertex_curvatures_cg(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    Real              smoothing_distance,
    int               curvtype,    
    Real              curvatures[] );
