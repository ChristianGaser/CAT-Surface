/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

typedef struct {
	long x;
	long y;
} Header;

typedef struct {
	float x;
	float y;
} Vector2D;


public  void  point_to_uv(
    Point            *point,
    double           *u,
    double           *v );

public  double point_to_uv_radius(
    Point            *point,
    double           *u,
    double           *v );

public void  uv_to_point(
    double           u,
    double           v,
    Point            *point );

void map_smoothed_curvature_to_sphere(
    polygons_struct  *polygons,
    double           *values,
    double           *data,
    double           fwhm,
    int              *size_map);

void map_sheet2d_to_sphere(
    double           *sheet2d,
    double           *values,
    polygons_struct  *polygons,
    int              interpolate,
    int              *size_xy);
