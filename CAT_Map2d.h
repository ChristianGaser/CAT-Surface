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
    Real             *u,
    Real             *v );

public  Real point_to_uv_radius(
    Point            *point,
    Real             *u,
    Real             *v );

public void  uv_to_point(
    Real             u,
    Real             v,
    Point            *point );
