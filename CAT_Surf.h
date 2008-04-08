/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

float* get_surface_ratio(
	float r,
	polygons_struct *polygons);

float get_area_of_points(
    polygons_struct     *polygons,
    float               *area_values);

void translate_to_center_of_mass(
    polygons_struct     *polygons);

float get_largest_dist(
    polygons_struct     *polygons);

void set_vector_length(
    Point   *p,
    float   newLength);

void get_radius_of_points(
    polygons_struct     *polygons,
    float               *radius);

void get_bounds(
    polygons_struct     *polygons,
    float bounds[6]);

int  count_edges(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[] );
    
int euler_characteristic(
    polygons_struct   *polygons);

void convert_ellipsoid_to_sphere_with_surface_area(
    polygons_struct     *polygons,
    float         desiredSurfaceArea);

void  linear_smoothing(
    polygons_struct     *polygons,
    float                strength,
    int                 iterations,
    int                 smoothEdgesEveryXIterations,
    int                 *smoothOnlyTheseNodes,
    int                 projectToSphereEveryXIterations);

void  areal_smoothing(
    polygons_struct     *polygons,
    float                strength,
    int                 iterations,
    int                 smoothEdgesEveryXIterations,
    int                 *smoothOnlyTheseNodes,
    int                 projectToSphereEveryXIterations);

void  distance_smoothing(
    polygons_struct     *polygons,
    float               strength,
    int                 iterations,
    int                 smoothEdgesEveryXIterations,
    int                 *smoothOnlyTheseNodes,
    int                 projectToSphereEveryXIterations);

void inflate_surface_and_smooth_fingers(
    polygons_struct     *polygonsIn,
    const int numberSmoothingCycles,
    const float regularSmoothingStrength,
    const int regularSmoothingIterations,
    const float inflationFactorIn,
    const float compressStretchThreshold,
    const float fingerSmoothingStrength,
    const int fingerSmoothingIterations);
