#ifndef  DEF_deform_prototypes
#define  DEF_deform_prototypes

void       deform_polygons(polygons_struct *, deform_struct *);
void       deform_lines(lines_struct *, deform_struct *);
void       deform_lines_one_iter(lines_struct *, deform_struct *, int);

void       get_line_equilibrium_point(lines_struct *, int, int, int [],
                                      Real [], Real, int, Volume, Volume,
                                      bound_def_struct *, deform_model_struct *,
                                      Point *, Point *);
int        find_axial_plane(lines_struct *);

void       deform_polygons(polygons_struct *, deform_struct *);
void       deform_polygons_one_iter(polygons_struct *, deform_struct *, int);

BOOLEAN    find_boundary_in_direction(Volume, Volume, voxel_coeff_struct *,
                                      bitlist_3d_struct *, bitlist_3d_struct *,
                                      Real, Point *, Vector *, Vector *, Real,
                                      Real, int, bound_def_struct *, Real *);
int        find_voxel_line_polynomial(Real [], int, int, int, int,
                                      Real [], Real [], Real []);
int        find_voxel_line_value_intersection(Real [], int, int, int, int,
                                              Real [], Real [], Real, Real,
                                              Real, Real [3]);

void       initialize_deformation_model(deform_model_struct *);
void       print_deformation_model(deform_model_struct *);
Status     add_deformation_model(deform_model_struct *, int, Real,
                                 char [], Real, Real);
void       delete_deformation_model(deform_model_struct *);

Status     input_original_positions(deform_model_struct *, char [], Real, int);

BOOLEAN    check_correct_deformation_polygons(polygons_struct *,
                                              deform_model_struct *);
BOOLEAN    check_correct_deformation_lines(lines_struct *,
                                           deform_model_struct *);

void       get_model_shape_point(Point *, Vector *, Vector *, Real, Point *);
void       compute_equilibrium_point(int, BOOLEAN, Real, Real, Real, Vector *,
                                     Vector *, Point *, deform_model_struct *,
                                     Point *);
void       compute_model_dirs(Point *, Vector *, Real, Point *, Real *,
                              Point *, Vector *, Vector *);
void       get_model_point(deform_model_struct *, Point [], int, int, int [],
                           Real [], Point *, Vector *, Real, Point *);

void       get_neighbours_of_line_vertex(lines_struct *, int, int [2]);
BOOLEAN    deformation_model_includes_average(deform_model_struct *);
Real       compute_line_curvature(lines_struct *, int, int, int, int);
Real       deform_point(int, Point [], Point *, Real, Real, BOOLEAN,
                        Real, Point [], Point *);

void       compute_line_centroid_and_normal(lines_struct *, int, int, int,
                                            Point *, Vector *, Real *);
int        get_subsampled_neighbours_of_point(deform_model_struct *,
                                              polygons_struct *, int, int,
                                              int [], int, BOOLEAN *);
BOOLEAN    is_point_inside_surface(Volume, Volume, int, Real [],
                                   Vector *, bound_def_struct *);

void       get_centre_of_cube(Point *, int [3], Point *);
BOOLEAN    contains_value(Real [2][2][2], int [3]);
BOOLEAN    cube_is_small_enough(Point [2], int [3], Real);

void       initialize_deform_stats(deform_stats  *);
void       record_error_in_deform_stats(deform_stats *, Real);
void       print_deform_stats(deform_stats *, int);
BOOLEAN    get_max_point_cube_dist(Point [2], int [3], Point *, Real *);
void       init_deform_params(deform_struct *);
void       delete_deform_params(deform_struct *);

void       set_boundary_definition(bound_def_struct *, Real, Real,
                                   Real, Real, char, Real);

void       initialize_lookup_volume_coeffs(voxel_coeff_struct *);
void       lookup_volume_coeffs(voxel_coeff_struct *, Volume, int,
                                int, int, int, Real []);
void       delete_lookup_volume_coeffs(voxel_coeff_struct *);

deform_model_struct * find_relevent_model(deform_model_struct *, int);

#endif
