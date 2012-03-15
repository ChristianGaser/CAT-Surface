/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_DeformPolygons.h 239 2012-03-14 23:57:40Z gaser $
 *
 */

#define  MAX_NEIGHBOURS  2000

void perturb_points(polygons_struct *, Point [], Real, Real, Real, int,
                    deform_data_struct *, boundary_definition_struct *,
                    deformation_model_struct *, Real, float [], deform_stats *);
void perturb_points_points(polygons_struct *, Point [], Real, Real, Real, int,
                           deform_data_struct *, boundary_definition_struct *,
                           deformation_model_struct *, Real, float [],
                           deform_stats *, int *);
double one_iter_polygons(polygons_struct *, deform_struct *, int);
double one_iter_polygons_points(polygons_struct *, deform_struct *, int, int *);
void check_polygons_shape_integrity(polygons_struct *, Point []);
void check_shape_integrity_points(polygons_struct *, Point [], int *);
void deform_polygons_points(polygons_struct *polygons, deform_struct *deform_parms,
                       int *flag);
