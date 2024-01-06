/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl/deform.h>

#include "CAT_DeformPolygons.h"
#include "CAT_SurfaceIO.h"
#include "CAT_NiftiIO.h"
#include "CAT_Surf.h"

void
usage(char *executable)
{
        fprintf(stderr, "%s volume_file\n", executable);
        fprintf(stderr, "   activity_file|none   nx ny nz\n");
        fprintf(stderr, "   input_surface_file output_surface_file\n");
        fprintf(stderr, "   original_positions|none max_distance\n");
        fprintf(stderr, "   n_models\n");
        fprintf(stderr, "   up_to_n_points model_weight model_file|avg|flat|parametric\n");
        fprintf(stderr, "   min_curvature max_curvature\n");
        fprintf(stderr, "   [up_to_n_points model_weight model_file|avg|none\n");
        fprintf(stderr, "   min_curvature max_curvature]\n");
        fprintf(stderr, "   fract_step max_step\n");
        fprintf(stderr, "   max_search_distance degrees_continuity\n");
        fprintf(stderr, "   min_isovalue max_isovalue +/-/n\n");
        fprintf(stderr, "   gradient_threshold angle tolerance\n");
        fprintf(stderr, "   max_iterations  movement_threshold stop_threshold\n");
        fprintf(stderr, "   force_no_selfintersections\n");
        fprintf(stderr, "   correct_mesh\n");
}

int
main(int argc, char *argv[])
{
        Status            status;
        double            start_time, end_time;
        char              *volume_file = NULL, *activity_file = NULL;
        char              *input_file = NULL, *output_surface_file = NULL;
        char              *model_file = NULL, *original_file = NULL;
        char              *normal_direction = NULL;
        double            min_isovalue, max_isovalue, gradient_thresh;
        double            model_weight, min_curv_offset, max_curv_offset;
        double            angle, tolerance, max_dist = 0;
        double            separations[N_DIMENSIONS];
        double            xfilt_width = 0, yfilt_width = 0, zfilt_width = 0;
        int               i, n_models = 0, up_to_n_points, correct_mesh;
        deform_struct     deform;
        File_formats      file_format;
        int               n_objects, check_every_iteration, force_no_selfintersections;
        object_struct     **object_list;
        Volume            volume, label_volume, tmp;
        polygons_struct   *polygons;

        set_alloc_checking(FALSE);

        initialize_argument_processing(argc, argv);

        initialize_deformation_parameters(&deform);

        if (get_string_argument("", &volume_file) == 0 ||
            get_string_argument("", &activity_file) == 0 ||
            get_real_argument(0.0, &xfilt_width) == 0 ||
            get_real_argument(0.0, &yfilt_width) == 0 ||
            get_real_argument(0.0, &zfilt_width) == 0 ||
            get_string_argument("", &input_file) == 0 ||
            get_string_argument("", &output_surface_file) == 0 ||
            get_string_argument("", &original_file) == 0 ||
            get_real_argument(0.0, &max_dist) == 0 ||
            get_int_argument(1, &n_models) == 0) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        for (i = 0; i < n_models; i++) {
                if (get_int_argument(0, &up_to_n_points) == 0 ||
                    get_real_argument(0.0, &model_weight) == 0 ||
                    get_string_argument("", &model_file) == 0 ||
                    get_real_argument(0.0, &min_curv_offset) == 0 ||
                    get_real_argument(0.0, &max_curv_offset) == 0) {
                        usage(argv[0]);
                        exit(EXIT_FAILURE);
                }

                if (add_deformation_model(&deform.deformation_model,
                                          up_to_n_points,
                                          model_weight, model_file,
                                          min_curv_offset,
                                          max_curv_offset) != OK)
                        exit(EXIT_FAILURE);
        }

        if (get_real_argument(0.0, &deform.fractional_step) == 0 ||
            get_real_argument(0.0, &deform.max_step) == 0 ||
            get_real_argument(0.0, &deform.max_search_distance) == 0 ||
            get_int_argument(0.0, &deform.degrees_continuity) == 0 ||
            get_real_argument(0.0, &min_isovalue) == 0 ||
            get_real_argument(0.0, &max_isovalue) == 0 ||
            get_string_argument("", &normal_direction) == 0 ||
            get_real_argument(0.0, &gradient_thresh) == 0 ||
            get_real_argument(0.0, &angle) == 0 ||
            get_real_argument(0.0, &tolerance) == 0 ||
            get_int_argument(0, &deform.max_iterations) == 0 ||
            get_real_argument(0.0, &deform.movement_threshold) == 0 ||
            get_real_argument(0.0, &deform.stop_threshold) == 0) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(0, &force_no_selfintersections);
        get_int_argument(0, &correct_mesh);
        
        set_boundary_definition(&deform.boundary_definition, min_isovalue,
                                max_isovalue, gradient_thresh, angle,
                                normal_direction[0], tolerance);

        deform.deform_data.type = VOLUME_DATA;

        if (input_volume_all(volume_file, 3, XYZ_dimension_names,
                             NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE,
                             &volume, (minc_input_options *) NULL) == ERROR)
                exit(EXIT_FAILURE);

        label_volume = (Volume) NULL;

        if (xfilt_width > 0.0 && yfilt_width > 0.0 && zfilt_width > 0.0) {
                get_volume_separations(volume, separations);

                xfilt_width /= fabs(separations[X]);
                yfilt_width /= fabs(separations[Y]);
                zfilt_width /= fabs(separations[Z]);

                tmp = create_box_filtered_volume(volume, NC_BYTE, FALSE,
                                                 0.0, 0.0, xfilt_width,
                                                 yfilt_width, zfilt_width);

                delete_volume(volume);
                volume = tmp;
        }

        deform.deform_data.volume = volume;
        deform.deform_data.label_volume = label_volume;

        if (input_graphics_any_format(input_file, &file_format,
                                      &n_objects, &object_list) == ERROR ||
            n_objects != 1 || object_list[0]->object_type != POLYGONS) {
                fprintf(stderr, "File must contain 1 polygons struct.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (!equal_strings(original_file, "none")) {
                if (input_original_positions(&deform.deformation_model,
                                             original_file, max_dist,
                                             polygons->n_points) == ERROR)
                        exit(EXIT_FAILURE);
        }

        start_time = current_cpu_seconds();
        
        /* check every 100 iterations for self interactions */
        if (force_no_selfintersections) {
                check_every_iteration = 25;
                deform_polygons_check_selfintersection(polygons, &deform, check_every_iteration, force_no_selfintersections);
        }
        else {
                check_every_iteration = 100;
                deform_polygons_check_selfintersection_old(polygons, &deform, check_every_iteration, force_no_selfintersections);
        }

        /* optionally correct mesh w.r.t. folding to achieve better fit to isovolue */
        if (correct_mesh)
                correct_mesh_folding(polygons, NULL, volume, min_isovalue);

        compute_polygon_normals(polygons);

        end_time = current_cpu_seconds();

        if (output_graphics_any_format(output_surface_file, ASCII_FORMAT,
                                       n_objects, object_list, NULL) == ERROR)
                exit(EXIT_FAILURE);

        if (deform.prev_movements != NULL)
                free(deform.prev_movements);
        delete_volume(deform.deform_data.volume);
        return(EXIT_SUCCESS);
}
