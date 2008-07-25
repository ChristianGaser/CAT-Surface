/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>

#include "CAT_Deform.h"

void
usage(char *executable)
{
        fprintf(stderr, "%s  volume_filename\n", executable);
        fprintf(stderr, "   activity_filename|none   nx ny nz\n");
        fprintf(stderr, "   input_polygons output_polygons\n");
        fprintf(stderr, "   original_positions|none max_distance\n");
        fprintf(stderr, "   n_models\n");
        fprintf(stderr, "   up_to_n_points model_weight model_file|avg|none\n");
        fprintf(stderr, "   min_curvature max_curvature\n");
        fprintf(stderr, "   [up_to_n_points model_weight model_file|avg|none\n");
        fprintf(stderr, "   min_curvature max_curvature]\n");
        fprintf(stderr, "   fract_step max_step\n");
        fprintf(stderr, "   max_search_distance degrees_continuity\n");
        fprintf(stderr, "   min_isovalue max_isovalue +/-/n\n");
        fprintf(stderr, "   gradient_threshold angle tolerance\n");
        fprintf(stderr, "   max_iterations  movement_threshold stop_threshold\n");
}

int
main(int argc, char *argv[])
{
        Status            status;
        Real              start_time, end_time;
        char              *volume_file, *activity_file;
        char              *input_file, *output_file;
        char              *model_file, *normal_direction, *original_file;
        Real              min_isovalue, max_isovalue, gradient_thresh;
        Real              model_weight, min_curv_offset, max_curv_offset;
        Real              angle, tolerance, max_dist;
        Real              separations[N_DIMENSIONS];
        Real              x_filter_width, y_filter_width, z_filter_width;
        int               i, n_models, up_to_n_points;
        deform_struct     deform;
        File_formats      file_format;
        int               n_objects;
        object_struct     **object_list;
        Volume            volume, label_volume, tmp;
        polygons_struct   *polygons;

        set_alloc_checking(FALSE);

        initialize_argument_processing(argc, argv);

        initialize_deformation_parameters(&deform);

        if (!get_string_argument("", &volume_file) ||
            !get_string_argument("", &activity_file) ||
            !get_real_argument(0.0, &x_filter_width) ||
            !get_real_argument(0.0, &y_filter_width) ||
            !get_real_argument(0.0, &z_filter_width) ||
            !get_string_argument("", &input_file) ||
            !get_string_argument("", &output_file) ||
            !get_string_argument("", &original_file) ||
            !get_real_argument(0.0, &max_dist) ||
            !get_int_argument(1, &n_models)) {
                usage(argv[0]);
                return(1);
        }

        for (i = 0; i < n_models; i++) {
                if (!get_int_argument(0, &up_to_n_points) ||
                    !get_real_argument(0.0, &model_weight) ||
                    !get_string_argument("", &model_file) ||
                    !get_real_argument(0.0, &min_curv_offset) ||
                    !get_real_argument(0.0, &max_curv_offset)) {
                        usage(argv[0]);
                        return(1);
                }

                if (add_deformation_model(&deform.deform_model, up_to_n_points,
                                          model_weight, model_file,
                                          min_curv_offset,
                                          max_curv_offset) != OK)
                        return(1);
        }

        if (!get_real_argument(0.0, &deform.fractional_step) ||
            !get_real_argument(0.0, &deform.max_step) ||
            !get_real_argument(0.0, &deform.max_search_dist) ||
            !get_int_argument(0.0, &deform.degrees_continuity) ||
            !get_real_argument(0.0, &min_isovalue) ||
            !get_real_argument(0.0, &max_isovalue) ||
            !get_string_argument("", &normal_direction) ||
            !get_real_argument(0.0, &gradient_thresh) ||
            !get_real_argument(0.0, &angle) ||
            !get_real_argument(0.0, &tolerance) ||
            !get_int_argument(0, &deform.max_iters) ||
            !get_real_argument(0.0, &deform.movement_thresh) ||
            !get_real_argument(0.0, &deform.stop_thresh)) {
                usage(argv[0]);
                return(1);
        }

        set_boundary_definition(&deform.bound_def, min_isovalue, max_isovalue,
                                gradient_thresh, angle, normal_direction[0],
                                tolerance);

        deform.deform_data.type = VOLUME_DATA;

        if (input_volume(volume_file, 3, XYZ_dimension_names,
                         NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE,
                         &volume, (minc_input_options *) NULL) == ERROR)
                return(1);

        label_volume = (Volume) NULL;

        if (x_filter_width > 0.0 && y_filter_width > 0.0 &&
                                    z_filter_width > 0.0) {
                get_volume_separations(volume, separations);

                x_filter_width /= fabs(separations[X]);
                y_filter_width /= fabs(separations[Y]);
                z_filter_width /= fabs(separations[Z]);

                tmp = create_box_filtered_volume(volume, NC_BYTE, FALSE,
                                                 0.0, 0.0,
                                                 x_filter_width,
                                                 y_filter_width,
                                                 z_filter_width);

                delete_volume(volume);
                volume = tmp;
        }

        deform.deform_data.volume = volume;
        deform.deform_data.label_volume = label_volume;

        if (input_graphics_any_format(input_file, &file_format,
                                      &n_objects, &object_list) == ERROR ||
            n_objects != 1 || object_list[0]->object_type != POLYGONS) {
                fprintf(stderr, "File must contain 1 polygons struct.\n");
                return(1);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (!equal_strings(original_file, "none")) {
                if (input_original_positions(&deform.deform_model,
                                             original_file, max_dist,
                                             polygons->n_points) == ERROR)
                        return(1);
        }

        start_time = current_cpu_seconds();
        deform_polygons(polygons, &deform);

        compute_polygon_normals(polygons);

        end_time = current_cpu_seconds();

        if (output_graphics_any_format(output_file, ASCII_FORMAT,
                                       n_objects, object_list) == ERROR)
                return(1);

        return(0);
}
