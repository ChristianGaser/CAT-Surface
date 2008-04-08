/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

Status   input_graphics_any_format(
    STRING         filename,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list );

Status   output_graphics_any_format(
    STRING         filename,
    File_formats   format,
    int            n_objects,
    object_struct  **object_list );

int input_oogl(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);

int output_oogl(
    char *fname,
    File_formats  format,
    int           n_objects,
    object_struct *object_list[]);

int output_freesurfer(
    char *fname,
    File_formats   format,
    int            n_objects,
    object_struct  *object_list[]);

int input_freesurfer(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);

int input_freesurfer_curv(
    char *fname,
    int	 *n_values,
    Real  *input_values[]);

int input_dx(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);
    
int input_dfs(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);

