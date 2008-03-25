/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

public int output_freesurfer(
    char *fname,
    File_formats   format,
    int            n_objects,
    object_struct  *object_list[]);

public int input_freesurfer(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);

public int input_freesurfer_curv(
    char *fname,
    int	 *n_values,
    Real  *input_values[]);

public int input_dx(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);
    
public int input_dfs(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list);

