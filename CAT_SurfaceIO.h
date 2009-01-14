/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

Status input_values_any_format(char *, int *, Real **);
Status output_values_any_format(char *, int, Real *);
Status input_graphics_any_format(char *, File_formats *, int *,
                                 object_struct  ***); 
Status output_graphics_any_format(char *, File_formats, int, object_struct  **);
int    input_oogl(char *, File_formats *, int *, object_struct  ***);
int    output_oogl(char *, File_formats, int, object_struct * []);
int    output_freesurfer(char *, File_formats, int, object_struct * []);
int    input_freesurfer(char *, File_formats *, int *, object_struct  ***);
int    input_freesurfer_curv(char *, int *, Real **);
int    input_dx(char *, File_formats *, int *, object_struct  ***);
int    input_dfs(char *, File_formats *, int *, object_struct  ***);
