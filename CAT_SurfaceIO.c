/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_SurfaceIO.h"

#define TRIANGLE_FILE_MAGIC_NUMBER  16777214
#define QUAD_FILE_MAGIC_NUMBER  16777215
#define NEW_VERSION_MAGIC_NUMBER 16777215

Status
input_values(char *file, int n_values, Real *values)
{
        FILE     *fp;
        int      i;

        if (filename_extension_matches(file,"txt")) {
                if (open_file(file, READ_FILE, ASCII_FORMAT, &fp) != OK)
                        return(1);

                for (i = 0; i < n_values; i++) {
                        if (input_real(fp, &values[i]) != OK)
                                return( 1 );
                }
                close_file(fp);
        } else {
                return(input_freesurfer_curv(file, &n_values, &values));
        }
}

Status
input_graphics_any_format(char *file, File_formats *format, int *n_objects,
                          object_struct ***object_list)
{
        Status status;

        if (filename_extension_matches(file, "obj")) {
                status = input_graphics_file(file, format,
                                             n_objects, object_list);
        } else if (filename_extension_matches(file, "off")) {
                status = input_oogl(file, format,
                                    n_objects, object_list);
        } else if (filename_extension_matches(file, "dfs")) {
                status = input_dfs(file, format,
                                   n_objects, object_list);
        } else if (filename_extension_matches(file, "dx")) {
                status = input_dx(file, format,
                                  n_objects, object_list);
        } else {
                status = input_freesurfer(file, format,
                                          n_objects, object_list);
        }

        return(status);
}

Status
output_graphics_any_format(char *file, File_formats format, int n_objects,
                           object_struct **object_list)
{
        Status     status;

        if (filename_extension_matches(file, "obj")) {
                status = output_graphics_file(file, format,
                                              n_objects, object_list);
        } else if (filename_extension_matches(file, "off")) {
                status = output_oogl(file, format,
                                     n_objects, object_list);
        } else {
                status = output_freesurfer(file, format,
                                           n_objects, object_list);
        }

        return(status);
}

void
swapDouble(double *n)
{
        char *by = (char *) n;
        char sw[4] = {by[3], by[2], by[1], by[0]};
  
        *n =* (double *) sw;
}

void
swapInt(int *n)
{
        char *by = (char *) n;
        char sw[4] = {by[3], by[2], by[1], by[0]};
  
        *n =* (int *) sw;
}

int
fwriteDouble(double f, FILE *fp)
{

#if (BYTE_ORDER == LITTLE_ENDIAN)
       swapDouble(&f);  
#endif
       return(fwrite(&f, sizeof(double), 1, fp));
}

int
fwriteInt(int v, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
        swapInt(&v);
#endif
        return(fwrite(&v, sizeof(int), 1, fp));
}

int
freadInt(FILE *fp)
{
        int temp;
        int count;
  
        count = fread(&temp, 4, 1, fp);
#if __LITTLE_ENDIAN__
        swapInt(&temp);
#endif
        return(temp);
}

double
freadDouble(FILE *fp)
{
        double temp;
        int count;
  
        count = fread(&temp, 4, 1, fp);
#if __LITTLE_ENDIAN__
        swapDouble(&temp);
#endif
        return(temp);  
}

int
fread3(int *v, FILE *fp)
{
        int i = 0;
        int  ret ;

        ret = fread(&i, 3, 1, fp);
#if __LITTLE_ENDIAN__
        swapInt(&i) ;
#endif
        *v = ((i >> 8) & 0xffffff);
        return(ret);
}

int
input_oogl(char *file, File_formats *format, int *n_objects,
           object_struct ***object_list)
{
        FILE              *fp;
        int               i, f, n_edges, dummy;
        char              line[1000];
        polygons_struct   *polygons;
        Point             point;
        object_struct     *object;

        /* prepare object and polygons */
        *n_objects = 0;
        object = create_object(POLYGONS);
        add_object_to_list(n_objects, object_list, object);
        polygons = get_polygons_ptr(object);
        initialize_polygons(polygons, WHITE, NULL);
        *format = ASCII_FORMAT;

        if ((fp = fopen(file, "rb")) == 0) {
                fprintf(stderr, "input_oogl: Couldn't open file %s.\n", file);
                return(0);
        }

        fscanf(fp, "%s", line);

        if (!strcmp(line,"OFF")) {
                fscanf(fp, "%d %d %d", &polygons->n_points,
                       &polygons->n_items, &dummy);
                ALLOC(polygons->points, polygons->n_points);
                ALLOC(polygons->normals, polygons->n_points);
                ALLOC(polygons->end_indices, polygons->n_items);

                polygons->bintree = (bintree_struct_ptr) NULL;

                for (i = 0; i < polygons->n_points; i++) {
                        fscanf(fp, "%f %f %f\n", &Point_x(point),
                               &Point_y(point), &Point_z(point)); 
                        polygons->points[i] = point;
                }
    
                for (i = 0; i < (polygons->n_items); i++)
                        polygons->end_indices[i] = (i + 1) * 3;

                ALLOC(polygons->indices,
                      polygons->end_indices[polygons->n_items-1]);
                for (f = 0; f < polygons->n_items; f++) {
                        fscanf(fp, "%d %d %d %d\n", &n_edges, 
                               &polygons->indices[POINT_INDEX(
                                                 polygons->end_indices, f, 0)],
                               &polygons->indices[POINT_INDEX(
                                                 polygons->end_indices, f, 1)],
                               &polygons->indices[POINT_INDEX(
                                                 polygons->end_indices, f, 2)]);
                        if (n_edges != 3) {
                                fprintf(stderr,
                                        "Only 3 elements per item allowed.\n");
                                return(0);
                        }
                }
    
                /* compute normals */
                compute_polygon_normals(polygons);
        } else {
                fprintf(stderr, "Wrong oogl format\n");
                return(0);
        }
  
        fclose(fp);
        return(OK);
}

int
output_oogl(char *file, File_formats format, int n_objects,
            object_struct *object_list[])
{
        int i, e, f, index, n_face_edges;
        FILE  *fp;
        polygons_struct   *polygons;

        if ((fp = fopen(file, "w")) == 0) {
                fprintf(stderr, "output_oogl: Couldn't open file %s.\n", file);
                return(0);
        }

        polygons = get_polygons_ptr(object_list[0]);

        fprintf(fp, "OFF\n");
        fprintf(fp, "%d %d %d\n", polygons->n_points, polygons->n_items, -1);
        for (i = 0; i < polygons->n_points; i++) {
                fprintf(fp, "%f %f %f\n", 
                        Point_x(polygons->points[i]),
                        Point_y(polygons->points[i]), 
                        Point_z(polygons->points[i]));
        }
        for (f = 0; f < polygons->n_items; f++) {
                n_face_edges = GET_OBJECT_SIZE(*polygons, f);

                fprintf(fp, "%d ", n_face_edges);
                for (e = 0; e < n_face_edges; e++) {
                        index = polygons->indices[POINT_INDEX(
                                                  polygons->end_indices, f, e)];
                        fprintf(fp, "%d ", index);
                }
                fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
        fclose(fp);
        return(1);
}

int
output_freesurfer(char *file, File_formats format, int n_objects,
                  object_struct *object_list[])
{
        FILE              *fp;
        int               i;
        unsigned char     buffer;
        polygons_struct   *polygons;

        if ((fp = fopen(file, "w")) == 0) {
                fprintf(stderr, "output_freesurfer: Couldn't open file %s.\n",
                        file);
                return(0);
        }

        polygons = get_polygons_ptr(object_list[0]);

        /* write 3 byte magic number for triangle */
        buffer = 255;
        fwrite(&buffer, 1, 1, fp);
        fwrite(&buffer, 1, 1, fp);
        buffer = 254;
        fwrite(&buffer, 1, 1, fp);

        output_newline(fp); 
        output_newline(fp); 

        /* # of vertices and faces */
        fwriteInt(polygons->n_points, fp);
        fwriteInt(polygons->n_items, fp);

        /* write points */
        for (i = 0; i < (polygons->n_points); i++) {
                fwriteDouble(Point_x(polygons->points[i]), fp);
                fwriteDouble(Point_y(polygons->points[i]), fp);
                fwriteDouble(Point_z(polygons->points[i]), fp);
        }
  
        /* write indices */
        for (i = 0; i < 3*(polygons->n_items); i++)
                fwriteInt(polygons->indices[i], fp);

        fclose(fp);
        return(1);
}

int
input_freesurfer(char *file, File_formats *format, int *n_objects,
                 object_struct  ***object_list)
{
        FILE              *fp;
        int               i, magic;
        char              line[1000];
        polygons_struct   *polygons;
        Point             point;
        object_struct     *object;
      
        /* prepare object and polygons */
        *n_objects = 0;
        object = create_object(POLYGONS);
        add_object_to_list(n_objects, object_list, object);
        polygons = get_polygons_ptr(object);
        initialize_polygons(polygons, WHITE, NULL);
        *format = ASCII_FORMAT;

        if ((fp = fopen(file, "rb")) == 0) {
                fprintf(stderr, "input_freesurfer: Couldn't open file %s.\n",
                        file);
                return(0);
        }

        /* read magic number for checking filetype */
        fread3(&magic, fp);
        if (magic == QUAD_FILE_MAGIC_NUMBER) {
                fprintf(stderr, "QUAD_FILE_MAGIC_NUMBER not yet prepared.\n");
                return(0);
        } else if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
                fgets(line, 1024, fp);
                fscanf(fp, "\n");
                /* read # of vertices and faces */
                polygons->n_points = freadInt(fp);
                polygons->n_items = freadInt(fp);
                ALLOC(polygons->points, polygons->n_points);
                ALLOC(polygons->normals, polygons->n_points);
                ALLOC(polygons->end_indices, polygons->n_items);
                polygons->bintree = (bintree_struct_ptr) NULL;
                for (i = 0; i < polygons->n_items; i++)
                        polygons->end_indices[i] = (i + 1) * 3;
                ALLOC(polygons->indices,
                      polygons->end_indices[polygons->n_items-1]);
                for (i = 0; i < polygons->n_points; i++) {
                        Point_x(point) = freadDouble(fp);
                        Point_y(point) = freadDouble(fp);
                        Point_z(point) = freadDouble(fp);
                        polygons->points[i] = point;
                }
                for (i = 0; i < 3*(polygons->n_items); i++) {
                        polygons->indices[i] = freadInt(fp);
                }
                /* compute normals */
                compute_polygon_normals(polygons);
        } else {
                fprintf(stderr, "input_freesurfer: Unknown magic identifier: %d.\n", magic);
                return(0);
        }

        fclose(fp);
        return(OK);
}

int
input_freesurfer_curv(char *file, int *vnum, Real *input_values[])
{
        FILE  *fp;
        int   i, magic, fnum, vals_per_vertex;
        Real  value;
      
        *vnum = 0;
  
        if ((fp = fopen(file, "rb")) == 0) {
                fprintf(stderr,
                        "input_freesurfer_curv: Couldn't open file %s.\n",
                        file);
                return(0);
        }

        /* read magic number for checking filetype */
        fread3(&magic,fp);
  
        if (magic != NEW_VERSION_MAGIC_NUMBER) {
                fprintf(stderr, "MAGIC_NUMBER %d not yet prepared.\n", magic);
                return(0);
        } else {
                /* read # of vertices and faces */
                *vnum = freadInt(fp);
                fnum = freadInt(fp);
                vals_per_vertex = freadInt(fp);
                if (vals_per_vertex != 1) {
                        fprintf(stderr, "Only one value per vertex allowed.\n");
                        return(0);
                }
                ALLOC(*input_values, *vnum);
                for (i = 0; i < *vnum; i++) {
                        value = freadDouble(fp);
                        (*input_values)[i]= value;
                }
        }

        fclose(fp);
        return(OK);
}

int
input_dx(char *file, File_formats *format, int *n_objects,
         object_struct  ***object_list)
{
        FILE              *fp;
        int               i, pos;
        char              line[1000];
        char              ch;
        polygons_struct   *polygons;
        Point             point;
        object_struct     *object;
        Status            status;
      
        /* prepare object and polygons */
        *n_objects = 0;
        object = create_object(POLYGONS);
        add_object_to_list(n_objects, object_list, object);

        polygons = get_polygons_ptr(object);
        initialize_polygons(polygons, WHITE, NULL);

        *format = ASCII_FORMAT;

        if ((fp = fopen(file, "r")) == 0) {
                fprintf(stderr, "input_dx: Couldn't open file %s.\n", file);
                return(0);
        }

        fgets(line, 54, fp);
        if (equal_strings(line,
            "object 1 class array type float rank 1 shape 3 items " )) {
                pos = 0;
                status = OK;

                do {
                        status = input_character(fp, &ch);
                        pos++;
                } while (status == OK && ch != ' ');

                fseek(fp, -pos, SEEK_CUR);
                fgets(line, pos, fp);
                sscanf(line,"%d", &polygons->n_points);

                polygons->bintree = (bintree_struct_ptr) NULL;
                ALLOC(polygons->points, polygons->n_points);
                ALLOC(polygons->normals, polygons->n_points);
                for (i = 0; i <= polygons->n_points; i++) {
                        fgets(line, 1000, fp);
                        sscanf(line, "%f %f %f", &Point_x(point),
                               &Point_y(point), &Point_z(point));
                        polygons->points[i] = point;
                }
                fgets(line, 1000, fp);
                fgets(line, 1000, fp);
                fgets(line, 52, fp);

                if (equal_strings(line,
                    "object 2 class array type int rank 1 shape 3 items ")) {
                        pos = 0;
                        status = OK;

                        do {
                                status = input_character(fp, &ch);
                                pos++;
                        } while (status == OK && ch != ' ');

                        fseek(fp, -pos, SEEK_CUR);
                        fgets(line, pos, fp);
                        sscanf(line,"%d", &polygons->n_items);
                } else {
                        fprintf(stderr, "input_dx: Error reading %s\n", file);
                        return(0);
                }

                ALLOC(polygons->end_indices, polygons->n_items);
                for (i = 0; i < polygons->n_items; i++)
                        polygons->end_indices[i] = (i + 1) * 3;

                ALLOC(polygons->indices,
                      polygons->end_indices[polygons->n_items-1]);

                for (i = 0; i <= polygons->n_items; i++) {
                        fgets(line, 1000, fp);
                        sscanf(line, "%d %d %d", &polygons->indices[3*i],
                               &polygons->indices[3*i + 1],
                               &polygons->indices[3*i + 2]);
                }
                compute_polygon_normals(polygons);
        } else {
                fprintf(stderr, "input_dx: Unknown dx format..\n");
                return(0);
        }

        fclose(fp);
        return(OK);
}

int
input_dfs(char *file, File_formats *format, int *n_objects,
          object_struct ***object_list)
{
        FILE              *fp;
        int               i, hdr_size;
        char              dummy[256];
        polygons_struct   *polygons;
        Point             *point;
        object_struct     *object;
      
        /* prepare object and polygons */
        *n_objects = 0;
        object = create_object(POLYGONS);
        add_object_to_list(n_objects, object_list, object);

        polygons = get_polygons_ptr(object);
        initialize_polygons(polygons, WHITE, NULL);

        *format = ASCII_FORMAT;

        if ((fp = fopen(file, "r")) == 0) {
                fprintf(stderr, "input_dfs: Couldn't open file %s.\n", file);
                return(0);
        }

        fread(&dummy, sizeof(char), 12, fp);
        fprintf(stderr, "%s\n", dummy);
        fread(&hdr_size, 4, 1, fp);
        fread(&dummy, sizeof(char), 8, fp);

        fseek(fp, 24, 0);
        fread(&polygons->n_points, sizeof(int), 1, fp);
        fread(&polygons->n_items, sizeof(int), 1, fp);
        fprintf(stderr, "%s %d %d %d\n", dummy, hdr_size,
                polygons->n_points, polygons->n_items);
        fseek(fp, hdr_size, -1); 

        polygons->bintree = (bintree_struct_ptr) NULL;

        ALLOC(polygons->normals, polygons->n_points);
        ALLOC(polygons->end_indices, polygons->n_items);
        for (i = 0; i < polygons->n_items; i++)
                polygons->end_indices[i] = (i + 1) * 3;

        ALLOC(polygons->indices, polygons->end_indices[polygons->n_items-1]);
/* -------------------
        ALLOC(polygons->points, polygons->n_points);
        for (i = 0; i < polygons->n_points; i++) {
                fread(Point_x(point), 4, 1, fp);
                fread(&Point_y(point), 4, 1, fp);
                fread(&Point_z(point), 4, 1, fp);
                polygons->points[i] = point;
        }
----------------- */
        for (i = 0; i < 3*(polygons->n_items); i++)
                fread(&polygons->indices[i], 4, 1, fp);

        /* compute normals */
        compute_polygon_normals(polygons);
  
        fclose(fp);
        return(OK);
}
