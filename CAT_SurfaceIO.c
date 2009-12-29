/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_SurfaceIO.h"

#define TRIANGLE_FILE_MAGIC_NUMBER  16777214
#define QUAD_FILE_MAGIC_NUMBER      16777215
#define NEW_VERSION_MAGIC_NUMBER    16777215

Status
input_values_any_format(char *file, int *n_values, Real **values)
{
        Status status;

        if (filename_extension_matches(file,"txt"))
                status = input_texture_values(file, n_values, values);
        else
                status = input_freesurfer_curv(file, n_values, values);

        return(status);
}

Status
output_values_any_format(char *file, int n_values, void *values, int flag)
{
        Status status;
        Real *buffer;
        double *d;
        int i, *r;

        /* a workaround for Real stupidity */
        if (flag == TYPE_REAL) {
                buffer = (Real *) values;
        } else {
                buffer = (Real *) malloc(sizeof(Real) * n_values);

                if (flag == TYPE_DOUBLE) {
                        d = (double *) values;
                        for (i = 0; i < n_values; i++)
                                buffer[i] = (Real) d[i];
                } else if (flag == TYPE_INTEGER) {
                        r = (int *) values;
                        for (i = 0; i < n_values; i++)
                                buffer[i] = (Real) r[i];
                }
        }

        if (filename_extension_matches(file, "txt"))
                status = output_texture_values(file, ASCII_FORMAT, n_values,
                                               buffer);
        else
                status = output_freesurfer_curv(file, n_values, buffer);

        if (flag != TYPE_REAL)
                free(buffer);

        return(status);
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
swapFloat(float *n)
{
        char *by = (char *) n;
        char sw[4] = {by[3], by[2], by[1], by[0]};
  
        *n =* (float *) sw;
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

void
swapShort(short *n)
{
        char *by = (char *) n;
        char sw[2] = {by[1], by[0]};
  
        *n =* (short *) sw;
}

int
fwriteFloat(float f, FILE *fp)
{

#if (BYTE_ORDER == LITTLE_ENDIAN)
       swapFloat(&f);  
#endif
       return(fwrite(&f, sizeof(float), 1, fp));
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
fwrite2(int v, FILE *fp)
{
        short s ;

        if (v > 0x7fff)    /* don't let it overflow */
                v = 0x7fff ;
        else if (v < -0x7fff)
                v = -0x7fff ;
        s = (short)v;
#if (BYTE_ORDER == LITTLE_ENDIAN)
        swapShort(&s) ;
#endif
        return(fwrite(&s, 2, 1, fp));
}

int
fwrite3(int v, FILE *fp)
{
        int i = (v << 8);

#if (BYTE_ORDER == LITTLE_ENDIAN)
        swapInt(&i) ;
#endif
        return(fwrite(&i, 3, 1, fp));
}

int
fwrite4(float v, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
        swapFloat(&v);
#endif
        return(fwrite(&v, 4, 1, fp));
}

int
freadInt(FILE *fp)
{
        int temp;
        int count;
  
        count = fread(&temp, 4, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
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
#if (BYTE_ORDER == LITTLE_ENDIAN)
        swapDouble(&temp);
#endif
        return(temp);  
}

float
freadFloat(FILE *fp)
{
        float temp;
        int count;
  
        count = fread(&temp, 4, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
        swapFloat(&temp);
#endif
        return(temp);  
}

int
fread3(int *v, FILE *fp)
{
        int i = 0;
        int  ret ;

        ret = fread(&i, 3, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
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
                return(-1);
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
                                return(-1);
                        }
                }
    
                /* compute normals */
                compute_polygon_normals(polygons);
        } else {
                fprintf(stderr, "Wrong oogl format\n");
                return(-1);
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
                return(-1);
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
        return(OK);
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
                return(-1);
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
                fwriteFloat(Point_x(polygons->points[i]), fp);
                fwriteFloat(Point_y(polygons->points[i]), fp);
                fwriteFloat(Point_z(polygons->points[i]), fp);
        }
  
        /* write indices */
        for (i = 0; i < 3*(polygons->n_items); i++)
                fwriteInt(polygons->indices[i], fp);

        fclose(fp);
        return(OK);
}

int
output_freesurfer_curv(char *fname, int nvertices, Real *data)
{
        FILE *fp;
        Real f;
        int k;

        fp = fopen(fname,"wb");
        if(fp == NULL){
                fprintf(stderr, "output_freesurfer_curv: Couldn't open file %s.\n",
                        fname);
                return(-1);
        }
        fwrite3(NEW_VERSION_MAGIC_NUMBER, fp); 
        fwriteInt(nvertices, fp);
        fwriteInt(0, fp);
        fwriteInt(1, fp);
        for (k = 0; k < nvertices; k++) {
                f = data[k];
                fwriteFloat(f, fp);
        }
        fclose(fp);

        return(OK);
}

int
input_freesurfer(char *file, File_formats *format, int *n_objects,
                 object_struct ***object_list)
{
        FILE              *fp;
        int               i, magic;
        char              line[1024];
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
                return(-1);
        }

        /* read magic number for checking filetype */
        fread3(&magic, fp);
        if (magic == QUAD_FILE_MAGIC_NUMBER) {
                fprintf(stderr, "QUAD_FILE_MAGIC_NUMBER not yet prepared.\n");
                return(-1);
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
                        Point_x(point) = freadFloat(fp);
                        Point_y(point) = freadFloat(fp);
                        Point_z(point) = freadFloat(fp);
                        polygons->points[i] = point;
                }
                for (i = 0; i < 3*(polygons->n_items); i++) {
                        polygons->indices[i] = freadInt(fp);
                }
                /* compute normals */
                compute_polygon_normals(polygons);
        } else {
                fprintf(stderr, "input_freesurfer: Unknown magic identifier: %d.\n", magic);
                return(-1);
        }

        fclose(fp);
        return(OK);
}

int
input_freesurfer_curv(char *file, int *vnum, Real **input_values)
{
        FILE  *fp;
        int   i, magic, fnum, vals_per_vertex;
        Real  value;
      
        *vnum = 0;
  
        if ((fp = fopen(file, "rb")) == 0) {
                fprintf(stderr,
                        "input_freesurfer_curv: Couldn't open file %s.\n",
                        file);
                return(-1);
        }

        /* read magic number for checking filetype */
        fread3(&magic,fp);
  
        if (magic != NEW_VERSION_MAGIC_NUMBER) {
                fprintf(stderr, "MAGIC_NUMBER %d not yet prepared.\n", magic);
                return(-1);
        } else {
                /* read # of vertices and faces */
                *vnum = freadInt(fp);
                fnum = freadInt(fp);
                vals_per_vertex = freadInt(fp);
                if (vals_per_vertex != 1) {
                        fprintf(stderr, "Only one value per vertex allowed.\n");
                        return(-1);
                }
                ALLOC(*input_values, *vnum);
                for (i = 0; i < *vnum; i++) {
                        value = freadFloat(fp);
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
                return(-1);
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
                        return(-1);
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
                return(-1);
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
                return(-1);
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
