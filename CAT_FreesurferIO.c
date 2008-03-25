/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl/objects.h>

#define TRIANGLE_FILE_MAGIC_NUMBER  16777214
#define QUAD_FILE_MAGIC_NUMBER  16777215

public int output_freesurfer_bin(
    char *fname,
    File_formats   format,
    int            n_objects,
    object_struct  *object_list[])
{
    FILE    *fp;
    int     i;
    unsigned char   buffer;
    polygons_struct     *polygons;

    if((fp = fopen(fname, "w")) == 0) {
        fprintf(stderr, "output_freesurfer_bin: Couldn't open file %s.\n", fname);
        return(0);
    }

    polygons = get_polygons_ptr(object_list[0]);

    /* write 3 byte magic number for triangle */
    buffer = 255;
    fwrite(&buffer, 1, 1,fp);
    fwrite(&buffer, 1, 1,fp);
    buffer = 254;
    fwrite(&buffer, 1, 1,fp);

    output_newline( fp ); 
    output_newline( fp ); 

    /* # of vertices and faces */
    fwrite(&polygons->n_points, 4, 1, fp);
    fwrite(&polygons->n_items, 4, 1, fp);

    /* write points */
    for (i = 0; i < (polygons->n_points); i++) {
        fwrite(&Point_x(polygons->points[i]), 4, 1, fp);
        fwrite(&Point_y(polygons->points[i]), 4, 1, fp);
        fwrite(&Point_z(polygons->points[i]), 4, 1, fp);
    }
    
    /* write indices */
    for (i = 0; i < 3*(polygons->n_items); i++)
        fwrite(&polygons->indices[i], 4, 1, fp);

    fclose(fp);
}

public int input_freesurfer_bin(
    char *fname,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list)
{
    FILE    *fp;
    int     i, magic, x, y, z;
    char    line[256];
    Real *vertices, *faces;
    polygons_struct     *polygons;
    Point   point;
    object_struct  *object;
            
    /* prepare object and polygons */
    *n_objects = 0;
    object = create_object( POLYGONS );
    add_object_to_list( n_objects, object_list, object );
    polygons = get_polygons_ptr(object);
    initialize_polygons( polygons, WHITE, NULL );
    *format = ASCII_FORMAT;

    if((fp = fopen(fname, "r")) == 0) {
        fprintf(stderr, "output_freesurfer_bin: Couldn't open file %s.\n", fname);
        return(0);
    }

    /* read magic number for checking filetype */
    magic = fread3(fp);
    
    if( magic == QUAD_FILE_MAGIC_NUMBER) {
        fprintf(stderr, "QUAD_FILE_MAGIC_NUMBER not yet prepared %s.\n");
        return(0);
    } else if( magic == TRIANGLE_FILE_MAGIC_NUMBER) {
        fgets(line, 1024, fp);
        fgets(line, 1024, fp);
        /* read # of vertices and faces */
        fread(&polygons->n_points, 4, 1, fp);
        fread(&polygons->n_items, 4, 1, fp);
        ALLOC( polygons->points, polygons->n_points );
        ALLOC( polygons->normals, polygons->n_points );
        ALLOC( polygons->end_indices, polygons->n_items );
        polygons->bintree = (bintree_struct_ptr) NULL;
        for (i = 0; i < (polygons->n_items); i++)
            polygons->end_indices[i] = (i+1) * 3;
        ALLOC( polygons->indices, polygons->end_indices[polygons->n_items-1] );
        for (i = 0; i < (polygons->n_points); i++) {
            fread(&Point_x(point), 4, 1, fp);
            fread(&Point_y(point), 4, 1, fp);
            fread(&Point_z(point), 4, 1, fp);
            polygons->points[i] = point;
        }
        for (i = 0; i < 3*(polygons->n_items); i++)
            fread(&polygons->indices[i], 4, 1, fp);
        /* compute normals */
        compute_polygon_normals( polygons );
    } else {
        fprintf(stderr, "read_freesurfer_bin: Unknown magic identifier: %d.\n", magic);
        return(0);
    }
    
    fclose(fp);
    return(OK);
}

int fread3(FILE *fp)
{
    unsigned char b1,b2,b3;
    int count;
    
    count = fread (&b1, 1, 1, fp);
    count = fread (&b2, 1, 1, fp);
    count = fread (&b3, 1, 1, fp);
    return((b1 << 16) + (b2 << 8) + b3) ;

}