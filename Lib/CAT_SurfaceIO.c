/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 * Most of the basic input functions are from fio.c, machine.c,
 * annotations.c and gifti_local.c from freesurfer
 *
 */

#include "CAT_SurfaceIO.h"

Status
bicpl_to_facevertexdata(polygons_struct *polygons, double **faces, double **vertices)
{
    int i, j;
    Status status;

    status = OK;    
    for (i = 0; i < polygons->n_points; i++) {
         (*vertices)[i]            = Point_x(polygons->points[i]);
         (*vertices)[i+polygons->n_points]   = Point_y(polygons->points[i]);
         (*vertices)[i+2*polygons->n_points] = Point_z(polygons->points[i]);
    }
    
    for (i = 0; i < polygons->n_items; i++) {
        if (GET_OBJECT_SIZE(*polygons, i) != 3) {
            status = ERROR;
            return(status);
        }
           
        /* add "1" to faces (for matlab arrays) */
        for (j = 0; j < 3; j++) {
            (*faces)[i + j*polygons->n_items] = polygons->indices[POINT_INDEX(
                          polygons->end_indices, i, j)] + 1;
        }
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
fwriteInt(int v, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
    swapInt(&v);
#endif
    return(fwrite(&v, sizeof(int), 1, fp));
}

int
fwrite3(int v, FILE *fp)
{
    int i = (v << 8);

#if (BYTE_ORDER == LITTLE_ENDIAN)
    swapInt(&i);
#endif
    return(fwrite(&i, 3, 1, fp));
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
    int  ret;

    ret = fread(&i, 3, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
    swapInt(&i);
#endif
    *v = ((i >> 8) & 0xffffff);
    return(ret);
}

static giiDataArray* 
gifti_alloc_and_add_darray (gifti_image* image)
{
    if (!image) {
        fprintf (stderr,"** gifti_alloc_and_add_darray: NULL image\n");
        return NULL;
    }

    /* Try to add an empty array. */
    if (gifti_add_empty_darray(image,1)) {
        fprintf (stderr,"** gifti_alloc_and_add_darray: gifti_add_empty_darray "
                         "failed\n");
        return NULL;
    }

    /* Return the array we just allocated. */
    return image->darray[image->numDA-1];
}

static double 
gifti_get_DA_value_2D (giiDataArray* da, int row, int col)
{
    int dim0_index, dim1_index;
    int dims_0=0, dims_1=0;

    if (!da || !da->data) {
        fprintf (stderr,"** gifti_get_DA_value_2D, invalid params: data=%p\n", da);
        exit(1);
    }

    if (da->num_dim == 1) {
        // support for using this routine to read 1D data, under one condition...
        if (col != 0)
        {
            fprintf (stderr,"** gifti_get_DA_value_2D, array dim is 1 "
                "but trying to access 2D data element (col=%d)\n",col);
            exit(1);
        }
        dims_0 = da->dims[0];
        dims_1 = 1; // 1D data
    }
    else if (da->num_dim != 2) {
        fprintf (stderr,"** gifti_get_DA_value_2D, array dim is %d\n", da->num_dim);
        exit(1);
    }
    else {
        dims_0 = da->dims[0];
        dims_1 = da->dims[1];
    }

    /* Get the dim0 and dims[1] indices based on our order. */
    if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord) {
        dim0_index = row;
        dim1_index = col;
    }
    else if (GIFTI_IND_ORD_COL_MAJOR == da->ind_ord) {
        // NJS NOTE: notice that order is treated as row/col, so that the
        // calling sequence can just assume row major
        dim0_index = row;//col;
        dim1_index = col;//row;
    }
    else {
        fprintf (stderr,"** gifti_get_DA_value_2D, unknown ind_ord: %d\n", da->ind_ord);
        exit(1);
    }
    if (da->num_dim == 1) /* support for using this routine to read 1D data */ {
        dim0_index = row;
        dim1_index = col;
    }

    /* Check the indices. */
    if (dim0_index < 0 || dim0_index >= dims_0 || dim1_index < 0 || dim1_index >= dims_1) {
        fprintf(stderr,"** gifti_get_DA_value_2D, invalid params: "
            "dim0_index=%d (max=%d), dim1_index=%d (max=%d)\n",
            dim0_index, dims_0, dim1_index, dims_1);
        exit(1);
    }

    /* Switch on the data type and return the appropriate
         element. Indexing depends on the data order. */
    switch (da->datatype) {
    default :
        fprintf(stderr,"** gifti_get_DA_value_2D, unsupported type %d-"
            "unknown, or can't convert to double\n",da->datatype);
        exit(1);
    case NIFTI_TYPE_UINT8: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((unsigned char*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((unsigned char*)
            (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_INT16: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((short*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((short*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_INT32: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((int*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((int*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_FLOAT32: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((float*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((float*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_FLOAT64: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((double*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((double*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_INT8: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((char*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((char*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_UINT16: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((unsigned short*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((unsigned short*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    case NIFTI_TYPE_UINT32: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            return (double)*((unsigned int*) (da->data) + (dim0_index*dims_1) + dim1_index);
        else
            return (double)*((unsigned int*) (da->data) + dim0_index + (dim1_index*dims_0));
        break;
    }
    }

    exit(1);
}


/*
 *
 */
static void 
gifti_set_DA_value_2D (giiDataArray* da, int row, int col, double value)
{
    int dim0_index, dim1_index;
    int dims_0=0, dims_1=0;

    if (!da || !da->data) {
        fprintf (stderr,"** gifti_set_DA_value_2D, invalid params: data=%p\n", da);
        exit(1);
    }

    if (da->num_dim == 1) {
        // support for using this routine to write 1D data, under one condition...
        if (col != 0) {
            fprintf (stderr,"** gifti_set_DA_value_2D, array dim is 1 "
                "but trying to access 2D data element (col=%d)\n",col);
            exit(1);
        }
        dims_0 = da->dims[0];
        dims_1 = 1; // 1D data
    }
    else if (da->num_dim != 2) {
        fprintf (stderr,"** gifti_set_DA_value_2D, array dim is %d\n", da->num_dim);
        exit(1);
    }
    else {
        dims_0 = da->dims[0];
        dims_1 = da->dims[1];
    }

    /* Get the dim0 and dims[1] indices based on our order. */
    if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord) {
        dim0_index = row;
        dim1_index = col;
    }
    else {
        dim0_index = col;
        dim1_index = row;
    }
    if (da->num_dim == 1) /* support for using this routine to read 1D data */ {
        dim0_index = row;
        dim1_index = col;
    }

    /* Check the indices. */
    if (dim0_index < 0 || dim0_index >= dims_0 || dim1_index < 0 || dim1_index >= dims_1) {
        fprintf(stderr,"** gifti_set_DA_value_2D, invalid params: "
            "dim0_index=%d (max=%d), dim1_index=%d (max=%d)\n",
            dim0_index, dims_0, dim1_index, dims_1);
        return;
    }

    /* Switch on the data type and write the appropriate
         element. Indexing depends on the data order. */
    switch (da->datatype) {
    default :
        fprintf(stderr,"** gifti_set_DA_value_2D, unsupported type %d-"
            "unknown, or can't convert to double\n",da->datatype);
        return;
    case NIFTI_TYPE_UINT8: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((unsigned char*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (unsigned char)value;
        else
            *((unsigned char*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (unsigned char)value;
        break;
    }
    case NIFTI_TYPE_INT16: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((short*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (short)value;
        else
            *((short*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (short)value;
        break;
    }
    case NIFTI_TYPE_INT32: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((int*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (int)value;
        else
            *((int*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (int)value;
        break;
    }
    case NIFTI_TYPE_FLOAT32: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((float*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (float)value;
        else
            *((float*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (float)value;
        break;
    }
    case NIFTI_TYPE_INT8: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((char*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (char)value;
        else
            *((char*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (char)value;
        break;
    }
    case NIFTI_TYPE_UINT16: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((unsigned short*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (unsigned short)value;
        else
            *((unsigned short*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (unsigned short)value;
        break;
    }
    case NIFTI_TYPE_UINT32: {
        if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
            *((unsigned int*)(da->data) + (dim0_index*dims_1) + dim1_index) =
                (unsigned int)value;
        else
            *((unsigned int*)(da->data) + dim0_index + (dim1_index*dims_0)) =
                (unsigned int)value;
        break;
    }
    }

    return;
}

int
input_oogl(char *file, File_formats *format, int *n_objects,
       object_struct ***object_list)
{
    FILE        *fp;
    int         i, f, n_edges, dummy;
    char        line[1000];
    polygons_struct   *polygons;
    Point       point;
    object_struct   *object;

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
output_gifti(char *fname, File_formats format, int n_objects,
          object_struct *object_list[], double *values)
{

    int j,k;
    polygons_struct   *polygons;

/* handle external binary .dat output */
    int use_dat = filename_extension_matches(fname, "dat");
    char header_fname[1024];
    char *out_name = fname;
    if (use_dat) {
        strncpy(header_fname, fname, sizeof(header_fname)-1);
        header_fname[sizeof(header_fname)-1] = '\0';
        char *dot = strrchr(header_fname, '.');
        if (dot) strcpy(dot, ".gii");
        out_name = header_fname;
        remove(fname);
    }

    gifti_image* image = (gifti_image *)calloc(1,sizeof(gifti_image));
    if (NULL == image) {
        fprintf (stderr,"output_gifti: couldn't allocate image\n");
        return(-1);
    }

    image->version = (char *) calloc(strlen(GIFTI_XML_VERSION)+1,sizeof(char));;
    strcpy(image->version,GIFTI_XML_VERSION);

    gifti_add_to_meta( &image->meta, "Name", out_name, 1 );

    giiDataArray* coords = gifti_alloc_and_add_darray (image);
    if (NULL == coords) {
        fprintf (stderr,"output_gifti: couldn't allocate giiDataArray\n");
        gifti_free_image (image);
        return(-1);
    }

    polygons = get_polygons_ptr(object_list[0]);

    /* Set its attributes. */
    coords->intent   = NIFTI_INTENT_POINTSET;
    coords->datatype = NIFTI_TYPE_FLOAT32;
    coords->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;
    coords->num_dim  = 2;
    coords->dims[0]  = polygons->n_points; /* In highest first, dim0 = rows */
    coords->dims[1]  = 3;       /* In highest first, dim1 = cols */
    coords->encoding = GIFTI_ENCODING_B64BIN;
#if (BYTE_ORDER == LITTLE_ENDIAN)
    coords->endian   = GIFTI_ENDIAN_LITTLE;
#else
    coords->endian   = GIFTI_ENDIAN_BIG;
#endif

    gifti_add_empty_CS( coords );
    coords->coordsys[0]->dataspace  = (char *) calloc(strlen("NIFTI_XFORM_UNKNOWN")  +1,sizeof(char));;
    coords->coordsys[0]->xformspace = (char *) calloc(strlen("NIFTI_XFORM_TALAIRACH")+1,sizeof(char));;

    strcpy(coords->coordsys[0]->dataspace, "NIFTI_XFORM_UNKNOWN");
    strcpy(coords->coordsys[0]->xformspace,"NIFTI_XFORM_TALAIRACH");
    int r,c;
    for (r=1; r <= 4; r++)
        for (c=1; c <= 4; c++)
            if (r==c) coords->coordsys[0]->xform[r-1][c-1] = 1.0;
            else    coords->coordsys[0]->xform[r-1][c-1] = 0.0;

    coords->nvals = gifti_darray_nvals (coords);
    gifti_datatype_sizes (coords->datatype, &coords->nbyper, NULL);

    /* Allocate the data array. */
    coords->data = NULL;
    coords->data = (void*) calloc (coords->nvals, coords->nbyper);
    if (NULL == coords->data) {
        fprintf (stderr,"output_gifti: couldn't allocate coords data of "
            "length %d, element size %d\n",
        (int)coords->nvals, coords->nbyper);
    gifti_free_image (image);
    return -1;
    }

    /* Copy in all our data. */
    int vertex_index;
    for (vertex_index = 0; vertex_index < polygons->n_points; vertex_index++) {
        gifti_set_DA_value_2D (coords, vertex_index, 0,
               Point_x(polygons->points[vertex_index]));
        gifti_set_DA_value_2D (coords, vertex_index, 1,
               Point_y(polygons->points[vertex_index]));
        gifti_set_DA_value_2D (coords, vertex_index, 2,
               Point_z(polygons->points[vertex_index]));
    }

    /* 
     * Faces
     */
    giiDataArray* faces = gifti_alloc_and_add_darray (image);
    if (NULL == faces) {
        fprintf (stderr,"output_gifti: couldn't allocate giiDataArray\n");
        gifti_free_image (image);
        return -1;
    }

    int numFaces = polygons->n_items;

    /* Set its attributes. */
    faces->intent = NIFTI_INTENT_TRIANGLE;
    faces->datatype = NIFTI_TYPE_INT32;
    faces->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;
    faces->num_dim  = 2;
    faces->dims[0]  = numFaces;    /* In highest first, dim0 = rows */
    faces->dims[1]  = 3;       /* In highest first, dim1 = cols */
    faces->encoding = GIFTI_ENCODING_B64GZ;
#if (BYTE_ORDER == LITTLE_ENDIAN)
    faces->endian = GIFTI_ENDIAN_LITTLE;
#else
    faces->endian = GIFTI_ENDIAN_BIG;
#endif
    faces->coordsys = NULL;
    faces->nvals  = gifti_darray_nvals (faces);
    gifti_datatype_sizes (faces->datatype, &faces->nbyper, NULL);

    /* Allocate the data array. */
    faces->data = NULL;
    faces->data = (void*) calloc (faces->nvals, faces->nbyper);
    if (NULL == faces->data) {
        fprintf (stderr,"output_gifti: couldn't allocate faces data of "
            "length %d, element size %d\n",
        (int)faces->nvals, faces->nbyper);
        gifti_free_image (image);
        return -1;
    }


    /* Copy in all our face data */
    int face_index;
    for (face_index = 0; face_index < polygons->n_items; face_index++)
        for (j = 0; j < 3; j++)
            gifti_set_DA_value_2D (faces, face_index, j,
                 polygons->indices[POINT_INDEX(polygons->end_indices, face_index, j)]);

    /* standard meta data for surfaces */
    if (fname) {
        const char *primary=NULL, *secondary=NULL, *geotype=NULL;
        char *name = fname;
        char *topotype="Closed";
        if (strstr(name, "lh.")) primary = "CortexLeft";
        if (strstr(name, "rh.")) primary = "CortexRight";
        if (strstr(name, ".orig"))     secondary = "GrayWhite";
        if (strstr(name, ".smoothwm")) secondary = "GrayWhite";
        if (strstr(name, ".white"))    secondary = "GrayWhite";
        if (strstr(name, ".central"))  secondary = "Central (Layer 4)";
        if (strstr(name, ".graymid"))  secondary = "MidThickness";
        if (strstr(name, ".gray"))     secondary = "Pial";
        if (strstr(name, ".pial"))     secondary = "Pial";
        if (strstr(name, ".orig"))     geotype = "Reconstruction";
        if (strstr(name, ".smoothwm")) geotype = "Anatomical";
        if (strstr(name, ".white"))    geotype = "Anatomical";
        if (strstr(name, ".central"))  geotype = "Anatomical";
        if (strstr(name, ".gray"))     geotype = "Anatomical";
        if (strstr(name, ".graymid"))  geotype = "Anatomical";
        if (strstr(name, ".pial"))     geotype = "Anatomical";
        if (strstr(name, ".inflated")) geotype = "Inflated";
        if (strstr(name, ".sphere"))   geotype = "Sphere";
        if (strstr(name, ".qsphere"))  geotype = "Sphere";
        if (strstr(name,"pial-outer")) geotype = "Hull";
    
        if (primary) gifti_add_to_meta( &coords->meta,
                    "AnatomicalStructurePrimary",
                    primary,
                    1 );
        if (secondary) gifti_add_to_meta( &coords->meta,
                    "AnatomicalStructureSecondary",
                    secondary,
                    1 );
        if (geotype) gifti_add_to_meta( &coords->meta,
                    "GeometricType",
                    geotype,
                    1 );
        gifti_add_to_meta( &faces->meta, "TopologicalType", topotype, 1 );
        gifti_add_to_meta( &coords->meta, "Name", out_name, 1 );
        gifti_add_to_meta( &faces->meta, "Name", out_name, 1 );
    }



    /* 
     * Shape (textures)
     */
    if (values != NULL) {

        giiDataArray* shape = gifti_alloc_and_add_darray (image);
        if (NULL == shape) {
            fprintf (stderr,"output_gifti_curv: couldn't allocate giiDataArray\n");
            gifti_free_image (image);
            return(-1);
        }
        
        /* Set its attributes. */
        shape->intent = NIFTI_INTENT_SHAPE;
        shape->datatype = NIFTI_TYPE_FLOAT32;
        shape->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
        shape->num_dim = 1;
        shape->dims[0] = polygons->n_points;
        shape->dims[1] = 0;
        shape->encoding = GIFTI_ENCODING_B64GZ; 
#if (BYTE_ORDER == LITTLE_ENDIAN)
        shape->endian = GIFTI_ENDIAN_LITTLE;
#else
        shape->endian = GIFTI_ENDIAN_BIG;
#endif
        shape->coordsys = NULL;
        shape->nvals = gifti_darray_nvals (shape);
        gifti_datatype_sizes (shape->datatype, &shape->nbyper, NULL);

        /* Allocate the data array. */
        shape->data = NULL;
        shape->data = (void*) calloc (shape->nvals, shape->nbyper);
        if (NULL == shape->data) {
            fprintf (stderr,"output_gifti_curv: couldn't allocate shape data of "
                "length %d, element size %d\n", (int)shape->nvals,shape->nbyper);
            gifti_free_image (image);
            return(-1);
        }

        /* Copy in all our data. */
        for (k = 0; k < polygons->n_points; k++)
            gifti_set_DA_value_2D (shape, k, 0, values[k]);
    }


    if (use_dat)
        gifti_set_extern_filelist(image, 1, &fname);

    /* check for compliance */
    int valid = gifti_valid_gifti_image (image, 1);
    if (valid == 0) {
        fprintf (stderr,"output_gifti_curv: GIFTI file %s is invalid!\n", fname);
        gifti_free_image (image);
        return(-1);
    }

    /* Write the file. */
    if (gifti_write_image (image, out_name, 1)) {
        fprintf (stderr,"output_gifti_curv: couldn't write image\n");
        gifti_free_image (image);
        return(-1);
    }

    gifti_free_image (image);

    return(OK);
}

int
output_gifti_curv(char *fname, int nvertices, double *data)
{

    int k;

    /* handle external binary .dat output */
    int use_dat = filename_extension_matches(fname, "dat");
    char header_fname[1024];
    char *out_name = fname;
    if (use_dat) {
        strncpy(header_fname, fname, sizeof(header_fname)-1);
        header_fname[sizeof(header_fname)-1] = '\0';
        char *dot = strrchr(header_fname, '.');
        if (dot) strcpy(dot, ".gii");
        out_name = header_fname;
        remove(fname);
    }

    gifti_image* image = (gifti_image *)calloc(1,sizeof(gifti_image));
    if (NULL == image) {
        fprintf (stderr,"output_gifti_curv: couldn't allocate image\n");
        return(-1);
    }

    image->version = (char *) calloc(strlen(GIFTI_XML_VERSION)+1,sizeof(char));;
    strcpy(image->version,GIFTI_XML_VERSION);

    gifti_add_to_meta( &image->meta, "Name", out_name, 1 );

    giiDataArray* shape = gifti_alloc_and_add_darray (image);
    if (NULL == shape) {
        fprintf (stderr,"output_gifti_curv: couldn't allocate giiDataArray\n");
        gifti_free_image (image);
        return(-1);
    }

    /* Set its attributes. */
    shape->intent = NIFTI_INTENT_SHAPE;
    shape->datatype = NIFTI_TYPE_FLOAT32;
    shape->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;
    shape->num_dim  = 1;
    shape->dims[0]  = nvertices;
    shape->dims[1]  = 0;
    shape->encoding = GIFTI_ENCODING_B64BIN; 
#if (BYTE_ORDER == LITTLE_ENDIAN)
    shape->endian = GIFTI_ENDIAN_LITTLE;
#else
    shape->endian = GIFTI_ENDIAN_BIG;
#endif
    shape->coordsys = NULL;
    shape->nvals  = gifti_darray_nvals (shape);
    gifti_datatype_sizes (shape->datatype, &shape->nbyper, NULL);

    /* include some metadata describing this shape */
    gifti_add_to_meta( &shape->meta, "Name", out_name, 1 );
    char *meta=NULL;
    if (strstr(fname, ".thickness")) meta = "Thickness";
    if (strstr(fname, ".curv"))    meta = "CurvatureRadial";
    if (strstr(fname, ".sulc"))    meta = "SulcalDepth";
    if (strstr(fname, ".area"))    meta = "Area";
    if (strstr(fname, ".volume"))  meta = "Volume";
    if (strstr(fname, ".jacobian"))  meta = "Jacobian";
    if (meta) gifti_add_to_meta( &shape->meta, "ShapeDataType", meta, 1 );

    /* Allocate the data array. */
    shape->data = NULL;
    shape->data = (void*) calloc (shape->nvals, shape->nbyper);
    if (NULL == shape->data) {
        fprintf (stderr,"output_gifti_curv: couldn't allocate shape data of "
            "length %d, element size %d\n", (int)shape->nvals,shape->nbyper);
        gifti_free_image (image);
        return(-1);
    }

    /* Copy in all our data. */
    for (k = 0; k < nvertices; k++)
        gifti_set_DA_value_2D (shape, k, 0, data[k]);

    if (use_dat)
        gifti_set_extern_filelist(image, 1, &fname);

    /* check for compliance */
    int valid = gifti_valid_gifti_image (image, 1);
    if (valid == 0) {
        fprintf (stderr,"output_gifti_curv: GIFTI file %s is invalid!\n", fname);
        gifti_free_image (image);
        return(-1);
    }

    /* Write the file. */
    if (gifti_write_image (image, out_name, 1)) {
        fprintf (stderr,"output_gifti_curv: couldn't write image\n");
        gifti_free_image (image);
        return(-1);
    }

    gifti_free_image (image);

    return(OK);

}

int
input_gifti(char *file, File_formats *format, int *n_objects,
         object_struct ***object_list)
{

    int         i, j, k, valid, numDA;
    polygons_struct   *polygons;
    Point       point;
    object_struct   *object;
  
    gifti_image* image = gifti_read_image (file, 1);
    if (NULL == image) {
        fprintf (stderr,"input_gifti: cannot read image\n");
        return(-1);
    }

    valid = gifti_valid_gifti_image (image, 1);
    if (valid == 0) {
        fprintf (stderr,"input_gifti: GIFTI file %s is invalid!\n", file);
        gifti_free_image (image);
        return(-1);
    }

    /* prepare object and polygons */
    *n_objects = 0;
    object = create_object(POLYGONS);
    add_object_to_list(n_objects, object_list, object);
    polygons = get_polygons_ptr(object);
    initialize_polygons(polygons, WHITE, NULL);
    *format = ASCII_FORMAT;

    giiDataArray* coords = NULL;
    giiDataArray* faces  = NULL;

    for (numDA = 0; numDA < image->numDA; numDA++) {
        if (image->darray[numDA]->intent == NIFTI_INTENT_POINTSET) {
            coords = image->darray[numDA];
        }
        else if (image->darray[numDA]->intent == NIFTI_INTENT_TRIANGLE) {
            faces = image->darray[numDA];
        }
    }

    if (coords && faces) {
    
        /* Check the number of vertices and faces. */
        long long num_vertices = 0;
        long long num_cols = 0;
        if (coords->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) // RowMajorOrder
            gifti_DA_rows_cols (coords, &num_vertices, &num_cols);
        else // ColumnMajorOrder
            gifti_DA_rows_cols (coords, &num_cols, &num_vertices);

        if (num_vertices <= 0 || num_cols != 3) {
            fprintf (stderr,"input_gifti: malformed coords data array in file "
                "%s: num_vertices=%d num_cols=%d\n",
                file, (int)num_vertices, (int)num_cols);
            gifti_free_image (image);
            return(-1);
        }

        long long num_faces = 0;
        num_cols = 0;
        if (faces->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) // RowMajorOrder
            gifti_DA_rows_cols (faces, &num_faces, &num_cols);
        else // ColumnMajorOrder
            gifti_DA_rows_cols (faces, &num_cols, &num_faces);

        if (num_faces <= 0 || num_cols != 3) {
            fprintf (stderr,"mrisReadGIFTIfile: malformed faces data array in file "
                "%s: num_faces=%d num_cols=%d\n",
                file, (int)num_faces, (int)num_cols);
            gifti_free_image (image);
            return(-1);
        }
        
        polygons->n_points = num_vertices;
        polygons->n_items = num_faces;

        ALLOC(polygons->points, polygons->n_points);
        ALLOC(polygons->normals, polygons->n_points);
        ALLOC(polygons->end_indices, polygons->n_items);
        polygons->bintree = (bintree_struct_ptr) NULL;
        for (i = 0; i < polygons->n_items; i++)
            polygons->end_indices[i] = (i + 1) * 3;
        ALLOC(polygons->indices,
            polygons->end_indices[polygons->n_items-1]);

        int vertex_index;
        for (vertex_index = 0; vertex_index < polygons->n_points; vertex_index++) {
            Point_x(polygons->points[vertex_index]) = (float) gifti_get_DA_value_2D (coords, vertex_index, 0);
            Point_y(polygons->points[vertex_index]) = (float) gifti_get_DA_value_2D (coords, vertex_index, 1);
            Point_z(polygons->points[vertex_index]) = (float) gifti_get_DA_value_2D (coords, vertex_index, 2);
        }

        int face_index;
        for (face_index = 0; face_index < polygons->n_items; face_index++)
            for (j = 0; j < 3; j++)
                polygons->indices[POINT_INDEX(polygons->end_indices, face_index, j)] = gifti_get_DA_value_2D (faces, face_index, j);

        /* compute normals */
        compute_polygon_normals(polygons);
    } else {
        fprintf (stderr,"input_gifti: GIFTI file %s does not contain vertices and faces!\n", file);
        gifti_free_image (image);
        return(-1);
    }

    for (numDA = 0; numDA < image->numDA; numDA++) {
        giiDataArray* darray = image->darray[numDA];
        
        /* did these already */
        if ((darray->intent == NIFTI_INTENT_POINTSET) ||
          (darray->intent == NIFTI_INTENT_TRIANGLE)) continue;

        if (darray->intent == NIFTI_INTENT_SHAPE) 
            fprintf(stderr, "input_gifti: Reading of shape data not yet implemented.\n");
    }


    return(OK);

}

int
input_gifti_curv(char *file, int *vnum, double **input_values)
{
    int k, valid, numDA;
  
    gifti_image* image = gifti_read_image (file, 1);
    if (NULL == image) {
        fprintf (stderr,"input_gifti_curv: cannot read image\n");
        return(-1);
    }

    valid = gifti_valid_gifti_image (image, 1);
    if (valid == 0) {
        fprintf (stderr,"input_gifti_curv: GIFTI file %s is invalid!\n", file);
        gifti_free_image (image);
        return(-1);
    }

    *vnum = image->darray[0]->dims[0];
    ALLOC(*input_values, image->darray[0]->dims[0]);

    for (k = 0; k < image->darray[0]->dims[0]; k++)
        (*input_values)[k] = (double) gifti_get_DA_value_2D (image->darray[0], k, 0);

    return(OK);
}

int
output_freesurfer(char *file, File_formats format, int n_objects,
          object_struct *object_list[])
{
    FILE        *fp;
    int         i,n;
    unsigned char   buffer;
    polygons_struct   *polygons;

    if ((fp = fopen(file, "w")) == 0) {
        fprintf(stderr, "output_freesurfer: Couldn't open file %s.\n",
            file);
        return(-1);
    }
    
    polygons = get_polygons_ptr(object_list[0]);

    /* write 3 byte magic number for triangle */
    fwrite3(TRIANGLE_FILE_MAGIC_NUMBER,fp);
    fprintf(fp, "created with CAT\n");

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
    for (i = 0; i < polygons->n_items;i++) 
        for (n = 0; n < 3; n++) 
            fwriteInt(polygons->indices[POINT_INDEX(polygons->end_indices, i, n)],fp);

    fclose(fp);
    return(OK);
}

int
output_freesurfer_curv(char *fname, int nvertices, double *data)
{
    FILE *fp;
    double f;
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
    FILE        *fp;
    int         i, magic;
    char        line[1024];
    polygons_struct   *polygons;
    Point       point;
    object_struct   *object;
    
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
        for (i = 0; i < 3*(polygons->n_items); i++) 
            polygons->indices[i] = freadInt(fp);
        
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
input_freesurfer_curv(char *file, int *vnum, double **input_values)
{
    FILE  *fp;
    int   i, magic, fnum, vals_per_vertex;
    double  value;
    
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
    FILE        *fp;
    int         i, pos;
    char        line[1000];
    char        ch;
    polygons_struct   *polygons;
    Point       point;
    object_struct   *object;
    Status        status;
    
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
    FILE        *fp;
    int         i, hdr_size, mdoffset, pdoffset, nStrips, stripSize, normals;
    int         uvStart, vcoffset, labelOffset, vertexAttributes, tmp;
    char        dummy[256], pad2[124];
    float       tmp2;
    polygons_struct   *polygons;
    Point       point;
    object_struct   *object;
    
    fprintf(stderr,"Input of DFS data not working.\n");      
    return(-1);
    
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
    fread(&hdr_size, 4, 1, fp);
    fread(&mdoffset, 4, 1, fp);
    fread(&pdoffset, 4, 1, fp);
    fread(&polygons->n_points, 4, 1, fp); /* Number of triangles */
    fread(&polygons->n_items, 4, 1, fp);  /* Number of vertices */
    fprintf(stderr, "%s %d %d %d\n", dummy, hdr_size,
        polygons->n_points, polygons->n_items);
    fread(&nStrips, 4, 1, fp);
    fread(&stripSize, 4, 1, fp);
    fread(&normals, 4, 1, fp);
    fread(&uvStart, 4, 1, fp);
    fread(&vcoffset, 4, 1, fp);
    fread(&labelOffset, 4, 1, fp);
    fread(&vertexAttributes, 4, 1, fp);

    fread(&pad2, sizeof(char), 124, fp);

    fseek(fp, hdr_size, -1); 

    polygons->bintree = (bintree_struct_ptr) NULL;

    ALLOC(polygons->points, polygons->n_points);
    ALLOC(polygons->normals, polygons->n_points);
    ALLOC(polygons->end_indices, polygons->n_items);

    for (i = 0; i < polygons->n_items; i++)
        polygons->end_indices[i] = (i + 1) * 3;

    ALLOC(polygons->indices, polygons->end_indices[polygons->n_items-1]);

    for (i = 0; i < 3*(polygons->n_items); i++) {
        fread(&tmp, 4, 1, fp);
        polygons->indices[i] = tmp;
    }
    
    for (i = 0; i < polygons->n_points; i++) {
        fread(&tmp2, 4, 1, fp); Point_x(point) = tmp2;
        fread(&tmp2, 4, 1, fp); Point_y(point) = tmp2;
        fread(&tmp2, 4, 1, fp); Point_z(point) = tmp2;
        polygons->points[i] = point;
    }

    /* compute normals */
    compute_polygon_normals(polygons);
  
    fclose(fp);
    return(OK);
}

double *
read_pgm(char *file, int *nx, int *ny)
{
    FILE      *fp;
    int       i, size, max_val;
    char      line[1000];
    unsigned char *data_char;
    double      *data;
      
    if ((fp = fopen(file, "r")) == 0) {
        fprintf(stderr, "read_pgm: Couldn't open file %s.\n", file);
        *nx = *ny = 0;
        return(NULL);
    }
  
    /* read PGM header */
    fgets(line, 256, fp);
    if (strncmp(line, "P5", 2)) {
        fprintf(stderr, "read_pgm: Can only read PGM binary files.\n");
        *nx = *ny = 0;
        return(NULL);
    }

    fgets(line, 256, fp);
    while (line[0] == '#')
        fgets(line, 256, fp);
    while (line[0] == '\n')
        fgets(line, 256, fp);
    sscanf(line, "%d %d", nx, ny);
    fgets(line, 256, fp);
    sscanf(line, "%d", &max_val);

    size = (*nx) * (*ny);
    ALLOC(data_char, size);

    if (fread(data_char, sizeof(unsigned char), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        *nx = *ny = 0;
        return(NULL);
    }

    ALLOC(data, size);   
    for (i = 0; i < size; i++) {
        data[i] = (double) data_char[i];
    }

    FREE(data_char);
    fclose(fp);
    return(data);
}

int
write_pgm(char *file, double *data, int nx, int ny)
{
    FILE  *fp;
    double    min_value, max_value, scale, offset;
    int   i;
    char  *data_char;
  
    ALLOC(data_char, nx*ny);
  
    min_value = 1e15;
    max_value = -1e15;

    /* calculate min/max */
    for (i = 0; i < nx*ny; i++) {
        if (data[i] < min_value)
            min_value = data[i];
        if (data[i] > max_value)
            max_value = data[i];
    }

    offset = min_value;
    scale = 255 / (max_value - min_value);
    /* scale image */
    for (i = 0; i < nx*ny; i++)
        data_char[i] = (char) ROUND((data[i] - offset) * scale);
  
    if ((fp = fopen(file, "w")) == 0) {
        fprintf(stderr, "write_pgm: Couldn't open file %s.\n", file);
        return(1);
    }
  
    /* write PGM header */
    fprintf(fp, "P5\n");
    /* modification to include scaling and offset */
    fprintf(fp, "# Scale: %3.5f Offset: %3.5f\n", scale, offset);
    fprintf(fp, "%d %d\n",nx,ny);
    fprintf(fp, "255\n");
  
    if (fwrite(data_char, sizeof(char), nx * ny, fp) != nx * ny) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
  
    fclose(fp);
    return(0);
}

Status
input_txt_values(
    STRING       filename,
    int        *n_values,
    double     *values[] )
{
    FILE       *fp;
    Status       status;
    double     value;
  
  
    if ((fp = fopen(filename, "r")) == 0) {
        fprintf(stderr, "input_txt_values: Couldn't open file %s.\n", filename);
        return(-1);
    }
  
    *n_values = 0;
    *values = NULL;
  
    while( input_real( fp, &value ) == OK )
        ADD_ELEMENT_TO_ARRAY( *values, *n_values, value, DEFAULT_CHUNK_SIZE );
  
    (void) close_file( fp );
  
  
    return( OK );
}


Status
input_values_any_format(char *file, int *n_values, double **values)
{
    Status status;
    FILE   *fp;

    if (filename_extension_matches(file,"txt"))
        status = input_txt_values(file, n_values, values);
    else if (filename_extension_matches(file,"gii"))
        status = input_gifti_curv(file, n_values, values);
    else
        status = input_freesurfer_curv(file, n_values, values);

    return(status);
}


Status
output_values_any_format(char *file, int n_values, void *values, int flag)
{
    Status status;
    double *buffer;
    double *d;
    int    i, *r;
    FILE   *fp;

    buffer = (double *) malloc(sizeof(double) * n_values);

    if (flag == TYPE_DOUBLE) {
        d = (double *) values;
        for (i = 0; i < n_values; i++)
            buffer[i] = (double) d[i];
    } else if (flag == TYPE_INTEGER) {
        r = (int *) values;
        for (i = 0; i < n_values; i++)
            buffer[i] = (double) r[i];
    }

    if (filename_extension_matches(file, "txt")) {
        if ((fp = fopen(file, "w")) == 0) {
            fprintf(stderr, "write_txt: Couldn't open file %s.\n", file);
            return(-1);
        }
        for (i = 0; i < n_values; i++) {
            if (flag == TYPE_DOUBLE)
                fprintf(fp, "%f\n",d[i]);
            else  fprintf(fp, "%d\n",r[i]);
        }
            
        fclose(fp);
    } else if (filename_extension_matches(file, "gii") ||
               filename_extension_matches(file, "dat"))
        status = output_gifti_curv(file, n_values, buffer);
    else
        status = output_freesurfer_curv(file, n_values, buffer);

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
    } else if (filename_extension_matches(file, "gii")) {
        status = input_gifti(file, format,
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
               object_struct **object_list, double *values)
{
    Status     status;

    if (filename_extension_matches(file, "obj")) {
        status = output_graphics_file(file, format,
                n_objects, object_list);
    } else if (filename_extension_matches(file, "off")) {
        status = output_oogl(file, format,
                n_objects, object_list);
    } else if (filename_extension_matches(file, "gii") ||
               filename_extension_matches(file, "dat")) {
        status = output_gifti(file, format,
                n_objects, object_list, values);
    } else {
        status = output_freesurfer(file, format,
                n_objects, object_list);
    }

    return(status);
}

int
read_annotation_table(char *file, int *n_array, int **out_array, int *n_labels, ATABLE **out_atable)
{
    FILE  *fp;
    char  annot_name[1000], **struct_names;
    int   i, flag, label, len, vno, tag, *array;
    int   version, structure;
    ATABLE *atable;

    fp = fopen(file, "r");
    if (!fp) {
        fprintf(stderr, "Could not open annotation file %s\n", file);
        return(-1);
    }

    *n_array = freadInt(fp);
    array = (int*) calloc (*n_array, sizeof(int));

    for (i = 0; i < *n_array; i++) {
        vno = freadInt(fp);
        label = freadInt(fp);
        if (vno >= *n_array || vno < 0) {
            fprintf(stderr, "Vertex index out of range: "
              "%d i=%8.8X, array_size=%d\n",
              vno,label,*n_array);
        } else  array[vno] = label;
    }

    *out_array = array;
    
    tag = freadInt (fp);
    if (feof(fp)) {
        fclose(fp);
        fprintf(stderr, "No colortable found\n");
        return(OK);
    }
    
    *n_labels = freadInt(fp);
    
    if (*n_labels > 0) {
        atable = (ATABLE*) malloc(*n_labels * sizeof(ATABLE));
        len = freadInt(fp);
        fread(&annot_name, sizeof(char), len, fp);
        for (i = 0; i < *n_labels; i++) {
            len = freadInt(fp);
            fread(&atable[i].name, sizeof(char), len, fp);
            atable[i].r = freadInt(fp);
            atable[i].g = freadInt(fp);
            atable[i].b = freadInt(fp);
            flag = freadInt(fp);
            atable[i].annotation = atable[i].r+atable[i].g*256+atable[i].b*65536;
        }
    } else {
        version = -(*n_labels);
        if (feof(fp)) {
            fclose(fp);
            fprintf(stderr, "Does not handle version %d\n",version);
            return(OK);
        }
        *n_labels = freadInt(fp);
        len = freadInt(fp);
        fread(&annot_name, sizeof(char), len, fp);
        *n_labels = freadInt(fp);
        atable = (ATABLE*) malloc(*n_labels * sizeof(ATABLE));
        for (i = 0; i < *n_labels; i++) {
            structure = freadInt(fp) + 1;
            len = freadInt(fp);
            fread(&atable[i].name, sizeof(char), len, fp);
            atable[i].r = freadInt(fp);
            atable[i].g = freadInt(fp);
            atable[i].b = freadInt(fp);
            flag = freadInt(fp);
            atable[i].annotation = atable[i].r+atable[i].g*256+atable[i].b*65536;
        }
    }
    *out_atable = atable;
    
    fclose(fp);

    return(OK);
}

int
write_annotation_table(char *file, int n_array, int *array, int n_labels, ATABLE *atable)
{
    FILE  *fp;
    char  **struct_names;
    int   i, flag, label, len, vno, tag;
    int   version, structure;

    fp = fopen(file, "w");
    if (!fp) {
        fprintf(stderr, "Could not open annotation file %s\n", file);
        return(-1);
    }

    fwriteInt(n_array, fp);

    for (i = 0; i < n_array; i++) {
        fwriteInt(i, fp);
        fwriteInt(array[i], fp);
    }
    
    fwriteInt (1, fp);      
    fwriteInt(-2, fp);      
    fwriteInt(n_labels, fp);    
    fwriteInt(22, fp);      
    fwrite("write_annotation_table", sizeof(char), 22, fp);
    fwriteInt(n_labels, fp);    
    
    for (i = 0; i < n_labels; i++) {
        fwriteInt(i, fp);
        fwriteInt(strlen(atable[i].name)+1, fp);
        fwrite(strcat(atable[i].name," "), sizeof(char), strlen(atable[i].name), fp);
        fwriteInt(atable[i].r, fp);
        fwriteInt(atable[i].g, fp);
        fwriteInt(atable[i].b, fp);
        fwriteInt(0, fp);
    }

    fclose(fp);

    return(OK);
}
