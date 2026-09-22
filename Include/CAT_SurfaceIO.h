/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#ifndef CAT_SURFACEIO_H
#define CAT_SURFACEIO_H

#include <bicpl.h>
#include <stdio.h>
#include <string.h>
#include <gifti_io.h>
#include <nifti1.h>
#include <nifti1_io.h>
#include <quadric.h>

#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)
#define NEW_VERSION_MAGIC_NUMBER    16777215

#define TYPE_DOUBLE 1
#define TYPE_INTEGER 2
#define TYPE_CHAR 3

/* allow override */
#ifndef BYTE_ORDER

/////////////Linux////////////////////////////
#ifdef linux
#include <endian.h>

#ifndef BYTE_ORDER
#define BYTE_ORDER __BYTE_ORDER
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN __LITTLE_ENDIAN
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN __BIG_ENDIAN
#endif

#endif

/////////////Windows Cygwin////////////////////////////
#ifdef WIN32

#define BIG_ENDIAN  4321
#define LITTLE_ENDIAN 1234
#define BYTE_ORDER  LITTLE_ENDIAN

#endif

////////////MacOS X and BSD ////////////////////////////
#if defined(__APPLE__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <machine/endian.h>
#endif

////////////Solaris 2.5.1//////////////////////
#ifdef sun

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN    4321
#endif

#include <sys/isa_defs.h>
/* only defines one of _LITTLE_ENDIAN or _BIG_ENDIAN */
#ifdef _LITTLE_ENDIAN
#define BYTE_ORDER LITTLE_ENDIAN
#endif

#ifdef _BIG_ENDIAN
#define BYTE_ORDER BIG_ENDIAN
#endif

#endif

/////////////IRIX  ////////////////////////////
#if defined(__sgi) || defined(Mips)
#include <sys/endian.h>
#endif

///////////////////////////////////////////////////
#endif /* BYTE_ORDER */
///////////////////////////////////////////////////

typedef struct
{
  int    r, g, b ;
  int    annotation ;
  char   name[1000] ;
}
ATABLE ;

/**
 * \brief Convert a BICPL polygons_struct into flat arrays (V,F) of vertices and
 *        triangles for quadric.c.
 *
 * \param poly             (in)  input mesh (may contain n-gons)
 * \param[out] V           (out) newly malloc'd array of vec3d of size nv
 * \param[out] F           (out) newly malloc'd array of vec3i of size nf
 * \param[out] nv          (out) number of vertices copied
 * \param[out] nf          (out) number of triangles written
 * \param[out] fan_used    (out) non-zero if any face required fan-triangulation
 * \return 0 on success, -1 on allocation/argument errors or if no triangles produced
 */
int polygons_to_tri_arrays(const polygons_struct *poly, vec3d **V, vec3i **F,
                           int *nv, int *nf, int *fan_used);
/**
 * \brief Write (V,F) triangle arrays from quadric.c back into a BICPL polygons_struct.
 *
 * \param poly (in/out) target mesh to overwrite (triangles only on exit)
 * \param V    (in)     vertex positions (size nv)
 * \param F    (in)     triangle index triplets (size nf)
 * \param nv   (in)     number of vertices
 * \param nf   (in)     number of triangles
 * \return 0 on success, -1 on allocation errors
 */
int tri_arrays_to_polygons(polygons_struct *poly, const vec3d *V,
                           const vec3i *F, int nv, int nf);
/**
 * \brief Utility to clamp/derive a valid QEM target (triangle) count.
 *
 * \param nf_total (int) number of input triangles
 * \param target   (int) requested target; if <=0 uses half; clamped to [1, nf_total]
 */
int qem_target(int nf_total, int target);
/**
 * \brief Read 1D scalar values from file in auto-detected format.
 *
 * \param file     (in)  input file path
 * \param n_values (out) number of values read
 * \param values   (out) allocated array of scalar values
 * \return OK on success; ERROR if format is unrecognized or read fails
 */
Status input_values_any_format(char *file, int *n_values, double **values);
/**
 * \brief Write 1D scalar values to file in auto-detected format.
 *
 * \param file     (in) output file path
 * \param n_values (in) number of scalar values
 * \param values   (in) array of values (type determined by flag)
 * \param flag     (in) data type: TYPE_DOUBLE, TYPE_INTEGER, or TYPE_CHAR
 * \return OK on success; ERROR on write failure
 */
Status output_values_any_format(const char *file, int n_values, void *values,
                                int flag);
/**
 * \brief Read polygon mesh in auto-detected geometry format.
 *
 * \param file         (in)  mesh file path
 * \param format       (out) detected file format
 * \param n_objects    (out) number of mesh objects (typically 1)
 * \param object_list  (out) allocated polygon mesh objects
 * \return OK on success; ERROR if format unrecognized or read fails
 */
Status input_graphics_any_format(char *file, File_formats *format,
                                 int *n_objects, object_struct ***object_list); 
/**
 * \brief Write polygon mesh to auto-detected output format.
 *
 * \param file         (in) output mesh file path
 * \param format       (in) format specification
 * \param n_objects    (in) number of mesh objects (typically 1)
 * \param object_list  (in) array of polygon mesh objects
 * \param values       (in) optional per-vertex scalar data; NULL to skip
 * \return OK on success; ERROR if format unrecognized or write fails
 */
Status output_graphics_any_format(char *file, File_formats format,
                                  int n_objects, object_struct **object_list,
                                  double *values);
/**
 * \brief Read both mesh and texture data from a single GIFTI file.
 *
 * \param file         (in)  path to GIFTI file
 * \param format       (out) file format identifier
 * \param n_objects    (out) number of objects (typically 1)
 * \param object_list  (out) allocated mesh objects
 * \param n_values     (out) number of per-vertex scalar values
 * \param values       (out) per-vertex data array
 * \return OK on success; ERROR if file format is invalid
 */
int input_gifti_mesh_and_texture(char *file, File_formats *format,
                                 int *n_objects, object_struct ***object_list,
                                 int *n_values, double **values);
/**
 * \brief Read OFF format mesh file (OOGL/Geomview format).
 *
 * \param file         (in)  path to OFF file
 * \param format       (out) set to ASCII_FORMAT
 * \param n_objects    (out) set to 1 (single polygon object)
 * \param object_list  (out) allocated array of objects (single POLYGONS object)
 * \return OK on success; ERROR if file format is invalid or non-triangular faces found
 */
int input_oogl(char *file, File_formats *format, int *n_objects,
               object_struct ***object_list);
/**
 * \brief Write polygons to OFF format file (OOGL/Geomview).
 *
 * \param file        (in) output file path
 * \param format      (in) format specification (typically ASCII_FORMAT)
 * \param n_objects   (in) number of objects to write (typically 1)
 * \param object_list (in) array of objects; first must be a mesh
 * \return OK on success; ERROR on file write failure
 */
int output_oogl(char *file, File_formats format, int n_objects,
                object_struct *object_list[]);
/**
 * \brief Write the first polygons object of a list as a FreeSurfer triangle surface.
 *
 * \param file        (in) output file name
 * \param format      (in) unused; FreeSurfer surfaces are always binary
 * \param n_objects   (in) number of objects in object_list (only the first is written)
 * \param object_list (in) objects; the first must be POLYGONS
 * \return OK on success
 */
int output_freesurfer(char *file, File_formats format, int n_objects,
                      object_struct *object_list[]);
/**
 * \brief Write per-vertex values as a FreeSurfer curvature file (new format).
 *
 * \param fname     (in) output file name
 * \param nvertices (in) number of values
 * \param data      (in) values, stored as float
 * \return OK on success
 */
int output_freesurfer_curv(char *fname, int nvertices, double *data);
/**
 * \brief Read a FreeSurfer triangle surface into a new polygons object.
 *
 * \param file        (in)     input file name
 * \param format      (out)    set to ASCII_FORMAT
 * \param n_objects   (out)    number of objects read (1)
 * \param object_list (in/out) list the new POLYGONS object is appended to
 * \return OK on success, ERROR for quad files or an unknown magic number
 */
int input_freesurfer(char *file, File_formats *format, int *n_objects,
                     object_struct ***object_list);
/**
 * \brief Read per-vertex scalar data from FreeSurfer curv format.
 *
 * \param file         (in)  path to FreeSurfer curv file
 * \param vnum         (out) number of vertices
 * \param input_values (out) allocated array of per-vertex data
 * \return OK on success; ERROR on file read failure
 */
int input_freesurfer_curv(char *file, int *vnum, double **input_values);
/**
 * \brief Write polygon mesh and optional scalar data to GIFTI format.
 *
 * \param fname       (in) output GIFTI filename (.gii); if .dat path given, triggers external binary storage
 * \param format      (in) format specification
 * \param n_objects   (in) number of mesh objects (typically 1)
 * \param object_list (in) array of objects; first should be a POLYGONS mesh
 * \param values      (in) optional scalar data array (per-vertex); NULL to skip
 * \return 0 on success; -1 on allocation or write error
 */
int output_gifti(char *fname, File_formats format, int n_objects,
                 object_struct *object_list[], double *values);
/**
 * \brief Write per-vertex scalar data to GIFTI format (shape data).
 *
 * \param fname     (in) output GIFTI filename (.gii); or .dat to trigger external binary pair
 * \param nvertices (in) number of vertices (length of data array)
 * \param data      (in) per-vertex scalar values (size nvertices)
 * \return 0 on success; -1 on allocation or write error
 */
int output_gifti_curv(char *fname, int nvertices, double *data);
/**
 * \brief Read polygon mesh and optional texture data from GIFTI format.
 *
 * \param file         (in)  path to GIFTI file (.gii)
 * \param format       (out) set to ASCII_FORMAT
 * \param n_objects    (out) number of objects loaded (typically 1)
 * \param object_list  (out) allocated array of objects containing the mesh
 * \param n_values     (out) number of per-vertex values (0 if no texture)
 * \param values       (out) allocated shape data array (NULL if not present)
 * \return OK on success; ERROR if file is invalid or corrupted
 */
int input_gifti(char *file, File_formats *format, int *n_objects,
                object_struct ***object_list, int *n_values, double **values);

/**
 * \brief Description of one DataArray stored in a GIFTI file.
 *
 * A .gii file holds the mesh as two DataArrays (NIFTI_INTENT_POINTSET and
 * NIFTI_INTENT_TRIANGLE) but may carry any number of further arrays next to
 * them -- thickness, curvature, labels, time series.  input_gifti() reads only
 * the mesh and the first shape array, so this structure exists to report what
 * else the file contains.  All strings are NUL-terminated and never NULL.
 */
typedef struct {
    int    intent;             /**< NIFTI_INTENT_* code                      */
    char   intent_name[64];    /**< intent as text, e.g. NIFTI_INTENT_SHAPE  */
    char   datatype_name[32];  /**< value type, e.g. NIFTI_TYPE_FLOAT32      */
    char   encoding_name[32];  /**< ASCII, Base64Binary, ExternalFileBinary  */
    char   ext_fname[256];     /**< external data file, empty when embedded  */
    char   name[128];          /**< Name metadata entry, empty when absent   */
    int    num_dim;            /**< number of dimensions                     */
    int    dims[6];            /**< dimension lengths, first num_dim set     */
    long long n_values;        /**< total number of values                   */
    long long n_nonfinite;     /**< NaN and infinite values, excluded below  */
    int    has_range;          /**< 1 when min, mean and max are valid       */
    double min;                /**< smallest finite value                    */
    double mean;               /**< mean of the finite values                */
    double max;                /**< largest finite value                     */
} gifti_darray_info;

/**
 * \brief List the DataArrays a GIFTI file contains.
 *
 * Reads the file and describes every DataArray in it, including those
 * input_gifti() ignores.  Value statistics are filled in for one- and
 * two-dimensional arrays of a numeric type and skipped otherwise, in which
 * case has_range is 0.  NaN and infinite entries are counted in n_nonfinite
 * and left out of the statistics, so one masked-out vertex cannot turn the
 * whole summary into nan.
 *
 * \param file     (in)  path to the GIFTI file
 * \param n_arrays (out) number of entries written to arrays
 * \param arrays   (out) allocated array of descriptions, free with free()
 * \return OK on success, ERROR if the file cannot be read or is not GIFTI
 */
Status input_gifti_darrays(char *file, int *n_arrays, gifti_darray_info **arrays);
/**
 * \brief Read per-vertex scalar data from GIFTI format.
 *
 * \param file     (in)  path to GIFTI file
 * \param vnum     (out) number of vertices (length of data array)
 * \param input_values (out) allocated array of per-vertex scalar values
 * \return OK on success; ERROR if file format is invalid or corrupted
 */
int input_gifti_curv(char *file, int *vnum, double **input_values);
/**
 * \brief Read polygon mesh from OpenDX (Data Explorer) format.
 *
 * \param file         (in)  path to OpenDX file
 * \param format       (out) set to ASCII_FORMAT
 * \param n_objects    (out) number of objects (typically 1)
 * \param object_list  (out) allocated polygon mesh object
 * \return OK on success; ERROR if file format is invalid
 */
int input_dx(char *file, File_formats *format, int *n_objects,
             object_struct ***object_list);
/**
 * \brief Read polygon mesh from DFS (Freesurfer-specific) format.
 *
 * \param file         (in)  path to DFS file
 * \param format       (out) format identifier
 * \param n_objects    (out) number of objects (typically 1)
 * \param object_list  (out) allocated polygon mesh object
 * \return OK on success; ERROR if format is not supported or invalid\n
 */
int input_dfs(char *file, File_formats *format, int *n_objects,
              object_struct ***object_list);
/**
 * \brief Read 2D image data from PGM (Portable Graymap) format.
 *
 * \param file (in)  path to PGM file
 * \param nx   (out) image width (number of columns)
 * \param ny   (out) image height (number of rows)
 * \return allocated double array of size nx*ny (row-major order); NULL on error
 */
double * read_pgm(char *file, int *nx, int *ny);
/**
 * \brief Write 2D image data to PGM (Portable Graymap) format.
 *
 * \param file (in) output PGM file path
 * \param data (in) 2D image array (size nx*ny, row-major order)
 * \param nx   (in) image width (number of columns)
 * \param ny   (in) image height (number of rows)
 * \return 0 on success; -1 on file write error
 */
int write_pgm(char *file, double *data, int nx, int ny);
/**
 * \brief Read FreeSurfer annotation table (ROI labels) from file.
 *
 * \param file       (in)  path to .annot file
 * \param n_array    (out) number of vertices in annotation
 * \param out_array  (out) allocated array of per-vertex labels (indices into label table)
 * \param n_labels   (out) number of unique labels in table
 * \param out_atable (out) allocated label table with names and colors
 * \return number of entries read on success; ERROR on read failure
 */
int read_annotation_table(char *file, int *n_array, int **out_array,
                          int *n_labels, ATABLE **out_atable);
/**
 * \brief Write FreeSurfer annotation table (ROI labels) to file.
 *
 * \param file     (in) output .annot file path
 * \param n_array  (in) number of vertices to write
 * \param array    (in) per-vertex label indices (size n_array)
 * \param n_labels (in) number of labels in table
 * \param atable   (in) label table with names and RGB colors
 * \return OK on success; ERROR on write failure
 */
int write_annotation_table(char *file, int n_array, int *array, int n_labels,
                           ATABLE *atable);

#endif
