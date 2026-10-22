/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Clusters of a surface map: the vertices with a value > 0 that are connected
 * along mesh edges, numbered by size.  This is a thin CLI wrapper around
 * label_value_clusters() and extract_label_polygons() in CAT_Separate.
 *
 * Input  : a surface and one value per vertex.
 * Output : whatever is asked for with -label and -surface.
 */

#include <bicpl.h>
#include <string.h>

#include "ParseArgv.h"
#include "CAT_Separate.h"
#include "CAT_SurfaceIO.h"
#include "CAT_ArgvCheck.h"

char *label_file = NULL;
char *surface_prefix = NULL;
int which = 0;

static ArgvInfo argTable[] = {
    {"-label", ARGV_STRING, (char *) 1, (char *) &label_file,
          "Write the cluster number of every vertex, 0 outside the clusters."},
    {"-surface", ARGV_STRING, (char *) 1, (char *) &surface_prefix,
          "Write each cluster as its own surface, prefix_1.gii for the largest."},
    {"-which", ARGV_INT, (char *) 1, (char *) &which,
          "Restrict both outputs to cluster k (1 = largest); 0 keeps all."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};

static void usage(const char *executable)
{
    (void) fprintf(stderr,
    "\nUsage: %s [-label labels] [-surface prefix] [-which k]\n"
    "          surface_file values_file\n\n"
    "Finds the clusters of a surface map: the vertices with a value > 0 that are\n"
    "connected along mesh edges.  They are numbered by size, 1 being the largest;\n"
    "clusters of equal size keep the order of their lowest vertex index.\n\n"
    "-label writes the cluster number of every vertex, 0 outside the clusters, so\n"
    "the result stays in register with the input surface.  -surface writes the\n"
    "polygons whose vertices all belong to a cluster as a surface of their own; a\n"
    "cluster without such a polygon has no surface and is skipped.\n\n"
    "All positional arguments are inputs, for example\n"
    "     %s -label lh.clusters.txt -surface lh.cluster lh.central.gii lh.tmap.txt\n"
    "     %s -help\n\n", executable, executable, executable);
}

int main(int argc, char *argv[])
{
    char *surface_file, *values_file, out_file[1024];
    int i, k, first, last, base_len, n_objects, n_values, n_clusters;
    int n_items, n_skipped;
    int *labels, *sizes;
    double *values;
    File_formats format;
    object_struct **object_list, *cluster_object;
    polygons_struct *polygons;

    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* a leftover that looks like an option is a typo or a removed flag;
       ParseArgv() keeps it, so it would be read as a file name below */
    if (cat_check_unknown_options(argc, argv, argv[0]) || (argc != 3)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    surface_file = argv[1];
    values_file = argv[2];

    if ((label_file == NULL) && (surface_prefix == NULL)) {
        (void) fprintf(stderr,
        "Error: No output requested. Write the cluster numbers with -label and the\n"
        "cluster surfaces with -surface.\n");
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (which < 0) {
        (void) fprintf(stderr, "Error: -which must be 0 (all clusters) or a cluster number.\n");
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(surface_file, &format, &n_objects,
                                  &object_list) != OK ||
        n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
        (void) fprintf(stderr, "Error: %s must contain one polygons object.\n", surface_file);
        exit(EXIT_FAILURE);
    }
    polygons = get_polygons_ptr(object_list[0]);

    if (input_values_any_format(values_file, &n_values, &values) != OK) {
        (void) fprintf(stderr, "Error: Cannot read values from %s.\n", values_file);
        exit(EXIT_FAILURE);
    }
    if (n_values != polygons->n_points) {
        (void) fprintf(stderr, "Error: %s has %d values for %d vertices.\n",
                       values_file, n_values, polygons->n_points);
        exit(EXIT_FAILURE);
    }

    labels = (int *) malloc(sizeof(int) * polygons->n_points);
    n_clusters = label_value_clusters(polygons, values, labels);

    if (which > n_clusters) {
        (void) fprintf(stderr, "Error: There is no cluster %d; the map has %d.\n",
                       which, n_clusters);
        exit(EXIT_FAILURE);
    }

    sizes = (int *) calloc(n_clusters + 1, sizeof(int));
    for (i = 0; i < polygons->n_points; i++)
        sizes[labels[i]]++;

    if (n_clusters > 0)
        fprintf(stdout, "%d cluster%s, the largest with %d vertices\n", n_clusters,
                (n_clusters == 1) ? "" : "s", sizes[1]);
    else
        fprintf(stdout, "No vertex has a value > 0\n");

    first = (which > 0) ? which : 1;
    last = (which > 0) ? which : n_clusters;

    if (label_file != NULL) {
        if (which > 0) {
            for (i = 0; i < polygons->n_points; i++) {
                if (labels[i] != which)
                    labels[i] = 0;
            }
        }
        if (output_values_any_format(label_file, polygons->n_points, labels,
                                     TYPE_INTEGER) != OK)
            exit(EXIT_FAILURE);
    }

    if (surface_prefix != NULL) {
        /* a prefix given with the extension must not end up as name.gii_1.gii */
        base_len = (int) strlen(surface_prefix);
        if (base_len > 4 && strcmp(surface_prefix + base_len - 4, ".gii") == 0)
            base_len -= 4;

        n_skipped = 0;
        for (k = first; k <= last; k++) {
            (void) snprintf(out_file, sizeof(out_file), "%.*s_%d.gii",
                            base_len, surface_prefix, k);

            cluster_object = create_object(POLYGONS);
            n_items = extract_label_polygons(polygons, labels, k,
                                             get_polygons_ptr(cluster_object));
            if (n_items == 0) {
                n_skipped++;
            } else {
                if (output_graphics_any_format(out_file, format, 1, &cluster_object,
                                               NULL) != OK)
                    exit(EXIT_FAILURE);
                fprintf(stdout, "%s: cluster %d, %d vertices, %d polygons\n",
                        out_file, k, sizes[k], n_items);
            }
            delete_object(cluster_object);
        }

        if (which > 0 && n_skipped > 0) {
            (void) fprintf(stderr, "Error: Cluster %d has no polygon whose vertices all "
                           "belong to it, so it has no surface.\n", which);
            exit(EXIT_FAILURE);
        }
        if (n_skipped > 0)
            fprintf(stdout, "Skipped %d cluster%s without a polygon of their own\n",
                    n_skipped, (n_skipped == 1) ? "" : "s");
    }

    free(labels);
    free(sizes);
    free(values);
    delete_object_list(n_objects, object_list);

    return (EXIT_SUCCESS);
}
