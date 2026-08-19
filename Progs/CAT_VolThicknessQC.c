/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_ThicknessQC.h"

double thresh = 4.5;
double plate_radius = 2.5;
double min_volume = 20.0;
int conn = 26;
char *label_file = NULL;
char *classmap_file = NULL;
int verbose = 0;

static ArgvInfo argTable[] = {
    {"-thresh", ARGV_FLOAT, (char *)1, (char *)&thresh,
     "Thickness above this is implausible, in mm (default: 4.5). Human cortex\n\
         reaches about 4.5 mm at the thickest."},
    {"-plate-radius", ARGV_FLOAT, (char *)1, (char *)&plate_radius,
     "Largest inscribed radius still reported as a plate, in mm (default: 2.5).\n\
         Two cortical banks back to back form a band of at most about 5 mm, so\n\
         a glued sulcus cannot exceed half that; anything thicker is a mass."},
    {"-min-volume", ARGV_FLOAT, (char *)1, (char *)&min_volume,
     "Discard components below this volume, in mm^3 (default: 20)."},
    {"-conn", ARGV_INT, (char *)1, (char *)&conn,
     "Voxel connectivity: 6, 18 or 26 (default: 26)."},
    {"-label", ARGV_STRING, (char *)1, (char *)&label_file,
     "PVE label map restricting the search to the cortical band. Recommended:\n\
         without it the thickness map is used as it stands, including whatever\n\
         it holds outside the ribbon."},
    {"-classmap", ARGV_STRING, (char *)1, (char *)&classmap_file,
     "Write a map tagging each flagged voxel with 1 = plate, 2 = solid, for\n\
         overlaying on the thickness map."},
    {"-v", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <gmt.nii>\n\
\n\
    Triage implausibly thick cortex by shape.\n\
\n\
    A thickness map normally carries a population of voxels above any\n\
    defensible value -- cortex does not exceed roughly 4.5 mm. Two different\n\
    faults produce them, and the repair for one is harmful applied to the other:\n\
\n\
      plate   two sulcal banks with no CSF between them, measured as one band.\n\
              There is a sulcus to recover; CAT_VolSulcusRepair is the tool.\n\
      solid   cortex merged with subcortical grey matter, or a genuinely thick\n\
              region such as a temporal or orbitofrontal pole. There is no\n\
              sulcus, and carving one invents anatomy.\n\
\n\
    Thickness alone cannot separate them, which makes it easy to attribute the\n\
    whole population to gluing and then build a separator for a fault that is\n\
    not there. Shape separates them with one number. A glued sulcus is a band\n\
    at most about 5 mm across however far it runs along the sulcus, so the\n\
    largest sphere fitting inside it has a radius of about 2.5 mm. A solid mass\n\
    has no such bound. The maximum inscribed radius of a connected component --\n\
    the maximum of the Euclidean distance transform taken inside it -- is\n\
    therefore the discriminator, and it depends neither on where the component\n\
    sits nor on how large it is.\n\
\n\
    For scale: a normal 2.5 mm ribbon gives about 1.25 mm, a glued pair of\n\
    banks stays under 2.5 mm, and a subcortical mass or a pole runs past\n\
    3.5 mm.\n\
\n\
Options:\n\
    -thresh       <float>  Implausible thickness, mm (default: 4.5).\n\
    -plate-radius <float>  Largest plate inscribed radius, mm (default: 2.5).\n\
    -min-volume   <float>  Smallest reported component, mm^3 (default: 20).\n\
    -conn         <int>    Connectivity 6, 18 or 26 (default: 26).\n\
    -label        <file>   PVE label map confining the search to the ribbon.\n\
    -classmap     <file>   Write per-voxel classes (1 = plate, 2 = solid).\n\
    -v                     Be verbose.\n\
\n\
Example:\n\
    %s -label p0.nii -classmap qc.nii gmt.nii\n\
    CAT_VolThicknessPbt p0.nii gmt.nii ppm.nii && %s -label p0.nii gmt.nii\n\n";

    fprintf(stderr, usage_str, executable, executable, executable);
}

/* Order the report worst-first, so the component most in need of attention is
   the one at the top rather than whichever happened to be found first. */
static int cmp_volume(const void *a, const void *b)
{
    const CAT_ThicknessComponent *ca = (const CAT_ThicknessComponent *)a;
    const CAT_ThicknessComponent *cb = (const CAT_ThicknessComponent *)b;
    if (ca->volume_mm3 < cb->volume_mm3)
        return 1;
    if (ca->volume_mm3 > cb->volume_mm3)
        return -1;
    return 0;
}

int main(int argc, char *argv[])
{
    char *gmt_file;
    int dims[3], nvox, rc, i;
    float *gmt = NULL, *label = NULL, *classmap = NULL;
    double voxelsize[3];
    nifti_image *gmt_ptr, *label_ptr = NULL;
    CAT_ThicknessQCOpts opts;
    CAT_ThicknessComponent *comps = NULL;
    int n_comps = 0, n_plate = 0, n_solid = 0;
    double vol_plate = 0.0, vol_solid = 0.0;

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    gmt_file = argv[1];

    gmt_ptr = read_nifti_float(gmt_file, &gmt, 0);
    if (!gmt_ptr)
    {
        fprintf(stderr, "Error reading %s.\n", gmt_file);
        return EXIT_FAILURE;
    }

    dims[0] = gmt_ptr->nx;
    dims[1] = gmt_ptr->ny;
    dims[2] = gmt_ptr->nz;
    voxelsize[0] = gmt_ptr->dx;
    voxelsize[1] = gmt_ptr->dy;
    voxelsize[2] = gmt_ptr->dz;
    nvox = dims[0] * dims[1] * dims[2];

    if (label_file)
    {
        label_ptr = read_nifti_float(label_file, &label, 0);
        if (!label_ptr)
        {
            fprintf(stderr, "Error reading %s.\n", label_file);
            return EXIT_FAILURE;
        }
        if (label_ptr->nx != dims[0] || label_ptr->ny != dims[1] ||
            label_ptr->nz != dims[2])
        {
            fprintf(stderr, "Label map must have the same dimensions as the thickness map.\n");
            return EXIT_FAILURE;
        }
    }

    if (classmap_file)
    {
        classmap = (float *)malloc(sizeof(float) * (size_t)nvox);
        if (!classmap)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }

    CAT_ThicknessQCOptionsInit(&opts);
    opts.thresh = thresh;
    opts.plate_radius = plate_radius;
    opts.min_volume = min_volume;
    opts.conn = conn;
    opts.verbose = verbose;

    rc = CAT_VolThicknessQC(gmt, label, dims, voxelsize, &opts, &comps, &n_comps,
                            classmap);
    if (rc != 0)
    {
        fprintf(stderr, "Thickness QC failed (%d).\n", rc);
        exit(EXIT_FAILURE);
    }

    qsort(comps, (size_t)n_comps, sizeof(CAT_ThicknessComponent), cmp_volume);

    printf("\n%s: %d component(s) above %.2f mm\n\n", gmt_file, n_comps, thresh);
    if (n_comps > 0)
    {
        printf("  %5s %10s %9s %9s %8s %8s   %-24s\n",
               "rank", "voxels", "vol/mm^3", "radius/mm", "mean/mm", "max/mm",
               "centroid (voxels)");
        for (i = 0; i < n_comps; i++)
        {
            printf("  %5d %10ld %9.1f %9.2f %8.2f %8.2f   %5.0f %5.0f %5.0f   %s\n",
                   i + 1, comps[i].n_voxels, comps[i].volume_mm3,
                   comps[i].max_radius, comps[i].gmt_mean, comps[i].gmt_max,
                   comps[i].centroid[0], comps[i].centroid[1], comps[i].centroid[2],
                   comps[i].shape == CAT_QC_PLATE ? "PLATE (glued sulcus?)"
                                                  : "SOLID (not a sulcus)");
        }
    }

    for (i = 0; i < n_comps; i++)
    {
        if (comps[i].shape == CAT_QC_PLATE)
        {
            n_plate++;
            vol_plate += comps[i].volume_mm3;
        }
        else
        {
            n_solid++;
            vol_solid += comps[i].volume_mm3;
        }
    }

    printf("\n  plate-like : %4d component(s), %9.1f mm^3  -- recoverable sulci\n",
           n_plate, vol_plate);
    printf("  solid      : %4d component(s), %9.1f mm^3  -- not sulci; masking, not repair\n",
           n_solid, vol_solid);
    if (vol_plate + vol_solid > 0.0)
        printf("  plate share: %.1f%% of the flagged volume\n\n",
               100.0 * vol_plate / (vol_plate + vol_solid));
    else
        printf("\n");

    if (classmap_file)
    {
        if (!write_nifti_float(classmap_file, classmap, DT_FLOAT32, 1.0, dims,
                               voxelsize, gmt_ptr))
            exit(EXIT_FAILURE);
    }

    free(gmt);
    free(label);
    free(classmap);
    free(comps);

    return EXIT_SUCCESS;
}
