/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_VolPbt.h"

int verbose = 0;

static ArgvInfo argTable[] = {
    {"-v", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

int main(int argc, char *argv[])
{
    char *infile;
    char outfile[1024];
    int dims[3];
    float *input;
    double voxelsize[3];
    nifti_image *nii_ptr;

#if !defined(_WIN32) && !defined(_WIN64)
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2))
#else
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3))
#endif
    {
        (void)fprintf(stderr,
                      "\nUsage: %s [options] in.nii [out.nii]\n", argv[0]);
        (void)fprintf(stderr,
                      "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    infile = argv[1];

    if (verbose)
        fprintf(stdout, "Blood-vessel correction for %s.\n", infile);

    nii_ptr = read_nifti_float(infile, &input, 0);
    if (nii_ptr == NULL)
    {
        fprintf(stderr, "Error reading %s.\n", infile);
        return EXIT_FAILURE;
    }

    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;

    blood_vessel_correction_pve_float(input, dims, voxelsize);

    if (argc == 3)
        (void)sprintf(outfile, "%s", argv[2]);
    else
    {
#if !defined(_WIN32) && !defined(_WIN64)
        (void)sprintf(outfile, "%s/bvc_%s", dirname(infile), basename(infile));
#endif
    }

    if (!write_nifti_float(outfile, input, nii_ptr->datatype, 0.0, dims, voxelsize, nii_ptr))
        exit(EXIT_FAILURE);

    return EXIT_SUCCESS;
}
