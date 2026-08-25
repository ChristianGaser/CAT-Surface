/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <float.h>
#include <stdlib.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"

double fwhm = 8.0;
int use_mask = 0;
int verbose = 0;

static
ArgvInfo argTable[] = {
    {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm, 
             "FWHM in mm."},
    {"-mask", ARGV_CONSTANT, (char *) 1, (char *) &use_mask,
             "Use masked smoothing (default no masking)."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
             "Be verbose."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};


static void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <in.nii> [<out.nii>]\n\
\n\
    Smooth a volume.  By default this is an isotropic Gaussian of -fwhm mm; the\n\
    output is written next to the input with an s<fwhm> prefix when no output\n\
    name is given.\n\
\n\
    With -oriented the filter instead diffuses *along* thin sheets and not\n\
    across them.  Each 3x3x3 neighbour is weighted by how well its direction\n\
    lies in the plane of the local sheet, which is found by a Hessian plate\n\
    filter (see CAT_VolSheetness), so a sulcal CSF sheet or a gyral white-matter\n\
    blade is smoothed along its length and never bridged.  Where no sheet is\n\
    detected the weights reduce to plain inverse distance and the operator is\n\
    ordinary local averaging.  This is an iterated local kernel rather than a\n\
    Gaussian: -fwhm does not apply and the amount of smoothing is set by -iter.\n\
\n\
    Pass the intensity image with -guide when smoothing a label or probability\n\
    map, so the orientation comes from the image that still shows the\n\
    structure.\n\
\n\
    Every option is listed with its default under 'Command-specific options'\n\
    above.\n\
\n\
Examples:\n\
    %s -fwhm 8 in.nii out.nii\n\
    %s -oriented -iter 3 -guide t1_corr.nii label.nii label_smoothed.nii\n\n";

    fprintf(stderr, usage_str, executable, executable, executable);
}

/* Main program */
int
main(int argc, char *argv[])
{
    char *infile, outfile[1024];
    int i, dims[3];
    float *input;
    double voxelsize[3], s[3];
    nifti_image *nii_ptr;

    /* Get arguments */
#if !defined(_WIN32) && !defined(_WIN64)
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
#else
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 3)) {
#endif
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    infile = argv[1];
    
    if (verbose)
        fprintf(stdout,"Filtering %s with FWHM of %gmm.\n", infile, fwhm);
    
    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(infile, &input, 0);
    if (nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }

    /* only allow isotropic filtering */
    for (i=0; i<3; i++) s[i] = fwhm;
    
    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    
    smooth3(input, dims, voxelsize, s, use_mask, DT_FLOAT32);
    
    /* if not defined use original name as basename for output */
    if (argc == 3)
            (void) sprintf(outfile, "%s", argv[2]); 
    else {
        #if !defined(_WIN32) && !defined(_WIN64)
            (void) sprintf(outfile, "%s/s%g%s", dirname(infile), fwhm, basename(infile)); 
        #endif
    }

    /* write data using same data type and rescale */
    if (!write_nifti_float(outfile, input, nii_ptr->datatype, 0.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    return(EXIT_SUCCESS);

}
