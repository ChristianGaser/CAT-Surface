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

#if !defined(_WIN32)
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


/* Main program */
int
main(int argc, char *argv[])
{
        char *infile, outfile[1024];
        int i, dims[3];
        float *input;
        double separations[3], s[3];
        nifti_image *nii_ptr;

        /* Get arguments */
#if !defined(_WIN32)
        if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
#else
        if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 3)) {
#endif
                (void) fprintf(stderr, 
                "\nUsage: %s [options] in.nii [out.nii]\n", argv[0]);
                (void) fprintf(stderr, 
                "         %s -help\n\n", argv[0]);
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
        
        separations[0] = nii_ptr->dx;
        separations[1] = nii_ptr->dy;
        separations[2] = nii_ptr->dz;
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;
        
        smooth3(input, dims, separations, s, use_mask, DT_FLOAT32);
        
        /* if not defined use original name as basename for output */
        if (argc == 3)
                        (void) sprintf(outfile, "%s", argv[2]); 
        else {
                #if !defined(_WIN32)
                        (void) sprintf(outfile, "%s/s%g%s", dirname(infile), fwhm, basename(infile)); 
                #endif
        }

        /* write data using same data type and rescale */
        if (!write_nifti_float(outfile, input, nii_ptr->datatype, 0.0, dims, separations, nii_ptr)) 
                exit(EXIT_FAILURE);
        
        return(EXIT_SUCCESS);

}
