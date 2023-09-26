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

int verbose = 0;

extern int ParseArgv(int *argcPtr, char **argv, ArgvInfo *argTable, int flags);

static ArgvInfo argTable[] = {
        {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
                  "Be verbose."},
         {NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[])
{
        char *infile, outfile[1024];
        int i, j, dims[3];
        float F, *input;
        double separations[3];
        nifti_image *nii_ptr;
        
        if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
                 (void) fprintf(stderr, "\nUsage: %s [options] in.nii [out.nii]\nProjection-based thickness estimation\n\n", argv[0]);
                 (void) fprintf(stderr, "         %s -help\n\n", argv[0]);
         exit(EXIT_FAILURE);
        }
        
        infile  = argv[1];

        /* if not defined use original name as basename for output */
        if(argc == 3)
                (void) sprintf(outfile, "%s", argv[2]); 
        else {
                #if !defined(_WIN32)
                        (void) sprintf(outfile, "%s/n%s", dirname(infile), basename(infile)); 
                #else
                        fprintf(stderr,"\n\Usage: %s input.nii output.nii\n\n\
                Projection-based thickness estimation.\n\n", argv[0]);
                        return( 1 );
                #endif
        }
        
        /* read first image to get image parameters */
        nii_ptr = read_nifti_float(infile, &input, 0);
        if(nii_ptr == NULL) {
                fprintf(stderr,"Error reading %s.\n", infile);
                return(EXIT_FAILURE);
        }

        separations[0] = nii_ptr->dx;
        separations[1] = nii_ptr->dy;
        separations[2] = nii_ptr->dz;
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;

        /* 

    % estimate thickness with PBT approach
    Ygmt1 = cat_vol_pbtp(round(Ymf),Ywmd,Ycsfd);  
    Ygmt2 = cat_vol_pbtp(round(4-Ymf),Ycsfd,Ywmd);

    % avoid meninges !
    Ygmt1 = min(Ygmt1,Ycsfd+Ywmd);
    Ygmt2 = min(Ygmt2,Ycsfd+Ywmd); 

    Ygmt  = min(cat(4,Ygmt1,Ygmt2),[],4); %clear Ygmt1 Ygmt2; 
            
        */        
        
        if (!write_nifti_float( outfile, input, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
                exit(EXIT_FAILURE);

        free(input);
        
        return(EXIT_SUCCESS);

}