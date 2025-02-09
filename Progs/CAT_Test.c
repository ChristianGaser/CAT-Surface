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

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"


#define IDX(x, y, z, nx, ny) ((z) * (nx) * (ny) + (y) * (nx) + (x))


void correct_topology(float *volume, float thresh, int dims[3], int conn_arr[2], int n_loops) {
    float *vol_euler;
    unsigned short *vol_bin;
    int loop, i, x, y, z, nx, ny, nz, iter, sign, chi_error;
    int V, E, F, conn, chi_local, n_errors, cube[8];
    int edges18[18][2] = {
        {0,1}, {1,3}, {3,2}, {2,0}, // Bottom face edges
        {4,5}, {5,7}, {7,6}, {6,4}, // Top face edges
        {0,4}, {1,5}, {2,6}, {3,7}, // Vertical edges
        {0,5}, {1,4}, {2,7}, {3,6}, // Diagonal edges across faces
        {0,7}, {3,4}  // Extended diagonal connections
    };
    int edges26[26][2] = {
        {0,1}, {1,3}, {3,2}, {2,0}, // Bottom face edges
        {4,5}, {5,7}, {7,6}, {6,4}, // Top face edges
        {0,4}, {1,5}, {2,6}, {3,7}, // Vertical edges
        {0,5}, {1,4}, {2,7}, {3,6}, // Diagonal edges across faces
        {0,7}, {3,4}, {2,5}, {1,6}, // Extended diagonal connections
        {0,6}, {1,7}, {2,4}, {3,5}  // Corner-to-corner vertex edges
    };
    // Count faces (6 faces in a cube)
    int faces[6][4] = {
        {0,1,3,2}, {4,5,7,6}, // XY faces
        {0,2,6,4}, {1,3,7,5}, // YZ faces
        {0,1,5,4}, {2,3,7,6}  // XZ faces
    };
    
    vol_euler = (float *)malloc(sizeof(float)*dims[0]*dims[1]*dims[2]);
    vol_bin = (unsigned short *)malloc(sizeof(unsigned short)*dims[0]*dims[1]*dims[2]);

    nx = dims[0], ny = dims[1], nz = dims[2];

    // Initialize Euler map
    for (i = 0; i < nx * ny * nz; i++)
        vol_euler[i] = (volume[i] > thresh) ? 1.0 : 0.0;

    for (iter = 0; iter < 50; iter++) {
        n_errors = 0;
        for (loop = 0; loop < n_loops; loop++) {
    
            // Threshold volume
            for (i = 0; i < nx * ny * nz; i++)
                vol_bin[i] = (volume[i] > thresh) ? 1.0 - (float)loop : 0.0 + (float)loop;
                
            conn = conn_arr[((loop + iter) % 2)];
            sign = (loop == 0) ? 0.0 : 1.0;
            chi_error = (conn == 18) ? 2 : -6;

            for (z = 0; z < nz - 1; z++) {
                for (y = 0; y < ny - 1; y++) {
                    for (x = 0; x < nx - 1; x++) {
                            
                        // Extract 2x2x2 voxel cube
                        cube[0] = vol_bin[IDX(x, y, z, nx, ny)];
                        cube[1] = vol_bin[IDX(x+1, y, z, nx, ny)];
                        cube[2] = vol_bin[IDX(x, y+1, z, nx, ny)];
                        cube[3] = vol_bin[IDX(x+1, y+1, z, nx, ny)];
                        cube[4] = vol_bin[IDX(x, y, z+1, nx, ny)];
                        cube[5] = vol_bin[IDX(x+1, y, z+1, nx, ny)];
                        cube[6] = vol_bin[IDX(x, y+1, z+1, nx, ny)];
                        cube[7] = vol_bin[IDX(x+1, y+1, z+1, nx, ny)];
                        
                        // Count vertices, edges, faces
                        V = 0; E = 0; F = 0;
        
                        // Count vertices (fully occupied corners)
                        for (i = 0; i < 8; i++) V += cube[i];
                        
                        if (conn == 18) {
                            for (i = 0; i < 18; i++) 
                                if (cube[edges18[i][0]] && cube[edges18[i][1]]) E++;
                        } else if (conn == 26) {
                            for (i = 0; i < 26; i++) 
                                if (cube[edges26[i][0]] && cube[edges26[i][1]]) E++;
                        }
        
                        for (i = 0; i < 6; i++)
                            if (cube[faces[i][0]] && cube[faces[i][1]] && cube[faces[i][2]] && cube[faces[i][3]]) F++;
        
                        // Compute local Euler number
                        chi_local = V - E + F;
        
                        if (vol_bin[IDX(x, y, z, nx, ny)] > 0) {
                            if (chi_local == chi_error) {
                                vol_euler[IDX(x, y, z, nx, ny)] = sign;
                                vol_euler[IDX(x+1, y, z, nx, ny)] = sign;
                                vol_euler[IDX(x, y+1, z, nx, ny)] = sign;
                                vol_euler[IDX(x+1, y+1, z, nx, ny)] = sign;
                                vol_euler[IDX(x, y, z+1, nx, ny)] = sign;
                                vol_euler[IDX(x+1, y, z+1, nx, ny)] = sign;
                                vol_euler[IDX(x, y+1, z+1, nx, ny)] = sign;
                                vol_euler[IDX(x+1, y+1, z+1, nx, ny)] = sign;
                                n_errors++;
                            }
                        }
                    }
                }
            }
        }
        printf("%d\n", n_errors);
        for (i = 0; i < dims[0]*dims[1]*dims[2]; i++)
            volume[i] = vol_euler[i];
            
        if (n_errors == 0) break;
    }
    
    free(vol_euler);
    free(vol_bin);

}



/* Main program */
int main(int argc, char *argv[])
{
        char *infile, *outfile;
        int i, j, dims[3], conn_arr[2], n_loops;
        float *input;
        double separations[3];
        nifti_image *nii_ptr;
        
        if (argc < 2) {
                 (void) fprintf(stderr, "\nUsage: %s [options] in.nii [out.nii]\nSpatial adaptive non-local means denoising filter.\n\n", argv[0]);
                 (void) fprintf(stderr, "         %s -help\n\n", argv[0]);
         exit(EXIT_FAILURE);
        }
        
        infile  = argv[1];
        outfile = argv[2];
        
        /* read first image to get image parameters */
        nii_ptr = read_nifti_float(infile, &input, 0);
        if (nii_ptr == NULL) {
                fprintf(stderr,"Error reading %s.\n", infile);
                return(EXIT_FAILURE);
        }

        separations[0] = nii_ptr->dx;
        separations[1] = nii_ptr->dy;
        separations[2] = nii_ptr->dz;
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;

        conn_arr[0] = 18; conn_arr[1] = 26;
        n_loops = 1;
        correct_topology(input, 0.5, dims, conn_arr, n_loops);

        if (!write_nifti_float(outfile, input, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
                exit(EXIT_FAILURE);

        free(input);
        
        return(EXIT_SUCCESS);

}