/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

public Real* read_pgm(
    char *fname,
    int *nx,
    int *ny);

public int write_pgm(
    char *fname,
    Real *data,
    int nx,
    int ny);
