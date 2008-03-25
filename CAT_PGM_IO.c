/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>

public Real* read_pgm(
    char *fname,
    int *nx,
    int *ny)
{
    FILE    *fp;
    int     i, size, max_val;
    char    line[256];
    unsigned char *data_char;
    Real *data;
            
    if((fp = fopen(fname, "r")) == 0) {
        fprintf(stderr, "read_pgm: Couldn't open file %s.\n", fname);
        *nx = *ny = 0;
        return(NULL);
    }
    
    /* read PGM header */
    fgets(line, 256, fp);
    if (strncmp(line,"P5", 2)) {
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

    size = (*nx)*(*ny);
    ALLOC(data_char, size);

    if (fread(data_char, sizeof(unsigned char), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        *nx = *ny = 0;
        return(NULL);
    }
         
    ALLOC(data, size);   
    for ( i=0; i < size; i++ ) {
        data[i] = (Real) data_char[i];
    }

    FREE(data_char);
    fclose(fp);
    return(data);
}

public int write_pgm(
    char *fname,
    Real *data,
    int nx,
    int ny)
{
    FILE    *fp;
    Real    min_value, max_value, scale, offset;
    int     i;
    char    *data_char;
    
    ALLOC(data_char, nx*ny);
    
    min_value = 1e15;
    max_value = -1e15;

    /* calculate min/max */
    for ( i=0; i < nx*ny; i++ ) {
        if( data[i] < min_value ) min_value = data[i];
        if( data[i] > max_value ) max_value = data[i];
    }

    offset = min_value;
    scale = 255 / (max_value - min_value);
    /* scale image */
    for ( i=0; i < nx*ny; i++ ) {
        data_char[i] = (char) ROUND((data[i] - offset) * scale);
    }
    
    if((fp = fopen(fname, "w")) == 0) {
        fprintf(stderr, "write_pgm: Couldn't open file %s.\n", fname);
        return(1);
    }
    
    /* write PGM header */
    fprintf(fp, "P5\n");
    /* modification to include scaling and offset */
    fprintf(fp, "# Scale: %3.5f Offset: %3.5f\n",scale, offset);
    fprintf(fp, "%d %d\n",nx,ny);
    fprintf(fp, "255\n");
    
    if (fwrite(data_char, sizeof(char), nx*ny, fp) != nx*ny) {
        fprintf(stderr, "Error writing data.\n");
        return(1);
    }
    
    fclose(fp);
    return(0);
}
