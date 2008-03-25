/* CAT_GLM_Estimate.c                                                        */
/*                                                                           */
/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <CAT_Pinv.h>
#include <time_stamp.h>

#define VERBOSE 0
#define EPS 1e-15
#define DEBUG 1

/* Definitions for accessing information on each dimension */
#define NUMBER_OF_DIMENSIONS 3		/* Number of dimensions of volume */
#define MAX_NUMBER_OF_FILES 1500	/* Number of possible files */

/* Types */
typedef struct {
   long nslices;             		/* Number of slices */
   long nrows;               		/* Number of rows in inputdata */
   long ncolumns;            		/* Number of columns in inputdata */
   double maximum;           		/* Volume maximum */
   double minimum;           		/* Volume minimum */
   double step[NUMBER_OF_DIMENSIONS];	/* Step sizes for dimensions */
   double start[NUMBER_OF_DIMENSIONS];	/* Start positions for dimensions */
   char dimension_names[NUMBER_OF_DIMENSIONS][MAX_NC_NAME]; /* Dimension names */
} Volume_Info;

void  usage(
    STRING   executable )
{    
    STRING   usage_str = "\n\
NAME\n\
	CAT_GLM_Estimate - estimation of a General Linear Model (GLM)\n\n\
SYNOPSIS\n\
	CAT_GLM_Estimate file1_group1 file2_group1 ... + file1_group2 file2_group2 ... : covariate_file \n\n\
DESCRIPTION\n\
	CAT_GLM_Estimate estimates the beta parameters for the GLM and writes the resulting parameters\n\
	and the t-values for each column of the design matrix to the current working directory.\n\
	The files are named beta_xxxx and T_xxxx, where xxxx are numbered according to the corresponding\n\
	column of the design matrix. Additionally the residual standard deviation is saved in a file\n\
	ResSD_xxxx.\n\
	The input files can be either in minc-format (using the extension .mnc) or in ascii-format\n\
	(using .txt as extension) and the respective result files are saved in the same format and\n\
	with the same extension.\n\
	The basic model at each element of the input files is of the form Y = X*B + e, for data Y,\n\
	design matrix X, (unknown) parameters B, and residual errors e. The errors are assumed to have\n\
	normal distribution.\n\
	Covariates for regression (correlation) or AnCova models can be defined using ':' as delimiter.\n\
	Each covariate should be defined as seperate file.\n\n\
	The following files are written:\n\
	ResMS       - residual mean square\n\
	beta_xxxx	- parameter estimates, which are numbered according to the corresponding\n\
				  column of the design matrix.\n\
	T_xxxx		- T-values, which are numbered according to the corresponding column of\n\
				  the design matrix. These values are the parameter estimates divided by\n\
				  residual standard deviation.\n\
	ResSD_xxxx	- estimated residual standard deviation, which are numbered according to\n\
				  the corresponding column of the design matrix.\n\
	These files can be used with glm_mat to calculate different contrasts of factor levels.\n";

    print_error( usage_str, executable );
}

public void normalize_vector(
    Real *data,
    int length )
{
    Real mean;
    int i;
    
    mean = 0.0;
    for (i=0; i<length; i++)
        mean += data[i];
    mean = mean/length;

    for (i=0; i<length; i++)
        data[i] -= mean;
        
}

/* --------------------------------------------------------------------------------------------------------- */
/*  estimation of GLM */
/* --------------------------------------------------------------------------------------------------------- */
int  estimate(
    char **infiles,
    char *arg_string,
    int argc )
{
    char                 *outfile;
    char                 buffer[EXTREMELY_LARGE_STRING_SIZE];
    int                  i, j, k, n_subj, n_values, n_tmp, prev_n_values;
    int                  n_beta, rank, erdf, counter;
    int                  input_icvid[MAX_NUMBER_OF_FILES], output_beta_icvid[MAX_NUMBER_OF_FILES];
    int                  output_T_icvid[MAX_NUMBER_OF_FILES], output_residual_icvid[MAX_NUMBER_OF_FILES];
    int                  output_ResMS_icvid, *index;
    int                  islice, n_dims, n_cov;
    float                *outputdata;
    Real                 **values, *values_tmp, **G, *v, **inv_G, **transp_G, **pinv_GG;
    Real                 *data, **inputdata, **beta, **estimates, **ResSD, sum;
    STRING               input_filename, output_filename;
    FILE                 *file;
    progress_struct      progress;
    Volume_Info          volume_info[MAX_NUMBER_OF_FILES];
      
    if(     filename_extension_matches( infiles[0], "txt" ) ) n_dims = 1;
    else if(filename_extension_matches( infiles[0], "mnc" ) ) n_dims = 3;
    else {
        print_error( "\nError: unsupported file format\n");
        return( 1 );
    }    

    counter = 0;
    n_beta = 1;
    n_cov = 0;
    
    ALLOC(index, argc-1);

    for (i=0; i<argc-1; i++) {
        if( equal_strings( infiles[i], "+" ) ) {
            n_beta++;
        }
        else if( equal_strings( infiles[i], ":" ) ) {
            n_beta++;
            n_cov++;
            /* skip the file containing covariates by increasing i */
            i++;
        }
        else {
            index[counter] = i;
            /* check for files */
            if( !file_exists( infiles[i] ) ) {
                print_error( "\nFile %s not found.\n", infiles[i]);
                return( 1 );                
            }
            /* read input icvid for minc files */
            if( n_dims == 3 )
                input_icvid[counter] = get_volume_info(infiles[i], &volume_info[i], NC_DOUBLE);

            counter++;
            }
    }
    
    n_subj = argc - n_cov - n_beta;
    
    ALLOC2D(G, n_subj, n_beta);
    ALLOC2D(inv_G, n_beta, n_subj);

    /* initialize design matrix G to zero */
    for (j=0; j < n_beta; j++)
        for (i=0; i < n_subj; i++)
	       G[i][j] = 0;     

    j = 0;
    counter = 0;

    /* read data and build design matrix */
    if( n_dims == 1 ) print("Read files:"); fflush(stdout);
    
    for (i=0; i < argc-1; i++) {
        /* define covariates */
        if( equal_strings( infiles[i], ":" ) ) {
            /* count columns and files */
            j++; i++;
            if( input_texture_values( infiles[i], &n_tmp, &values_tmp ) != OK ) {
                print_error( "\nError reading file %s\n", infiles[i]);
                return( 1 );
            }
                if( n_subj != n_tmp) {
                print_error( "\nError: Number of files differs from number of rows in design matrix.\n");
	            return( 1 );
	       }
	       /* normalize covariates to zero mean */
	       normalize_vector(values_tmp, n_tmp);
	       
	       /* build design matrix G */
	       for (k=0; k < n_subj; k++) {
	           G[k][j] = values_tmp[k];
	       }        
	    }
	    /* define groups */
        else if( equal_strings( infiles[i], "+" ) ) {
            j++ ;
        }
        else if( n_dims == 1 ) {
	       if( input_texture_values( infiles[i], &n_values, &values_tmp ) != OK ) {
	           print_error( "\nError reading file %s\n", infiles[i]);
	           return( 1 );
	       } /* if */
	       if( i == 0) ALLOC2D(values, n_subj, n_values);
	       else {
	           if( prev_n_values != n_values) {
		          print_error("\nError: Wrong number of values in %s\n",infiles[i]);
		          return( 1 );
	           } /* if */
	       } /* else */
	       for (k=0; k < n_values; k++)
	           values[counter][k] = values_tmp[k];
	       prev_n_values = n_values;
	       G[counter][j] = 1;
	       counter++;
        } /* else if */
        else if( n_dims == 3 ) {
	       G[counter][j] = 1;
	       counter++;
        } /* else if */
        if( n_dims == 1 ) { print( "."); fflush(stdout); }
    } /* i++ */

    if( n_dims == 1 ) print( "\n");

    if( n_dims == 3)
    	n_values = volume_info[0].nrows * volume_info[0].ncolumns;
	
    /* compute pseudo inverse from design matrix */        
    rank = pinv(n_subj, n_beta, G, inv_G);
    
    /* effective residual d.f. */        
    erdf = n_subj - rank;

    print("Effective residual d.f.: %d\n",erdf);
    
    /* check estimability */
    if( erdf < 0 ) {
        print_error( "This design is unestimable! (df=%d).\n",erdf );
        return( 1 );
    }
    
    if( erdf == 0 ) {
        print_error( "This design has no res! (df=0).\n" );
        return( 1 );
    }

    /* open log file */
    if((file = fopen("glm.log", "w")) == 0) {
        fprintf(stderr, "Couldn't open file glm.log.\n");
        return(1);
    }
    
    fprintf(file, "[df]\n%d\n", erdf);
    fprintf(file, "\n[design matrix]\n");
    for (i=0; i < (n_subj); i++) {
        if (VERBOSE) print("%s\t", infiles[index[i]]);
        fprintf(file,"%s\t", infiles[index[i]]);          
        for (j=0; j < n_beta; j++) {
            if (VERBOSE) print("%6.3f ", G[i][j]);
            fprintf(file, "%6.3f ", G[i][j]);
        }
        if (VERBOSE) print("\n");
        fprintf(file, "\n");
    }
    
    fprintf(file, "\n[history]\n%s\n", arg_string);
    fclose(file);

    ALLOC(data, n_values);
    ALLOC(outputdata, n_values);
    ALLOC(v, n_values);
    ALLOC2D(ResSD, n_beta, n_values);
    ALLOC2D(inputdata, n_subj, n_values);
    ALLOC2D(estimates, n_subj, n_values);
    ALLOC2D(beta, n_beta, n_values);
    ALLOC2D(pinv_GG, n_beta, n_subj);
    ALLOC2D(transp_G, n_beta, n_subj);

    /* calculate pinv(G'*G) */
    transpose(n_subj, n_beta, G, transp_G);
    matrix_multiply(n_beta, n_subj, n_beta, transp_G, G, pinv_GG);
    (void) pinv(n_beta, n_beta, pinv_GG, pinv_GG);

    /* --------------------------------------------------------------------------------------------------------- */
    /*  estimation for minc-files*/
    /* --------------------------------------------------------------------------------------------------------- */
    if( n_dims == 3) {
        /* Check volume information about dimensions, voxel size and origin */
        for (i=1; i < n_subj; i++) {
            if( volume_info[index[i]].nslices != volume_info[0].nslices ||
        	volume_info[index[i]].nrows != volume_info[0].nrows ||
        	volume_info[index[i]].ncolumns != volume_info[0].ncolumns ) {
                    print_error( "Dimensions of file %s differ.\n", infiles[index[i]] );		
                    return( 1 );
        	}
            if( volume_info[index[i]].step[0] != volume_info[0].step[0] ||
        	volume_info[index[i]].step[1] != volume_info[0].step[1] ||
        	volume_info[index[i]].step[2] != volume_info[0].step[2] ) {
                    print_error( "Voxelsize (step value) of file %s differs.\n", infiles[index[i]] );		
                    return( 1 );
            }
            if( volume_info[index[i]].start[0] != volume_info[0].start[0] ||
        	volume_info[index[i]].start[1] != volume_info[0].start[1] ||
        	volume_info[index[i]].start[2] != volume_info[0].start[2] ) {
                    print_error( "Origin (start value) of file %s differs.\n", infiles[index[i]] );		
                    return( 1 );
            }
        } /* i++ */
        
        /* prepare output for ResMS */
        outfile = create_string("ResMS.mnc");
        output_ResMS_icvid = save_volume_info(input_icvid[0], outfile, arg_string, 
                                   &volume_info[0], NC_FLOAT);

        /* prepare output files for beta-, T-, and ResSD-images */
        for (j=0; j < n_beta; j++) {
            outfile = create_string("beta");
            (void) sprintf( buffer, "_%04d.mnc", j+1 );
            concat_to_string( &outfile, buffer );
            output_beta_icvid[j] = save_volume_info(input_icvid[0], outfile, arg_string, 
                                   &volume_info[0], NC_FLOAT);

            outfile = create_string("T");
            (void) sprintf( buffer, "_%04d.mnc", j+1 );
            concat_to_string( &outfile, buffer );
            output_T_icvid[j] = save_volume_info(input_icvid[0], outfile, arg_string, 
                                   &volume_info[0], NC_FLOAT);

            outfile = create_string("ResSD");
            (void) sprintf( buffer, "_%04d.mnc", j+1 );
            concat_to_string( &outfile, buffer );
            output_residual_icvid[j] = save_volume_info(input_icvid[0], outfile, arg_string, 
                                   &volume_info[0], NC_FLOAT);
        }  /* j++ */

        initialize_progress_report(&progress, FALSE, volume_info[0].nslices, "Estimate GLM");

        /* Loop through slices */
        for (islice=0; islice < volume_info[0].nslices; islice++) {
            /* get input data */	
            for (i=0; i < n_subj; i++) {
                get_volume_slice_double(input_icvid[i], &volume_info[index[i]], islice, data);
                for (k=0; k < n_values; k++)
                     inputdata[i][k] = data[k];
            } /* i++ */
	
           /* calculate beta by multiplying pseudo inverse from design matrix with transposed data */	
            matrix_multiply(n_beta, n_subj, n_values, inv_G, inputdata, beta);

           /* calculate fitted data: estimates = G*beta */
            matrix_multiply(n_subj, n_beta, n_values, G, beta, estimates);

            /* calculate estimated residual standard deviation: ResSD = sqrt(sum(res.^2)/erdf) */
            for (k=0; k < n_values; k++) {
                sum = 0.0;
                for (i=0; i < n_subj; i++) {
        	       estimates[i][k] = inputdata[i][k] - estimates[i][k];
        	       sum += estimates[i][k] * estimates[i][k];
                }  /* i++ */
        	    v[k] = sum/ (Real) erdf;
            }  /* k++ */

            for (k=0; k < n_values; k++)
                outputdata[k] = (float) v[k];
            save_volume_slice_float(output_ResMS_icvid, &volume_info[0], islice, outputdata);

            for (j=0; j < n_beta; j++)
                for (k=0; k < n_values; k++)
	               ResSD[j][k] = sqrt(v[k]*pinv_GG[j][j]);

            /* save slices */
            for (j=0; j < n_beta; j++) {
                for (k=0; k < n_values; k++)
        	         outputdata[k] = (float) ResSD[j][k];
                save_volume_slice_float(output_residual_icvid[j], &volume_info[0], islice, outputdata);

                for (k=0; k < n_values; k++)
        	         outputdata[k] = (float) beta[j][k];
                save_volume_slice_float(output_beta_icvid[j], &volume_info[0], islice, outputdata);

                for (k=0; k < n_values; k++)
                    outputdata[k] = (float) beta[j][k]/(ResSD[j][k] + EPS);
                save_volume_slice_float(output_T_icvid[j], &volume_info[0], islice, outputdata);
           }  /* j++ */
            update_progress_report(&progress, islice + 1);	
        }  // islice++
        terminate_progress_report(&progress);

        /* Free up data and close files */
        close_volume(output_ResMS_icvid);
        for (i=0; i < n_subj; i++)
            close_volume(input_icvid[i]);
        for (j=0; j < n_beta; j++) {
            close_volume(output_beta_icvid[j]);
            close_volume(output_T_icvid[j]);
        } /* j++ */
    }  // n_dims == 3
    
    /* --------------------------------------------------------------------------------------------------------- */
    /*  estimation for ascii-files*/
    /* --------------------------------------------------------------------------------------------------------- */

    if( n_dims == 1 ) {
    /* multiply pseudo inverse from design matrix with transposed data */	
        for (k=0; k < n_values; k++)
        {
            for (j=0; j < n_beta; j++)
	    {
	        sum = 0.0;
	        for (i=0; i < n_subj; i++)
                    sum += inv_G[j][i] * values[i][k];		
                beta[j][k] = sum;
	    } /* j++ */
        } /* k++ */

        /* write betas for each column of design matrix */
        for (j=0; j < n_beta; j++) {
            output_filename = create_string("beta");
            (void) sprintf( buffer, "_%04d.txt", j+1 );
            concat_to_string( &output_filename, buffer );
	
            if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
                return( 1 );
            for (k=0; k < n_values; k++) {
                if( output_real( file, beta[j][k] ) != OK ||
                    output_newline( file ) != OK )
                    break;
            } /* k++ */
            (void) close_file( file );	    
        } /* j++ */    

        /* calculate fitted data: estimates = G*beta */
        matrix_multiply(n_subj, n_beta, n_values, G, beta, estimates);

        /* calculate estimated squared residual standard deviation: v = sum((values - estimates).^2)/erdf */
        for (k=0; k < n_values; k++) {
            sum = 0.0;
            for (i=0; i < n_subj; i++) {
	           estimates[i][k] = values[i][k] - estimates[i][k];
	           sum += estimates[i][k]*estimates[i][k];
	       } /* i++ */
	       v[k] = sum/erdf;
        } /* k++ */

        output_filename = create_string("ResMS.txt");	
        if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
            return( 1 );
        for (k=0; k < n_values; k++) {
            if( output_real( file, v[k]) != OK ||
                output_newline( file ) != OK )
                break;
        } /* k++ */
        (void) close_file( file );	    

        /* write beta and beta/ResSD for each column of design matrix */
        for (j=0; j < n_beta; j++) {
	        /* ResSD values*/
	        output_filename = create_string("ResSD");
            (void) sprintf( buffer, "_%04d.txt", j+1 );
            concat_to_string( &output_filename, buffer );
	
            if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
                return( 1 );
            for (k=0; k < n_values; k++) {
                if( output_real( file, sqrt(v[k]*pinv_GG[j][j])) != OK ||
                    output_newline( file ) != OK )
                    break;
            } /* k++ */
            (void) close_file( file );	    

	    /* T-values*/
            output_filename = create_string("T");
            (void) sprintf( buffer, "_%04d.txt", j+1 );
            concat_to_string( &output_filename, buffer );
	
            if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
                return( 1 );
            for (k=0; k < n_values; k++) {
                if( output_real( file, beta[j][k]/((sqrt(v[k]*pinv_GG[j][j])) + EPS)) != OK ||
                    output_newline( file ) != OK )
                    break;
            } /* k++ */
            (void) close_file( file );	
        } /* j++ */
    } // n_dims == 1
    
    FREE(data);
    FREE(outputdata);
    FREE(v);
    FREE(index);
    FREE2D(ResSD);
    FREE2D(inputdata);
    FREE2D(estimates);
    FREE2D(transp_G);
    FREE2D(pinv_GG);
    FREE2D(beta);
   
    return( 0 );

}

int  main(
    int   argc,
    char  *argv[] )
{
    char  *arg_string, **infiles;

    /* Save time stamp and args */
    arg_string = time_stamp(argc, argv);
    
    initialize_argument_processing( argc, argv );

    infiles = &argv[1];

    /* get filenames */
    if( argc < 2 ) {
        usage( argv[0] );
        return( 1 );
    }
    estimate(infiles, arg_string, argc);
}
