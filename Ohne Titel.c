/**
 * apply correction that removes high-intensity values that were found in the 
 * outer CSF rim. This rim is estimated by checking those areas that are removed
 * after eroding label mask. This helps in removing meninges/dura
 */
void
correct_outer_rim(float *src, unsigned char *label, int *dims, double *voxelsize)
{
    int i, j, nvol, n[MAX_NC], n_classes = 3, replace = 0;
    float *outer_rim;
    double max_val, CSF_ratio;
    unsigned char *mask;

    nvol = dims[0]*dims[1]*dims[2];

    outer_rim = (float *)malloc(sizeof(float)*nvol);
    mask = (unsigned char *)malloc(sizeof(unsigned char)*nvol);
    if (!outer_rim || !mask) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* estimate mean for each label */
    for (j = 0; j < n_classes; j++)
        n[j] = 0;

    for (i = 0; i < nvol; i++) {
        if (label[i] == 0) continue;
        n[label[i]-1]++;
    }
    CSF_ratio = (double)n[0]/(double)(n[0] + n[1] + n[2]);


    for (i = 0; i < nvol; i++) {
        outer_rim[i] = (src[i] >= GWM) ? 1.0f : 0.0f;
        mask[i]  = (src[i] > CSF) ? 1 : 0;
    }    

    /* obtain WM distance map */
    vbdist(outer_rim, mask, dims, voxelsize, replace);
    
    /* create initial outer rim using label */
    for (i = 0; i < nvol; i++) outer_rim[i] = (label[i] > 0) ? 0.0 : 1.0;

    vbdist(outer_rim, NULL, dims, voxelsize, 0);
    for (i = 0; i < nvol; i++)
        outer_rim[i] = ((outer_rim[i] < (27.0*CSF_ratio)) && (label[i] == CSF)) ? sqrt(outer_rim[i]) : 0.0;

    for (i=0; i < nvol; i++)
        max_val = MAX((double)outer_rim[i], max_val);
        fprintf(stderr,"%g\n",CSF_ratio);
    
    /* assume that larger values in outer rim should be rather CSF and don't have that high
     * intensities and set these therefore to zero */
    for (i = 0; i < nvol; i++) src[i] = (outer_rim[i] == 0.0) ? src[i] : (src[i]*outer_rim[i]/(float)max_val);

    free(outer_rim);
    free(mask);
}

/**
 * apply correction that removes high-intensity values that were found in the 
 * outer CSF rim. This rim is estimated by checking those areas that are removed
 * after eroding label mask. This helps in removing meninges/dura
 */
void
correct_outer_rim_old(float *src, unsigned char *label, int *dims, double *voxelsize)
{
    int i, j, nvol, n[MAX_NC], n_classes = 3;
    float *outer_rim;
    double mu[MAX_NC], threshold;

    nvol = dims[0]*dims[1]*dims[2];
    int erosion_steps = 3;

    outer_rim = (float *)malloc(sizeof(float)*nvol);
    if (!outer_rim) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* estimate mean for each label */
    for (j = 0; j < n_classes; j++) {
        n[j] = 0;
        mu[j] = 0.0;
    }
    for (i = 0; i < nvol; i++) {
        if (label[i] == 0) continue;
        n[label[i]-1]++;
        mu[label[i]-1] += (double)src[i];
    }
    for (j = 0; j < n_classes; j++) mu[j] /= n[j];

    /* create initial outer rim using label */
    for (i = 0; i < nvol; i++) outer_rim[i] = (label[i] > 0) ? 1.0 : 0.0;

    /* estimate outer rim of CSF that changed by erosion steps and fill it with
       input values inside that mask */
    morph_erode(outer_rim, dims, erosion_steps, 0, DT_FLOAT32);
    
    /* fill with original values if it's CSF and was removed by erosion */
    for (i = 0; i < nvol; i++)
        outer_rim[i] = ((label[i] == CSF) && (outer_rim[i] == 0)) ? src[i] : 0.0;
    
    fprintf(stderr,"%g %g %g\n",mu[0],mu[1],mu[2]);
    threshold = (mu[0] + mu[1])/2.0;

    /* assume that larger values in outer rim should be rather CSF and don't have that high
     * intensities and set these therefore to zero */
    for (i = 0; i < nvol; i++) src[i] = (outer_rim[i] > threshold) ? 0.9*mu[0] : src[i];

    free(outer_rim);

}
