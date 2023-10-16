/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Vol.h"

/* estimate minimum of A and its index in A */
void
pmin(float *A, int sA, float *minimum, int *index)
{
    int i; 
    
    *minimum = FLT_MAX;
    *index = 0;
    
    for (i=0; i<sA; i++) {
        if ((A[i] > 0.0) && (*minimum > A[i])) { 
            *minimum = A[i]; 
            *index   = i;
        }
    }
}

/* subfunction for SBT to get all values of the voxels which are in WMD-range (children of this voxel) */
float
pmax(const float GMT[], const float PPM[], const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA) {
    float n=0.0, maximum=WMD;
    int i;

    /* the pure maximum */
    for (i=0; i<=sA; i++) {
        if ((GMT[i] < FLT_MAX) && (maximum < GMT[i]) &&              /* thickness/WMD of neighbors should be larger */
                (SEG[i] >= 1.0) && (SEGI>1.2 && SEGI<=2.75) &&       /* projection range */
                (((PPM[i] - ND[i] * 1.2) <= WMD)) &&                 /* upper boundary - maximum distance */
                (((PPM[i] - ND[i] * 0.5) >  WMD) || (SEG[i]<1.5)) && /* lower boundary - minimum distance - corrected values outside */
                ((((SEGI * MAX(1.0,MIN(1.2,SEGI-1.5))) >= SEG[i])) || (SEG[i]<1.5))) /* for high values will project data over sulcal gaps */
        {
            maximum = GMT[i];
        }
    }

    
    /* the mean of the highest values */
    float maximum2=maximum; float m2n=0.0; 
    for (i=0; i<=sA; i++) {
        if ((GMT[i] < FLT_MAX) && ((maximum - 1) < GMT[i]) && 
                 (SEG[i] >= 1.0) && (SEGI>1.2 && SEGI<=2.75) && 
                 (((PPM[i] - ND[i] * 1.2) <= WMD)) && 
                 (((PPM[i] - ND[i] * 0.5) >  WMD) || (SEG[i]<1.5)) &&
                 ((((SEGI * MAX(1.0,MIN(1.2,SEGI-1.5))) >= SEG[i])) || (SEG[i]<1.5))) {
            maximum2 += GMT[i]; 
            m2n++;
        } 
    }
    if (m2n > 0.0)
        maximum = (maximum2 - maximum) / m2n;

    return maximum;
}

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void
ind2sub(int i, int *x, int *y, int *z, int sxy, int sx)
{
    int j = i % sxy;

    *z = (int)floor((double)i / (double)sxy); 
    *y = (int)floor((double)j / (double)sx);   
    *x = j % sx;
}

/* 
 * Estimate index i of a voxel x,y,z in an array size s.
 * See also for ind2sub.
 */
int
sub2ind(int x, int y, int z, int s[]) {
    /* handling on boundaries */
    if (x<0) x=0; if (x>s[0]-1) x=s[0]-1; 
    if (y<0) y=0; if (y>s[1]-1) y=s[1]-1; 
    if (z<0) z=0; if (z>s[2]-1) z=s[2]-1; 
  
    /*   z * (number of voxels within a slice) 
       + y * (number of voxels in a column)  
       + x   (position within the column)   
    */   
    return z*s[0]*s[1] + y*s[0] + x;
}


/* 
 * Read out the linear interpolated value of a volume SEG with the size
 * s on the position x,y,z (c-notation). See also ind2sub for details of
 * the c-notation.
 */
float
isoval(float vol[], float x, float y, float z, int s[]){

    int i;
    float seg=0.0, n=0.0;
    float fx = floor(x),   fy = floor(y),   fz = floor(z);
    float cx = floor(x+1), cy = floor(y+1), cz = floor(z+1);
    
    float wfx = cx-x, wfy = cy-y, wfz = cz-z;
    float wcx = x-fx, wcy = y-fy, wcz = z-fz;
    float N[8], W[8];  
    
    /* value of the 8 neighbors and there distance weight */
    N[0] = vol[sub2ind((int)fx,(int)fy,(int)fz,s)];  W[0] = wfx * wfy * wfz; 
    N[1] = vol[sub2ind((int)cx,(int)fy,(int)fz,s)];  W[1] = wcx * wfy * wfz;
    N[2] = vol[sub2ind((int)fx,(int)cy,(int)fz,s)];  W[2] = wfx * wcy * wfz;
    N[3] = vol[sub2ind((int)cx,(int)cy,(int)fz,s)];  W[3] = wcx * wcy * wfz;
    N[4] = vol[sub2ind((int)fx,(int)fy,(int)cz,s)];  W[4] = wfx * wfy * wcz;
    N[5] = vol[sub2ind((int)cx,(int)fy,(int)cz,s)];  W[5] = wcx * wfy * wcz;
    N[6] = vol[sub2ind((int)fx,(int)cy,(int)cz,s)];  W[6] = wfx * wcy * wcz; 
    N[7] = vol[sub2ind((int)cx,(int)cy,(int)cz,s)];  W[7] = wcx * wcy * wcz;
    
    for (i=0; i<8; i++) {
        if (!isnan(N[i]) && isfinite(N[i])) {
            seg = seg + (N[i] * W[i]);
            n+= W[i];
        }
    }
     
    if (n>0.0) return seg/n; 
    else return FNAN;
}

void 
get_prctile(float *src, int *dims, double threshold[2], double prctile[2], int exclude_zeros)
{
    double mn_thresh, mx_thresh;
    double min_src = FLT_MAX, max_src = -FLT_MAX;
    int *cumsum, *histo;
    int i, sz_histo = 10000, nvol = dims[0]*dims[1]*dims[2];
    
    cumsum  = (int *)malloc(sizeof(int)*sz_histo);
    histo   = (int *)malloc(sizeof(int)*sz_histo);
    
    for (i=0; i<nvol; i++) {
        min_src = MIN((double)src[i], min_src);
        max_src = MAX((double)src[i], max_src);
    }

    /* build histogram */
    for (i = 0; i < sz_histo; i++) histo[i] = 0;
    for (i = 0; i < nvol; i++) {
        if ((exclude_zeros > 0) && (src[i] == 0)) continue; /* exclude zeros from histogram */
        histo[(int)round((double)sz_histo*((double)src[i]-min_src)/(max_src-min_src))]++;
    }

    /* find values defined prctile */
    cumsum[0] = histo[0];
    for (i = 1; i < sz_histo; i++) cumsum[i] = cumsum[i-1] + histo[i];
    for (i = 0; i < sz_histo; i++) cumsum[i] = (int) round(100000.0*(double)cumsum[i]/(double)cumsum[sz_histo-1]);
    for (i = 0; i < sz_histo; i++) if (cumsum[i] >= (int)round(prctile[0]*1000.0)) break;
    threshold[0] = (double)i/(double)sz_histo*(max_src-min_src) + min_src;
    for (i = sz_histo-1; i > 0; i--) if (cumsum[i] <= (int)round(prctile[1]*1000.0)) break;
    threshold[1] = (double)i/(double)sz_histo*(max_src-min_src) + min_src;
 
    free(cumsum);
    free(histo);
 
}

static void 
convxy(double out[], int xdim, int ydim, double filtx[], double filty[], int fxdim, int fydim, int xoff, int yoff, double buff[])
{
    int x,y,k;

    for (y=0; y<ydim; y++) {
        for (x=0; x<xdim; x++) {
            buff[x] = out[x+y*xdim];
            if (!isfinite(buff[x]))
                buff[x] = 0.0;
        }
        for (x=0; x<xdim; x++) {
            double sum1 = 0.0;
            int fstart, fend;
            fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
            fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

            for (k=fstart; k<fend; k++)
                sum1 += buff[x-xoff-k]*filtx[k];
            out[x+y*xdim] = sum1;
        }
    }
    for (x=0; x<xdim; x++) {
        for (y=0; y<ydim; y++)
            buff[y] = out[x+y*xdim];

        for (y=0; y<ydim; y++) {
            double sum1 = 0.0;
            int fstart, fend;
            fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
            fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

            for (k=fstart; k<fend; k++)
                sum1 += buff[y-yoff-k]*filty[k];
            out[y*xdim+x] = sum1;
        }
    }
}

static void 
convxy_float(float out[], int xdim, int ydim, double filtx[], double filty[], int fxdim, int fydim, int xoff, int yoff, float buff[])
{
    int x,y,k;

    for (y=0; y<ydim; y++) {
        for (x=0; x<xdim; x++) {
            buff[x] = out[x+y*xdim];
            if (!isfinite(buff[x]))
                buff[x] = 0.0;
        }
        for (x=0; x<xdim; x++) {
            double sum1 = 0.0;
            int fstart, fend;
            fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
            fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

            for (k=fstart; k<fend; k++)
                sum1 += (double)buff[x-xoff-k]*filtx[k];
            out[x+y*xdim] = (float)sum1;
        }
    }
    for (x=0; x<xdim; x++) {
        for (y=0; y<ydim; y++)
            buff[y] = out[x+y*xdim];

        for (y=0; y<ydim; y++) {
            double sum1 = 0.0;
            int fstart, fend;
            fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
            fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

            for (k=fstart; k<fend; k++)
                sum1 += (double)buff[y-yoff-k]*filty[k];
            out[y*xdim+x] = (float)sum1;
        }
    }
}


int 
convxyz_double(double *iVol, double filtx[], double filty[], double filtz[],
    int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
    double *oVol, int dims[3])
{
    double *tmp = NULL, *buff = NULL, **sortedv = NULL, *obuf;
    int xy, z, y, x, k, fstart, fend, startz, endz;
    int xdim, ydim, zdim;

    xdim = dims[0];
    ydim = dims[1];
    zdim = dims[2];

    tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
    buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
    sortedv = (double **)malloc(sizeof(double *)*fzdim);

    if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
    endz     = zdim+fzdim+zoff-1;

    for (z=startz; z<endz; z++) {
        double sum2 = 0.0;

        if (z >= 0 && z<zdim) {
            for (y=0; y<ydim; y++) for (x=0; x<xdim; x++)
                tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];  
            convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
                filtx, filty, fxdim, fydim, xoff, yoff, buff);
        }
        if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim) {
            fstart = ((z >= zdim) ? z-zdim+1 : 0);
            fend = ((z-fzdim < 0) ? z+1 : fzdim);

            for (k=0; k<fzdim; k++) {
                int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
                sortedv[k] = &(tmp[z1*xdim*ydim]);
            }

            for (k=fstart, sum2=0.0; k<fend; k++)
                sum2 += filtz[k];

            obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
            if (sum2) {
                for (xy=0; xy<xdim*ydim; xy++) {
                    double sum1=0.0;
                    for (k=fstart; k<fend; k++)
                        sum1 += filtz[k]*sortedv[k][xy];

                    obuf[xy] = sum1/sum2;
                }
            }
            else
                for (xy=0; xy<xdim*ydim; xy++)
                    obuf[xy] = 0.0;
        }
    }
    free(tmp);
    free(buff);
    free(sortedv);
    return(0);
}

int 
convxyz_float(float *iVol, double filtx[], double filty[], double filtz[],
    int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
    float *oVol, int dims[3])
{
    float *tmp = NULL, *buff = NULL, **sortedv = NULL, *obuf;
    int xy, z, y, x, k, fstart, fend, startz, endz;
    int xdim, ydim, zdim;

    xdim = dims[0];
    ydim = dims[1];
    zdim = dims[2];

    tmp = (float *)malloc(sizeof(float)*xdim*ydim*fzdim);
    buff = (float *)malloc(sizeof(float)*((ydim>xdim) ? ydim : xdim));
    sortedv = (float **)malloc(sizeof(float *)*fzdim);

    if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
    endz     = zdim+fzdim+zoff-1;

    for (z=startz; z<endz; z++) {
        double sum2 = 0.0;

        if (z >= 0 && z<zdim) {
            for (y=0; y<ydim; y++) for (x=0; x<xdim; x++)
                tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];  
                convxy_float(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
                        filtx, filty, fxdim, fydim, xoff, yoff, buff);
        }
        if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim) {
            fstart = ((z >= zdim) ? z-zdim+1 : 0);
            fend = ((z-fzdim < 0) ? z+1 : fzdim);

            for (k=0; k<fzdim; k++) {
                int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
                sortedv[k] = &(tmp[z1*xdim*ydim]);
            }

            for (k=fstart, sum2=0.0; k<fend; k++)
                sum2 += filtz[k];

            obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
            if (sum2) {
                for (xy=0; xy<xdim*ydim; xy++) {
                    double sum1=0.0;
                    for (k=fstart; k<fend; k++)
                        sum1 += filtz[k]*(double)sortedv[k][xy];

                    obuf[xy] = (float)sum1/sum2;
                }
            }
            else
                for (xy=0; xy<xdim*ydim; xy++)
                    obuf[xy] = 0.0;
        }
    }
    free(tmp);
    free(buff);
    free(sortedv);
    return(0);
}

int
convxyz_uint8(unsigned char *iVol, double filtx[], double filty[], double filtz[],
    int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
    unsigned char *oVol, int dims[3])
{
    double *tmp = NULL, *buff = NULL, **sortedv = NULL;
    int xy, z, y, x, k, fstart, fend, startz, endz;
    int xdim, ydim, zdim;
    double tmp2;
    unsigned char *obuf;

    xdim = dims[0];
    ydim = dims[1];
    zdim = dims[2];

    tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
    buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
    sortedv = (double **)malloc(sizeof(double *)*fzdim);

    if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
    endz     = zdim+fzdim+zoff-1;

    for (z=startz; z<endz; z++) {
        double sum2 = 0.0;

        if (z >= 0 && z<zdim) {
            for (y=0; y<ydim; y++) for (x=0; x<xdim; x++)
                tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];  
            convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
                filtx, filty, fxdim, fydim, xoff, yoff, buff);
        }
        if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim) {
            fstart = ((z >= zdim) ? z-zdim+1 : 0);
            fend = ((z-fzdim < 0) ? z+1 : fzdim);

            for (k=0; k<fzdim; k++) {
                int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
                sortedv[k] = &(tmp[z1*xdim*ydim]);
            }

            for (k=fstart, sum2=0.0; k<fend; k++)
                sum2 += filtz[k];

            obuf = oVol;
            obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
            if (sum2) {
                for (xy=0; xy<xdim*ydim; xy++)
                {
                    double sum1=0.0;
                    for (k=fstart; k<fend; k++)
                        sum1 += filtz[k]*sortedv[k][xy];
                    tmp2 = sum1/sum2;
                    if (tmp2<0.0) tmp2 = 0.0;
                    else if (tmp2>255.0) tmp2 = 255.0;
                    obuf[xy] = (unsigned char)RINT(tmp2);
                }
            }
            else
                for (xy=0; xy<xdim*ydim; xy++)
                    obuf[xy] = 0;
        }
    }
    free(tmp);
    free(buff);
    free(sortedv);
    return(0);
}

/* 
    [Ygmt,Ypp] = projection_based_thickness(Yp0, dist_WM, dist_CSF); 
 */
void
projection_based_thickness(float *SEG, float *WMD, float *CSFD, float *GMT, int dims[3], double *voxelsize) 
{     
    /* main information about input data (size, dimensions, ...) */
    const int   nvol = dims[0]*dims[1]*dims[2];
    const int   x  = dims[0];
    const int   y  = dims[1];
    const int   xy = x*y;
    
    const float s2 = sqrt(2.0);
    const float s3 = sqrt(3.0);
    
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    const int   NI[] = {0, -1, -x+1,- x, -x-1, -xy+1, -xy, -xy-1, -xy+x+1, -xy+x, -xy+x-1, -xy-x+1, -xy-x, -xy-x-1};    
    const float ND[] = {0.0, 1.0, s2, 1.0, s2, s2, 1.0, s2, s3, s2, s3, s3, s2, s3};
    
    const int   sN = sizeof(NI)/sizeof(NI[0]); /* division by 4 to get from the number of bytes to the number of elements */ 
    float DN[sN], DI[sN], GMTN[sN], WMDN[sN], SEGN[sN], DNm;
    
    int   i,n,ni,u,v,w,nu,nv,nw, count_WM=0, count_CSF=0;
        
    /* initialisiation */
    for (i=0; i<nvol; i++) {
        GMT[i] = WMD[i] + 0.0;
        
        /* proof distance input */
        if (SEG[i] >= GWM) count_WM++;
        if (SEG[i] <= CGM) count_CSF++;
    }

    if (count_WM == 0) {
        fprintf(stderr,"ERROR: no WM voxels\n");
        exit(EXIT_FAILURE);
    }    
    if (count_CSF == 0) {
        fprintf(stderr,"ERROR: no CSF voxels\n");
        exit(EXIT_FAILURE);
    }    
    
    /* thickness calculation
     ======================================================================= */
    for (i = 0; i < nvol; i++) {
        if (SEG[i] > CSF && SEG[i] < WM) {
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++) {
                ni = i + NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);

                if ((ni < 0) || (ni >= nvol) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;

                GMTN[n] = GMT[ni];
                WMDN[n] = WMD[ni];
                SEGN[n] = SEG[ni];
            }

            /* find minimum distance within the neighborhood */
            DNm = pmax(GMTN, WMDN, SEGN, ND, WMD[i], SEG[i], sN);
            GMT[i] = DNm;
        }
    }

    for (i = nvol - 1; i >= 0; i--) {
        if (SEG[i] > CSF && SEG[i] < WM) {
            ind2sub(i, &u, &v, &w, xy, x);

            float GMTN[sN], WMDN[sN], SEGN[sN];

            /* read neighbour values */
            for (n = 0; n < sN; n++) {
                ni = i - NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);

                if ((ni < 0) || (ni >= nvol) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;

                GMTN[n] = GMT[ni];
                WMDN[n] = WMD[ni];
                SEGN[n] = SEG[ni];
            }

            /* find minimum distance within the neighborhood */
            DNm = pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN);
            if ((GMT[i] < DNm) && (DNm > 0)) GMT[i] = DNm;
        }
    }

    for (i = 0; i < nvol; i++) {
        if (SEG[i] < CGM || SEG[i] > GWM)
            GMT[i] = 0.0;
    }

    float GMTi, CSFDi;
    for (i = 0; i < nvol; i++) {
        if (SEG[i] >= CGM && SEG[i] <= GWM) {
            GMTi = CSFD[i] + WMD[i];
            CSFDi = GMT[i] - WMD[i];

            if (CSFD[i] <= CSFDi)
                GMT[i] = GMTi;
        }
    }
}

void
vbdist(float *V, unsigned int *M, int dims[3], double *voxelsize) 
{
    
    /* main information about input data (size, dimensions, ...) */
    const int nvol = dims[0]*dims[1]*dims[2];
    const int x    = dims[0];
    const int y    = dims[1];
    const int xy   = x*y;
    
    float s1 = (float)fabs(voxelsize[0]);
    float s2 = (float)fabs(voxelsize[1]);
    float s3 = (float)fabs(voxelsize[2]);
    const float s12  = (float) sqrt((double)s1*s1   + s2*s2); /* xy - voxel size */
    const float s13  = (float) sqrt((double)s1*s1   + s3*s3); /* xz - voxel size */
    const float s23  = (float) sqrt((double)s2*s2   + s3*s3); /* yz - voxel size */
    const float s123 = (float) sqrt((double)s12*s12 + s3*s3); /* xyz - voxel size */
    
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    const int NI[] = { 0, -1,-x+1, -x,-x-1, -xy+1,-xy,-xy-1, -xy+x+1,-xy+x,-xy+x-1, -xy-x+1,-xy-x,-xy-x-1}; 
    const float  ND[] = {0.0, s1, s12, s2, s12, s13, s3,    s13, s123, s23, s123, s123, s23, s123};
    const int sN = sizeof(NI)/4; /* division by 4 to get from the number of bytes to the number of elements */ 
    float DN[sN];
    float DNm = FLT_MAX;
    int  i, n, ni, DNi;
    int  u,v,w,nu,nv,nw; 
    
    /* data */
    float        *D = NULL;
    unsigned int *I = NULL;
    
    D = (float *)malloc(sizeof(float)*nvol);
    I = (unsigned int *)malloc(sizeof(unsigned int)*nvol);
    
    if ((D == NULL) || (I == NULL)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Initiaize mask with ones if not defined */
    if (M == NULL) {
        M = (unsigned int *)malloc(sizeof(unsigned int)*nvol);
        for (i=0; i<nvol; i++)
            M[i] = 1.0;
    }
    
    /* initialisation of D and I */
    for (i=0; i<nvol; i++) {
        if ((round(V[i])<0.5) || isnan(V[i])) D[i] = FLT_MAX; else D[i] = 0.0; 
        I[i] = (unsigned int)i;
    }
    
    /* forward direction that consider all points smaller than i */
    for (i=0; i<nvol; i++) {
        if ((D[i]>0) && (M[i]>0)) {
            ind2sub(i,&u,&v,&w,xy,x);
            
            /* read neighbor values */
            for (n=0; n<sN; n++) {
                ni = i + NI[n];
                ind2sub(ni,&nu,&nv,&nw,xy,x);
                if ((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
                DN[n] = D[ni] + ND[n];
            }
    
            /* find minimum distance within the neighborhood */
            pmin(DN,sN,&DNm,&DNi);
    
            /* update values */
            if (DNi>0) {
                I[i] = (unsigned int) I[i+NI[DNi]];
                D[i] = DNm; 
                ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
                D[i] = (float)sqrt(pow((double)(u-nu),2) + pow((double)(v-nv),2) + pow((double)(w-nw),2));
            }
         }
    }
    
    /* backward direction that consider all points larger than i */
    for (i=nvol-1;i>=0;i--) {
        if ((D[i]>0) && (M[i]>0)) {
            ind2sub(i,&u,&v,&w,xy,x);
        
            /* read neighbour values */
            for (n=0; n<sN; n++) {
                ni = i - NI[n];
                ind2sub(ni,&nu,&nv,&nw,xy,x);
                if ((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
                DN[n] = D[ni] + ND[n];
            }
        
            /* find minimum distance within the neighborhood */
            pmin(DN,sN,&DNm,&DNi);
        
            /* update values */
            if (DNi>0) {
                I[i] = (unsigned int) I[i-NI[DNi]];
                D[i] = DNm; 
                ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
                D[i] = (float)sqrt(pow((double)(u-nu),2) + pow((double)(v-nv),2) + pow((double)(w-nw),2));
            }
        }
    }

    /* finally return output to original variables V and M */
    for (i=0; i<nvol; i++) {
        if ((M[i]==0) || (D[i] == FLT_MAX))
            V[i] = 0.0; else V[i] = D[i];
        M[i] = I[i];
    }
        
    free(D);
    free(I);
}


/* laplace calculation
 * ________________________________________________________________________
 * Filter SEG within the intensity range of low and high until the changes
 * are below TH. 
 *
 * L = laplace3(SEG,TH)
 *
 * SEG  = 3d single input matrix
 * TH       = threshold to control the number of iterations
 *              maximum change of an element after iteration
 *
 */
void 
laplace3(float *SEG, int dims[3], int maxiter)
{
    /* main information about input data (size, dimensions, ...) */
    const int x = dims[0];
    const int y = dims[1];
    const int z = dims[2];
    const int xy = x*y;
    const int nvol = x*y*z;
        
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    const int NI[6] = { -1, 1, -x, x, -xy, xy}; 
    const int sN = sizeof(NI)/4;   
    
    float *L1 = (float *)malloc(sizeof(float)*nvol);
    float *L2 = (float *)malloc(sizeof(float)*nvol);
    unsigned char *LN = (unsigned char *)malloc(sizeof(unsigned char)*nvol);
    int i, n;
    
    /* initialisiation */
    for (i=0; i<nvol; i++) {
        if (isnan(SEG[i]))
            L1[i] = FLT_MAX; else L1[i] = SEG[i];
        L2[i] = L1[i];
        if (SEG[i] == 0)
            LN[i] = 1; else LN[i] = 0;
    }

    int u,v,w,nu,nv,nw,ni,iter=0;
    float Nn;
    while (iter < maxiter) {
        iter++;
        for (i=0; i<nvol; i++) {
            if ((SEG[i] == 0) && LN[i]) {
                ind2sub(i,&u,&v,&w,xy,x);

                /* read neighbor values */
                L2[i]=0.0; Nn=0.0;
                for (n=0;n<sN;n++) {
                    ni = i + NI[n];
                    ind2sub(ni,&nu,&nv,&nw,xy,x);
                    if (((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX))==0) {
                        L2[i] += L1[ni];
                        Nn++;
                    }
                }
                if (Nn>0) {L2[i]/=(float)Nn;} else {L2[i]=L1[i];}
                
                for (n=0;n<sN;n++) {
                    ni = i + NI[n];
                    ind2sub(ni,&nu,&nv,&nw,xy,x);
                    if (((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX))==0) 
                        LN[ni] = 1; /* if I change his neigbors it has to be recalculated */
                    }
                LN[i] = 0;
            }
        }
        
        /* update of L1 */
        for (i=0; i<nvol; i++) 
            L1[i] = L2[i];
        
    }
    
    for (i=0; i<nvol; i++) 
        SEG[i] = L1[i];

    free(L1);
    free(L2);
    free(LN);
}

/* laplace calculation
 * ________________________________________________________________________
 * Filter SEG within the intensity range of low and high until the changes
 * are below TH. 
 *
 * L = laplace3(SEG,TH)
 *
 * SEG  = 3d single input matrix
 * TH       = threshold to control the number of iterations
 *              maximum change of an element after iteration
 *
 */
void 
laplace3R(float *SEG, unsigned char *M, int dims[3], double TH)
{
    /* main information about input data (size, dimensions, ...) */
    const int x = dims[0];
    const int y = dims[1];
    const int z = dims[2];
    const int xy = x*y;
    const int nvol = x*y*z;
        
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    const int NI[] = { -1, 1, -x, x, -xy, xy}; 
    const int sN = sizeof(NI)/4;   

    float *L1 = (float *)malloc(sizeof(float)*nvol);
    float *L2 = (float *)malloc(sizeof(float)*nvol);
    unsigned char *LN = (unsigned char *)malloc(sizeof(unsigned char)*nvol);
    int i,n;
    
    /* initialisiation */
    for (i=0; i<nvol; i++) {
        if (isnan(SEG[i]))
            L1[i] = FLT_MAX; else L1[i] = SEG[i];
        L2[i] = L1[i];
        LN[i] = M[i];
    }

    int u,v,w,nu,nv,nw,ni,iter=0,maxiter=2000;
    float Nn, diff, maxdiffi, maxdiff=1.0;
    while (maxdiff > TH && iter < maxiter) {
        maxdiffi=0; iter++;
        for (i=0; i<nvol; i++) {
            if (M[i] && LN[i]) {  
                ind2sub(i,&u,&v,&w,xy,x);

                /* read neighbor values */
                L2[i]=0.0; Nn=0.0;
                for (n=0;n<sN;n++) {
                    ni = i + NI[n];
                    ind2sub(ni,&nu,&nv,&nw,xy,x);
                    if (((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX))==0) {
                        L2[i] += L1[ni];
                        Nn++;
                    }
                }
                if (Nn>0) L2[i] /= (float)Nn; else L2[i] = L1[i];
                
                diff    = fabs(L1[i] - L2[i]); 
                if (diff>(TH/10.0)) { 
                    for (n=0;n<sN;n++) {
                        ni = i + NI[n];
                        ind2sub(ni,&nu,&nv,&nw,xy,x);
                        if (((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX))==0) 
                            LN[ni] = 1; /* if I change his neigbors it has to be recalculated */
                    }
                }
                LN[i] = 0;
                if (maxdiffi<diff) maxdiffi = diff; 
            }
        }
        maxdiff = maxdiffi;
        
        /* update of L1 */
        for (i=0; i<nvol; i++) 
            L1[i] = L2[i];
        
    }
    
    for (i=0; i<nvol; i++) 
        SEG[i] = L1[i];

    free(L1);
    free(L2);
    free(LN);
}

void
distclose_uint8(unsigned char *vol, int dims[3], double voxelsize[3], int niter, double th)
{
    float *buffer = NULL;
    int i,x,y,z,j,band,dims2[3];
    unsigned char max_vol;
    int nvol2,nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;
    nvol2 = dims2[0]*dims2[1]*dims2[2];

    buffer = (float *)malloc(sizeof(float)*dims2[0]*dims2[1]*dims2[2]);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(float)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
        buffer[index(x+band,y+band,z+band,dims2)] = (float)((double)vol[index(x,y,z,dims)]>th);
                
    vbdist(buffer, NULL, dims2, voxelsize);
    for (i=0;i<nvol2;i++)
        buffer[i] = buffer[i] > (float)niter;

    vbdist(buffer, NULL, dims2, voxelsize);
    for (i=0;i<nvol2;i++)
        buffer[i] = buffer[i] > (float)niter;

    /* return image */
    for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
        vol[index(x,y,z,dims)] = (unsigned char)buffer[index(x+band,y+band,z+band,dims2)];
        
    free(buffer);
}

void
distclose_float(float *vol, int dims[3], double voxelsize[3], int niter, double th)
{
    float *buffer = NULL;
    int i,x,y,z,j,band,dims2[3];
    float max_vol;
    int nvol2,nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;
    nvol2 = dims2[0]*dims2[1]*dims2[2];

    buffer = (float *)malloc(sizeof(float)*dims2[0]*dims2[1]*dims2[2]);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(float)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
        buffer[index(x+band,y+band,z+band,dims2)] = (vol[index(x,y,z,dims)]>(float)th);
                
    vbdist(buffer, NULL, dims2, voxelsize);
    for (i=0;i<nvol2;i++)
        buffer[i] = buffer[i] > (float)niter;

    vbdist(buffer, NULL, dims2, voxelsize);
    for (i=0;i<nvol2;i++)
        buffer[i] = buffer[i] > (float)niter;

    /* return image */
    for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
        vol[index(x,y,z,dims)] = buffer[index(x+band,y+band,z+band,dims2)];
        
    free(buffer);
}

void
distopen_uint8(unsigned char *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    float *buffer = NULL;
    int i,j;
    unsigned char max_vol;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (dist == 0.0) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;
    
    buffer = (float *)malloc(sizeof(float)*nvol);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* threshold input */
    for (i=0; i<nvol; i++)
        buffer[i] = (float)((double)vol[i] <= th);
                
    vbdist(buffer, NULL, dims, voxelsize);
    for (i=0; i<nvol; i++)
        buffer[i] = buffer[i] > (float)dist;

    vbdist(buffer, NULL, dims, voxelsize);
    for (i=0; i<nvol; i++)
        buffer[i] = buffer[i] <= (float)dist;

    /* return image */
    for (i=0; i<nvol; i++)
        vol[i] = (unsigned char)buffer[i];

    free(buffer);
}

void
distopen_float(float *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    float *buffer = NULL;
    int i,j;
    float max_vol;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (dist == 0.0) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;
    
    buffer = (float *)malloc(sizeof(float)*nvol);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* threshold input */
    for (i=0; i<nvol; i++)
        buffer[i] = 1.0 - ((float)vol[i]>th);
                
    vbdist(buffer, NULL, dims, voxelsize);
    for (i=0; i<nvol; i++)
        buffer[i] = buffer[i] > (float)dist;

    vbdist(buffer, NULL, dims, voxelsize);
    for (i=0; i<nvol; i++)
        buffer[i] = buffer[i] <= (float)dist;

    /* return image */
    for (i=0; i<nvol; i++)
        vol[i] = buffer[i];

    free(buffer);
}

void
morph_erode_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i,j;
    unsigned char max_vol;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;
    
    /* threshold input */
    for (j=0;j<nvol;j++)
        vol[j] = (unsigned char)((double)vol[j]>th);

    for (i=0;i<niter;i++) {
        convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
        for (j=0;j<nvol;j++)
            vol[j] = (vol[j]>=9);
    }
}

void
morph_erode_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i,j;
    float max_vol;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* threshold input */
    for (j=0;j<nvol;j++)
        vol[j] = vol[j]>(float)th;

    for (i=0;i<niter;i++) {
        convxyz_float(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
        for (j=0;j<nvol;j++)
            vol[j] = vol[j]>=9.0;
    }
}

void
morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i,x,y,z,j,band,dims2[3];
    unsigned char max_vol;
    unsigned char *buffer = NULL;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);

    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = buffer[j]>0;
    }

    /* return image */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        vol[index(x,y,z,dims)] = buffer[index(x+band,y+band,z+band,dims2)];
        
    free(buffer);
    
}

void
morph_dilate_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i,x,y,z,j,band,dims2[3];
    float max_vol;
    unsigned char *buffer = NULL;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);

    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = buffer[j]>0;
    }

    /* return image */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        vol[index(x,y,z,dims)] = (float)buffer[index(x+band,y+band,z+band,dims2)];
        
    free(buffer);
    
}

void
morph_close_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    unsigned char *buffer = NULL;
    int i,x,y,z,j,band,dims2[3];
    unsigned char max_vol;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);

    /* dilate */
    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = (buffer[j]>0);
    }

    /* erode */
    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = (buffer[j]>=9);
    }

    /* return image */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        vol[index(x,y,z,dims)] = buffer[index(x+band,y+band,z+band,dims2)];
        
    free(buffer);
}

void
morph_close_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    unsigned char *buffer = NULL;
    int i,x,y,z,j,band,dims2[3];
    float max_vol;
    int nvol = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        buffer[index(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[index(x,y,z,dims)]>th);
                
    /* dilate */
    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = (buffer[j]>0);
    }

    /* erode */
    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = (buffer[j]>=9);
    }

    /* return image */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        vol[index(x,y,z,dims)] = (float)buffer[index(x+band,y+band,z+band,dims2)];
        
    free(buffer);
}

void
morph_open_uint8(unsigned char *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i, j, nvol;
    unsigned char max_vol;
    
    if (niter < 1) return;

    nvol = dims[0]*dims[1]*dims[2];
    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* threshold input */
    for (j=0;j<nvol;j++)
        vol[j] = (unsigned char)(vol[j]>th);

    for (i=0;i<niter;i++) {
        convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
        for (j=0;j<nvol;j++)
            vol[j] = (vol[j]>=9);
    }
    for (i=0;i<niter;i++) {
        convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
        for (j=0; j<nvol; j++) 
            vol[j] = (vol[j]>0);
    }

    for (j=0;j<nvol;j++)
        vol[j] = 255*vol[j];

}

void
morph_open_float(float *vol, int dims[3], int niter, double th)
{
    unsigned char *buffer = NULL;
    double filt[3]={1,1,1};
    int i, j, nvol;
    float max_vol;
    
    if (niter < 1) return;

    nvol = dims[0]*dims[1]*dims[2];
    for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*nvol);

    if (buffer == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* threshold input */
    for (j=0;j<nvol;j++)
        buffer[j] = (unsigned char)((double)vol[j]>th);

    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
        for (j=0;j<nvol;j++)
            buffer[j] = (buffer[j]>=9);
    }

    for (i=0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
        for (j=0; j<nvol; j++) 
            buffer[j] = (buffer[j]>0);
    }

    for (i=0; i<nvol; i++)
        vol[i] = (float)buffer[i];
        
    free(buffer);
}

/* First order hold resampling - trilinear interpolation */
void 
subsample_double(double *in, double *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out)
{
    int i, x, y, z;
    double k111,k112,k121,k122,k211,k212,k221,k222;
    double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
    int off1, off2, xcoord, ycoord, zcoord;

    for (i=0; i<3; i++) {
        if (dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
        else samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
    }
    
    for (z=0; z<dim_out[2]; z++) {
        zi = 1.0+(double)z/samp[2];
        for (y=0; y<dim_out[1]; y++) {
            yi = 1.0+(double)y/samp[1];
            for (x=0; x<dim_out[0]; x++) {
                xi = 1.0+(double)x/samp[0];
                i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x + offset_out;

                if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0])    {
                    xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
                    ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
                    zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

                    off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1)) + offset_in;
                    k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
                    k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
                    k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
                    k211 = (double)in[off2]; k111 = (double)in[off2+1];

                    out[i] = ((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                        + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                                 
                } else out[i] = 0;
            }
        }
    }
}

/* First order hold resampling - trilinear interpolation */
void subsample_uint8(unsigned char *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out)
{
    int i, x, y, z;
    double k111,k112,k121,k122,k211,k212,k221,k222;
    double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
    int off1, off2, xcoord, ycoord, zcoord;

    for (i=0; i<3; i++) {
        if (dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
        else samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
    }
    
    for (z=0; z<dim_out[2]; z++) {
        zi = 1.0+(double)z/samp[2];
        for (y=0; y<dim_out[1]; y++) {
            yi = 1.0+(double)y/samp[1];
            for (x=0; x<dim_out[0]; x++) {
                xi = 1.0+(double)x/samp[0];
                i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x + offset_out;

                if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0]) {
                    xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
                    ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
                    zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

                    off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1)) + offset_in;
                    k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
                    k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
                    k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
                    k211 = (double)in[off2]; k111 = (double)in[off2+1];

                    out[i] = (float)((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                        + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                                 
                } else out[i] = 0;
            }
        }
    }
}

/* First order hold resampling - trilinear interpolation */
void subsample_float(float *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out)
{
    int i, x, y, z;
    double k111,k112,k121,k122,k211,k212,k221,k222;
    double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
    int off1, off2, xcoord, ycoord, zcoord;
        
    for (i=0; i<3; i++) {
        if (dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
        else samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
    }
    
    for (z=0; z<dim_out[2]; z++) {
        zi = 1.0+(double)z/samp[2];
        for (y=0; y<dim_out[1]; y++) {
            yi = 1.0+(double)y/samp[1];
            for (x=0; x<dim_out[0]; x++) {
                xi = 1.0+(double)x/samp[0];
                i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x + offset_out;

                if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0]) {
                    xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
                    ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
                    zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

                    off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1)) + offset_in;
                    k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
                    k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
                    k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
                    k211 = (double)in[off2]; k111 = (double)in[off2+1];

                    out[i] = (float)((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                        + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                                 
                } else out[i] = 0;
            }
        }
    }
}

void
smooth_double(double *vol, int dims[3], double voxelsize[3], double s0[3], int use_mask)
{
    int i;
    double xsum, ysum, zsum;
    double *x, *y, *z, s[3];
    int xyz[3], nvol, sum_mask;
    double *mask;
    unsigned char *mask2;
    
    nvol = dims[0]*dims[1]*dims[2];

    for (i=0; i<3; i++) {
        s[i] = s0[i]/voxelsize[i];
        if(s[i] < 1.0) s[i] = 1.0;
        s[i] /= sqrt(8.0*log(2.0));
        xyz[i] = (int) RINT(6.0*s[i]);
    }
    
    x = (double *) malloc(sizeof(double)*((2*xyz[0])+1));
    y = (double *) malloc(sizeof(double)*((2*xyz[1])+1));
    z = (double *) malloc(sizeof(double)*((2*xyz[2])+1));
    
    /* build mask for masked smoothing */
    if(use_mask) {
        mask  = (double *) malloc(sizeof(double)*nvol);
        mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvol);
        sum_mask = 0;
        for (i=0; i<nvol; i++) {
            if(vol[i] == 0.0) {
                mask[i]  = 0.0;
                mask2[i] = 0;
            } else {
                mask[i]  = 1.0;
                mask2[i] = 1;
                sum_mask++;
            }
        }
    }
    
    for (i=-xyz[0]; i <= xyz[0]; i++) x[i+xyz[0]] = (double)i;
    for (i=-xyz[1]; i <= xyz[1]; i++) y[i+xyz[1]] = (double)i;
    for (i=-xyz[2]; i <= xyz[2]; i++) z[i+xyz[2]] = (double)i;
    
    xsum = 0.0; ysum = 0.0; zsum = 0.0;
    for (i=0; i < ((2*xyz[0])+1); i++) {
        x[i] = exp(-pow(x[i],2) / (2.0*pow(s[0],2)));
        xsum += x[i];
    }
    for (i=0; i < ((2*xyz[1])+1); i++) {
        y[i] = exp(-pow(y[i],2) / (2.0*pow(s[1],2)));
        ysum += y[i];
    }
    for (i=0; i < ((2*xyz[2])+1); i++) {
        z[i] = exp(-pow(z[i],2) / (2.0*pow(s[2],2)));
        zsum += z[i];
    }
    
    for (i=0; i < ((2*xyz[0])+1); i++) x[i] /= xsum;
    for (i=0; i < ((2*xyz[1])+1); i++) y[i] /= ysum;
    for (i=0; i < ((2*xyz[2])+1); i++) z[i] /= zsum;
    
    convxyz_double(vol,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],vol,dims);
    if(use_mask) {
        /* only smooth mask if mask has values > 0 */
        if(sum_mask>0) {
            convxyz_double(mask,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],mask,dims);
            for (i=0; i<nvol; i++) {
                if(mask2[i]>0) vol[i] /= (double)mask[i];    
                else vol[i] = 0.0; 
            }
        }
        free(mask);
        free(mask2);
    }
    
    free(x);
    free(y);
    free(z);

}

void
smooth_float(float *vol, int dims[3], double voxelsize[3], double s0[3], int use_mask)
{
    int i;
    double xsum, ysum, zsum;
    double *x, *y, *z, s[3];
    int xyz[3], nvol, sum_mask;
    float *mask;
    unsigned char *mask2;
    
    nvol = dims[0]*dims[1]*dims[2];

    for (i=0; i<3; i++) {
        s[i] = s0[i]/voxelsize[i];
        if(s[i] < 1.0) s[i] = 1.0;
        s[i] /= sqrt(8.0*log(2.0));
        xyz[i] = (int) RINT(6.0*s[i]);
    }
    
    x = (double *) malloc(sizeof(double)*((2*xyz[0])+1));
    y = (double *) malloc(sizeof(double)*((2*xyz[1])+1));
    z = (double *) malloc(sizeof(double)*((2*xyz[2])+1));
    
    /* build mask for masked smoothing */
    if(use_mask) {
        mask  = (float *) malloc(sizeof(float)*nvol);
        mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvol);
        sum_mask = 0;
        for (i=0; i<nvol; i++) {
            if(vol[i] == 0.0) {
                mask[i]  = 0.0;
                mask2[i] = 0;
            } else {
                mask[i]  = 1.0;
                mask2[i] = 1;
                sum_mask++;
            }
        }
    }
    
    for (i=-xyz[0]; i <= xyz[0]; i++) x[i+xyz[0]] = (double)i;
    for (i=-xyz[1]; i <= xyz[1]; i++) y[i+xyz[1]] = (double)i;
    for (i=-xyz[2]; i <= xyz[2]; i++) z[i+xyz[2]] = (double)i;
    
    xsum = 0.0; ysum = 0.0; zsum = 0.0;
    for (i=0; i < ((2*xyz[0])+1); i++) {
        x[i] = exp(-pow(x[i],2) / (2.0*pow(s[0],2)));
        xsum += x[i];
    }
    for (i=0; i < ((2*xyz[1])+1); i++) {
        y[i] = exp(-pow(y[i],2) / (2.0*pow(s[1],2)));
        ysum += y[i];
    }
    for (i=0; i < ((2*xyz[2])+1); i++) {
        z[i] = exp(-pow(z[i],2) / (2.0*pow(s[2],2)));
        zsum += z[i];
    }
    
    for (i=0; i < ((2*xyz[0])+1); i++) x[i] /= xsum;
    for (i=0; i < ((2*xyz[1])+1); i++) y[i] /= ysum;
    for (i=0; i < ((2*xyz[2])+1); i++) z[i] /= zsum;
    
    convxyz_float(vol,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],vol,dims);
    if(use_mask) {
        /* only smooth mask if mask has values > 0 */
        if(sum_mask>0) {
            convxyz_float(mask,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],mask,dims);
            for (i=0; i<nvol; i++) {
                if(mask2[i]>0) vol[i] /= (float)mask[i];    
                else vol[i] = 0.0; 
            }
        }
        free(mask);
        free(mask2);
    }
    
    free(x);
    free(y);
    free(z);

}

void
smooth_subsample_double(double *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp)
{
    int i, nvol_samp, nvol;
    int dims_samp[3];
    double *vol_samp, voxelsize_samp[3];
    
    /* define grid dimensions */
    for (i=0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
    for (i=0; i<3; i++) voxelsize_samp[i] = voxelsize[i]*((double)dims[i]/(double)dims_samp[i]);

    nvol      = dims[0]*dims[1]*dims[2];
    nvol_samp = dims_samp[0]*dims_samp[1]*dims_samp[2];
    vol_samp  = (double *)malloc(sizeof(double)*nvol_samp);

    subsample_double(vol, vol_samp, dims, dims_samp, 0, 0);     
    smooth_double(vol_samp, dims_samp, voxelsize_samp, s, use_mask);
    subsample_double(vol_samp, vol, dims_samp, dims, 0, 0);     

    free(vol_samp);
}

void
smooth_subsample_float(float *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp)
{
    int i, nvol_samp, nvol;
    int dims_samp[3];
    float *vol_samp;
    double voxelsize_samp[3];
    
    /* define grid dimensions */
    for (i=0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
    for (i=0; i<3; i++) voxelsize_samp[i] = voxelsize[i]*((double)dims[i]/(double)dims_samp[i]);

    nvol      = dims[0]*dims[1]*dims[2];
    nvol_samp = dims_samp[0]*dims_samp[1]*dims_samp[2];
    vol_samp  = (float *)malloc(sizeof(float)*nvol_samp);

    subsample_float(vol, vol_samp, dims, dims_samp, 0, 0);   
    smooth_float(vol_samp, dims_samp, voxelsize_samp, s, use_mask);
    subsample_float(vol_samp, vol, dims_samp, dims, 0, 0);   

    free(vol_samp);
}

void
vol_approx(float *vol, int dims[3], double voxelsize[3], int samp)
{
    int i, nvolr, nvol;
    int dimsr[3];
    float *volr, *buffer, *TAr;
    double voxelsizer[3];
    float min_vol = FLT_MAX, max_vol = -FLT_MAX;
    unsigned int *MIr;
    unsigned char *BMr, *BMr2;
    double threshold[2], prctile[2] = {5,95};
        
    /* define grid dimensions */
    for (i=0; i<3; i++) {
        dimsr[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
        voxelsizer[i] = voxelsize[i]*samp;
    }

    nvol   = dims[0]*dims[1]*dims[2];
    nvolr  = dimsr[0]*dimsr[1]*dimsr[2];
    volr   = (float *)malloc(sizeof(float)*nvolr);
    buffer = (float *)malloc(sizeof(float)*nvolr);
    TAr    = (float *)malloc(sizeof(float)*nvolr);
    MIr    = (unsigned int *)malloc(sizeof(unsigned int)*nvolr);
    BMr    = (unsigned char *)malloc(sizeof(unsigned char)*nvolr);
    BMr2   = (unsigned char *)malloc(sizeof(unsigned char)*nvolr);

    /* find values between 0.1% and 99.9% percentile */
    for (i = 0; i < nvol; ++i)
        vol[i] -= 0;

    for (i = 0; i < nvol; ++i) {
        min_vol = MIN(vol[i], min_vol);
        max_vol = MAX(vol[i], max_vol);
    }
    
    /* only keep values between 5..95% percentiles 
       to remove extremes that occur at edges */
    get_prctile(vol, dims, threshold, prctile, 1);  
    for (i = 0; i < nvol; ++i)
        if ((vol[i] < threshold[0]) || (vol[i] > threshold[1])) vol[i] = 0;;

    /* scale vol to 0..1 but keep zero-background */
    for (i = 0; i < nvol; ++i)
        if (vol[i] !=0) min_vol = MIN(vol[i], min_vol);
    for (i = 0; i < nvol; ++i) {
        if (vol[i] != 0) {
            vol[i] -= min_vol;   
            max_vol = MAX(vol[i], max_vol);
        }
        vol[i] /= max_vol;
    }
    
    /* downsample to lower resolution */
    subsample_float(vol, volr, dims, dimsr, 0, 0);

    /* create mask by closing holes */ 
    for (i = 0; i < nvol; ++i) BMr[i] = volr[i] > 0;
    morph_close_uint8(BMr, dimsr, 20, 0);

    /* vbdist to fill values in background with neighbours */
    memcpy(buffer,volr,nvolr*sizeof(float));    
    vbdist(buffer, MIr, dimsr, voxelsizer);
    for (i = 0; i < nvolr; ++i) {
        if (volr[i] != 0) buffer[i] = volr[i];
        else buffer[i] = volr[MIr[i]];
        TAr[i] = buffer[i];
    }
    
    /* smooth values outside mask */
    double s[3] = {20,20,20};
    smooth_float(buffer, dimsr, voxelsizer, s, 0);
    for (i = 0; i < nvolr; ++i)
        if (BMr[i] == 0) TAr[i] = buffer[i];

    /* rescue mask and create new mask that is only defined inside */
    for (i = 0; i < nvolr; ++i) {
        BMr2[i] = BMr[i];
        BMr[i] = (BMr[i] > 0) && (volr[i] == 0);
    }
    
    laplace3R(TAr, BMr, dimsr, 0.4);
    median3_float(TAr, dimsr);
    laplace3R(TAr, BMr, dimsr, 0.4);

    /* only keep TAr inside (closed) mask */
    for (i = 0; i < nvolr; ++i)
        if ((BMr2[i] == 0) && (vol[i] == 0)) TAr[i] = 0.0;

    /* again apply vbdist to fill values in background with neighbours */
    memcpy(buffer,TAr,nvolr*sizeof(float));  
    vbdist(buffer, MIr, dimsr, voxelsizer);
    for (i = 0; i < nvolr; ++i)
        buffer[i] = TAr[MIr[i]];
    for (i = 0; i < nvolr; ++i)
        TAr[i] = buffer[i];

    for (i = 0; i < 3; ++i) s[i] *= 2.0;
    smooth_float(buffer, dimsr, voxelsizer, s, 0);
    for (i = 0; i < nvolr; ++i)
        if (BMr[i] == 0) TAr[i] = buffer[i];

    for (i = 0; i < nvolr; ++i)
        BMr[i] = (BMr2[i] == 0);
    laplace3R(TAr, BMr, dimsr, 0.4);

    median3_float(TAr, dimsr);

    for (i = 0; i < nvolr; ++i)
        BMr[i] = (volr[i] == 0);
    laplace3R(TAr, BMr, dimsr, 0.4);

    memcpy(volr,TAr,nvolr*sizeof(float));
    
    subsample_float(volr, vol, dimsr, dims, 0, 0);  

    /* get old range back */
    for (i = 0; i < nvol; ++i)
        vol[i] *= max_vol;

    free(volr);
    free(buffer);
    free(MIr);
    free(BMr);
    free(BMr2);
    free(TAr);
}

void
initial_cleanup(unsigned char *probs, unsigned char *label, int dims[3], double *voxelsize, int strength, int remove_sinus)
{
    
    double scale = 3.0/(voxelsize[0] + voxelsize[1] + voxelsize[2]);
    int nvol, th, i;
    int n_initial_openings = MAX(1,round(scale*strength));
    float *sum;
    double filt[3] = {0.75, 1.0, 0.75};
    
    nvol = dims[0]*dims[1]*dims[2];

    sum = (float *)malloc(sizeof(float)*nvol);

    for (i = 0; i < nvol; ++i)
        sum[i] = (float)probs[i];

    morph_open_float(sum, dims, 1, 0.1);
    morph_dilate_float(sum, dims, 1, 0.5);

    /* build a first rough mask to remove noisy parts */
    for (i = 0; i < nvol; ++i)
        sum[i] = (float)probs[i]*sum[i] + (float)probs[i + (int)WM*nvol];

    morph_open_float(sum, dims, n_initial_openings, 0.1);
    morph_dilate_float(sum, dims, round(scale*1), 0.5);
    distclose_float(sum, dims, voxelsize, round(scale*20), 0.5);

    if(remove_sinus) {
        /* remove sinus sagittalis */
        for (i = 0; i < nvol; i++)
            sum[i] = sum[i] && (label[i] < 4);
    }

    distclose_float(sum, dims, voxelsize, round(scale*2), 0.5);

    for (i = 0; i < nvol; ++i)
        if(remove_sinus) {
            label[i] = (unsigned char)(label[i] < 4)*label[i]*sum[i];
            probs[i] = (label[i] < 4) * (sum[i] > 0) * probs[i];
        } else
            label[i] = label[i]*(unsigned char)sum[i];
    
    free(sum);
}

void
cleanup_orig(unsigned char *probs, unsigned char *mask, int dims[3], double *voxelsize, int strength)
{
    
    double scale = 3.0/(voxelsize[0] + voxelsize[1] + voxelsize[2]);
    int niter, iter, nvol, th, th_erode, th_dilate, i, j;
    int n_initial_openings = MAX(1,round(scale*strength));
    float *sum;
    double filt[3] = {0.75, 1.0, 0.75};
    
    niter     =    45;
    th_erode  = 153;  /* initial threshold for erosion 0.6*255.0 */
    th_dilate = (5*strength + 1)*16; /* threshold for dilation */
    
    nvol = dims[0]*dims[1]*dims[2];

    sum = (float *)malloc(sizeof(float)*nvol);

    /* build a first rough mask to remove noisy parts */
    for (i = 0; i < nvol; ++i)
        sum[i] = (float)probs[i] + (float)probs[i + (int)WM*nvol];

    morph_open_float(sum, dims, n_initial_openings, 0.25);
    
    /* init mask with WM values that are larger than GM and CSF and threshold for erosion */
    for (i = 0; i < nvol; ++i)
        if ((probs[i + (int)WM*nvol] > probs[i]) && (probs[i + (int)WM*nvol] > probs[i + (int)CSF*nvol]) && (probs[i + (int)WM*nvol] > th_erode) && (sum[i] > 0))
            mask[i] = probs[i + (int)WM*nvol];
        else mask[i] = 0;

    /* use masked WM image for all subsequent operations */
    for (i = 0; i < nvol; ++i) probs[i + (int)WM*nvol] = mask[i];
    
    /* mask computed from gm and wm */
    /* erosions and conditional dilations */
    for (iter=0; iter < niter; iter++) {

        /*  start with 2 iterations of erosions*/
        if(iter < 2) th = th_erode;
        else th = th_dilate;
        
        /* mask = (mask>th).*(white+gray) */
        for (i = 0; i < nvol; ++i) {
            if(mask[i] > th) {
                sum[i] = (float)probs[i] + (float)probs[i + (int)WM*nvol];
                mask[i] = (unsigned char)MIN(sum[i], 255.0);
            } else  mask[i] = 0;             
        }

        /* convolve mask with filter width of 3 voxel */
        convxyz_uint8(mask,filt,filt,filt,3,3,3,-1,-1,-1,mask,dims);
    }
    
    for (i = 0; i < nvol; ++i)
        sum[i] = (float)mask[i];

    /* use copy of mask to erode and fill holes */
    morph_erode_float(sum, dims, round(scale*4), 0.5);
    morph_close_float(sum, dims, round(scale*20), 0.5);


    /* use either original mask or new eroded and filled mask */
    for (i = 0; i < nvol; ++i)
        sum[i] = (sum[i] > 0) || (mask[i] > 0);

    /* fill remaining CSF spaces */
    distclose_float(sum, dims, voxelsize, round(scale*2), 0.5);
    distclose_uint8(mask, dims, voxelsize, round(scale*4), 0.5);

if (0) {
    for (i = 0; i < nvol; ++i) {
        for (j = 0; j < 6; ++j)
            probs[i + j*nvol] = probs[i + j*nvol] * (sum[i] > 0);
        mask[i] = probs[i] + probs[i + nvol] + probs[i + 2*nvol]; 
    }
    }
    free(sum);
    
}

void
cleanup(unsigned char *probs, unsigned char *mask, int dims[3], double *voxelsize, int strength, int gmwm_only)
{
    
    double scale = 3.0/(voxelsize[0] + voxelsize[1] + voxelsize[2]);
    int niter, iter, nvol, th, th_erode, th_dilate, i, j;
    int n_initial_openings = MAX(1,round(scale*strength));
    unsigned char *b, *c;
    float bp, tot;
    double filt[3] = {0.75, 1.0, 0.75};
    
    niter     = 32;
    th_erode  = 153;  /* initial threshold for erosion 0.6*255.0 */
    th_dilate = 38 + (strength*13);     /* threshold for dilation (0.15 + strength*0.05)*255 */
    
    /* ensure that sum of filt is 1 */
    tot = 0.0;
    for (i = 0; i < 3; ++i)
        tot += filt[i];
    for (i = 0; i < 3; ++i)
        filt[i] /= tot;

    nvol = dims[0]*dims[1]*dims[2];

    b = (unsigned char *)malloc(sizeof(unsigned char)*nvol);
    c = (unsigned char *)malloc(sizeof(unsigned char)*nvol);

    /* init with WM */
    for (i = 0; i < nvol; ++i)
        b[i] = probs[i + (int)WM*nvol];
    
    /* mask computed from gm and wm */
    /* erosions and conditional dilations */
    for (iter=0; iter < niter; iter++) {
        
        fprintf(stderr,"."); 
        /*  start with 2 iterations of erosions */
        if(iter < 2) th = th_erode;
        else th = th_dilate;
        
        /* b = (b>th).*(white+gray) */
        for (i = 0; i < nvol; ++i) {
            if(b[i] > th) {
                bp = (float)probs[i] + (float)probs[i + (int)WM*nvol];
                b[i] = (unsigned char)MIN(round(bp), 255.0);
            } else b[i] = 0;               
        }

        /* convolve mask with filter width of 3 voxel */
        convxyz_uint8(b,filt,filt,filt,3,3,3,-1,-1,-1,b,dims);           
    }
    fprintf(stderr,"\n");
    
    morph_open_uint8(b, dims, n_initial_openings, 0.05);

    if (gmwm_only == 0) {
        for (i = 0; i < nvol; ++i)
            c[i] = b[i];
            
        /* use copy of mask to fill holes */
        morph_close_uint8(c, dims, round(scale*20), 0.5);
    
        /* use either original mask or new eroded and filled mask */
        for (i = 0; i < nvol; ++i)
            c[i] = (c[i] > 0) ||    (b[i] > 0);
    
        /* fill remaining CSF spaces */
        distclose_uint8(c, dims, voxelsize, round(scale*4), 0.5);
    }
    
    th = 13; /* 0.05*255 */
    for (i = 0; i < nvol; ++i) {
        if((b[i] < th) || (((float)probs[i] + (float)probs[i + (int)WM*nvol]) < th)) {
            probs[i] = 0;
            probs[i + (int)WM*nvol] = 0;
        }        
    }
    
    if (gmwm_only == 0) {
        for (i = 0; i < nvol; ++i) {
            if((c[i] < th) || (((float)probs[i] + (float)probs[i + (int)WM*nvol] + (float)probs[i + (int)CSF*nvol]) < th)) {
                probs[i + (int)CSF*nvol] = 0;
            }        
        }
    }
    
    for (i = 0; i < nvol; ++i) {
        tot = 0.0;
        for (j = 0; j < 3; ++j)
            tot += (float)probs[i + j*nvol];
        for (j = 0; j < 3; ++j)
            probs[i + j*nvol] = (unsigned char)(round((float)probs[i + j*nvol]/tot*255.0));
        mask[i] = probs[i] + probs[i + (int)WM*nvol];
    }

    free(b);
    free(c);
    
}

/* qicksort */
void
swap_float(float *a, float *b)
{
    float t=*a;
    *a=*b;
    *b=t;
}

void
sort_float(float arr[], int start, int end)
{
    if (end > start + 1)
    {
        float piv = arr[start];
        int l = start + 1, r = end;
        while (l < r)
        {
            if (arr[l] <= piv) l++;
            else swap_float(&arr[l], &arr[--r]);
        }
        swap_float(&arr[--l], &arr[start]);
        sort_float(arr, start, l);
        sort_float(arr, r, end);
    }
}

/* simple median function for uint8 */
void 
median3_uint8(unsigned char *D, int dims[3])
{
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    float NV[27];
    int i,j,k,ind,ni,x,y,z,n;
    unsigned char *M;
    int nvol = dims[0]*dims[1]*dims[2];
                
    /* output */
    M = (unsigned char *)malloc(sizeof(unsigned char)*nvol);

    /* filter process */
    for (z=0; z<dims[2]; z++) for (y=0; y<dims[1]; y++) for (x=0; x<dims[0]; x++) {
        ind = index(x,y,z,dims);
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (i=-1; i<=1; i++) for (j=-1; j<=1; j++) for (k=-1; k<=1; k++) {
            /* check borders */ 
            if (((x+i)>=0) && ((x+i)<dims[0]) && ((y+j)>=0) && ((y+j)<dims[1]) && ((z+k)>=0) && ((z+k)<dims[2])) {
                ni = index(x+i,y+j,z+k,dims);
                /* check masks and NaN or Infinities */
                if (isnan(D[ni]) || D[ni]==FLT_MAX || D[ind]==-FLT_MAX) ni = ind;
                NV[n] = (float)D[ni];
                n++;
            }
        }
        /* get correct n */
        n--;
        /* sort and get the median by finding the element in the middle of the sorting */
        sort_float(NV,0,n);
        M[ind] = (unsigned char)NV[(int)(n/2)];
    }
     
    for (i=0; i<nvol; i++) D[i] = M[i];
    
    free(M);
    
}

/* simple median function for unsigned short */
void 
median3_short(unsigned short *D, int dims[3])
{
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    float NV[27];
    int i,j,k,ind,ni,x,y,z,n;
    unsigned short *M;
    int nvol = dims[0]*dims[1]*dims[2];
                
    /* output */
    M = (unsigned short *)malloc(sizeof(unsigned short)*nvol);

    /* filter process */
    for (z=0; z<dims[2]; z++) for (y=0; y<dims[1]; y++) for (x=0; x<dims[0]; x++) {
        ind = index(x,y,z,dims);
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (i=-1; i<=1; i++) for (j=-1; j<=1; j++) for (k=-1; k<=1; k++) {
            /* check borders */ 
            if (((x+i)>=0) && ((x+i)<dims[0]) && ((y+j)>=0) && ((y+j)<dims[1]) && ((z+k)>=0) && ((z+k)<dims[2])) {
                ni = index(x+i,y+j,z+k,dims);
                /* check masks and NaN or Infinities */
                if (isnan(D[ni]) || D[ni]==FLT_MAX || D[ind]==-FLT_MAX) ni = ind;
                NV[n] = (float)D[ni];
                n++;
            }
        }
        /* get correct n */
        n--;
        /* sort and get the median by finding the element in the middle of the sorting */
        sort_float(NV,0,n);
        M[ind] = (unsigned short)NV[(int)(n/2)];
    }
     
    for (i=0; i<nvol; i++) D[i] = M[i];
    
    free(M);
    
}

/* simple median function for float */
void 
median3_float(float *D, int dims[3])
{
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    float NV[27];
    int i,j,k,ind,ni,x,y,z,n;
    float *M;
    int nvol = dims[0]*dims[1]*dims[2];
                
    /* output */
    M = (float *)malloc(sizeof(float)*nvol);

    /* filter process */
    for (z=0; z<dims[2]; z++) for (y=0; y<dims[1]; y++) for (x=0; x<dims[0]; x++) {
        ind = index(x,y,z,dims);
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (i=-1; i<=1; i++) for (j=-1; j<=1; j++) for (k=-1; k<=1; k++) {
            /* check borders */ 
            if (((x+i)>=0) && ((x+i)<dims[0]) && ((y+j)>=0) && ((y+j)<dims[1]) && ((z+k)>=0) && ((z+k)<dims[2])) {
                ni = index(x+i,y+j,z+k,dims);
                /* check masks and NaN or Infinities */
                if (isnan(D[ni]) || D[ni]==FLT_MAX || D[ind]==-FLT_MAX) ni = ind;
                NV[n] = D[ni];
                n++;
            }
        }
        /* get correct n */
        n--;
        /* sort and get the median by finding the element in the middle of the sorting */
        sort_float(NV,0,n);
        M[ind] = NV[(int)(n/2)];
    }
     
    for (i=0; i<nvol; i++) D[i] = M[i];
    
    free(M);
    
}
