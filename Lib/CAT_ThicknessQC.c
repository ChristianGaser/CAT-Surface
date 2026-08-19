/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "CAT_ThicknessQC.h"
#include "CAT_Vol.h"

/**
 * \brief Fill a CAT_ThicknessQCOpts with defaults.
 *
 * The threshold is the point above which cortex stops being defensible: human
 * cortex reaches about 4.5 mm at the thickest. The plate radius follows from
 * the geometry of the fault being looked for -- two cortical banks back to back
 * form a band of at most about 5 mm, so nothing thicker than a 2.5 mm inscribed
 * radius can be a glued sulcus.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_ThicknessQCOptionsInit(CAT_ThicknessQCOpts *opts)
{
    if (!opts)
        return;
    opts->thresh = 4.5;
    opts->plate_radius = 2.5;
    opts->min_volume = 20.0;
    opts->conn = 26;
    opts->verbose = 0;
}

/* Offsets of the 26 neighbours, ordered so that the first 6 are face
   neighbours and the first 18 are face plus edge. That lets one table serve
   all three connectivities by taking a prefix of it. */
static const int NB[26][3] = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
    {-1, -1, 0}, {-1, 1, 0}, {1, -1, 0}, {1, 1, 0},
    {-1, 0, -1}, {-1, 0, 1}, {1, 0, -1}, {1, 0, 1},
    {0, -1, -1}, {0, -1, 1}, {0, 1, -1}, {0, 1, 1},
    {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},
    {1, -1, -1}, {1, -1, 1}, {1, 1, -1}, {1, 1, 1}};

/**
 * \brief Find implausibly thick components and classify them by shape.
 *
 * The inscribed radius is read off a Euclidean distance transform seeded on the
 * *complement* of the thick set, so that at a thick voxel it gives the distance
 * to the nearest voxel outside it -- which is exactly the radius of the largest
 * sphere centred there that still fits inside the component. The maximum of
 * that over a component is the quantity that distinguishes a band from a mass,
 * and unlike the component's volume or its bounding box it does not grow with
 * extent along the sulcus.
 *
 * Components are grown iteratively rather than by recursion: one can hold tens
 * of thousands of voxels, far past a comfortable call depth.
 *
 * \param gmt       (in)  thickness map in mm
 * \param label     (in)  PVE label map in [0..3], or NULL for no mask
 * \param dims      (in)  {nx, ny, nz}
 * \param voxelsize (in)  voxel spacing in mm
 * \param opts      (in)  parameters; NULL selects the defaults
 * \param comps     (out) newly allocated array of components, or NULL to skip
 * \param n_comps   (out) number of components returned, or NULL to skip
 * \param classmap  (out) float[nvox] class tags, or NULL to skip
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolThicknessQC(const float *gmt, const float *label, int dims[3],
                       double voxelsize[3], const CAT_ThicknessQCOpts *opts,
                       CAT_ThicknessComponent **comps, int *n_comps,
                       float *classmap)
{
    CAT_ThicknessQCOpts defaults;
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    unsigned char *thick = NULL, *seen = NULL;
    float *edt = NULL;
    int *stack = NULL;
    CAT_ThicknessComponent *out = NULL;
    int n_out = 0, cap = 0;
    double vox_mm3;
    int n_nb;
    int i, rc = 0;

    if (!gmt || !dims || !voxelsize || nvox <= 0)
        return -1;
    if (nx < 1 || ny < 1 || nz < 1)
        return -1;

    if (!opts)
    {
        CAT_ThicknessQCOptionsInit(&defaults);
        opts = &defaults;
    }

    n_nb = (opts->conn == 6) ? 6 : ((opts->conn == 18) ? 18 : 26);
    vox_mm3 = fabs(voxelsize[0]) * fabs(voxelsize[1]) * fabs(voxelsize[2]);

    thick = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    seen = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    edt = (float *)malloc(sizeof(float) * (size_t)nvox);
    stack = (int *)malloc(sizeof(int) * (size_t)nvox);
    if (!thick || !seen || !edt || !stack)
    {
        rc = -2;
        goto cleanup;
    }

    /* the implausibly thick set, optionally confined to the cortical band */
    for (i = 0; i < nvox; i++)
    {
        int ok = ((double)gmt[i] > opts->thresh);
        if (ok && label)
            ok = (label[i] > CGM && label[i] < GWM);
        thick[i] = ok ? 1 : 0;
        seen[i] = 0;
    }

    /* Distance to the nearest voxel *outside* the thick set. Seeding on the
       complement is what turns euclidean_distance() into an inscribed radius:
       at a thick voxel it then reports how far the component extends around
       it in every direction. */
    for (i = 0; i < nvox; i++)
        edt[i] = thick[i] ? 0.0f : 1.0f;
    euclidean_distance(edt, NULL, dims, voxelsize, 0);

    if (opts->verbose)
    {
        long n_thick = 0;
        for (i = 0; i < nvox; i++)
            n_thick += thick[i];
        fprintf(stderr, "Thickness QC: %ld voxels above %.2f mm.\n",
                n_thick, opts->thresh);
    }

    if (classmap)
        for (i = 0; i < nvox; i++)
            classmap[i] = 0.0f;

    for (i = 0; i < nvox; i++)
    {
        int top = 0, head = 0;
        long n_vox = 0;
        double sum_gmt = 0.0, max_gmt = 0.0, max_r = 0.0;
        double cx = 0.0, cy = 0.0, cz = 0.0;
        int shape, j;

        if (!thick[i] || seen[i])
            continue;

        /* Grow the component as a queue rather than a stack, so that when the
           fill is done stack[0 .. top-1] still holds every voxel of it. The
           class is only known once the whole component has been measured, and
           keeping the list is what lets it be tagged without a second sweep. */
        stack[top++] = i;
        seen[i] = 1;
        while (head < top)
        {
            const int idx = stack[head++];
            const int z = idx / xy;
            const int y = (idx - z * xy) / nx;
            const int x = idx - z * xy - y * nx;
            int k;

            n_vox++;
            sum_gmt += (double)gmt[idx];
            if ((double)gmt[idx] > max_gmt)
                max_gmt = (double)gmt[idx];
            if ((double)edt[idx] > max_r)
                max_r = (double)edt[idx];
            cx += x;
            cy += y;
            cz += z;

            for (k = 0; k < n_nb; k++)
            {
                const int px = x + NB[k][0];
                const int py = y + NB[k][1];
                const int pz = z + NB[k][2];
                int nidx;

                if (px < 0 || px >= nx || py < 0 || py >= ny || pz < 0 || pz >= nz)
                    continue;
                nidx = px + py * nx + pz * xy;
                if (!thick[nidx] || seen[nidx])
                    continue;
                seen[nidx] = 1;
                stack[top++] = nidx;
            }
        }

        if ((double)n_vox * vox_mm3 < opts->min_volume)
            continue;

        shape = (max_r <= opts->plate_radius) ? CAT_QC_PLATE : CAT_QC_SOLID;

        if (classmap)
            for (j = 0; j < top; j++)
                classmap[stack[j]] = (float)shape;

        if (comps && n_comps)
        {
            if (n_out == cap)
            {
                CAT_ThicknessComponent *grown;
                int newcap = cap ? cap * 2 : 64;
                grown = (CAT_ThicknessComponent *)realloc(
                    out, sizeof(CAT_ThicknessComponent) * (size_t)newcap);
                if (!grown)
                {
                    rc = -2;
                    goto cleanup;
                }
                out = grown;
                cap = newcap;
            }
            out[n_out].n_voxels = n_vox;
            out[n_out].volume_mm3 = (double)n_vox * vox_mm3;
            out[n_out].max_radius = max_r;
            out[n_out].centroid[0] = cx / (double)n_vox;
            out[n_out].centroid[1] = cy / (double)n_vox;
            out[n_out].centroid[2] = cz / (double)n_vox;
            out[n_out].gmt_mean = sum_gmt / (double)n_vox;
            out[n_out].gmt_max = max_gmt;
            out[n_out].shape = shape;
        }
        n_out++;
    }

    if (comps && n_comps)
    {
        *comps = out;
        *n_comps = n_out;
        out = NULL;
    }
    else if (n_comps)
        *n_comps = n_out;

cleanup:
    free(thick);
    free(seen);
    free(edt);
    free(stack);
    free(out);
    return rc;
}
