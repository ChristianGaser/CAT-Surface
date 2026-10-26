/*
 * John Ashburner
 * Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging
 * Modified for CAT-Surface: new header; SPM declares invdef() in shoot_diffeo3d.h (see README.md).
 */

#include "shoot_namespace.h"

/* Spelled as in shoot_diffeo3d.h, so that both can be included together
   (shoot_invdef.c, which spells the matrices M1[][3], does not include it) */
extern void invdef(mwSize dim_y[3], float y[], mwSize dim_iy[3], float iy[], float M1[4][3], float M2[4][3]);
