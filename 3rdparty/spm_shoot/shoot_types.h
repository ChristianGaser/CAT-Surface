/*
 * Replacement for SPM's "spm_mex.h" (and MATLAB's "mex.h") in the
 * CAT-Surface port of the SPM26 geodesic-shooting primitives.
 *
 * Only the ported .c files include this header.  It provides the MATLAB
 * index types and the few mx and mex calls the computational code uses:
 *  - mwSize/mwIndex/mwSignedIndex as with "mex -largeArrayDims";
 *  - mxIsFinite/mxGetNaN from <math.h>;
 *  - mexErrMsgTxt records the message (spm_shoot_last_error()) and returns
 *    instead of unwinding out of the mex call; the callers were adapted
 *    where continuing would be unsafe.
 *
 * It also switches floating-point contraction off for these files.  The
 * SPM26 mex files (spm_diffeo, spm_bsplinc, spm_bsplins) contain no fused
 * multiply-add instructions, so contraction has to stay off for the serial
 * code path to reproduce them bit for bit.  GCC ignores the pragma; build
 * with -ffp-contract=off there (GCC defaults to "fast" in GNU C modes).
 * Define SPM_SHOOT_ALLOW_FP_CONTRACT to leave the compiler's setting alone.
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later, as the SPM code it accompanies.
 */

#ifndef SHOOT_TYPES_H
#define SHOOT_TYPES_H

#include <stddef.h>
#include <math.h>

#if defined(__clang__) && !defined(SPM_SHOOT_ALLOW_FP_CONTRACT)
#pragma STDC FP_CONTRACT OFF
#endif

#ifndef SPM_SHOOT_MWTYPES
#define SPM_SHOOT_MWTYPES
typedef size_t    mwSize;
typedef size_t    mwIndex;
typedef ptrdiff_t mwSignedIndex;
#endif

#define mxIsFinite(x)     (isfinite((double)(x)) != 0)
#define mxGetNaN()        ((double)NAN)
#define mexErrMsgTxt(msg) spm_shoot_error(msg)

extern void spm_shoot_error(const char *msg);

#endif /* SHOOT_TYPES_H */
