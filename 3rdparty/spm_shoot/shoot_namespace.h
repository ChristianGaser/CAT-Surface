/*
 * Symbol namespace for the CAT-Surface port of the SPM26 shooting primitives.
 *
 * Included first by every ported .c and .h file: it renames each symbol
 * with external linkage (and the two that SPM's headers declare extern but
 * its sources define static: sumsq, pushpull) to an spm_ prefixed name, so
 * generic names such as bound, norm, grad, kernel or wrap cannot clash with
 * libCAT, bicpl, nifti, s2kit or FFTW.  The sources keep SPM's names and stay
 * diffable against upstream.
 *
 * Code outside 3rdparty/spm_shoot should include "spm_shoot.h" instead,
 * which declares the prefixed names directly and defines none of these
 * macros.
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later, as the SPM code it accompanies.
 */

#ifndef SHOOT_NAMESPACE_H
#define SHOOT_NAMESPACE_H

/* shoot_boundary.c */
#define bound                   spm_bound
#define set_bound               spm_set_bound
#define get_bound               spm_get_bound

/* shoot_regularisers.c */
#define trapprox                spm_trapprox
#define vel2mom                 spm_vel2mom
#define relax                   spm_relax
#define Atimesp                 spm_Atimesp
#define sumsq                   spm_sumsq
#define kernel                  spm_kernel

/* shoot_optim3d.c */
#define fmg3                    spm_fmg3
#define fmg3_scratchsize        spm_fmg3_scratchsize
#define cgs3                    spm_cgs3
#define norm                    spm_norm

/* shoot_multiscale.c */
#define resize_vol              spm_resize_vol
#define restrict_vol            spm_restrict_vol

/* shoot_diffeo3d.c */
#define grad                    spm_grad
#define bracket                 spm_bracket
#define composition             spm_composition
#define composition_jacobian    spm_composition_jacobian
#define composition_jacdet      spm_composition_jacdet
#define def2det                 spm_def2det
#define def2jac                 spm_def2jac
#define pushpull                spm_pushpull
#define pullc                   spm_pullc
#define pull                    spm_pull
#define pushc                   spm_pushc
#define push                    spm_push
#define pushc_grads             spm_pushc_grads
#define smalldef                spm_smalldef
#define smalldef_jac            spm_smalldef_jac
#define smalldef_jac1           spm_smalldef_jac1
#define divergence              spm_divergence
#define minmax_div              spm_minmax_div
#define determinant             spm_determinant

/* shoot_expm3.c */
#define expm22                  spm_expm22
#define expm33                  spm_expm33

/* shoot_invdef.c */
#define invdef                  spm_invdef

/* bsplines.c */
#define splinc_wrap             spm_splinc_wrap
#define splinc_mirror           spm_splinc_mirror
#define get_poles               spm_get_poles
#define mirror                  spm_mirror
#define wrap                    spm_wrap
#define sample3                 spm_sample3
#define dsample3                spm_dsample3

#endif /* SHOOT_NAMESPACE_H */
