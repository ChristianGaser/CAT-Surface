# spm_shoot: SPM26 geodesic-shooting primitives

The C code behind SPM's `spm_diffeo`, `spm_bsplinc` and `spm_bsplins` mex
files, as CAT12's geodesic shooting (`cat_main_registration.m` →
`run_Shooting` → `spm_shoot_update`, `spm_shoot3d`) uses them, without
MATLAB.

## Origin and licence

- SPM26, version 26.10.01 (13-Oct-2026), `spm12/src`.
- John Ashburner, Copyright (C) 2001-2022 Wellcome Centre for Human
  Neuroimaging; `bsplines.c` is based on code by Philippe Thévenaz.
- GNU General Public License, version 2 or later (as SPM). See
  `../README.md` for what that means for the CAT-Surface binaries.

Every ported file keeps its copyright header, followed by one
`Modified for CAT-Surface:` line.

## Files

From SPM (`spm12/src`), modified as listed below:

| File | Content |
|------|---------|
| `shoot_boundary.c/.h` | boundary condition (circulant, Neumann) |
| `shoot_regularisers.c/.h` | `vel2mom`, `kernel`, `relax`, `Atimesp`, `trapprox` |
| `shoot_optim3d.c/.h` | `fmg3`, `fmg3_scratchsize`, `cgs3`, `norm` |
| `shoot_multiscale.c/.h` | `restrict_vol`, `resize_vol` |
| `shoot_diffeo3d.c/.h` | `pull(c)`, `push(c)`, `pushc_grads`, `composition*`, `smalldef*`, `determinant`, `def2det`, `def2jac`, `grad`, `divergence`, `minmax_div`, `bracket` |
| `shoot_expm3.c/.h` | `expm22`, `expm33` |
| `shoot_invdef.c` | `invdef` |
| `bsplines.c/.h` | B-spline deconvolution and sampling, double precision |

New:

| File | Content |
|------|---------|
| `spm_shoot.h` | **the public header**: `spm_` prefixed prototypes and conventions |
| `shoot_types.h` | replaces `spm_mex.h`/`mex.h`: index types, `mxIsFinite`, `mxGetNaN`, `mexErrMsgTxt`; floating-point contraction off |
| `shoot_namespace.h` | `#define`s every exported SPM name to `spm_<name>` |
| `shoot_parallel.h/.c` | thread control, parallel for, sweep teams (pthreads via libCAT's `cat_parallel_run`) |
| `shoot_invdef.h` | `invdef` prototype |
| `spm_shoot.c` | error reporting, `spm_diffeo` gateway equivalents (`spm_shoot_param`, `spm_shoot_fmg_param`, `spm_shoot_fmg`, `spm_shoot_mat`) |
| `spm_bsplines3.c` | `spm_bsplinc3`, `spm_bsplins3`: `spm_bsplinc.c`/`spm_bsplins.c` for volumes in memory |

Not ported: the mex gateways (`spm_diffeo.c`, `spm_bsplinc.c`,
`spm_bsplins.c`, `shoot_bsplines.c`), `shoot_dartel.c`, `shoot_optimN.c`
and the single-precision B-splines of `spm_diffeo('bsplinc'/'bsplins')`.

## Modifications

1. **No mex.** `#include "spm_mex.h"` became `#include "shoot_types.h"`.
   `mexErrMsgTxt` records the first message (`spm_shoot_get_error()`) and
   returns instead of unwinding the call. Where returning would have been
   unsafe the caller was changed: `scan_line()` (`shoot_invdef.c`) stops at
   `MAXV` voxels instead of writing past its array (SPM's non-mex
   `spm_mex.h` only prints the message and overruns).
2. **Namespace.** `shoot_namespace.h`, included first by every ported file,
   renames all exported symbols (`bound`, `norm`, `grad`, `kernel`, `wrap`,
   `push`, `pull`, `relax`, ...) to `spm_<name>`. The sources keep SPM's
   names, so they stay diffable against upstream.
3. **pthreads instead of OpenMP.** All `#pragma omp` lines are gone and
   `spm_openmp.h` is not used. The `#ifdef _OPENMP` blocks are left in place;
   without OpenMP they are inactive, and they only compute the same offsets
   elsewhere. Threads come from libCAT's `cat_parallel_run()` (prototype in
   `shoot_parallel.h`, which does not include `CAT_Vol.h`: its macros `GM`,
   `WM`, `CSF` and `index` clash). See *Threads* below.
4. **Range functions.** To split loops among threads, `trapprox`, `sumsq`,
   `vel2mom`, `Atimesp1`, `composition_stuff`, `def2jac_fun`, `pushpull`,
   `pushc_grads`, `smalldef*`, `determinant`, `restrict_vol`, `resize_vol`
   and `invert_it` got loop bounds as parameters; the exported functions
   call them over the whole range in serial mode, which is the SPM loop.
5. **In-memory B-splines.** `spm_bsplines3.c` holds `vol_coeffs` of
   `spm_bsplinc.c` and `fun`/`dfun` of `spm_bsplins.c` for a float volume
   and double coordinates: the copy reproduces `resample()` (`1.0*v + 0.0`,
   non-finite → 0), then deconvolution along x, y, z.
6. **Undefined behaviour removed** (results unchanged where defined):
   - no pointer arithmetic on null pointers (`Atimesp1` without Hessian,
     `composition_stuff` without Jacobians, `pushc_grads` without
     Jacobians);
   - `expm22`/`expm33` of a zero matrix (`smalldef_jac1` of a zero velocity,
     e.g. the first shooting iteration) cast `log(0) = -Inf` to `int`; the
     mex gets K ≤ 0 there, the port sets K = 0, i.e. the same Padé
     approximant without squaring;
   - `spm_bsplins3` returns 0 (derivatives 0) for non-finite or |x| > 2^30
     coordinates in periodic dimensions, where the mex casts them to `int`
     (on arm64 its weights then vanish and it returns 0 as well).
7. **Floating-point contraction off.** The SPM26 mex files contain no fused
   multiply-add instruction (`otool -tv spm_diffeo.mexmaca64`: 2563 `fmul`,
   2055 `fadd`, no `fmadd`), so `shoot_types.h` sets
   `#pragma STDC FP_CONTRACT OFF` for clang. **GCC ignores the pragma and
   fuses by default in its GNU dialects** (346 `fmadd` in
   `shoot_regularisers.o` with gcc-14 on arm64): add `-ffp-contract=off`
   for these files when building with GCC on FMA hardware (aarch64, or x86
   with `-mfma`/`-march=native`). `-DSPM_SHOOT_ALLOW_FP_CONTRACT` disables
   the pragma. With contraction the results differ from the mex at rounding
   level (relative ~1e-7 in float, ~1e-16 in double; `invdef` up to 0.008
   voxels, because float rounding moves tetrahedron borders).
8. `invert_it`'s 384 KB voxel scratch moved from its stack frame to the heap
   (one per thread); `spm_shoot_error` is protected by a mutex.

## Threads

`spm_shoot_set_threads(0)`, the default, runs SPM's serial code, which is bit
for bit what the macOS mex computes (it was built without OpenMP). With
`n >= 1`:

- **Per-voxel loops** (`vel2mom`, `Atimesp`, `pull(c)`, `composition*`,
  `smalldef*`, `determinant`, `def2det/jac`, B-splines) are split into
  blocks of slices or voxels.
- **Relaxation** (`relax_me/be/le/all`, the core of `fmg3`): SPM's sweeps are
  multi-colour Gauss-Seidel (2, 8 or 27 colours, i.e. 2 or 3 per dimension).
  SPM's OpenMP branches use the same colour order as its serial code and only
  run the `(k,j)` rows of a colour concurrently. That is a race whenever the
  second or third dimension is not a multiple of 2 or 3, as with
  113x137x113: across the periodic boundary the last slice (or row) of a
  colour reads the first one. Here a team of
  threads claims the slices of each colour, the thread that takes the first
  slice also does the last one after it, and the team meets at the end of
  every sweep. Every slice thus reads exactly what it reads in the serial
  loop, so the result is the serial one. The team claims dynamically, so it
  also finishes when fewer threads run than requested.
- **Scatter** (`pushc_grads`, `invdef`): each thread owns a slab of output
  slices and walks all inputs in the serial order, writing only into its
  slab. Every output voxel receives the same additions (or assignments and
  averages, `invdef`) in the same order as in the serial loop.
- **Multigrid transfer**: each thread restricts or resizes a block of output
  slices with its own plane buffer.
- `push`/`pushc` (not used by shooting) stay serial.

Every function therefore gives the serial result for any number of threads.
The exceptions are `trapprox` and `sumsq`, whose sums are formed per slice
and added in slice order with threads: identical for every `n >= 1`, but not
bit-identical to `n = 0` (relative difference ≤ 6e-14).
`spm_shoot_set_grain()` sets the minimum work per thread; results do not
depend on it. Without pthreads (Windows, where `cat_parallel_run()` runs its
workers one after another) everything runs in the calling thread, with the
same results.

## Validation

`tests/spm_shoot_reference/`:

- `make_reference.m` runs the SPM26 mex files on seeded smooth random inputs
  (31x27x23 at 1.5x1.2x1.8 mm and 113x137x113 at 1.5 mm, CAT12's
  regularisation with schedules 32, 4, 1, all four relaxation variants,
  NaN/Inf inputs, B-splines `[2 2 2 1 1 1]`, `[3 3 3 0 0 0]`,
  `[4 1 2 1 0 1]`) and writes inputs and outputs as raw arrays.
- `compare_reference.c` (built by `build.sh`, no autotools) compares every
  primitive with them for several thread settings; `--time N` times them.

With contraction off, all 116 outputs (71 small, 45 real-size) are
bit-identical to the mex with 0, 1 and 8 threads, and with 8 threads forced
onto even the coarsest multigrid levels. The only differences are the
threaded `trapprox` sums above. The same holds for gcc-14 with
`-ffp-contract=off`, and repeated runs with 2-24 threads never differed.

`tests/test_spm_shoot.c` (minunit, no MATLAB) checks the properties that
hold without a reference: identity pull, adjointness of push and pull, zero
velocity → identity, composition with the identity, `fmg3` inverting
`vel2mom`, the kernel as the impulse response of `vel2mom`, B-spline
interpolation, bit-identical results for 1 and 4 threads, and the error
paths. It passes under `-fsanitize=undefined`.

## Upstream issues found

- `spm_diffeo('invdef', y)` with one argument uses two uninitialised 4x4
  matrices (the gateway sets them only when more arguments are given).
- `def2det`/`def2jac` with a slice number write that slice's full-volume
  offset into a one-slice output (out of bounds for slices > 1).
- `relax_*` read `dm[3]` (`dm[0]==1 && dm[2]==1 && dm[3]==1`, probably meant
  `dm[1]`).
- SPM's OpenMP build (off by default) races in `relax_*` across the periodic
  boundary, shares the loop index `i` among threads in `trapprox`, and makes
  `push`/`pushc` order-dependent through atomic float additions.
