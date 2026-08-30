# CLAUDE.md — Claude Code guide for CAT-Surface

> Full agent context is in [`Agents.md`](Agents.md). This file is the concise
> quick-reference loaded automatically by Claude Code.

## Sub-Agent Routing Rules
- **Always sequential:** All tasks (security, performance, style, refactoring) must be processed sequentially.
- **No parallelization:** Only one sub-agent or one check may be active at a time.
- **Workflow:** First execute `security`, then `performance`, then `style`. Wait for each step to complete.
- **Dependencies:** B tasks must wait for the output of A tasks.

## Background Execution Rules
 
Run in background automatically:
 
- Web research and documentation lookups
- Codebase exploration and analysis
- Security audits and performance profiling
- Any task where results aren't immediately needed
- Research or analysis tasks (not file modifications)
- Results aren't blocking your current work

## Project at a glance

C library + 70+ CLI tools for cortical surface mesh processing (neuroimaging).
Backend for the [CAT12](https://github.com/ChristianGaser/cat12) SPM toolbox.

- **Language:** C99
- **Build system:** GNU Autotools (Autoconf / Automake / Libtool)
- **Only external dependency:** FFTW3

## Build & test

```bash
# First time / after changing configure.ac or Makefile.am
./autogen.sh
./configure        # add --prefix=... if desired
make -j$(nproc)

# Run unit tests
make check

# Build the Python bindings against the tree and smoke-test them.
# Needed because cat_surf/_*.c are untracked build artifacts, so only an
# actual build proves the .pyx sources still cythonize, compile and link.
CAT_SURFACE_ROOT=$PWD CAT_BUILD_DIR=$PWD pip install ./cat_surface_cython
python cat_surface_cython/tests/smoke_test.py

# Spell-check before every PR
codespell --config .codespellrc
```

CI runs the build, `make check` and the binding smoke test on Ubuntu 22.04 via
`.github/workflows/ci.yml`. It does **not** run codespell — that stays a local
pre-PR step.

## Directory layout

```
Include/    Public headers   CAT_<Module>.h
Lib/        Library source   CAT_<Module>.c  → compiled into libCAT
Progs/      CLI entry points CAT_<Tool>.c    → thin wrappers over libCAT
3rdparty/   Vendored deps    DO NOT modify unless unavoidable
tests/      Unit tests       minunit framework; run with `make check`
docs/       Generated Doxygen output
cat_surface_cython/  Cython bindings (`cat-surf` package) exposing libCAT to Python
```

## Python bindings (`cat-surf`)

`cat_surface_cython/` builds the `cat-surf` PyPI package. Keep three layers consistent:

- **Low-level array API** (`cat_surf.<fn>`): numpy in/out, one function per algorithm.
- **File-based CLI mirror** (`cat_surf.cli.<fn>`): one function per `CAT_<Tool>` binary,
  same positional order/defaults; reads files → calls the array API → writes files.
- **Downstream callers**: T1Prep's `surface_estimation.py` and the Nipype
  `nipype.interfaces.t1prep.cat_surf` interfaces call `cat_surf.cli.*`.

The two spherical-registration back-ends must stay interface-compatible — both take
`(source, source_sphere, target, target_sphere)` and return/write the warped source
sphere, so they are drop-in interchangeable:

| Algorithm | Binary | Array API | CLI mirror |
| --- | --- | --- | --- |
| DARTEL | `CAT_SurfWarp` | `cat_surf.surf_warp` | `cat_surf.cli.surf_warp` |
| Spherical Demons | `CAT_SurfSphericalDemon` | `cat_surf.spherical_demon` | `cat_surf.cli.surf_spherical_demon` |

Defaults live in the C source of truth (`Include/CAT_WarpDemons.h` +
`CAT_WarpDemonsDefaults`); the Cython signature/docstring and the docs must match it.

`cat_surf/_*.c` are Cython build artifacts, gitignored and regenerated from the `.pyx` at
build time — never commit them. Building needs Cython (`pip install cython`, or just build
through pip, which installs it from `build-system.requires`).

## The sulcal barrier in PBT (`-sulcal-barrier`)

The defaults in `CAT_PbtOptionsInit()` are the values tuned on real data and both front-ends
defer to them: `n_avgs` 5, `n_median_filter` 0, `median_subsample` 2, `sulcal_width` 5.0,
`barrier_gmtfactor` 1.8, `barrier_q` 0.7. `sulcal_barrier` itself defaults to **off**, as
`CAT_PpmSulciOpts::strength` does, so the correction stays out of the way until asked for.
`correct_thickness` is deliberately left at its old default: it compensates the border shift
of whichever segmentation produced the label map, so it is a per-pipeline value (0.0 for
AMAP, -0.05 otherwise in T1Prep) that no single library default can express.


The buried sulcus is created in the **distance map**, not at marching cubes. Where the
classifier lost the CSF in a sulcus there is no boundary for the CSF distance to stop at,
so the front from one bank runs through the fused grey matter into the other; GMT follows
it and the PPM never drops below the isovalue.

`CAT_VolSulcalMedialSet()` recovers the midline the front should have stopped at **from
geometry alone** — it is where the fronts from the two banks collide, detected as
`||grad dist_WM|| <= q`. No intensity image, no sheetness, no per-subject threshold. Gyral
crowns cannot false-positive, and that is structural rather than a guard: a front leaving a
blade has nothing to collide with, because outside the crown is CSF and not more cortex.

Two properties make it safe to leave on, and both are asserted in the tests:

- **Applied as a minimum.** Where CSF *was* segmented it is nearer than the midline, so the
  result is bit-identical to not running it. A correctly segmented sulcus also produces no
  collisions at all, because the CSF splits the band.
- **One-sided.** The distance can only shrink, so an overestimated thickness can be fixed
  but a correct one can never be inflated.

It gates the projection too — capping the distance alone is not enough, because
`projection_based_thickness()` would hand the thickness straight across the fused banks
again. The medial voxels are marked CSF in the scratch copy the projection reads, never in
the segmentation.

Measured on a fused-bank phantom against the same phantom correctly segmented (GMT 2.19 mm,
surfaces at x = 25/38): fusing gives GMT 4.19 and surfaces at 26/36; the barrier restores
GMT 2.19 exactly and surfaces to 25/37. `barrier_q` is flat between 0.5 and 0.7 and only
fails below ~0.4; `barrier_halfwidth` does not move the surface at all and trades against
thickness accuracy only. That insensitivity is the point — it is what the pre-PBT repair
could not deliver.

### Why there is no gyral dual

The mirror construction — fronts inward from the pial boundary, bounding `dist_WM` — detects
lost blades correctly on phantoms but **must not be used**, and the asymmetry is structural,
not a tuning problem. `PPM = (GMT - dist_WM)/GMT` has `dist_WM` in the *numerator*, so a
false collision does not perturb the thickness, it flips the voxel to WM-like (PPM -> 1).
On real data that produces stripes of spurious white matter inside the cortex — breaks in
the surface — while also shrinking `sum_dist` and so the thickness. `dist_CSF` enters only
through GMT, a far gentler dependence, which is why the sulcal direction is usable and the
gyral one is not. It was implemented, tested on real data, and removed.

### `barrier_gmtfactor` — the gate follows the cortex

The collision test alone is far too permissive on real data. Measured on an ADNI subject
(`lh.seg`, 0.5 mm), with the thickness gate disabled it caps **20.5% of the cortex** at
`q = 0.8` and drops the mean thickness from 2.31 mm to 1.81 mm — the surface-sampled result
was 1.44 mm against ~2.0 mm for FreeSurfer 8, CAT12 and T1Prep. It was truncating the whole
distribution, not fixing sulci.

The discriminator is that **a glued sulcus does not merely look thick, it looks like two
cortices back to back** — 5-6 mm where 2-3 mm is normal. Gating on the implied thickness at
the voxel, `dist_WM + dist_CSF`, separates the two populations; gating on the CSF distance
alone (`barrier_dmin`) does not, because plenty of ordinary voxels sit far from CSF simply
by being near the white matter.

| gate | capped | mean GMT | p95 |
| --- | --- | --- | --- |
| off | 20.5% | 1.814 | 2.122 |
| 4 mm | 5.7% | 2.242 | 3.115 |
| 5 mm | 2.4% | 2.277 | 3.213 |
| 6 mm | 1.0% | 2.290 | 3.234 |
| *no barrier at all* | — | 2.311 | 3.246 |

**The threshold is derived from the data, not fixed.** A glued sulcus is two cortices back to
back, so the gate belongs at a multiple of *this brain's* thickness — `barrier_gmtfactor`,
default 2.0 — rather than at a millimetre value that is only right for the cortex it was
tuned on. The proxy is `dist_WM + dist_CSF` in the GM band: for a band of locally constant thickness
those are complementary and sum to it exactly. It is summarised by the **mean of the values
below `barrier_gmtpct` (default p90)**, not by a median — the glued sulci the gate exists to
find sit in that upper tail, and a median only limits their influence, it still sits inside a
distribution they have skewed. Measured against the GMT finally reported, on four hemispheres
from two datasets, the ratio spans 0.087 for the trimmed mean against 0.102 for the median,
so a factor calibrated against it transfers better between subjects. On the ADNI subject it derives 2.72 mm -> 5.45 mm;
on a deliberately thinned copy of the same brain, 2.48 mm -> 4.96 mm, with the capped
fraction holding at 1.6-1.9% either way.

The one assumption is that glued sulci are a **minority** — true on a brain, where the
capped fraction is under 2%, but false on a phantom built entirely out of fused banks, where
the median measures the fusion itself and the derived gate comes out too loose. Pin
`-barrier-gmtmax` in that case.

A second caveat: that median is the *pre-projection* proxy and runs about 15% high
against the GMT finally reported (2.72 vs 2.36 mm here), so a factor of 2.0 is nearer 2.3x
the reported thickness. `-verbose` prints the value actually in use. `-barrier-gmtmax`
overrides it with an absolute millimetre value.

With the gate the parameters stop mattering, which was the whole point. Across
`q` in 0.6-0.9 and `barrier_tmin` in 0.4-0.8 the mean spans 2.252-2.293 mm — a spread of
0.04 mm, against 0.34 mm without it. `-verbose` reports the capped percentage: a few percent
is glued sulci being fixed, anything approaching double digits means the gate is too loose.

Two things that do **not** work, both measured rather than assumed:

- **Smoothing the arrival time before differentiating.** Intuition says a derivative operator
  wants a smoothed input; in fact flattening the field lowers `||grad T||` everywhere and
  makes the false-positive rate *worse* (13.8% -> 19.8% at `q = 0.6`).
- **Raising `q` towards 1.** `||grad T|| = 1` is the regular value of an uncollided front, so
  `q = 1` admits everything: on an all-healthy phantom it caps 56448 voxels and costs
  0.57 mm of mean thickness, where 0.4-0.8 cap none at all.

## The myelination boundary offset (`CAT_VolBoundaryOffset`)

In the primary motor and somatosensory strip, and along the line of Gennari in V1, the deep
cortical layers carry enough myelin to approach white-matter intensity on T1w. The classifier
follows the intensity, so the GM/WM boundary is placed too far out and the cortex comes back
too thin.

**Correcting this by relabelling voxels cannot work, and that is structural.** PBT already
carries the boundary to sub-voxel precision — `dist_WM` is refined by the PVE value itself in
`CAT_VolPbt.c`. The myelination displacement is 0.3-1.0 mm, the *same order* as that
precision, but a label can only express it by crossing the 2.5 isovalue, a quantised
half-voxel-minimum step. A previous `CAT_VolCorrectMyelination` did exactly that and was
either inert or catastrophic at every setting; it was removed. The quantity is continuous, so
the lever has to be.

**Moving the border and correcting the distance are the same edit.** With the pial at depth 0,
a true boundary at `T` and a segmented one at `T' = T - d`, a GM voxel at depth `x` has
`dist_WM_true = dist_WM_seg + d` and `GMT_true = GMT_seg + d`. Substituting into
`PPM = (GMT - dist_WM)/GMT` gives `(GMT - dist_WM)/(GMT + d)`. So adding `d` to `dist_WM`
before the projection grows the thickness by `d` **and** moves the central surface inward by
`d/2`, consistently, with no second knob. Note `correct_thickness` is already this correction
in global form — a constant in mm for "the systematic border shift of whichever segmentation
produced the label map" — but it is applied after PPM, so it never moves the surface.

**The observable is the width of the intensity transition, not its position.** The label map
is derived from the intensity, so the label boundary and the intensity boundary agree by
construction and their difference measures nothing — the first estimator tried this and
returned zero. What distinguishes myelinated cortex is the *shape* of the profile: healthy
cortex steps from grey to white over about the partial-volume width, while a myelinated deep
layer turns that step into a ramp one to three millimetres wide. The excess over this brain's
own healthy majority (`width_pct`, default p25) is the signal, and `gain` (0.5) converts it to
a displacement — for a symmetric ramp the boundary sits about half the excess beyond the
intensity midpoint.

The profile is read along the normal from `grad(label)`, in the contrast-normalized coordinate
`t = (T1w - GM_local)/(WM_local - GM_local)` from `CAT_VolLocalTissueReference()`. Both ends
are local, so a multiplicative bias cancels out of the ratio: measured on a slab phantom, a 10%
bias field changes the recovered displacement by under 3%.

**The band must be geometric, not a range of label values.** Selecting voxels whose label lies
near 2.5 samples a wide transition far more densely than a sharp one — on the phantom it
selected *only* myelinated voxels, so the reference width was derived from the pathology itself
and the excess came out at zero. Taking the voxels adjacent to the isosurface, found by a sign
change across a face, samples proportionally to surface area instead.

**Reject profiles that do not cross the boundary here.** Where the normal comes out tangential
— at the edge of the mask, or wherever `grad(label)` answers to something else — the ray runs
along the surface rather than through it and returns a meaningless width. Requiring the label
to cross 2.5 within `BO_ANCHOR_MM` of the start took the residual on a uniformly sharp phantom
from 748 displaced voxels to **exactly zero**.

Measured on a slab phantom with a widened transition inside a disc, at stock defaults:

| ramp | noise 0 | sigma 2 | sigma 2, 10% bias | sigma 4, 10% bias |
| --- | --- | --- | --- | --- |
| none | 0.000 | 0.012 / 0.014 | 0.012 / 0.013 | 0.028 / 0.037 |
| 2 mm | 0.190 / 0.013 | 0.211 / 0.026 | 0.205 / 0.025 | 0.229 / 0.047 |
| 3 mm | 0.378 / 0.026 | 0.374 / 0.037 | 0.366 / 0.036 | 0.366 / 0.057 |
| 4 mm | 0.569 / 0.040 | 0.568 / 0.051 | 0.555 / 0.050 | 0.555 / 0.070 |

*(mean displacement in mm, inside the disc / outside it)*

The displacement recovers `gain * (width - reference)` exactly, is 8-14x larger inside the
patch than outside, and is zero when there is nothing to find. The residue outside is mostly
`smooth_fwhm` spilling across the patch edge.

**Smooth before clamping.** Myelination is centimetre-scale and stereotyped while the
per-voxel profile is noisy, so `smooth_fwhm` (8 mm, applied *within* the boundary sheet by
normalized convolution over the band) is where the estimate becomes usable. It is also what
keeps the non-negativity clamp from half-wave-rectifying noise into a systematic inflation.

**This is a diagnostic before it is a correction.** `-width-out` writes the raw transition
width, which is the map to look at first — the displacement is only a rescaling of it. Run it
on a subject and check that it picks out the central sulcus and V1 before wiring it into PBT.
Everything above is phantom evidence.

## The signed sheetness offset (`CAT_VolMarchingCubes -sheet-offset`)

A global isovalue shift cannot fix glued sulci without breaking thin gyri, because it moves
every voxel the same way regardless of what is there. Measured on an ADNI PPM, lowering the
isovalue by 0.05 produces **88402 downward crossings of 0.5 and zero upward** — every change
opens a sulcus or thins a blade, indiscriminately.

`CAT_SheetnessOpts::signed_response` makes a valley negative and a ridge positive, so adding
the map lowers the PPM along sulci and raises it along blades in one pass. Same subject:

| | down | up |
| --- | --- | --- |
| global isovalue -0.05 | 88402 | **0** |
| global isovalue -0.10 | 170511 | **0** |
| signed offset 0.05 | 5133 | 5141 |
| signed offset 0.10 | 9980 | 10308 |

Balanced, and two orders of magnitude more surgical.

The defaults in `CAT_PpmSulciOptionsInit()` are the values tuned on real data and are what
both front-ends use unless told otherwise: `strength` 15, `sheet_strength` 30, `offset` 0.6,
`sigma_factor` 0.9, `cutoff` 0.4, `sheet_skeleton` off. That last one is a deliberate
trade against the flank argument below — the unthinned map does push the flanks the wrong
way, but across subjects the stronger correction won. `strength` itself defaults to **0**,
so the correction stays off until it is asked for; the tuned values above are what it uses
once `strength_sulci` is raised.

**The map must be skeletonized before it is used as an offset**, and that is not a
preference. The flanks of a valley curve upward and therefore read as *ridges*, carrying the
opposite sign: on a profile through a one-voxel valley the response runs
`+0.19 +0.29 -0.86 +0.29 +0.19`. Added unthinned, the offset would push the surface the wrong
way immediately beside every structure. Non-maximum suppression along the normal removes the
flanks and leaves the structure's own value untouched; `tests/test_sheetness.c` asserts both
halves of that. Skeletonization is **off** by default (`-sulci-skeleton` enables it), and it
is a weaker correction at the same offset -- roughly a third of the crossings -- so raise
`-sheet-offset` by about 3x when comparing against the unthinned field.

Everything internal to the filter stays defined on the magnitude; the sign is applied last,
after the anchor, the gain and the skeleton. Never hand a signed map to the oriented filters
-- they clamp to [0,1] and would silently discard every sulcus.

### The two halves are not symmetric (`-sheet-offset-gyri`)

One offset scales both halves of the signed map, and that is not what the data wants. Only
the *raising* half can re-glue banks: lowering a valley opens a sulcus, while raising a ridge
protects a thin blade but lifts the sulcal floor beside it at the same time. On real PPMs the
raising half is also the larger one, so the "balanced" offset **adds** tissue on net --
measured at `-sheet-offset 0.6`:

| subject | down (sulci opened) | up | ratio |
| --- | --- | --- | --- |
| OASIS (2.1 mm cortex) | 87608 | 133322 | 0.66 |
| ADNI 014_S_0328 | 111507 | 172166 | 0.65 |
| IXI199 | 152404 | 157145 | 0.97 |
| Aarhus (3.1 mm cortex) | 260306 | 235039 | 1.11 |

`offset_gyri` scales the raising half alone; negative (the default) means "same as `offset`",
so nothing changes unless it is asked for. Lowering it does **not** reduce the number of sulci
opened -- that count is identical at every value, because the two halves act on disjoint
voxels -- it only removes the inflation. `-verbose` now reports both crossing counts, which is
the diagnostic that shows whether the offset is opening the surface or growing it.

### `sigma_factor` derives a scale that is too coarse

`sigma_factor` ties `sigma_max` to the PPM's own median thickness (0.9x). The derived value is
the **worst** of the tested range on all four subjects above, and by a wide margin on the two
thin-cortex ones -- it maximizes the raising half without opening more sulci. A single fine
scale (`n_scales = 1`, i.e. `sigma_min` = 0.3 mm alone) is better everywhere:

| subject | derived | `sigma_max` 0.9 | `n_scales` 1 |
| --- | --- | --- | --- |
| OASIS | 0.66 | 0.89 | **1.09** |
| ADNI | 0.65 | 0.84 | **1.20** |
| Aarhus | 1.11 | 1.07 | **1.47** |
| IXI | 0.97 | 1.04 | **1.19** |

On OASIS and ADNI it raises the absolute opening too (87608 -> 102329, 111507 -> 129095). The
reason is in the `-sigma-max` help text: a sulcal CSF sheet at 0.5 mm is one to three voxels,
and 1.9 mm of sigma is about four -- large scales start answering to the cortical ribbon
itself, which is a ridge, not the sulcus. The defaults are unchanged pending more subjects,
but `sigma_factor` should be regarded as unproven rather than tuned.

**`sigma_max` used to be silently ignored.** The derivation runs whenever `sigma_factor > 0`
and overwrote whatever the caller passed, so `-sulci-sigma-max` -- the one knob the docs point
at -- did nothing at all. Setting it explicitly now switches the derivation off unless
`sigma_factor` is passed too. Both front-ends were affected.

## The sheetness family (`Include/CAT_Sheetness.h`)## The sheetness family (`Include/CAT_Sheetness.h`)

One shared shape prior feeds four tools, so a change to `Lib/CAT_Sheetness.c` propagates to
all of them — treat them as one unit. It exists because every isotropic regularizer (local
median, Potts MRF, TV) penalizes boundary area and therefore deletes thin structures
whichever side of the label boundary they lie on: the same filter that opens a glued sulcus
closes a cerebellar fissure.

| Consumer | Option |
| --- | --- |
| `CAT_VolSheetness` | the tool itself — writes the response map for tuning |
| `CAT_VolLocalStat` | `-oriented` (with `-stat 7`) |
| `CAT_VolThicknessPbt` | `-oriented-filter` |
| `CAT_VolMarchingCubes` | `-strength-sulci` — on the PPM, no intensity image needed |

**The response is anchored, so thresholds are data-independent.**
`CAT_SheetnessOpts::normalize` (default `CAT_SHEETNESS_NORMALIZE` = 1.0) scales the
response so its p99.9 is 1. The automatic noise scale is half the largest Hessian norm in
the volume, so the raw level depends on whatever the strongest structure in that image is —
which is why fixed thresholds used to need a per-dataset gain of 20 or more before anything
happened. The anchor removes that, and every threshold is now read as a fraction of it:

- `csf_thresh` / `wm_thresh` / `CAT_PpmSulciOpts::thresh` = **0.3** — full effect at twice
  that, 0.60, which is where the p99 of the reference response sits (raw 0.20 / p99.9 0.33).
- `CAT_ORIENTED_MEDIAN_CUTOFF` = **0.10**, i.e. preservation from 0.20. Lowered from the
  0.30 that the same derivation gives, because on real data the median needs to protect
  well below the p99 level to be useful.

Scaling is a single positive factor, so ranking, winning scale and the zero set are
untouched and the s = 0 invariant still holds exactly.

**`skeletonize` tightens what the thresholds see.** The plate response is as wide as the
Gaussian that produced it, so a large `sigma_max` locates a structure well but answers
several voxels into the tissue on either side, and the per-voxel gates then correct that
tissue too. Non-maximum suppression along the sheet normal collapses the band onto its
ridge line — one voxel at any scale — leaving the ridge value exactly unchanged. Off by
default; `-skeleton` on `CAT_VolSheetness`, `-sulci-skeleton` on `CAT_VolMarchingCubes`.
It runs *before* the anchor, so p99.9 is then taken over ridge values.

`sheet_strength` / `gain` survives as a deliberate relative adjustment but should no longer
be needed for calibration. Pass `-sheet-normalize 0` (or `-normalize 0` on
`CAT_VolSheetness`) to get the raw response back — the one case that needs it is an image
that may contain no sheets at all, where a percentile anchor amplifies noise.

**Invariant:** where the sheetness is zero, every oriented operator must be numerically
identical to the isotropic one it replaces. `tests/test_sheetness.c` asserts this voxel by
voxel — do not break it.

**Two thresholds decide whether any of this is visible**, and both must match the response
the data actually produces — check it with `CAT_VolSheetness` before tuning anything else.

1. `CAT_VolOrientedMedian()` admits a neighbour when `s*(dhat.n)^2 < cutoff`. A
   one-voxel-thick sheet is preserved from `s >= 2*cutoff`, so the cutoff is *half* the
   sheetness at which a thin structure starts being protected. Default
   `CAT_ORIENTED_MEDIAN_CUTOFF = 0.10`; override with `-sheet-cutoff` /
   `-oriented-cutoff`. The `s = 0` invariant holds at every cutoff.
2. The response itself is usually far lower than expected, because the automatic noise
   scale `c` is half the largest Hessian norm in the volume. Raise `-sheet-strength` (a
   gain, `CAT_SheetnessOpts::gain`), lower `-c`, or widen the scale range to bracket the
   voxel size — the scale defaults are tuned for 0.5 mm data.

On a 0.5 mm MPRAGE the dark-sheet response has `p99 = 0.20` and max `0.56`, so the old
cutoff of 0.5 changed 0.000% of brain voxels while 0.10 changes ~2.6%; on the same scan the
repair went from 0.23% to 4.6% of the GM band.



## Architecture rules

1. **Library-first:** All non-trivial logic belongs in `Lib/`, not `Progs/`.
2. A CLI tool in `Progs/` should only: parse args → load data → call lib →
   write results.
3. New feature → create `Include/CAT_Feature.h` + `Lib/CAT_Feature.c` first,
   then write the slim CLI.

## Adding a new library module

1. Create `Include/CAT_Feature.h` and `Lib/CAT_Feature.c`.
2. Add the `.c` to `libCAT_la_SOURCES` in `Makefile.am`.
3. Add the header to `noinst_HEADERS` in `Makefile.am`.
4. Keep `Makefile.in` in sync (it is tracked in the repo).
5. Re-run `./autogen.sh && ./configure`.

## Documentation — mandatory

Every public function declared in `Include/CAT_*.h` that is called from
`Progs/` **must** have Doxygen docs at **both** the header declaration and
the `Lib/` definition. Missing docs = documentation bug.

```c
/**
 * \brief One-line description.
 *
 * Longer explanation of the algorithm / edge cases.
 *
 * \param name   (in)     Description
 * \param result (out)    Description
 * \param buf    (in/out) Description
 * \return Description of return value or error codes
 */
```

## Coding style

- 4-space indentation, no tabs.
- Function opening brace on its own line (BSD style).
- No trailing whitespace; files must end with a newline.
- Keep functions small.
- Portable C99 — no compiler-specific extensions.

## Commit conventions

```
type: short summary under 50 chars

Body wrapped at 80 chars. Reference issues with (#N).
```

Types: `feat` `fix` `docs` `chore` `refactor` `test`

Commits should be small and logically atomic.

## Scope / ignore rules

Treat anything matched by `.gitignore` as out of scope unless explicitly asked.
Do not edit `3rdparty/` except when strictly necessary.
