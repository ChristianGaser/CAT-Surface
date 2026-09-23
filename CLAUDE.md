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

Spherical registration is Spherical Demons:

| Algorithm | Binary | Array API | CLI mirror |
| --- | --- | --- | --- |
| Spherical Demons | `CAT_SurfSphericalDemon` | `cat_surf.spherical_demon` | `cat_surf.cli.surf_spherical_demon` |

The DARTEL back-end was deprecated: `CAT_SurfWarp`, `CAT_SurfApplyWarp`,
`CAT_SurfApplyWarpValues`, `CAT_SurfWarpDartel.[ch]` and the sphere solver
(`3rdparty/dartel/diffeosphere.c`, `optimizersphere.c`) now live in `deprecated/`, and
`cat_surf.surf_warp` / `cat_surf.cli.surf_warp` are gone. What remains of DARTEL is
`3rdparty/dartel/diffeopoly.c`: Spherical Demons uses its `init_dartel_poly()`, so four
small helpers (`pow2`, `dotprod`, `addscaled`, `norm`) were moved into it from the
deprecated files.

Defaults live in the C source of truth (`Include/CAT_WarpDemons.h` +
`CAT_WarpDemonsDefaults`); the Cython signature/docstring and the docs must match it.

For PBT, `cat_surf.cli.vol_thickness_pbt` reproduces `CAT_VolThicknessPbt` bit for bit: it
applies the blood-vessel correction first (`blood_vessel_correction=False` is `-no-bvc`). The
array API `cat_surf.vol_thickness_pbt` does **not** -- T1Prep runs its own vessel correction
before calling it, so a second one inside would apply it twice. Every other PBT keyword,
`oriented_filter` included, is a sentinel that keeps the `CAT_PbtOptionsInit()` default.

`cat_surf/_*.c` are Cython build artifacts, gitignored and regenerated from the `.pyx` at
build time — never commit them. Building needs Cython (`pip install cython`, or just build
through pip, which installs it from `build-system.requires`).

## The sulcal barrier in PBT (`-sulcal-barrier`)

The defaults in `CAT_PbtOptionsInit()` are the values tuned on real data and both front-ends
defer to them: `n_avgs` 5, `n_median_filter` 0, `median_subsample` 2, `sulcal_width` 5.0,
`barrier_gmtfactor` 1.5, `barrier_q` 0.7. `sulcal_barrier` itself defaults to **off**, as
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
default 1.5 — rather than at a millimetre value that is only right for the cortex it was
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

**The gate does not remove the parameter dependence -- it moves it into the factor.**
Measured on 5 subjects (T1Prep settings), the mean GM-band correction runs 0.215 / 0.151 /
0.117 / 0.094 mm at `barrier_gmtfactor` 1.3 / 1.5 / 1.7 / 2.0 and never levels off: about
0.027 mm per 10% change of the gate. No per-voxel or per-patch quantity separates fused from
normal cortex either (implied thickness, CSF distance and CSF/WM ratio at the medial voxels,
patch medians -- all unimodal). A stricter `q` cannot replace the gate: at `q = 0.2` without
it, 27-47% of the band is still capped, because `min(dist_CSF, dist_medial)` turns every
medial voxel into CSF for everything around it.

**Choosing the factor against the pial-white distance.** T1Prep uses 1.3, below the library
default of 1.5. Scored against Tfs on 8 hemispheres (precision: share of the removed thickness
lying where the unbarriered PBT exceeds Tfs by > 0.5 mm; recall: share of an excess > 1 mm that
is removed), 1.3 has the same precision as 1.5 (45-47% vs 45-46%), better recall (48-49% vs
35-38%) and puts the mean PBT within 0.02-0.03 mm of the mean Tfs (1.5: 0.04-0.05; no barrier:
0.13-0.14). Tfs is not independent of PBT -- the pial search starts at central + half the PBT
thickness -- so this was checked against two references, surfaces built from gate-1.3 PBT
(0.7.4) and from gate-1.5 PBT (0.7.1); the ranking is the same. Tfs is not ground truth either.

Lower factors do not improve on it (shared reference, same 8 hemispheres, precision / recall /
thickness removed where Tfs agrees / mean PBT - mean Tfs):

| factor | correction | capped | PPM pushed below 0.5 per 1000 | prec. | recall | removed @ agree | PBT - Tfs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.3 | 0.162 | 16% | 23 | 45-47% | 48-49% | 0.063 | +0.02/+0.03 |
| 1.2 | 0.198 | 20% | 28 | 43-46% | 55-56% | 0.074 | -0.01/+0.01 |
| 1.1 | 0.250 | 24% | 36 | 41-43% | 62-63% | 0.094 | -0.04/-0.03 |
| 1.0 | 0.319 | 27% | 47 | 39-41% | 70-71% | 0.127 | -0.09/-0.08 |
| 0.9 | 0.407 | 31% | 60 | 37-38% | 78-79% | 0.178 | -0.17/-0.15 |

Precision falls with every step and the false correction accelerates below 1.2, where the mean
also drops under Tfs (which itself runs short through shortcuts). 1.2 matches the Tfs mean best.

Through the full T1Prep surface pipeline (19 subjects, both hemispheres), 1.2 against 1.3 lowers
the mean PBT by 0.023 mm and the final thickness by 0.022 mm while Tfs moves only 0.009 mm, cuts
the vertices above 1.6x the median by 42% (0.93% -> 0.54%; BUSS02 rh 4.4% -> 2.2%) and the upper
thickness skewness by 0.11. The surfaces barely move (mean 0.043 mm, 0.58% of the vertices by
> 0.5 mm) and their quality does not change systematically: the glued fraction goes down in 20
and up in 18 hemispheres (it reacts strongly and in both directions to small changes), and one
Euler defect at 1.2 (002 rh) stands against one intersecting central surface at 1.3.

### Shared reference for both hemispheres (`barrier_gmtref`)

The reference estimate differs by up to 10% between the hemispheres of a subject (mean 4.2%
over 19 subjects), and the difference follows the amount of fusion, not the thickness: a
hemisphere with more fused sulci gets a *looser* gate and still *more* correction (r = 0.85
between gate and capped fraction). `CAT_VolPbtBarrierReference()` (`-barrier-ref-only`,
`cat_surf.vol_pbt_barrier_reference`) returns the reference a full run would derive -- same
preprocessing, bit-identical -- in ~13 s, and `barrier_gmtref` (`-barrier-gmtref`) hands it
back. The benefit is small and is about applying one criterion to both sides rather than about
the mean: at factor 1.3 on all 19 subjects the mean thickness is unchanged, a hemisphere's
band mean moves by 0.007 mm on average (at most 0.019 mm), and the mean |lh-rh| drops from
0.080 to 0.069 mm (smaller in 15/19; up to 0.04 mm where the references differ by 7-9%).
That matters for asymmetry analyses, where effects are of the same size. T1Prep does this by
default in the simplest way: each hemisphere process estimates both references itself (~20 s
extra) and passes the mean -- the estimate is deterministic, so both arrive at the same value
without any coordination -- and writes the references into its QA sidecar for the report
(`barrier_ref_lh/rh/shared`).

The reference was not what changed between T1Prep 0.7.1 and 0.7.3 (improved ventricle
filling): gate and correction moved by at most 0.7% and 0.011 mm. The hemisphere asymmetries
that appeared there came from `CAT_SurfCorrectThicknessFolding` -- see below.

### PBT is deterministic now

`projection_based_thickness()` used to vary between calls on the same input, by up to 1 mm
per voxel on the fused-bank phantom (most visible with the barrier on), for two reasons:
`pmax()` read a 15th neighbour past the end of its 14-element stack arrays, and `ornlm()`
added the contributions of neighbouring thread slabs to their shared slice in
timing-dependent order -- float sums depend on it, and the projection amplifies 1e-7 into
tenths of a millimetre. The slabs now run in two passes (even, then odd), which never share
a slice, so no lock is needed either. On real data this moves single voxels by up to 0.57 mm
and the mean by 1e-4 mm; the runtime is unchanged.

The CLI also ignored the library's `n_avgs` unless `-n-avgs` was given: it clamped the unset
sentinel to 1, and halved an explicit value twice under `-fast`. Results from
`CAT_VolThicknessPbt` without `-n-avgs` before this fix used a single distance level (the
reference on ADHD200 lh: 4.93 mm with one level, 4.85 mm with the default five).

Three things that do **not** work, all measured rather than assumed:

- **Smoothing the arrival time before differentiating.** Intuition says a derivative operator
  wants a smoothed input; in fact flattening the field lowers `||grad T||` everywhere and
  makes the false-positive rate *worse* (13.8% -> 19.8% at `q = 0.6`).
- **Raising `q` towards 1.** `||grad T|| = 1` is the regular value of an uncollided front, so
  `q = 1` admits everything: on an all-healthy phantom it caps 56448 voxels and costs
  0.57 mm of mean thickness, where 0.4-0.8 cap none at all.
- **Keeping only collisions with a label dip.** Across the ridge of `dist_WM` about 40% of the
  medial voxels show a dip of >= 0.1 in the label map (19% on the ADHD200 motion scan), the
  rest none -- two populations, but the dipping one has no gap to cut at. As the only gate it
  over-corrects: 0.39 mm mean correction at a dip of 0.1, 30% of the band capped, the mean PBT
  0.20 mm *below* Tfs and 0.29 mm removed where Tfs and the unbarriered PBT agree; the result
  still moves 0.43 -> 0.22 mm over dips of 0.05-0.3. Combined with the thickness gate it halves
  the correction without raising its precision (45% vs 47%). Dips of 0.1-0.3 are within the
  label noise of grey matter (p10-p90 1.85-2.20), so a dip does not prove a lost sulcus -- the
  information that makes the pial profile search useful is not available at the collision.

## Folding correction (`CAT_SurfCorrectThicknessFolding`)

The correction is applied where the smoothed mean curvature is positive. With the
outward-oriented surfaces CAT writes, that is **convex (gyral)** cortex -- the opposite of
FreeSurfer's `?h.curv` sign: on real central surfaces the mean curvature averages -0.12 in the
deepest quarter of the sulcal depth and +0.14 in the shallowest.

The sign used to be taken after the curvature was centred on its mean. A few hundred
degenerate vertices of the averaged central surface reach |H| ~ 1e4 against a p1-p99 range of
-0.9 to 0.3, so the mean -- and with it the selection -- jumped between hemispheres and
between versions: 0.1% of the vertices were corrected in one hemisphere, 99.6% in another.
That, not the sulcal barrier, is what made ADHD200 rh, BUSS02 rh, OASIS rh, HR075 lh and
yv98 lh change between T1Prep 0.7.1 and 0.7.3 (the PBT maps themselves agree to 0.015 mm).
With the sign taken before centring, 50-53% of every hemisphere is corrected, the
0.7.1 -> 0.7.3 change drops from 0.135 to 0.015 mm (max) and the lh-rh difference of the
correction from 0.144 to 0.043 mm (max). The same outliers make the Gaussian-curvature,
curvedness and mean-curvature columns of the regression nearly inert; in practice the shape
index drives it. Winsorizing them was measured and changes nothing worth having.

`pinv()` also left the part of `S` outside the leading rank block uninitialized, which only
matters for rank-deficient designs (`tests/test_folding.c`).

## Topology correction (`CAT_VolMarchingCubes`)

Three mechanisms decide whether a topological defect is opened or closed, and two
of them used to ignore anatomy entirely.

**The morphological step was the expensive one.** Before genus0 the binary volume is
opened or closed globally by `dist_morph`, and that distance used to be chosen by the
number of voxels *genus0* changed afterwards -- its own cost was never counted.
Measured on three hemispheres: it picked a 1.5 mm closing that changes 95448 voxels
(3.3% of the foreground, 94% of them inside sulcal sheets, i.e. every sulcus narrower
than 3 mm bridged) to save genus0 the 9 voxels it would otherwise have changed; on
another hemisphere a 1.0 mm opening, 16703 voxels, 89% on gyral ridges. The search now
counts the morphology's own voxels too, which picks `dist = 0` on all three, and Euler 2
is still reached in one or two iterations -- ADNI_014 lh gains 5.8% of surface area
(901 -> 954 cm2) that the opening had eaten.

**Filling is the safer default, and the order enforces it.** The pair of genus0 calls
resolves a defect by filling it if the first (filling, 6-connected) call can, and cuts
only what is left for the second. That asymmetry is deliberate: cutting a defect that sits
in a gyral blade severs the blade, and a severed blade is a hole in the surface, while an
unnecessary fill is a local thickening. `-topo-sheet` (default 0.3, 0 disables) does not
change the order. It runs genus0 twice on the same input to see both resolutions, and
takes away only the defects where filling is the wrong default -- a sulcus whose banks
nearly touch, which the filling closes into a bridge. Those are cut in the volume the pair
is then given; everything else reaches the pair untouched and is filled as before.

**The absolute sheetness of a defect says nothing, and reading it as if it did is what cut
gyri.** Filling adds *background* voxels, which are dark and therefore read as a valley;
cutting removes *foreground* voxels, which are bright and read as a ridge. Measured on four
hemispheres, all 15 contested regions had a negative fill mean (-0.05 to -0.57) and a
positive cut mean (+0.03 to +0.59) -- the sign is a property of the operation, not of the
anatomy. A rule that took the mean over both sets therefore cut 11 of those 15 regions,
i.e. it silently replaced fill-priority with cut-priority. What does carry information is
how strong the structure each action would damage is, so the two means are only ever
compared with each other: cut when the filling would run along a dark sheet of at least
`topo_sheet` **and** cutting damages less (`-mean_fill > mean_cut`). That selects 8 of the
15.

Two more things were measured rather than assumed:

- **Do not undo a decision afterwards.** Reverting genus0's cut only puts the handle back,
  a local closing rarely closes it, and the surface ends at Euler -6 instead of 2.
  Pre-resolving in the *input*, before genus0 runs, avoids that.
- **Do not run the local Euler pass on the result.** It only ever removed voxels, and the
  corner it picks is the one whose probability is closest to the isovalue -- exactly the
  voxels a fill has just added. It is only needed when two genus-0 results are mixed,
  which this no longer does.

**The rough pass before genus0 was removal-only for the same reason**, and there it is not
a mixing artefact but the rule itself: it resolved a 2x2x2 defect by flipping the most
ambiguous *foreground* corner, so a defect in a thin blade was always resolved by severing
the blade -- adding the one ambiguous hole voxel was never even a candidate, although the
apply step has always handled both directions. Both now compete on the same
|prob - thresh|. It fires rarely: over 36 hemispheres the surfaces move 0.002 mm on average
(0.03% of the vertices by more than 0.5 mm), the mean thickness by 0.001 mm and the area by
0.006%, while the glued fraction drops from 0.52% to 0.50% (lower in 10 hemispheres, higher
in 7) and one of the two central surfaces that still had self-intersections comes out clean.

Taking the cut run as the base also applied every removal it had made anywhere, not only
in the contested regions: 002 rh lost 1910 vertices that way on a hemisphere where genus0
found no defect at all.

`-verbose` prints one line per contested defect -- the two means and the decision -- which
is the diagnostic that separates the two failure modes.

At the T1Prep settings (`strength_sulci` 1.0, `sheet_offset` 0.2) the pre-cut lowers the
glued fraction of the marching-cubes surface -- 002 lh 0.37% -> 0.23%, BUSS02 rh 0.55% ->
0.25%, ADNI_014 lh 0.048% -> 0.042% -- and removes the intersections that the touching
banks caused (58 -> 0 and 37 -> 0 pairs), at Euler 2 everywhere. A hemisphere with no
contested defect is bit-identical to plain genus0.

Through the full T1Prep pipeline, on 19 subjects from the same segmentations (36
hemispheres; the 4397 scan is excluded, its GM/WM segmentation fails on huge ventricles):

| | genus0 alone | **pre-cut + symmetric pass** | mean over both sets, cut-priority |
| --- | --- | --- | --- |
| glued fraction | 0.574% | **0.502%** | 0.365% |
| mean thickness | 2.266 | 2.272 | 2.247 |
| area (mm2) | 99961 | 99987 | 99881 |
| white surfaces left with self-intersections | 4 (up to 121 pairs) | **0** | 1 (52) |
| pial | 2 (8, 24) | **0** | 0 |
| central | 1 (2) | 1 (11) | 0 |

Against plain genus0 the glued fraction falls in 17 of the 36 hemispheres and rises in 5,
every white and pial surface comes out clean, and the thickness and area barely move
(+0.006 mm, +0.03%). **16 hemispheres are bit-identical**: the rule only acts where a
defect is genuinely contested. The cut-priority column is the version this replaced -- it
reaches the lowest glued fraction of the three, and has the smallest area and the thinnest
cortex to go with it, which is what severed blades look like in the mean.

## The signed sheetness offset (`CAT_PpmSulciOpts::offset`)

`CAT_VolMarchingCubes` no longer exposes the offset knobs (removed in `7b8abd6`); they are set
through `CAT_PpmSulciOpts` or the `sheet_offset`, `sheet_offset_gyri` and `sulci_skeleton`
keywords of `cat_surf.vol_marching_cubes`. T1Prep passes those keywords, so they must stay
in the Python API. The option names below are the ones used when this was measured.

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

The defaults in `CAT_PpmSulciOptionsInit()` are what both front-ends use unless told
otherwise: `sheet_strength` 1.0, `offset` 0.2, `sigma_factor` 0.75, `cutoff` 0.2,
`sheet_skeleton` off -- the values T1Prep passes explicitly. (Before the response was
anchored they were `sheet_strength` 30, `offset` 0.6, `sigma_factor` 0.9, `cutoff` 0.4.)
`sheet_skeleton` off is a deliberate trade against the flank argument below — the unthinned
map does push the flanks the wrong way, but across subjects the stronger correction won.
`strength` itself defaults to **0**, so the correction stays off until it is asked for
(T1Prep uses 1.0); the values above are what it uses once `strength_sulci` is raised.

**The map must be skeletonized before it is used as an offset**, and that is not a
preference. The flanks of a valley curve upward and therefore read as *ridges*, carrying the
opposite sign: on a profile through a one-voxel valley the response runs
`+0.19 +0.29 -0.86 +0.29 +0.19`. Added unthinned, the offset would push the surface the wrong
way immediately beside every structure. Non-maximum suppression along the normal removes the
flanks and leaves the structure's own value untouched; `tests/test_sheetness.c` asserts both
halves of that. Skeletonization is **off** by default (`sulci_skeleton=True` enables it), and
it is a weaker correction at the same offset -- roughly a third of the crossings -- so raise
`offset` by about 3x when comparing against the unthinned field.

Everything internal to the filter stays defined on the magnitude; the sign is applied last,
after the anchor, the gain and the skeleton. Never hand a signed map to the oriented filters
-- they clamp to [0,1] and would silently discard every sulcus.

### The two halves are not symmetric (`offset_gyri`)

One offset scales both halves of the signed map, and that is not what the data wants. Only
the *raising* half can re-glue banks: lowering a valley opens a sulcus, while raising a ridge
protects a thin blade but lifts the sulcal floor beside it at the same time. On real PPMs the
raising half is also the larger one, so the "balanced" offset **adds** tissue on net --
measured at `offset` 0.6:

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

`sigma_factor` ties `sigma_max` to the PPM's own median thickness (0.9x when measured, 0.75x
since `7b8abd6`). The derived value is
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
itself, which is a ridge, not the sulcus. The factor has since been lowered to 0.75, but
`sigma_factor` should still be regarded as unproven rather than tuned.

**`sigma_max` used to be silently ignored.** The derivation runs whenever `sigma_factor > 0`
and overwrote whatever the caller passed, so `-sulci-sigma-max` -- the one knob the docs point
at -- did nothing at all. Setting it explicitly now switches the derivation off unless
`sigma_factor` is passed too. Both front-ends were affected.

## The sheetness family (`Include/CAT_Sheetness.h`)

One shared shape prior feeds four tools, so a change to `Lib/CAT_Sheetness.c` propagates to
all of them — treat them as one unit. It exists because every isotropic regularizer (local
median, Potts MRF, TV) penalizes boundary area and therefore deletes thin structures
whichever side of the label boundary they lie on: the same filter that opens a glued sulcus
closes a cerebellar fissure.

| Consumer | Option |
| --- | --- |
| `CAT_VolSheetness` | the tool itself — writes the response map for tuning |
| `CAT_VolLocalStat` | `-oriented` (with `-stat 7`) |
| `CAT_VolThicknessPbt` | always (`oriented_filter`, on by default; `-oriented-cutoff`) |
| `CAT_VolMarchingCubes` | `-strength-sulci` — on the PPM, no intensity image needed |

**The response is anchored, so thresholds are data-independent.**
`CAT_SheetnessOpts::normalize` (default `CAT_SHEETNESS_NORMALIZE` = 1.0) scales the
response so its p99.9 is 1. The automatic noise scale is half the largest Hessian norm in
the volume, so the raw level depends on whatever the strongest structure in that image is —
which is why fixed thresholds used to need a per-dataset gain of 20 or more before anything
happened. The anchor removes that, and every threshold is now read as a fraction of it:

- `CAT_PpmSulciOpts::thresh` = **0.3** — full effect at twice
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
default; `-skeleton` on `CAT_VolSheetness`, `sulci_skeleton=True` in
`cat_surf.vol_marching_cubes`.
It runs *before* the anchor, so p99.9 is then taken over ridge values.

`sheet_strength` / `gain` survives as a deliberate relative adjustment but should no longer
be needed for calibration. Pass `-normalize 0` to `CAT_VolSheetness` (`sulci_normalize=0`
in `cat_surf.vol_marching_cubes`) to get the raw response back — the one case that needs it is an image
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



## Pial placement by profile search (`CAT_Surf2PialWhite`, `Include/CAT_SurfPialProfile.h`)

The pial surface is placed per vertex from the label profile along its normal
(`CAT_PialWhiteOptions::pial_profile`, default on), starting from central surface + half
thickness; the white surface starts from ADE streamlines and keeps the balloon deformation
(`method` 2, the default -- see *Start surfaces* below). `-legacy-pial` / `pial_profile = 0`
restores the old pial path.

The balloon deformation stopped at label 1.7-1.8 instead of 1.5, for three separate reasons,
measured on HR075 -- tuning `w1`/`w2`/`w3`/`sigma` cannot fix any of them:

- **Smoothing the accumulated displacement** after the loop (with its `exp(-sigma)` damping)
  pulls convex crowns inwards: gyri are at 1.498 at the end of the loop and 1.68 after it.
  Smoothing the *update* inside the loop does not bias the fixed point, which is why annealing
  `sigma` changes nothing (1.693 vs 1.686).
- **`find_near_self_intersections` flags 23-37% of all vertices** per iteration, ~95% of them
  2-ring neighbours on the same sheet: the radius is 0.75x the mean edge on a decimated mesh with
  edges from 0.4 to 1.8 mm. Each flagged step is reverted 1.5x.
- **In glued sulci 1.5 does not exist**: 31% of vertices reach a valley first, 15% see a flat
  grey-matter plateau.

The profile search takes the isovalue crossing, or the valley bottom where the profile rises by
`valley_depth` above its running minimum first. Vertices without either in range are **held** --
letting them balloon ran them 3.8 mm out at p95. Facing sulcal walls (opposing normals, other
wall in front) clamp the outward step to half the gap, so both banks meet in the middle. The
normal Laplacian acts in concave regions only, because on a convex crown it pulls inwards.

**`search_out` is 2 mm, not 1.** With 1 mm, 60-69% of the vertices left more than 0.5 mm short
sat on a grey-matter plateau with the boundary 1.5-1.6 mm away. 2 mm halves that tail; 3 mm
adds nothing and more intersections. What remains (~5%) sits in concave fundi where walls touch
and where the central-surface normal used for evaluation diverges from the pial normal; a
smaller contact margin or a looser fold threshold recovers none of it.

Measured on 18 hemispheres (AD, controls, 7T from `T1Prep-0.7.1`, T1Prep settings `method 2
-remove_intersect`; label at the pial surface unless noted):

| | balloon | profile |
| --- | --- | --- |
| gyri | 1.78 (1.60-1.88) | 1.51 (1.509-1.517, median 1.500) |
| sulcal walls | 1.70 | 1.53 |
| offset from the 1.5 crossing | -0.44 mm | -0.02 mm |
| valley vertices: label - valley minimum | +0.155 | -0.009 |
| more than 0.5 mm short of the crossing | 37% | 5.1% |
| pial vertices in CSF (< 1.25) | 2.4% | 0.2% |

The pial-white distance moves from 1.89 to 2.29 mm against a PBT mean of 2.27; the exception is
the 3 mm Aarhus cortex (2.76 vs 2.93).

**The placed pial surface is smoothed with 2 HC Laplacian iterations** (`smooth_laplacian(.., 2,
0.1, 0.5)`, before the repair). Placement is per vertex and leaves the mesh visibly noisier than
the central surface; 2 iterations were the visual sweet spot. Measured on HR075 lh (umbrella =
mean distance of a vertex to its neighbour centroid; central surface 0.084):

| iterations | 0 | 1 | **2** | 3 | 5 | 10 |
| --- | --- | --- | --- | --- | --- | --- |
| umbrella | 0.102 | 0.086 | **0.082** | 0.080 | 0.078 | 0.076 |
| label mean | 1.658 | 1.658 | **1.659** | 1.659 | 1.660 | 1.662 |
| label MAE vs 1.5 | 0.173 | 0.184 | **0.189** | 0.192 | 0.196 | 0.201 |
| vertices in CSF (< 1.25) | 0.13% | 0.52% | **0.67%** | 0.78% | 0.92% | 1.10% |

The white surface needs nothing extra: `surf_deform_dual` already ends with 10 of these iterations.

**The label map, not the T1.** The bias- and LAS-corrected T1 is linear in the label scale
(3*T1: CSF 1.0, GM 2.0, WM 3.0) but opens only ~10% of the glued sulci, and its grey-matter spread
(p10-p90 1.85-2.20) exceeds most valley depths. Max-gradient targets were also worse than the
crossing on the label map (gyri 1.556, sulci 1.595).

**Gradient orientation.** `gradient3D()` differentiates along the voxel axes, but surface
normals and streamline positions live in world space. `surf_deform`, `surf_deform_dual` and the
ADE streamlines used it directly, so the gradient changed sign with every axis stored negatively
and was wrong by the rotation of oblique images -- 8 of the 9 test label maps.
`gradient3D_world_matrix()` now rotates it. On axis-aligned images stored with all axes negative
nothing changes (HR075 bit-identical); `tests/test_deform.c` asserts RAS and LAS storage give the
same surface for all three. The effect was small for the balloon terms (legacy pial gyri 1.86 ->
1.65 on RAS-stored images; `CAT_SurfDeform` on the PPM moves 0.002 mm, because `w3` = 1.0
dominates `w2` = 0.1) and large for ADE -- see below.

### Start surfaces: ADE for the white, thickness for the pial surface (`method` 2)

**The white intersections came from ADE.** Its streamlines stepped along the voxel-axis gradient
in world space, i.e. the wrong way: 361/143375 white streamlines converged on ADNI, 610/166444 on
HR075, and on mixed-sign storage they ran along corrupted paths (Aarhus: 282k intersecting pairs
at the start). `-remove_intersect` then left 0-440 per hemisphere. With world-space streamlines
99.8-100% converge.

Measured on the same 18 hemispheres; white error against the 2.5 crossing, CPU time per stage:

| white surface | thickness start + deformation | **ADE start + deformation** | ADE alone |
| --- | --- | --- | --- |
| label MAE | 0.121 | **0.091** | 0.247 |
| position MAE / bias | 0.174 / +0.124 mm | **0.135 / -0.047 mm** | 0.334 / -0.294 mm |
| self-intersections (repaired) | 135 (7T 365) | **0** | 0.5 |

| pial surface | **thickness start + profile** | ADE start + profile | ADE alone |
| --- | --- | --- | --- |
| gyri / sulcal walls | 1.512 / 1.529 | 1.502 / 1.482 | 1.168 / 1.226 |
| valley: label - valley minimum | **-0.009** | -0.172 | -0.431 |
| vertices in CSF | **0.2%** | 2.2% | 58% |
| self-intersections (repaired) | **0** | 600 | 724 |

ADE stops its streamlines at phi >= 0.999 / <= 0.001, the far end of the partial-volume ramp:
the white surface lands 0.29 mm inside WM, and in glued sulci the pial streamlines cross the
valley into the opposite bank, from where profile placement cannot bring them back. So ADE is
the better *start* for the white surface only, and never the final surface.

| CPU s (mean) | ADE | deform white | pial profile | repair white | repair pial | total |
| --- | --- | --- | --- | --- | --- | --- |
| `method` 0 | -- | 11.3 | 11.1 | 12.4 | 7.8 | 42.4 |
| `method` 1 | 6.6 | 9.4 | 13.1 | 5.5 | 25.5 | 60.1 |
| **`method` 2** | 6.6 | 9.4 | 11.1 | 5.4 | 7.6 | 40.1 |

In `method` 2 ADE traces only the white streamlines (`surf_ade_pial_white` takes `pial_out = NULL`),
which saves about 1 s of the 6-10 s; the solve dominates. The white fallback where |grad phi|
vanishes used to step outward (double negation) and now follows the inward normal.

ADE pays for itself through a shorter white deformation and a cheaper white repair; run alone,
`method` 2 took 43 s on HR075 and 40 s on yv98 against 35 / 48 s for `method` 0. With the ADE start the white target offset
is 0.1 (`GWM + 0.1`), not 0.2: 0.2 put the surface 0.10 mm inside WM, 0.0 is unbiased on average
but has the higher error (label MAE 0.097 vs 0.091).

**Why a few intersections always survived the repair.** `remove_intersections_iter` detected
defects once and then only re-tested the triangles it had labelled, but smoothing a defect can
make it cut unlabelled neighbours. Later passes kept smoothing the stale regions and the new
crossings were never seen: a pre-repair white surface went 1088 -> 114 pairs with "5 regions
left"; detecting afresh before every pass gives 0. On a marching-cubes mesh the old loop even
reported 0 regions while 1 pair remained. With the fix all 36 surfaces end at 0.

**What smoothing cannot repair at all: crossed sheets (`remove_intersections_ref`).** On the 38
hemispheres of the T1Prep test set, 2-4 white surfaces per run still ended with 86-269
intersecting pairs, and 1-3 pial surfaces with 3-126. They are not local folds: 99-100% of the
pairs join triangles more than 10 mesh steps apart, and on the central surface those partners
lie back to back 0.8-1.3 mm from each other. They are **thin gyral blades whose two sides were
driven through each other**: each side moves inward by 0.6-0.9 mm towards a white target the
label never reaches inside the blade (2.48 where 2.6 is asked for), the same way the pial
surface overshoots in a glued sulcus. The balloon deformation creates them -- ADNI_AD rh keeps
0 pairs from its ADE start surface and 86 after it -- and no amount of smoothing separates two
sheets that already cross: three times the passes or four times the iterations left 265 of 269.

The way back is the surface the deformation started from: the central surface for pial and
white, the start mesh for the central one. `remove_intersections_ref()` moves the vertices of
every defect the smoothing could not resolve, plus `CAT_RETREAT_RINGS` (2) rings of
neighbours, `CAT_RETREAT_FRACTION` (0.25) of the way back and repairs again, up to
`CAT_RETREAT_STEPS` (16) times. `CAT_Surf2PialWhite` and `CAT_SurfDeform` use it, and
`CAT_SurfFixSelfIntersect -reference` / `cat_surf.fix_self_intersect(reference=...)` expose it.
Measured on the five failing hemispheres: every one ends at 0 pairs, 0.13-0.55% of the vertices
move, the mean label error of the surface changes by at most 0.0004 and the pial-white distance
by at most 0.0004 mm. A reference is not needed everywhere, only where the defects are -- but
where the reference itself crosses (a central surface that came out of its own repair with
defects), the retreat cannot help either.

Profile placement of the white surface (the pial search on a mirrored label map, so a blade
ridge becomes a valley) was tried instead and rejected: more accurate (label MAE 0.03-0.05
against 0.07-0.10) but 200-600 pairs left after the repair, because the white surface has far
more thin structures for the facing-wall clamp to get wrong.

### `surf_deform` (`CAT_SurfDeform`, the central surface in T1Prep)

Two safeguards against self-intersections cost accuracy, measured on six central surfaces
(T1Prep steps 1-3 rebuilt, `-iter 75 -remove_intersect`; error = mean |PPM - 0.5| at the vertices):

- **The near-intersection check froze correct vertices.** `find_near_self_intersections` flags
  8-11% of a reduced central surface per iteration, 88-99.9% of them 2-ring neighbours of the
  same sheet. It now uses `find_near_facing_intersections` (opposing normals only): loop error
  0.0346 -> 0.0263, while the loop leaves 832 intersections for the repair (7410 without any
  check). A contact clamp or the fold revert of the pial placement do not help
  here -- these intersections are neither folds nor sulcal contacts.
- **The post-loop step moved outliers the wrong way and smoothed the result off the isovalue.**
  Displacements above the 95th percentile were replaced by one fixed vector (the per-axis
  percentiles, 6000-8400 vertices), and the total displacement was smoothed with a shrinking
  kernel. The outliers are now capped in length keeping their direction, and the smoothing is gone.

| | old | new |
| --- | --- | --- |
| PPM error | 0.0429 | **0.0275** |
| vertices within 0.1 of the isovalue | 91.5% | **95.5%** |
| self-intersections after `-remove_intersect` | 0 | 0 |

The surfaces move 0.08-0.15 mm on average. `surf_deform_dual` still uses the distance test.

The result is then smoothed with 2 HC Laplacian iterations, before `-remove_intersect`. Unlike
smoothing the accumulated displacement, this reduces the error: on HR075 lh the PPM error goes
0.0220 -> 0.0194 and the umbrella roughness 0.084 -> 0.072, with 0 intersections after the repair.

`-giter` (gradient refinement) was removed, from the Python binding (`gradient_iterations`) too:
it searched for the slope sample nearest to the vertex, which is the vertex itself, so it did
not move the surfaces.

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
4. Re-run `./autogen.sh && ./configure`.  `Makefile.in` is gitignored and
   regenerated by `autogen.sh`; `NOCONFIGURE=1 ./autogen.sh` regenerates it
   without configuring in-tree, which would block out-of-tree build directories.

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
