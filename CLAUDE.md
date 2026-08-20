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
| `CAT_VolThicknessPbt` | `-oriented-filter` |
| `CAT_VolSulcusRepair` | always (`Lib/CAT_SulcusRepair.c`) |
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
default; `-skeleton` on `CAT_VolSheetness`, `-sheet-skeleton` on `CAT_VolSulcusRepair`.
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
1b. `CAT_VolSulcusRepair` has the same structure with its own constants: it ignores a
   response below `csf_thresh` / `wm_thresh` (default 0.1) and ramps the blend weight to
   `csf_strength` / `wm_strength` at **twice** the threshold. Do not re-anchor that ramp
   at a sheetness of 1 — no real response reaches it, and that is what used to make the
   whole tool a no-op.
2. The response itself is usually far lower than expected, because the automatic noise
   scale `c` is half the largest Hessian norm in the volume. Raise `-sheet-strength` (a
   gain, `CAT_SheetnessOpts::gain`), lower `-c`, or widen the scale range to bracket the
   voxel size — the scale defaults are tuned for 0.5 mm data.

On a 0.5 mm MPRAGE the dark-sheet response has `p99 = 0.20` and max `0.56`, so the old
cutoff of 0.5 changed 0.000% of brain voxels while 0.10 changes ~2.6%; on the same scan the
repair went from 0.23% to 4.6% of the GM band.

Step 2 of the repair is `CAT_VolStrengthenWmBlades()` — it rescues thin WM fingers at the
gyral crowns. It replaced `CAT_VolReconnectGyri()`, whose two-sided gap test could not fire
at a blade tip (WM behind, GM in front) and so missed exactly the failure it was for.
Connectivity comes from the geodesic growth through the candidate set; sulci are excluded
by the sheetness polarity guard, not by the gap test. `-no-reconnect-gyri` and the Python
`reconnect_gyri` kwarg remain as aliases so downstream callers keep working.

The library module is `CAT_SulcusRepair` while the CLI is `CAT_VolSulcusRepair`: automake's
`subdir-objects` would otherwise collide on two objects of the same name.

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
