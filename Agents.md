# AGENTS

## Overview

CAT-Surface provides high-performance command-line tools for surface-based neuroimaging analysis, focusing on the processing, conversion, and analysis of cortical surface meshes.
These tools are integral to the [CAT12 toolbox](https://github.com/ChristianGaser/cat12) for SPM, and serve as the backend for most of CAT12’s surface-based processing, including extraction, smoothing, morphometric analysis, and format conversion. The tools are written in portable C and optimized for efficient batch processing of large cohorts in neuroimaging studies.

Most CAT-Surface functionality is implemented as a single static/shared library (`libCAT`) built from `Lib/` plus bundled third-party code in `3rdparty/`. Command-line tools in `Progs/` link against this library.

This file is written for coding agents working in this repository: it explains the layout, build system, and the most important dependencies between directories/files.

**2. Build the tools:**

```bash
./autogen.sh       # Generate the configure script
./configure        # Configure the build (add --prefix if needed)
make               # Build all tools
make install       # Optional: install system-wide
```

## Project Structure for OpenAI Codex Navigation

```bash
CAT-Surface/
├── 3rdparty/
│   ├── bicpl-surface/          # Core surface mesh utils (polygons, curvature, IO, args)
│   ├── volume_io/              # Low-level geometry + IO helpers used by bicpl/CAT
│   ├── nifti/                  # NIfTI IO
│   ├── gifticlib/              # GIFTI IO
│   ├── zlib/, expat/           # Bundled compression/XML deps
│   ├── dartel/, s2kit10/       # SPM dartel registration, spherical harmonics
│   ├── nii2mesh/               # Mesh simplification (quadric decimation)
│   ├── genus0/                 # Topology correction
├── Include/
│   ├── CAT_*.h                 # Public-ish headers for CAT library modules
│   └── ...
├── Lib/
│   ├── CAT_*.c                 # Library modules compiled into libCAT
│   └── ...
├── Progs/
│   ├── CAT_*.c                 # CLI entry points; typically thin wrappers over libCAT
│   └── ...
├── cat_surface_cython/         # Cython bindings exposing libCAT to Python
├── tests/                      # Small C unit tests (minunit)
├── ...
├── Makefile.am
├── configure.ac
├── autogen.sh
└── README.md
```

## Key Dependencies Between Files

### The "CAT" library boundary (Include/ + Lib/)

- `Include/CAT_<Module>.h` declares public functions/types.
- `Lib/CAT_<Module>.c` implements those functions.
- `libCAT` is built from `Lib/*.c` plus selected `3rdparty/*` sources.

Example: folding-related thickness correction

- Header: `Include/CAT_CorrectThicknessFolding.h`
- Implementation: `Lib/CAT_CorrectThicknessFolding.c`
- CLI tool (consumer): `Progs/CAT_SurfCorrectThicknessFolding.c`

If you move logic from a program into the library, the program should become a thin wrapper that:

1) parses args,
2) loads surface + data,
3) calls the library function,
4) writes results.

### Command-line tools (Progs/)

- **New CLI tools must delegate as much logic as possible to library functions in `Lib/`.** The program in `Progs/` should be a thin wrapper that:

 1) parses arguments,
 2) loads input data,
 3) calls one or more library functions to do the real work,
 4) writes results.

- When adding a new feature, first create a library module (`Include/CAT_<Feature>.h` + `Lib/CAT_<Feature>.c`), then write a slim CLI that calls it.
- Each `Progs/CAT_*.c` usually depends on:
  - `libCAT` (core algorithms in `Lib/`)
  - `bicpl-surface` (polygons structs, mesh utilities, argument parsing)
  - optional bundled IO libs (GIFTI/NIfTI)

### Bundled third-party code (3rdparty/)

- CAT-Surface intentionally vendors key dependencies so compilation does not require external downloads.
- Avoid modifying `3rdparty/` unless necessary; prefer fixing/adding functionality in `Lib/`.

### Python bindings (cat_surface_cython/ → `cat-surf`)

The Cython package exposes libCAT to Python in three layers that must stay consistent:

- **Array API** `cat_surf.<fn>` — numpy in/out, one function per algorithm.
- **CLI mirror** `cat_surf.cli.<fn>` — one function per `CAT_<Tool>` binary (snake_case,
  `CAT_` prefix dropped), same positional order and defaults; reads files → calls the
  array API → writes files. This is what downstream code (T1Prep `surface_estimation.py`,
  the Nipype `nipype.interfaces.t1prep.cat_surf` interfaces) calls.
- **Generated sources are not tracked.** `cat_surf/_*.c` are Cython build artifacts,
  ignored by `.gitignore` and produced from the `.pyx` at build time. Cython is in
  `pyproject.toml`'s `build-system.requires`, so every PEP 517 build (pip, cibuildwheel)
  has it; `setup.py` errors with an actionable message otherwise. Never commit them:
  they were ~25% of the repo's tracked bytes, and `setup.py build_ext` rewrites all nine
  whenever the local Cython version differs, producing huge unrelated diffs.
- **Source of truth for defaults** is the C side. For spherical registration the
  option struct and `CAT_WarpDemonsDefaults` in `Include/CAT_WarpDemons.h` /
  `Lib/CAT_WarpDemons.c` define the defaults that the Cython signature, its docstring,
  and the docs must all match (e.g. `n_steps` max is `CAT_WARP_DEMONS_MAX_STEPS` = 4).

The two spherical-registration back-ends are kept interface-compatible — both take
`(source, source_sphere, target, target_sphere)` and produce the warped source sphere:
`CAT_SurfWarp`/`surf_warp` (DARTEL) and `CAT_SurfSphericalDemon`/`surf_spherical_demon`
(Spherical Demons). When you change one, mirror the change in the other layer and the docs.

### The sheetness family (`Include/CAT_Sheetness.h`)

Several tools share one idea and therefore one implementation, so a change to
`Lib/CAT_Sheetness.c` propagates to all of them and they must be considered together.

The problem it solves: every isotropic regularizer — a local median, a Potts MRF, total
variation — penalizes boundary area, and a thin structure has an extreme area-to-volume
ratio, so deleting it is always the cheaper labelling (the *shrinking bias*). That is why
one and the same median filter opens glued sulci in one place and closes cerebellar
fissures in another. `CAT_VolSheetness()` replaces the smoothness prior with a *shape*
prior read off the Hessian eigenvalues, which keeps thin sheets, ignores blobs, and
shrinks nothing.

| Consumer | Option | What the field does there |
| --- | --- | --- |
| `Progs/CAT_VolSheetness` | (the tool itself) | writes the response map, for tuning scales and polarity |
| `Progs/CAT_VolLocalStat` | `-oriented` | median over a sheet-oriented neighbourhood (`-stat 7` only) |
| `Progs/CAT_VolSmooth` | `-oriented` | coherence-enhancing smoothing along the sheet |
| `Progs/CAT_VolAmap` | `-mrf-aniso 1\|2` | relaxes the Potts prior on sheets (local beta / direction weights) |
| `Progs/CAT_VolThicknessPbt` | `-oriented-filter` | replaces the three isotropic medians inside PBT |
| `Progs/CAT_VolSulcusRepair` | (always) | evidence term for the pre-PBT repair, `Lib/CAT_SulcusRepair.c` |

**Invariant every consumer relies on:** where the sheetness is zero the oriented operator
must be numerically identical to the isotropic one it replaces. That is what makes each of
these safe to enable by default-off. It is asserted voxel by voxel on both sides of the
binding boundary: `tests/test_sheetness.c` (C, via `make check`) and
`cat_surface_cython/tests/smoke_test.py` (Python, via CI). Do not break it.

Note the deliberate name split: the library module is `CAT_SulcusRepair` while the CLI is
`CAT_VolSulcusRepair`, because automake's `subdir-objects` would otherwise collide on two
objects named `CAT_VolSulcusRepair.o`.

## Build System Notes (Autotools/Libtool)

This project uses autotools.

Important files:

- `configure.ac` and `m4/`: configure checks (including FFTW detection).
- `Makefile.am`: authoritative source lists for the build.
- `Makefile.in`: tracked in-repo template; keep it in sync with `Makefile.am` changes.

When adding a new library module:

1) Add the `.c` file to `libCAT_la_SOURCES` in `Makefile.am`.
2) Add the public header to `noinst_HEADERS` in `Makefile.am`.
3) Keep `Makefile.in` in sync (this repo tracks it).
4) Re-run `./autogen.sh` and re-run `configure` in each build directory you care about.

When adding or changing a CLI tool:

- Update `Progs/Makefile.am` and ensure the program links `../libCAT.la`.

## Tests

- Unit tests live in `tests/` and currently focus on small, dependency-light C code.
- Typical flow: `make check` (after configure/make).
- The Python bindings have their own gate, `cat_surface_cython/tests/smoke_test.py`.
  It exists because the generated `cat_surf/_*.c` sources are untracked build
  artifacts, so nothing in the repository proves the `.pyx` sources still
  cythonize, compile and link — only an actual build does. It imports every
  extension module and re-checks the numeric contracts from the Python side;
  numpy is its only dependency, so it runs in the same environment the wheel is
  built in. CI runs it after `make check`:

  ```bash
  CAT_SURFACE_ROOT=$PWD CAT_BUILD_DIR=$PWD pip install ./cat_surface_cython
  python cat_surface_cython/tests/smoke_test.py
  ```

  Add a check here whenever you add a binding — a broken `.pyx` would otherwise
  first surface when a release wheel is built.

## Development Guidelines

- Use C99 for compatibility with older environments.
- Keep functions small and well documented.
- For public C APIs, use clear comment blocks in headers and at function definitions.
- **All public library functions used by CLI tools must be documented.**
  Any function declared in `Include/CAT_*.h` that is called (directly or
  indirectly) from `Progs/CAT_*.c` must have Doxygen-style documentation both:
  1) at the declaration in the header, and
  2) at the corresponding definition in `Lib/`.
  Missing docs for CLI-used public APIs are considered a documentation bug and
  should be fixed as part of the same change.
- **Documentation Style:** All public functions and significant internal utilities use **Doxygen-style documentation** with the following pattern:

  ```c
  /**
   * \brief One-line concise description of what the function does.
   *
   * Paragraph(s) explaining the algorithm, behavior, or special considerations.
   * Use \brief followed by multiple descriptive paragraphs for clarity.
   *
   * \param param1 (in)  type and description of first parameter
   * \param param2 (out) type and description of output parameter
   * \param param3 (in/out) bidirectional parameter description
   * \return Description of return value(s) or error codes
   */
  ```

  - **\brief** provides a concise one-line summary
  - **Paragraph descriptions** explain the algorithm and edge cases
  - **\param** tags document all parameters with direction indicators: (in), (out), or (in/out)
  - **\return** explains return values and error conditions
  - Examples: [CAT_Vol.c morphological operations](Lib/CAT_Vol.c#L1388-L1410), [CAT_Math.c type conversions](Lib/CAT_Math.c#L95-L108)
- **Library-first architecture:** place all non-trivial logic in `Lib/` library modules, not in `Progs/` CLI wrappers. CLI tools should only handle argument parsing, I/O, and calls to library functions.

## Ignore Rules

- Coding agents should treat any files or folders matched by `.gitignore` as out of scope (do not search, edit, or base decisions on them) unless the user explicitly asks to work with those ignored paths.

## Coding style

- Keep indentation at four spaces.
- Provide documentation for all public functions (header docs + brief implementation docs).
- Avoid large refactors unless required; keep changes surgical.

## Commit Guidance

- Break up work into small, logically separate commits.
- Commit messages use the form `type: short summary` where `type` can be `feat`, `fix`, `docs`, `chore`, etc.
- Reference issues when relevant, e.g. `fix: handle missing atlas path (#12)`.
- Use short (<50 char) summaries and include a blank line before the body.
- Wrap body lines at 80 characters.

## Pull Request Guidelines

When creating a PR:

1. Summarize the changes and link to relevant issues.
2. Mention any new dependencies or setup steps.

## Notes

These instructions are a starting point for contributors. Update this file as the project evolves.
