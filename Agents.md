# AGENTS

## Overview
CAT-Surface provides high-performance command-line tools for surface-based neuroimaging analysis, focusing on the processing, conversion, and analysis of cortical surface meshes.
These tools are integral to the [CAT12 toolbox](https://github.com/ChristianGaser/cat12) for SPM, and serve as the backend for most of CAT12’s surface-based processing, including extraction, smoothing, morphometric analysis, and format conversion. The tools are written in portable C and optimized for efficient batch processing of large cohorts in neuroimaging studies.

Most CAT-Surface functionality is implemented as a single static/shared library (`libCAT`) built from `Lib/` plus bundled third-party code in `3rdparty/`. Command-line tools in `Progs/` link against this library.

This file is written for coding agents working in this repository: it explains the layout, build system, and the most important dependencies between directories/files.

## External Packages
The only required external dependency is FFTW3.

**1. Install FFTW3:**

```
apt-get install libfftw3-dev
```

macOS:

```
brew install fftw
```

**2. Build the tools:**
```
./autogen.sh       # Generate the configure script
./configure        # Configure the build (add --prefix if needed)
make               # Build all tools
make install       # Optional: install system-wide
```

## Project Structure for OpenAI Codex Navigation

```
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
│   ├── MeshFix-V2.1/           # Mesh repair (self-intersections, GPL-3)
├── Include/
│   ├── CAT_*.h                 # Public-ish headers for CAT library modules
│   └── ...
├── Lib/
│   ├── CAT_*.c                 # Library modules compiled into libCAT
│   └── ...
├── Progs/
│   ├── CAT_*.c                 # CLI entry points; typically thin wrappers over libCAT
│   └── ...
├── cat_surface_python/         # Pure-Python subprocess wrapper around CLI tools
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
- Each `Progs/CAT_*.c` usually depends on:
	- `libCAT` (core algorithms in `Lib/`)
	- `bicpl-surface` (polygons structs, mesh utilities, argument parsing)
	- optional bundled IO libs (GIFTI/NIfTI)

### Bundled third-party code (3rdparty/)
- CAT-Surface intentionally vendors key dependencies so compilation does not require external downloads.
- Avoid modifying `3rdparty/` unless necessary; prefer fixing/adding functionality in `Lib/`.

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

## Development Guidelines
- Use C99 for compatibility with older environments.
- Keep functions small and well documented.
- For public C APIs, use clear comment blocks in headers and at function definitions.

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