# AGENTS

## Overview
CAT-Surface provides high-performance command-line tools for surface-based neuroimaging analysis, focusing on the processing, conversion, and analysis of cortical surface meshes.
These tools are integral to the [CAT12 toolbox](https://github.com/ChristianGaser/cat12) for SPM, and serve as the backend for most of CAT12’s surface-based processing, including extraction, smoothing, morphometric analysis, and format conversion. The tools are written in portable C and optimized for efficient batch processing of large cohorts in neuroimaging studies.

The library code lives in the `Lib` directory and code for the command line tools is in `Progs/`. The repository also contains additional third-party libraries in `3rdparty/`.

## External Packages
**1. Install FFTW3:**

```
apt-get install libfftw3-dev
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
T1Prep/
├── 3rdparty/
├── Include/
├── Lib/
├── Progs/
├── ...
├── Makefile.am
├── configure.ac
├── autogen.sh
└── README.md
```

## Development Guidelines
- Use C99 for compatibility to older environments.
- Keep functions small and well documented. Include docstrings for public functions.
- Check documentation of functions and add missing documentation.

## Coding style
- Keep indentation at four spaces.
- Provide docstrings for all public functions and classes.
- Use consistent code formatting tools
- Follow language-specific best practices
- Keep code clean and readable

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