# Contributing to CAT-Surface

Thanks for your interest in contributing to the CAT surface tools. This document provides a short overview of how to build the project, the expected coding style and the usual pull request workflow.

## Build process

The repository uses the GNU autotools build system. In order to compile the project you will need the FFTW library and the usual build tools. A typical build looks as follows:

```bash
./autogen.sh
./configure    # requires FFTW to be available
make
make install   # optional
```

The same commands are listed in the [`README.md`](README.md) file. Please verify that the project still builds before sending a patch.

## Coding style

Follow the style of the surrounding code. Source files use four spaces for indentation and function braces start on a new line. Trailing whitespace should be avoided and all files must end with a newline.

Typos are checked with [`codespell`](https://github.com/codespell-project/codespell). The repository contains a configuration file [`.codespellrc`](.codespellrc) with an ignore list. Run `codespell --config .codespellrc` before submitting changes to ensure that no new misspellings are introduced.

## Pull request workflow

1. Fork the repository and create a feature branch for your changes.
2. Make sure the project builds and run `codespell` as described above.
3. Commit your changes with clear commit messages.
4. Open a pull request against the `master` branch and describe the changes you made.

Maintainers will review the request and may ask for updates. Thank you for improving CAT-Surface!
