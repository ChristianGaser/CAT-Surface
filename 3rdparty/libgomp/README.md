# Bundled static libgomp

This folder contains architecture-specific static GNU OpenMP runtime archives
used for building portable CAT-Surface binaries.

## Expected files

- `lib/libgomp-x86_64.a`
- `lib/libgomp-arm64.a` (optional)

During `configure`, Linux builds prefer these bundled archives over system
runtime linkage. If no bundled archive is available for the target
architecture, `configure` falls back to the compiler-provided `libgomp.a`.
