# T1Prep `surface_estimation` port

Python port of the bash `surface_estimation()` function from
[T1Prep/scripts/T1Prep](https://github.com/ChristianGaser/T1Prep/blob/main/scripts/T1Prep),
using the `cat_surf` in-process wrappers instead of subprocesses to the
`CAT_*` command-line binaries.

This is a **draft** intended to be dropped into the T1Prep repo at
`src/t1prep/surface_estimation.py` and wired into the bash script with
the included `T1Prep.patch`.

## Files

| File | Purpose |
| --- | --- |
| `surface_estimation.py` | Python port of `surface_estimation()` |
| `T1Prep.patch` | Minimal patch for `scripts/T1Prep` that adds a feature-flagged Python path |

## Usage

After installing `cat-surf` into your T1Prep venv and dropping
`surface_estimation.py` into `src/t1prep/`, apply the patch:

```bash
cd /path/to/T1Prep
patch -p1 < /path/to/T1Prep.patch
```

Then opt in to the Python path with an environment variable:

```bash
T1PREP_USE_PYTHON_SURFACES=1 ./scripts/T1Prep sub-01_T1w.nii.gz
```

Unset the variable (or set it to `0`) to fall back to the original bash
binary-calling path, allowing A/B comparison.

## Coverage

13 of 14 CAT binary calls in `surface_estimation()` are fully replaced
in-process via `cat_surf` / `cat_surf.cli`. Two steps still fall through
to the binaries:

- **`CAT_SurfWarp -avg`** (rotated double-run) — the `-avg` mode of
  `cat_surf.surf_warp` isn't exposed yet. Adding it is ~15 lines in
  `_surf_warp.pyx`.
- **Atlas annot resampling** (`CAT_SurfResample -label … .annot`) —
  annot file I/O requires libCAT's `read_annotation_table` which isn't
  surfaced in `cat_surf._io` yet.
- **`CAT_SurfCurvature`** (fmriprep branch only) — no `cat_surf`
  wrapper exists. Skipped unless `--fmriprep 1`.

Set `T1PREP_DISABLE_FALLBACK=1` to make these raise instead of
silently invoking the binary (useful for CI).

## Behavioural differences from the bash version

- **`-downsample`** option of `CAT_VolThicknessPbt`: the post-PBT
  resampling step is not yet implemented in the Python port (the C
  library exposes only the core PBT call). Output is kept at native
  resolution; a warning is logged. Match the bash path by passing
  `--downsample 0`.
- **Topology-change map** (3rd positional argument of
  `CAT_VolMarchingCubes`): not currently surfaced through
  `cat_surf.vol_marching_cubes`. Only affects `--debug 1` runs.
- **Logging**: each step is logged as a one-line description plus
  `Execution time: <N>s`, mirroring the bash `run_cmd_log`. Live
  output goes to the same `report_log` file used by the bash script.
- **Progress bar**: the multi-process bar is preserved by sharing a
  count file between hemispheres. The bash and Python paths can be
  used interleaved (one hemisphere bash, one Python) without breaking
  the bar.

## Testing

```bash
# Smoke test (no real data needed)
python -m t1prep.surface_estimation --help

# Single-subject A/B comparison (requires a real T1Prep input)
T1PREP_USE_PYTHON_SURFACES=0 ./scripts/T1Prep sub-01_T1w.nii.gz  # baseline
T1PREP_USE_PYTHON_SURFACES=1 ./scripts/T1Prep sub-01_T1w.nii.gz  # ported
# diff the output surfaces / thickness values
```
