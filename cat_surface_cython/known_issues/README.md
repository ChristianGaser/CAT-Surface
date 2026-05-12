# Known issues — cat_surf Python bindings

Open bugs that don't block T1Prep but should be investigated. Each issue
has a self-contained repro script in this directory.

| Issue                                                        | Script                 | Status |
| ------------------------------------------------------------ | ---------------------- | ------ |
| `get_area_normalized` SIGABRT after long surface-op sequence | `smoke_test_repro.py`  | open   |

## get_area_normalized aborts after long surface-op sequence

### Symptom

`cat_surf.get_area_normalized(...)` aborts with `SIGABRT` (exit 134) when
called after a long sequence of other cat_surf surface operations in the
same Python process. The Python `faulthandler` reports the call site but
no C-level stack trace is emitted before the abort:

```text
Fatal Python error: Aborted

Current thread 0x00000001f38658c0 (most recent call first):
  File "smoke_test_repro.py", line NN in t_get_area_normalized
```

The abort happens reproducibly when the call follows the smoke-test
sequence in `smoke_test_repro.py` (≈20 prior surface ops
on a 640k-vertex MC mesh). Standalone calls — `reduce_mesh +
surf_to_sphere(stop_at=4) + get_area_normalized` in a fresh process —
do *not* abort. Removing arbitrary subsets of the prior ops also does
not abort, so the trigger is heap-layout-dependent rather than a
specific prior op.

### What we ruled out

- **macOS malloc-level heap corruption.** Running with
  `MallocScribble=1 MallocPreScribble=1 MallocGuardEdges=1
  MallocCheckHeapEach=1` does **not** emit any `malloc: ... corruption`
  / `GuardMalloc` warning before the abort. So it's not a vanilla
  buffer overrun the allocator can detect.
- **OpenMP thread-pool collision.** Repro uses `OMP_NUM_THREADS=1`.
- **Stop-at iteration count.** The actual aborting call is
  `get_area_normalized`, not `surf_to_sphere`. `surf_to_sphere(stop_at=6)`
  on the same mesh works fine standalone (78s on the full 640k-vertex
  mesh; 13s on the lh.central.gii from disk).
- **Mesh topology.** Both `polygons` and `sphere` arguments have
  matching topology (same V, same F).

### What this *isn't*

The bug **does not affect T1Prep's production path** — T1Prep does not
call `get_area_normalized` ([surface_estimation.py][1] only uses
`surf_to_sphere(stop_at=6)` and friends, all of which pass).

### Reproducer

```bash
PYTHONPATH=/path/to/cat_surface_cython \
  OMP_NUM_THREADS=1 \
  python -u -X faulthandler smoke_test_repro.py \
    --ppm /path/to/lh.ppm.<subject>.nii \
    --seg /path/to/lh.seg.<subject>.nii
```

Exit 0 → bug went away (heap layout shifted).
Exit 134 → still reproduces.

### Suggested next step

Rebuild libCAT and the cat_surf extensions with AddressSanitizer:

```bash
# libCAT
./configure CFLAGS="-fsanitize=address -fno-omit-frame-pointer -O1 -g"
make -j$(nproc)

# cat_surf extension (set CAT_BUILD_DIR to the ASan-built tree)
cd cat_surface_cython
CFLAGS="-fsanitize=address -fno-omit-frame-pointer -O1 -g" \
  LDFLAGS="-fsanitize=address" \
  python setup.py build_ext --inplace --force
```

Then run the repro under ASan — it should report the precise
heap-overflow / use-after-free / stack-overflow site with file:line.
That points directly at the culprit; bisecting allocation patterns is
not productive.

### Likely candidates inside `get_area_of_points_normalized_to_sphere`

[`Lib/CAT_Surf.c:219`][2] calls `resample_surface(polygons, sphere, 81920, NULL, NULL)`
and then `resample_values_sphere(...)`. Both routines do bintree
spatial searches and writes to caller-allocated buffers — natural
spots for a write-past-end that becomes lethal only when the heap
layout puts something important after the buffer.

[1]: ../../T1Prep/src/t1prep/surface_estimation.py
[2]: ../../Lib/CAT_Surf.c
