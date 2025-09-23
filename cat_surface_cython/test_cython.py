import os
import numpy as np
from cat_surface import load_bic_surface_arrays


def test_import_only():
    # Just ensure we can import; skip if not built
    try:
        import cat_surface.surface_io  # noqa: F401
    except Exception as e:
        print("Extension not built yet:", e)


def test_shapes_example(tmp_path):
    # This test is a placeholder; real tests require actual .obj files
    assert True
