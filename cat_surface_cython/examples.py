from cat_surface import load_bic_surface

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python examples.py /path/to/surface.obj")
        sys.exit(1)
    surf = load_bic_surface(sys.argv[1])
    print(surf)
    print("Vertices:", surf.vertices.shape, surf.vertices.dtype)
    print("Faces:", surf.triangles.shape, surf.triangles.dtype)
