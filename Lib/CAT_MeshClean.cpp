/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/**
 * @file CAT_MeshClean.cpp
 * @brief Surface mesh cleaning using MeshFix.
 *
 * C++ implementation wrapping MeshFix's meshclean functionality
 * with conversion between bicpl polygons_struct and TMesh.
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstdint>

/* CAT interface - extern "C" for C linkage */
/* NOTE: Include CAT headers BEFORE MeshFix to avoid name collisions */
extern "C" {
#include "CAT_MeshClean.h"
}

/* MeshFix headers - must come after CAT headers */
#include "tmesh.h"
#include "tin.h"
#include "vertex.h"

/* Typedefs to avoid namespace collision with Point/Vector */
typedef T_MESH::Basic_TMesh    TMesh;
typedef T_MESH::Node           TNode;
typedef T_MESH::Vertex         TVertex;
typedef T_MESH::Triangle       TTriangle;
typedef T_MESH::ExtVertex      TExtVertex;
typedef T_MESH::coord          TCoord;

/**
 * Convert bicpl polygons_struct to MeshFix TMesh.
 */
static int polygons_to_tmesh(const polygons_struct *polys, TMesh *mesh)
{
    if (!polys || !mesh || polys->n_points < 3 || polys->n_items < 1)
        return -1;

    int n_vertices = polys->n_points;
    int n_triangles = polys->n_items;

    /* Add vertices */
    for (int i = 0; i < n_vertices; i++) {
        TCoord x = (TCoord)Point_x(polys->points[i]);
        TCoord y = (TCoord)Point_y(polys->points[i]);
        TCoord z = (TCoord)Point_z(polys->points[i]);
        mesh->V.appendTail(mesh->newVertex(x, y, z));
    }

    /* Create ExtVertex array for indexed triangle creation */
    TExtVertex **var = (TExtVertex **)malloc(sizeof(TExtVertex *) * n_vertices);
    if (!var) return -1;

    int idx = 0;
    TNode *node;
    TVertex *v;
    for (node = mesh->V.head(); node != NULL; node = node->next()) {
        v = (TVertex *)node->data;
        var[idx++] = new TExtVertex(v);
    }

    /* Add triangles */
    int poly_idx = 0;
    for (int t = 0; t < n_triangles; t++) {
        int n_verts = polys->end_indices[t] - poly_idx;
        if (n_verts != 3) {
            /* Only triangles supported */
            fprintf(stderr, "CAT_MeshClean: Non-triangle face detected (polygon %d has %d vertices)\n", t, n_verts);
            poly_idx = polys->end_indices[t];
            continue;
        }

        int i1 = polys->indices[poly_idx];
        int i2 = polys->indices[poly_idx + 1];
        int i3 = polys->indices[poly_idx + 2];
        poly_idx = polys->end_indices[t];

        if (i1 < 0 || i1 >= n_vertices ||
            i2 < 0 || i2 >= n_vertices ||
            i3 < 0 || i3 >= n_vertices) {
            fprintf(stderr, "CAT_MeshClean: Invalid vertex index at triangle %d\n", t);
            continue;
        }

        if (i1 == i2 || i2 == i3 || i3 == i1) {
            /* Degenerate triangle, skip */
            continue;
        }

        mesh->CreateIndexedTriangle(var, i1, i2, i3);
    }

    /* Cleanup ExtVertex array */
    for (int i = 0; i < n_vertices; i++) {
        delete var[i];
    }
    free(var);

    /* Update Euler characteristics */
    mesh->eulerUpdate();

    return 0;
}

/**
 * Convert MeshFix TMesh back to bicpl polygons_struct.
 */
static int tmesh_to_polygons(TMesh *mesh, polygons_struct *polys)
{
    if (!mesh || !polys)
        return -1;

    int n_vertices = mesh->V.numels();
    int n_triangles = mesh->T.numels();

    if (n_vertices < 3 || n_triangles < 1)
        return -1;

    /* Free old polygon data if present */
    if (polys->n_points > 0 && polys->points) {
        free(polys->points);
        polys->points = NULL;
    }
    if (polys->n_items > 0) {
        if (polys->normals) { free(polys->normals); polys->normals = NULL; }
        if (polys->indices) { free(polys->indices); polys->indices = NULL; }
        if (polys->end_indices) { free(polys->end_indices); polys->end_indices = NULL; }
    }

    /* Initialize polygon structure */
    polys->n_points = n_vertices;
    polys->n_items = n_triangles;
    
    polys->points = (VIO_Point *)malloc(n_vertices * sizeof(VIO_Point));
    polys->normals = (VIO_Vector *)malloc(n_vertices * sizeof(VIO_Vector));
    polys->indices = (int *)malloc(n_triangles * 3 * sizeof(int));
    polys->end_indices = (int *)malloc(n_triangles * sizeof(int));
    
    if (!polys->points || !polys->normals || !polys->indices || !polys->end_indices) {
        return -1;
    }

    /* Copy vertices and assign indices for triangle lookup */
    TNode *node;
    TVertex *v;
    int idx = 0;
    for (node = mesh->V.head(); node != NULL; node = node->next()) {
        v = (TVertex *)node->data;
        Point_x(polys->points[idx]) = (Real)v->x;
        Point_y(polys->points[idx]) = (Real)v->y;
        Point_z(polys->points[idx]) = (Real)v->z;
        /* Store index in vertex info for later triangle lookup */
        v->info = (void *)(intptr_t)idx;
        idx++;
    }

    /* Copy triangles */
    TTriangle *t;
    idx = 0;
    for (node = mesh->T.head(); node != NULL; node = node->next()) {
        t = (TTriangle *)node->data;
        TVertex *v1 = t->v1();
        TVertex *v2 = t->v2();
        TVertex *v3 = t->v3();

        int i1 = (int)(intptr_t)v1->info;
        int i2 = (int)(intptr_t)v2->info;
        int i3 = (int)(intptr_t)v3->info;

        polys->indices[idx * 3] = i1;
        polys->indices[idx * 3 + 1] = i2;
        polys->indices[idx * 3 + 2] = i3;
        polys->end_indices[idx] = (idx + 1) * 3;
        idx++;
    }

    /* Compute normals */
    compute_polygon_normals(polys);

    return 0;
}

/* C interface implementation */

extern "C" {

void CAT_MeshCleanOptionsInit(CAT_MeshCleanOptions *opts)
{
    if (!opts) return;
    opts->max_iters = 10;
    opts->inner_loops = 3;
    opts->fill_holes = 1;
    opts->verbose = 0;
}

int CAT_SurfMeshClean(
    polygons_struct *polygons,
    const CAT_MeshCleanOptions *opts
)
{
    CAT_MeshCleanOptions default_opts;
    if (!opts) {
        CAT_MeshCleanOptionsInit(&default_opts);
        opts = &default_opts;
    }

    if (!polygons || polygons->n_points < 3 || polygons->n_items < 1) {
        fprintf(stderr, "CAT_SurfMeshClean: Invalid input surface\n");
        return -1;
    }

    /* Create TMesh from polygons */
    TMesh *mesh = new TMesh();
    
    if (polygons_to_tmesh(polygons, mesh) != 0) {
        fprintf(stderr, "CAT_SurfMeshClean: Failed to convert surface to TMesh\n");
        delete mesh;
        return -1;
    }

    if (opts->verbose) {
        printf("Input mesh: %d vertices, %d triangles\n",
               mesh->V.numels(), mesh->T.numels());
        printf("Boundaries: %d, Shells: %d\n",
               mesh->boundaries(), mesh->shells());
    }

    /* Remove small components */
    int sc = mesh->removeSmallestComponents();
    if (sc && opts->verbose) {
        printf("Removed %d small components\n", sc);
    }

    /* Fill holes if requested */
    if (opts->fill_holes && mesh->boundaries()) {
        if (opts->verbose) {
            printf("Filling %d boundary holes...\n", mesh->boundaries());
        }
        mesh->fillSmallBoundaries(0, true);
    }

    /* Run meshclean - the main intersection/degeneracy removal */
    int result = 0;
    if (!mesh->boundaries()) {
        bool success = mesh->meshclean(opts->max_iters, opts->inner_loops);
        if (!success) {
            if (opts->verbose) {
                fprintf(stderr, "MeshFix could not fix all issues\n");
            }
            result = 1;
        }
    } else {
        if (opts->verbose) {
            fprintf(stderr, "Mesh still has boundaries after hole filling\n");
        }
        result = 1;
    }

    if (opts->verbose) {
        printf("Output mesh: %d vertices, %d triangles\n",
               mesh->V.numels(), mesh->T.numels());
    }

    /* Convert back to polygons */
    if (tmesh_to_polygons(mesh, polygons) != 0) {
        fprintf(stderr, "CAT_SurfMeshClean: Failed to convert TMesh back to surface\n");
        delete mesh;
        return -1;
    }

    delete mesh;
    return result;
}

int CAT_SurfCountIntersections(polygons_struct *polygons)
{
    if (!polygons || polygons->n_points < 3 || polygons->n_items < 1) {
        return -1;
    }

    /* Create TMesh from polygons */
    TMesh *mesh = new TMesh();
    
    if (polygons_to_tmesh(polygons, mesh) != 0) {
        delete mesh;
        return -1;
    }

    /* Detect intersections */
    mesh->deselectTriangles();
    int count = mesh->selectIntersectingTriangles();

    delete mesh;
    return count;
}

} /* extern "C" */
