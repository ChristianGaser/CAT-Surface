/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"
#include  <assert.h>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/quadmesh.c,v 1.21 2005/08/17 22:28:27 bert Exp $";
#endif


/*! \defgroup grp_quadmesh Quadmeshes
 *  \ingroup grp_bicobj
 *
 * A quadmesh consists of \em m rows and \em n columns of points.
 * Neighbouring points along a row and along a column are joined with
 * an edge, as on a piece of graph paper.  The first and last
 * points of each row (or of each column) may optionally have an edge
 * between them; picture rolling the graph paper into a
 * cylinder.  If both the rows and the columns are closed in this
 * manner, you end up with a toroidal topology.
 *
 * Point locations in a quadmesh are addressed using integer pairs
 * (i,j) with 0 <= i < m, and 0 <= j < n.
 *
 * Each quadrilateral of the mesh, or \em patch, is defined by the
 * four points at its corners.  The set of points comprising the interior of
 * the patch is undefined.  Four points are not planar, in general.
 *
 * Patches are addressed using integer pairs (i,j) with 0 <= i < mp,
 * and 0 <= j < np.  Use get_quadmesh_n_objects() to compute \c mp and
 * \c np.
 * 
 * Routines dealing with quadmeshes generally take a pointer
 * to a quadmesh_struct as their first argument.  
 *
 * @{
 */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_quadmesh
@INPUT      : col
              spr
              m
              n
@OUTPUT     : quadmesh
@RETURNS    : 
@DESCRIPTION: Initializes a quadmesh to the given size.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Initialize a quadmesh to given size.
 *
 * The resulting quadmesh uses a single colour, and is
 * open along the row and column directions.
 * Space is allocated for points and normals, but not
 * initialized.
 *
 * \param quadmesh pointer to structure to fill in
 * \param col colour for mesh
 * \param spr surface properties of mesh; set to NULL
 *            for default properties
 * \param m number of rows in mesh
 * \param n number of columns in mesh
 */
BICAPI  void  initialize_quadmesh(
    quadmesh_struct   *quadmesh,
    Colour            col,
    Surfprop          *spr,
    int               m,
    int               n )
{
    ALLOC( quadmesh->colours, 1 );

    quadmesh->colour_flag = ONE_COLOUR;

    quadmesh->colours[0] = col;

    if( spr != (Surfprop *) 0 )
        quadmesh->surfprop = *spr;
    else
        get_default_surfprop( &quadmesh->surfprop );

    quadmesh->m = m;
    quadmesh->n = n;
    quadmesh->m_closed = FALSE;
    quadmesh->n_closed = FALSE;
    if( quadmesh->m > 0 && quadmesh->n > 0 )
    {
        ALLOC( quadmesh->points, m * n );
        ALLOC( quadmesh->normals, m * n );
    }

    quadmesh->bintree = (bintree_struct_ptr) NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_quadmesh
@INPUT      : quadmesh
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes a quadmesh structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Delete quadmesh structure.
 */
BICAPI  void  delete_quadmesh(
    quadmesh_struct *quadmesh )
{
    free_colours( quadmesh->colour_flag, quadmesh->colours,
                  quadmesh->m * quadmesh->n,
                  (quadmesh->m-1) * (quadmesh->n-1) );

    if( quadmesh->m > 0 && quadmesh->n > 0 )
    {
        FREE( quadmesh->points );

        if( quadmesh->normals != (Vector *) NULL )
            FREE( quadmesh->normals );
    }

    delete_bintree_if_any( &quadmesh->bintree );   
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_quadmesh_point
@INPUT      : quadmesh
              i
              j
              point
              normal
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the point and normal of the quadmesh.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Set one point and normal in quadmesh.
 *
 * \param i row number 
 * \param j column number
 * \param point the point assigned to location (i,j)
 * \param normal the normal assigned to location (i,j)
 *               use NULL to avoid setting the normal
 *
 * \note The normal is ignored if the quadmesh does not already
 * have normals for each point.
 */
BICAPI  void  set_quadmesh_point(
    quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    Point            *point,
    Vector           *normal )
{
    int  ind;

    assert( 0 <= i && i < quadmesh->m );
    assert( 0 <= j && j < quadmesh->n );

    ind = IJ( i, j, quadmesh->n );

    quadmesh->points[ind] = *point;
    if( normal != NULL && quadmesh->normals != NULL )
        quadmesh->normals[ind] = *normal;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_point
@INPUT      : quadmesh
              i
              j
@OUTPUT     : point
@RETURNS    : TRUE if i,j is valid
@DESCRIPTION: Passes back the specified point of the quadmesh.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Get coordinates of one point in quadmesh.
 *
 * \param i row number 
 * \param j column number
 * \param point to which coordinates are copied
 * \return TRUE if (i,j) is a valid location
 *
 * \note The coordinates are written into \c point
 * only if (i,j) is a valid location.
 */
BICAPI  BOOLEAN  get_quadmesh_point(
    const quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    Point            *point )
{
    if( i < 0 || i >= quadmesh->m ||
        j < 0 || j >= quadmesh->n )
        return( FALSE );

    *point = quadmesh->points[IJ(i,j,quadmesh->n)];

    return( TRUE );
}
    
/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_n_objects
@INPUT      : quadmesh
@OUTPUT     : mp
              np
@RETURNS    : 
@DESCRIPTION: Passes back the number of patches in the quadmesh.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! Compute the number of patches in the quadmesh.
 *
 * If the quadmesh has \c m rows and \c n columns of
 * points, there will be \c np rows and \c np columns
 * of patches in the quadmesh.  The value of \c mp will
 * be m or m-1, depending upon whether the mesh is open or
 * closed in the row direction.
 *
 * The mesh has mp*np patches, in total.
 *
 * \param mp location in which to store mp
 * \param np location in which to store np
 */
BICAPI  void  get_quadmesh_n_objects(
    const quadmesh_struct  *quadmesh,
    int              *mp,
    int              *np )
{
    if( quadmesh->m_closed )
        *mp = quadmesh->m;
    else
        *mp = quadmesh->m - 1;

    if( quadmesh->n_closed )
        *np = quadmesh->n;
    else
        *np = quadmesh->n - 1;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_quadmesh_normals
@INPUT      : quadmesh
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Computes the normals of a quadmesh.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! Compute all normal vectors in quadmesh.
 *
 * Use find_polygon_normal() to set normal from location of
 * the point's 4-neighbours.
 *
 * \bug Does not take into account the closed flags when 
 * computing the 4-neighbourhood.
 */
BICAPI  void  compute_quadmesh_normals(
    quadmesh_struct  *quadmesh )
{
    int                i, j, m, n, n_neighs;
    Point              neighbours[4];
    progress_struct    progress;

    m = quadmesh->m;
    n = quadmesh->n;

    if( quadmesh->normals == (Vector *) NULL )
        ALLOC( quadmesh->normals, m * n );

    for_less( i, 0, m * n )
        fill_Vector( quadmesh->normals[i], 0.0, 0.0, 0.0 );

    initialize_progress_report( &progress, FALSE, m, "Computing Normals" );

    for_less( i, 0, m )
    {
        for_less( j, 0, n )
        {
            n_neighs = 0;
            if( get_quadmesh_point( quadmesh, i,   j-1, &neighbours[n_neighs]) )
                ++n_neighs;
            if( get_quadmesh_point( quadmesh, i+1, j  , &neighbours[n_neighs]) )
                ++n_neighs;
            if( get_quadmesh_point( quadmesh, i,   j+1, &neighbours[n_neighs]) )
                ++n_neighs;
            if( get_quadmesh_point( quadmesh, i-1, j  , &neighbours[n_neighs]) )
                ++n_neighs;

            if( n_neighs < 2 )
            {
                handle_internal_error( "compute_quadmesh_normals" );
            }
            else if( n_neighs == 2 )
            {
                neighbours[n_neighs] = quadmesh->points[IJ(i,j,n)];
                ++n_neighs;
            }

            find_polygon_normal( n_neighs, neighbours,
                                 &quadmesh->normals[IJ(i,j,n)] );
            NORMALIZE_VECTOR( quadmesh->normals[IJ(i,j,n)],
                              quadmesh->normals[IJ(i,j,n)] );
        }

        update_progress_report( &progress, i + 1 );
    }

    terminate_progress_report( &progress );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_patch_indices
@INPUT      : quadmesh
              i
              j
@OUTPUT     : indices
@RETURNS    : 
@DESCRIPTION: Passes back the indices of the quadmesh points making up the
              i,j'th patch
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Compute the indices of points for patch.
 *
 * A patch is addressed by a pair of integers (i,j) for
 * 0 <= i < mp, and 0 <= j < np, where \c mp and \c np are the
 * numbers of patches in each row and column, respectively.
 * Use get_quadmesh_n_objects() to compute \c mp and \c np.
 *
 * \param i row number of the patch
 * \param j column number of the patch
 * \indices location in which to write the 4 point indices
 *
 * The four indices returned are in cyclic order around the 
 * patch.  The point indices are returned in the following
 * order: (i,j), (i+1,j), (i+1,j+1), (i,j+1).
 */
BICAPI  void  get_quadmesh_patch_indices(
    const quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    int              indices[] )
{
    indices[0] = IJ(                   i,                   j, quadmesh->n );
    indices[1] = IJ( (i+1) % quadmesh->m,                   j, quadmesh->n );
    indices[2] = IJ( (i+1) % quadmesh->m, (j+1) % quadmesh->n, quadmesh->n );
    indices[3] = IJ(                   i, (j+1) % quadmesh->n, quadmesh->n );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadmesh_patch
@INPUT      : quadmesh
              i
              j
@OUTPUT     : points
@RETURNS    : 
@DESCRIPTION: Passes back the 4 points defining the i,j'th patch.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/*! \brief Get coordinates of patch corners.
 *
 * A patch is addressed by a pair of integers (i,j) for
 * 0 <= i < mp, and 0 <= j < np, where \c mp and \c np are the
 * numbers of patches in each row and column, respectively.
 * Use get_quadmesh_n_objects() to compute \c mp and \c np.
 *
 * \param i row number of the patch
 * \param j column number of the patch
 * \points location in which to write coordinates for the 4 points
 */
BICAPI  void  get_quadmesh_patch(
    const quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    Point            points[] )
{
    int   p, indices[4];

    get_quadmesh_patch_indices( quadmesh, i, j, indices );

    for_less( p, 0, 4 )
        points[p] = quadmesh->points[indices[p]];
}


/*! \brief Flip columns around centre.
 * 
 * Flip the quadmesh about its middle column, swapping point and
 * normal information in column 0 with column n-1, etc.
 */
BICAPI  void  reverse_quadmesh_vertices(
    quadmesh_struct  *quadmesh )
{
    int     i, j, i1, i2;
    Point   tmp_point;
    Vector  tmp_normal;

    for_less( i, 0, quadmesh->m )
    {
        for_less( j, 0, quadmesh->n / 2 )
        {
            i1 = IJ( i, j, quadmesh->n );
            i2 = IJ( i, quadmesh->n - 1 - j, quadmesh->n );

            tmp_point = quadmesh->points[i1];
            quadmesh->points[i1] = quadmesh->points[i2];
            quadmesh->points[i2] = tmp_point;

            if( quadmesh->normals != NULL )
            {
                tmp_normal = quadmesh->normals[i1];
                quadmesh->normals[i1] = quadmesh->normals[i2];
                quadmesh->normals[i2] = tmp_normal;
            }
        }
    }
}


/*!
 * @}
 * End of grp_quadmesh.
 */
