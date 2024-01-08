#include "bicpl_internal.h"
#include "bicpl/marching.h"

#define  MAX_POLYGONS_PER_VOXEL    4
#define  MAX_INDICES_PER_VOXEL    12

typedef  struct
{
    int                n_polygons;
    int                sizes[MAX_POLYGONS_PER_VOXEL];
    voxel_point_type   indices[MAX_INDICES_PER_VOXEL];
} polygons_list;

typedef  struct
{
    int     offsets[4][3];
} face_struct;

typedef  struct
{
    int            n_ambiguities;
    face_struct    *ambiguity_faces;
    polygons_list  *polygons;

} case_struct;

typedef  enum  { PLUS_FLAG, MINUS_FLAG, MAX_CASES /*, ZERO_FLAG */ }
               Case_types;

static  case_struct   cases[2][2][2][2][2][2][2][2];

static   BOOLEAN  initialized = FALSE;

static  void  create_marching_cubes_lookup( void );
static  void  create_case(
    Case_types   case_flags[2][2][2],
    case_struct  *case_info );
static  void  get_face_axes(
    int          c,
    int          face,
    int          *a1,
    int          *a2 );
static  void  create_face_indices(
    int          c,
    int          face,
    face_struct  *face_indices );
static  void  get_face_flags(
    Case_types   case_flags[2][2][2],
    int          c,
    int          face,
    Case_types   face_flags[4] );
static  BOOLEAN  ambiguous_face_case(
    Case_types   face_flags[4] );
static  BOOLEAN  face_is_ambiguous(
    Case_types   case_flags[2][2][2],
    int          c,
    int          face );
static  BOOLEAN  surface_edge(
    Case_types     face_flags[4],
    int            edge );
static  void  get_edge_point(
    voxel_point_type  *edge_point,
    int               c,
    int               face,
    int               edge );
static  void  find_neighbour_face(
    int    c,
    int    face,
    int    edge,
    int    *new_c,
    int    *new_face,
    int    *new_edge );
static  void  delete_case(
    case_struct  *case_info );
static  void  create_case_polygons(
    Case_types     case_flags[2][2][2],
    BOOLEAN        face_ambiguity_flags[3][2],
    polygons_list  *polygons );

static  void  check_initialized( void )
{
    if( !initialized )
    {
        create_marching_cubes_lookup();

        initialized = TRUE;
    }
}

#ifdef DEBUG
static  void  test_case(
    case_struct  *c )
{
    int   p, n_sets, ind, poly, size, v, v1;

    n_sets = 1 << c->n_ambiguities;
    for_less( p, 0, n_sets )
    {
        ind = 0;
        for_less( poly, 0, c->polygons[p].n_polygons )
        {
            size = c->polygons[p].sizes[poly];
            for_less( v, 0, size-1 )
            {
                for_less( v1, v+1, size )
                {
                    if( c->polygons[p].indices[ind+v].coord[0] ==
                        c->polygons[p].indices[ind+v1].coord[0] &&
                        c->polygons[p].indices[ind+v].coord[1] ==
                        c->polygons[p].indices[ind+v1].coord[1] &&
                        c->polygons[p].indices[ind+v].coord[2] ==
                        c->polygons[p].indices[ind+v1].coord[2] &&
                        c->polygons[p].indices[ind+v].edge_intersected ==
                        c->polygons[p].indices[ind+v1].edge_intersected )
                        handle_internal_error( "test_case" );
                }
            }

            ind += size;
        }
    }
}
#endif

static  void  create_marching_cubes_lookup( void )
{
    Case_types   case_flags[2][2][2];

    for_enum( case_flags[0][0][0], MAX_CASES, Case_types )
    for_enum( case_flags[0][0][1], MAX_CASES, Case_types )
    for_enum( case_flags[0][1][0], MAX_CASES, Case_types )
    for_enum( case_flags[0][1][1], MAX_CASES, Case_types )
    for_enum( case_flags[1][0][0], MAX_CASES, Case_types )
    for_enum( case_flags[1][0][1], MAX_CASES, Case_types )
    for_enum( case_flags[1][1][0], MAX_CASES, Case_types )
    for_enum( case_flags[1][1][1], MAX_CASES, Case_types )
    {
        create_case( case_flags,
                     &cases[case_flags[0][0][0]]
                           [case_flags[0][0][1]]
                           [case_flags[0][1][0]]
                           [case_flags[0][1][1]]
                           [case_flags[1][0][0]]
                           [case_flags[1][0][1]]
                           [case_flags[1][1][0]]
                           [case_flags[1][1][1]] );

#ifdef DEBUG
        test_case( &cases[case_flags[0][0][0]]
                           [case_flags[0][0][1]]
                           [case_flags[0][1][0]]
                           [case_flags[0][1][1]]
                           [case_flags[1][0][0]]
                           [case_flags[1][0][1]]
                           [case_flags[1][1][0]]
                           [case_flags[1][1][1]] );

/*
{
    case_struct   c;
    int           p;

    c = cases[case_flags[0][0][0]]
             [case_flags[0][0][1]]
             [case_flags[0][1][0]]
             [case_flags[0][1][1]]
             [case_flags[1][0][0]]
             [case_flags[1][0][1]]
             [case_flags[1][1][0]]
             [case_flags[1][1][1]];

    print( "%d polygons, %d ambiguities:  ", c.polygons[0].n_polygons,
           c.n_ambiguities );

    for_less( p, 0, c.polygons[0].n_polygons )
    {
         print( " %d", c.polygons[0].sizes[p] );
    }

     print( "\n\n" );
}
*/
#endif
    }
}

static  void  create_case(
    Case_types   case_flags[2][2][2],
    case_struct  *case_info )
{
    int          i, amb, n_cases, c, face, n_ambiguities;
    int          ambiguous_faces[6][2];
    BOOLEAN      face_ambiguity_flags[3][2];
    face_struct  face_indices;

    n_ambiguities = 0;

    for_less( c, 0, N_DIMENSIONS )
    {
        for_less( face, 0, 2 )
        {
            if( face_is_ambiguous( case_flags, c, face ) )
            {
                ambiguous_faces[n_ambiguities][0] = c;
                ambiguous_faces[n_ambiguities][1] = face;
                ++n_ambiguities;
            }
        }
    }

    case_info->n_ambiguities = n_ambiguities;

    if( n_ambiguities > 0 )
    {
        ALLOC( case_info->ambiguity_faces, n_ambiguities );

        for_less( amb, 0, n_ambiguities )
        {
            c = ambiguous_faces[amb][0];
            face = ambiguous_faces[amb][1];

            create_face_indices( c, face, &face_indices );
   
            case_info->ambiguity_faces[amb] = face_indices;
        }
    }

    n_cases = 1 << n_ambiguities;

    ALLOC( case_info->polygons, n_cases );

    for_less( i, 0, n_cases )
    {
        for_less( c, 0, N_DIMENSIONS )
        {
            for_less( face, 0, 2 )
            {
                face_ambiguity_flags[c][face] = FALSE;
            }
        }

        for_less( amb, 0, n_ambiguities )
        {
            if( (i & (1 << amb)) != 0 )
            {
                c = ambiguous_faces[amb][0];
                face = ambiguous_faces[amb][1];

                face_ambiguity_flags[c][face] = TRUE;
            }
        }

        create_case_polygons( case_flags, face_ambiguity_flags,
                              &case_info->polygons[i] );
    }
}

static  void  get_face_axes(
    int          c,
    int          face,
    int          *a1,
    int          *a2 )
{
    if( face == 0 )
    {
        *a1 = (c + 2) % N_DIMENSIONS;
        *a2 = (c + 1) % N_DIMENSIONS;
    }
    else
    {
        *a1 = (c + 1) % N_DIMENSIONS;
        *a2 = (c + 2) % N_DIMENSIONS;
    }
}

static  void  create_face_indices(
    int          c,
    int          face,
    face_struct  *face_indices )
{
    int   a1, a2;

    get_face_axes( c, face, &a1, &a2 );

    face_indices->offsets[0][a1] = 0;
    face_indices->offsets[0][a2] = 0;
    face_indices->offsets[0][c]  = face;

    face_indices->offsets[1][a1] = 1;
    face_indices->offsets[1][a2] = 0;
    face_indices->offsets[1][c]  = face;

    face_indices->offsets[2][a1] = 1;
    face_indices->offsets[2][a2] = 1;
    face_indices->offsets[2][c]  = face;

    face_indices->offsets[3][a1] = 0;
    face_indices->offsets[3][a2] = 1;
    face_indices->offsets[3][c]  = face;
}

static  void  get_face_flags(
    Case_types   case_flags[2][2][2],
    int          c,
    int          face,
    Case_types   face_flags[4] )
{
    int          i;
    face_struct  face_indices;

    create_face_indices( c, face, &face_indices );

    for_less( i, 0, 4 )
    {
        face_flags[i] = case_flags[face_indices.offsets[i][0]]
                                  [face_indices.offsets[i][1]]
                                  [face_indices.offsets[i][2]];
    }
}

static  BOOLEAN  ambiguous_face_case(
    Case_types   face_flags[4] )
{
    BOOLEAN      ambiguous;

    ambiguous = (face_flags[0] == MINUS_FLAG && face_flags[2] == MINUS_FLAG &&
                 face_flags[1] == PLUS_FLAG  && face_flags[3] == PLUS_FLAG) ||
                (face_flags[0] == PLUS_FLAG  && face_flags[2] == PLUS_FLAG &&
                 face_flags[1] == MINUS_FLAG && face_flags[3] == MINUS_FLAG);

    return( ambiguous );
}

static  BOOLEAN  face_is_ambiguous(
    Case_types   case_flags[2][2][2],
    int          c,
    int          face )
{
    Case_types   corners[4];

    get_face_flags( case_flags, c, face, corners );

    return( ambiguous_face_case( corners ) );
}

typedef  struct
{
    int      n_edges;
    int      edge_points[4];
    BOOLEAN  edge_used[2];
} edges_struct;

static  void  follow_edge(
    polygons_list  *polygons,
    int            *ind,
    edges_struct   edges[3][2],
    int            c,
    int            face,
    int            edge_index );
static  void  find_voxel_edge_index(
    int            c,
    int            face,
    int            edge,
    edges_struct   edges[3][2],
    int            *edge_index );
static  void  create_edges(
    Case_types     case_flags[2][2][2],
    BOOLEAN        face_ambiguity_flags[3][2],
    edges_struct   edges[3][2] );
static  void  create_edge_for_face(
    Case_types     face_flags[4],
    BOOLEAN        face_ambiguity_flag,
    edges_struct   *edges );
static  void  follow_edge(
    polygons_list  *polygons,
    int            *ind,
    edges_struct   edges[3][2],
    int            c,
    int            face,
    int            edge_index );

static  void  create_case_polygons(
    Case_types     case_flags[2][2][2],
    BOOLEAN        face_ambiguity_flags[3][2],
    polygons_list  *polygons )
{
    int               c, face, edge_index, ind, prev_ind, size;
    edges_struct      edges[3][2];

    create_edges( case_flags, face_ambiguity_flags, edges );

    polygons->n_polygons = 0;

    ind = 0;

    for_less( c, 0, N_DIMENSIONS )
    {
        for_less( face, 0, 2 )
        {
            for_less( edge_index, 0, edges[c][face].n_edges )
            {
                if( !edges[c][face].edge_used[edge_index] )
                {
                    if( polygons->n_polygons >= MAX_POLYGONS_PER_VOXEL )
                    {
                        HANDLE_INTERNAL_ERROR( "n polygons" );
                    }

                    prev_ind = ind;

                    follow_edge( polygons, &ind, edges, c, face, edge_index );

                    size = ind - prev_ind;

                    polygons->sizes[polygons->n_polygons] = size;

                    if( (size < 3) ||
                        (size > MAX_INDICES_PER_VOXEL) )
                    {
                        HANDLE_INTERNAL_ERROR( "n ind" );
                    }

                    ++polygons->n_polygons;
                }
            }
        }
    }
}

static  void  create_edges(
    Case_types     case_flags[2][2][2],
    BOOLEAN        face_ambiguity_flags[3][2],
    edges_struct   edges[3][2] )
{
    int          c, face;
    Case_types   face_flags[4];

    for_less( c, 0, N_DIMENSIONS )
    {
        for_less( face, 0, 2 )
        {
            get_face_flags( case_flags, c, face, face_flags );

            create_edge_for_face( face_flags, face_ambiguity_flags[c][face],
                                  &edges[c][face] );
        }
    }
}

static  void  create_edge_for_face(
    Case_types     face_flags[4],
    BOOLEAN        face_ambiguity_flag,
    edges_struct   *edges )
{
    int   n_edges, edge, pair[2];

    edges->edge_used[0] = FALSE;
    edges->edge_used[1] = FALSE;

    if( ambiguous_face_case( face_flags ) )
    {
        if( face_flags[0] == PLUS_FLAG )
        {
            if( face_ambiguity_flag )
            {
                edges->edge_points[0] = 1;
                edges->edge_points[1] = 0;
                edges->edge_points[2] = 3;
                edges->edge_points[3] = 2;
            }
            else
            {
                edges->edge_points[0] = 3;
                edges->edge_points[1] = 0;
                edges->edge_points[2] = 1;
                edges->edge_points[3] = 2;
            }
        }
        else
        {
            if( face_ambiguity_flag )
            {
                edges->edge_points[0] = 0;
                edges->edge_points[1] = 3;
                edges->edge_points[2] = 2;
                edges->edge_points[3] = 1;
            }
            else
            {
                edges->edge_points[0] = 0;
                edges->edge_points[1] = 1;
                edges->edge_points[2] = 2;
                edges->edge_points[3] = 3;
            }
        }

        edges->n_edges = 2;
    }
    else
    {
        n_edges = 0;
        for_less( edge, 0, 4 )
        {
            if( surface_edge( face_flags, edge ) )
            {
                pair[n_edges] = edge;
                ++n_edges;
            }
        }

        if( n_edges != 0 && n_edges != 2 )
        {
            handle_internal_error( "n_edges in marching_cubes" );
        }

        edges->n_edges = n_edges / 2;

        if( n_edges == 2 )
        {
            if( face_flags[(pair[0]+1) % 4] == PLUS_FLAG )
            {
                edges->edge_points[0] = pair[0];
                edges->edge_points[1] = pair[1];
            }
            else
            {
                edges->edge_points[0] = pair[1];
                edges->edge_points[1] = pair[0];
            }
        }
    }
}

static  BOOLEAN  surface_edge(
    Case_types     face_flags[4],
    int            edge )
{
    Case_types   flag1, flag2;

    flag1 = face_flags[edge];
    flag2 = face_flags[(edge+1) % 4];

    return( (flag1 == PLUS_FLAG  && flag2 == MINUS_FLAG) ||
            (flag1 == MINUS_FLAG && flag2 == PLUS_FLAG) );
}

static  void  follow_edge(
    polygons_list  *polygons,
    int            *ind,
    edges_struct   edges[3][2],
    int            c,
    int            face,
    int            edge_index )
{
    int   edge;

    while( !edges[c][face].edge_used[edge_index] )
    {
        edge = edges[c][face].edge_points[2*edge_index+1];

        get_edge_point( &polygons->indices[(*ind)], c, face, edge );

        ++(*ind);

        edges[c][face].edge_used[edge_index] = TRUE;

        find_neighbour_face( c, face, edge, &c, &face, &edge );

        find_voxel_edge_index( c, face, edge, edges, &edge_index );
    }
}

static  void  get_edge_point(
    voxel_point_type  *edge_point,
    int               c,
    int               face,
    int               edge )
{
    int           a1, a2;
    face_struct   face_indices;

    create_face_indices( c, face, &face_indices );

    get_face_axes( c, face, &a1, &a2 );

    edge_point->coord[c] = face;
    edge_point->coord[a1] = 0;
    edge_point->coord[a2] = 0;

    switch( edge )
    {
    case 0:
        edge_point->edge_intersected = a1;
        break;

    case 1:
        ++edge_point->coord[a1];
        edge_point->edge_intersected = a2;
        break;

    case 2:
        ++edge_point->coord[a2];
        edge_point->edge_intersected = a1;
        break;

    case 3:
        edge_point->edge_intersected = a2;
        break;
    }
}

static  void  find_neighbour_face(
    int    c,
    int    face,
    int    edge,
    int    *new_c,
    int    *new_face,
    int    *new_edge )
{
    static  BOOLEAN  initialized = FALSE;
    static  struct
    {
        int   c, face, edge;
    }        neighbours[3][2][4];
    int          c1, face1, edge1, en1;
    int          c2, face2, edge2, en2, n_matches;
    face_struct  face1_ind, face2_ind;

    if( !initialized )
    {
        initialized = TRUE;

        for_less( c1, 0, N_DIMENSIONS )
        for_less( face1, 0, 2 )
        {
            create_face_indices( c1, face1, &face1_ind );
            for_less( edge1, 0, 4 )
            {
                en1 = (edge1 + 1) % 4;
                n_matches = 0;
                for_less( c2, 0, N_DIMENSIONS )
                for_less( face2, 0, 2 )
                {
                    create_face_indices( c2, face2, &face2_ind );
                    for_less( edge2, 0, 4 )
                    {
                        en2 = (edge2 + 1) % 4;
                        if( face1_ind.offsets[edge1][0] ==
                            face2_ind.offsets[en2][0] &&
                            face1_ind.offsets[edge1][1] ==
                            face2_ind.offsets[en2][1] &&
                            face1_ind.offsets[edge1][2] ==
                            face2_ind.offsets[en2][2] &&
                            face1_ind.offsets[en1][0] ==
                            face2_ind.offsets[edge2][0] &&
                            face1_ind.offsets[en1][1] ==
                            face2_ind.offsets[edge2][1] &&
                            face1_ind.offsets[en1][2] ==
                            face2_ind.offsets[edge2][2] )
                        {
                            neighbours[c1][face1][edge1].c = c2;
                            neighbours[c1][face1][edge1].face = face2;
                            neighbours[c1][face1][edge1].edge = edge2;
                            ++n_matches;
                        }
                    }
                }
                if( n_matches != 1 )
                    handle_internal_error( "n_matches" );
            }
        }
    }

    *new_c = neighbours[c][face][edge].c;
    *new_face = neighbours[c][face][edge].face;
    *new_edge = neighbours[c][face][edge].edge;

    if( neighbours[*new_c][*new_face][*new_edge].c != c ||
        neighbours[*new_c][*new_face][*new_edge].face != face ||
        neighbours[*new_c][*new_face][*new_edge].edge != edge )
    {
        HANDLE_INTERNAL_ERROR( "neighbours" );
    }
}

static  void  find_voxel_edge_index(
    int            c,
    int            face,
    int            edge,
    edges_struct   edges[3][2],
    int            *edge_index )
{
    BOOLEAN  found;

    found = FALSE;

    for_less( *edge_index, 0, edges[c][face].n_edges )
    {
        if( edges[c][face].edge_points[2*(*edge_index)] == edge )
        {
            found = TRUE;
            break;
        }
    }

    if( !found )
    {
        HANDLE_INTERNAL_ERROR( "find edge index" );
    }
}

BICAPI  int  get_holeless_isosurface_polygons(
    Real              corners[2][2][2],
    Real              isovalue,
    int               *sizes[],
    voxel_point_type  *points[] )
{
    Real         corner_values[2][2][2];
    int          amb, amb_index;
    face_struct  *face;
    int          c0, c1, c2, c3, c4, c5, c6, c7;
    case_struct  *voxel_case;
    Real         m1, m2, p1, p2;

    check_initialized();

    corner_values[0][0][0] = corners[0][0][0] - isovalue;
    corner_values[0][0][1] = corners[0][0][1] - isovalue;
    corner_values[0][1][0] = corners[0][1][0] - isovalue;
    corner_values[0][1][1] = corners[0][1][1] - isovalue;
    corner_values[1][0][0] = corners[1][0][0] - isovalue;
    corner_values[1][0][1] = corners[1][0][1] - isovalue;
    corner_values[1][1][0] = corners[1][1][0] - isovalue;
    corner_values[1][1][1] = corners[1][1][1] - isovalue;

    c0 = (corner_values[0][0][0] <= 0.0);
    c1 = (corner_values[0][0][1] <= 0.0);
    c2 = (corner_values[0][1][0] <= 0.0);
    c3 = (corner_values[0][1][1] <= 0.0);
    c4 = (corner_values[1][0][0] <= 0.0);
    c5 = (corner_values[1][0][1] <= 0.0);
    c6 = (corner_values[1][1][0] <= 0.0);
    c7 = (corner_values[1][1][1] <= 0.0);

    voxel_case = &cases[c0][c1][c2][c3][c4][c5][c6][c7];

    if( voxel_case->n_ambiguities > 0 )
    {
        amb_index = 0;
        for_less( amb, 0, voxel_case->n_ambiguities )
        {
            face = &voxel_case->ambiguity_faces[amb];

            if( corner_values[face->offsets[0][0]][face->offsets[0][1]]
                             [face->offsets[0][2]] <= 0.0 )
            {
                m1 = corner_values[face->offsets[0][0]]
                                  [face->offsets[0][1]]
                                  [face->offsets[0][2]];
                m2 = corner_values[face->offsets[2][0]]
                                  [face->offsets[2][1]]
                                  [face->offsets[2][2]];
                p1 = corner_values[face->offsets[1][0]]
                                  [face->offsets[1][1]]
                                  [face->offsets[1][2]];
                p2 = corner_values[face->offsets[3][0]]
                                  [face->offsets[3][1]]
                                  [face->offsets[3][2]];
            }
            else
            {
                p1 = corner_values[face->offsets[0][0]]
                                  [face->offsets[0][1]]
                                  [face->offsets[0][2]];
                p2 = corner_values[face->offsets[2][0]]
                                  [face->offsets[2][1]]
                                  [face->offsets[2][2]];
                m1 = corner_values[face->offsets[1][0]]
                                  [face->offsets[1][1]]
                                  [face->offsets[1][2]];
                m2 = corner_values[face->offsets[3][0]]
                                  [face->offsets[3][1]]
                                  [face->offsets[3][2]];
            }

            if( (p2 - m2) * m1 > (m1 - p1) * p2 )
                amb_index |= (1 << amb);
        }
    }
    else
    {
        amb_index = 0;
    }

    *sizes = voxel_case->polygons[amb_index].sizes;
    *points = voxel_case->polygons[amb_index].indices;

    return( voxel_case->polygons[amb_index].n_polygons );
}

BICAPI  void  delete_marching_cubes_table( void )
{
    Case_types   case_flags[2][2][2];

    if( initialized )
    {
        for_enum( case_flags[0][0][0], MAX_CASES, Case_types )
        for_enum( case_flags[0][0][1], MAX_CASES, Case_types )
        for_enum( case_flags[0][1][0], MAX_CASES, Case_types )
        for_enum( case_flags[0][1][1], MAX_CASES, Case_types )
        for_enum( case_flags[1][0][0], MAX_CASES, Case_types )
        for_enum( case_flags[1][0][1], MAX_CASES, Case_types )
        for_enum( case_flags[1][1][0], MAX_CASES, Case_types )
        for_enum( case_flags[1][1][1], MAX_CASES, Case_types )
        {
            delete_case( &cases[case_flags[0][0][0]]
                               [case_flags[0][0][1]]
                               [case_flags[0][1][0]]
                               [case_flags[0][1][1]]
                               [case_flags[1][0][0]]
                               [case_flags[1][0][1]]
                               [case_flags[1][1][0]]
                               [case_flags[1][1][1]] );
        }
    }
}

static  void  delete_case(
    case_struct  *case_info )
{
    if( case_info->n_ambiguities > 0 )
        FREE( case_info->ambiguity_faces );

    FREE( case_info->polygons );
}
