#include "bicpl_internal.h"
#include "bicpl/marching.h"

#define  MAX_POLYGONS_PER_VOXEL  5
#define  MAX_INDICES_PER_VOXEL   20

typedef  struct
{
    int                n_polygons;
    int                *sizes;
    voxel_point_type   *indices;
} case_struct;

typedef  enum  { PLUS_FLAG, MINUS_FLAG, MAX_CASES }
               Case_types;

static  case_struct   cases[2][2][2][2][2][2][2][2][2][2][2];

static   BOOLEAN  initialized = FALSE;

/*  ---------------- static prototypes ------------------------- */

static  void  create_marching_cubes_lookup( void );
static  void  create_case(
    int          x,
    int          y,
    int          z,
    Case_types   case_flags[2][2][2],
    case_struct  *case_info );

/*  ------------------------------------------------------------- */

static  void  check_initialized( void )
{
    if( !initialized )
    {
        create_marching_cubes_lookup();

        initialized = TRUE;
    }
}

static  void  create_marching_cubes_lookup( void )
{
    int          x, y, z;
    Case_types   case_flags[2][2][2];

    for_less( x, 0, 2 )
    for_less( y, 0, 2 )
    for_less( z, 0, 2 )
    for_enum( case_flags[0][0][0], MAX_CASES, Case_types )
    for_enum( case_flags[0][0][1], MAX_CASES, Case_types )
    for_enum( case_flags[0][1][0], MAX_CASES, Case_types )
    for_enum( case_flags[0][1][1], MAX_CASES, Case_types )
    for_enum( case_flags[1][0][0], MAX_CASES, Case_types )
    for_enum( case_flags[1][0][1], MAX_CASES, Case_types )
    for_enum( case_flags[1][1][0], MAX_CASES, Case_types )
    for_enum( case_flags[1][1][1], MAX_CASES, Case_types )
    {
        create_case( x, y, z, case_flags,
                     &cases[x][y][z][case_flags[0][0][0]]
                                    [case_flags[0][0][1]]
                                    [case_flags[0][1][0]]
                                    [case_flags[0][1][1]]
                                    [case_flags[1][0][0]]
                                    [case_flags[1][0][1]]
                                    [case_flags[1][1][0]]
                                    [case_flags[1][1][1]] );
    }
}

static  int  get_tetra_polygon(
    Case_types   tetra_flags[],
    int          tetra_indices[][2] )
{
    static  int   sizes[5] = { 0, 3, 4, 3, 0 };
    static  int   indices[5][4][2] = {
                    { {0,0} },
                    { {0,1}, {0,2}, {0,3} },
                    { {0,2}, {0,3}, {1,3}, {1,2} },
                    { {0,3}, {1,3}, {2,3} },
                    { {0,0} }
                                  };
    int  i, n_minus;

    n_minus = 0;
    for_less( i, 0, 4 )
    {
        if( tetra_flags[i] == MINUS_FLAG )
            ++n_minus;
    }

    for_less( i, 0, sizes[n_minus] )
    {
        tetra_indices[i][0] = indices[n_minus][i][0];
        tetra_indices[i][1] = indices[n_minus][i][1];
    }

    return( sizes[n_minus] );
}

static  int  get_tetra_case(
    Case_types   case_flags[2][2][2],
    int          indices[][3],
    int          edge_indices[][2][3] )
{
    int     ord, rank, best_ord, best_rank, vertex;
    int     tetra_indices[4][2], x, y, z, p, e, size;
    Case_types   tetra_flags[4];
    static  int  orders[12][4] = {
                 { 0, 1, 2, 3 },
                 { 0, 2, 3, 1 },
                 { 0, 3, 1, 2 },
                 { 1, 2, 0, 3 },
                 { 1, 0, 3, 2 },
                 { 1, 3, 2, 0 },
                 { 2, 0, 1, 3 },
                 { 2, 1, 3, 0 },
                 { 2, 3, 0, 1 },
                 { 3, 0, 2, 1 },
                 { 3, 2, 1, 0 },
                 { 3, 1, 0, 2 }
                                 };

    best_rank = 0;

    for_less( ord, 0, 12 )
    {
        rank = 0;
        for_less( vertex, 0, 4 )
        {
            x = indices[orders[ord][vertex]][0];
            y = indices[orders[ord][vertex]][1];
            z = indices[orders[ord][vertex]][2];
            if( case_flags[x][y][z] == PLUS_FLAG )
                rank |= (1 << (3 - vertex));
        }

        if( ord == 0 || rank < best_rank )
        {
            best_rank = rank;
            best_ord = ord;
        }
    }

    for_less( vertex, 0, 4 )
    {
        x = indices[orders[best_ord][vertex]][0];
        y = indices[orders[best_ord][vertex]][1];
        z = indices[orders[best_ord][vertex]][2];
        tetra_flags[vertex] = case_flags[x][y][z];
    }

    size = get_tetra_polygon( tetra_flags, tetra_indices );

    for_less( p, 0, size )
    {
        for_less( e, 0, 2 )
        {
            vertex = orders[best_ord][tetra_indices[p][e]];

            edge_indices[p][e][0] = indices[vertex][0];
            edge_indices[p][e][1] = indices[vertex][1];
            edge_indices[p][e][2] = indices[vertex][2];
        }
    }

    return( size );
}

static  int  lookup_case(
    Case_types   case_flags[2][2][2],
    int          sizes[],
    int          edge_indices[][2][3] )
{
    int     n_polygons, ind, size, tetra;
    static  int    indices[5][4][3] = {
                      { { 0, 0, 0 },
                        { 0, 0, 1 },
                        { 0, 1, 1 },
                        { 1, 0, 1 } },
                      { { 0, 0, 0 },
                        { 1, 0, 0 },
                        { 1, 0, 1 },
                        { 1, 1, 0 } },
                      { { 0, 1, 1 },
                        { 1, 1, 1 },
                        { 1, 1, 0 },
                        { 1, 0, 1 } },
                      { { 0, 1, 0 },
                        { 0, 1, 1 },
                        { 1, 1, 0 },
                        { 0, 0, 0 } },
                      { { 0, 0, 0 },
                        { 0, 1, 1 },
                        { 1, 1, 0 },
                        { 1, 0, 1 } }
                                     };
    
    n_polygons = 0;
    ind = 0;
    for_less( tetra, 0, 5 )
    {
        size = get_tetra_case( case_flags, indices[tetra], &edge_indices[ind] );
        if( size > 0 )
        {
            sizes[n_polygons] = size;
            ++n_polygons;
            ind += size;
        }
    }

    return( n_polygons );
}

static  int   offsets[N_MARCHING_TETRA_EDGES][3] = {
                                    { 1, 0, 0 },
                                    { 0, 1, 0 },
                                    { 0, 0, 1 },
                                    { 0, 1, 1 },
                                    { 1, 0, 1 },
                                    { 1, 1, 0 },
                                    { 1, 1, 1 },
                                    { 1, 0, -1 },
                                    { 0, 1, -1 },
                                    { 1, -1, 0 },
                                    { 0, -1, -1 }
                                 };

BICAPI  void  translate_to_edge_index(
    int                x1,
    int                y1,
    int                z1,
    int                x2,
    int                y2,
    int                z2,
    voxel_point_type   *edge_point )
{
    int     tmp, edge;

    if( x2 < x1 || (x1 == x2 && y2 < y1) || (x1 == x2 && y1 == y2 && z2 < z1) )
    {
        tmp = x1;
        x1 = x2;
        x2 = tmp;
        tmp = y1;
        y1 = y2;
        y2 = tmp;
        tmp = z1;
        z1 = z2;
        z2 = tmp;
    }

    for_less( edge, 0, SIZEOF_STATIC_ARRAY(offsets) )
    {
        if( x1 + offsets[edge][X] == x2 &&
            y1 + offsets[edge][Y] == y2 &&
            z1 + offsets[edge][Z] == z2 )
            break;
    }

    if( edge >= SIZEOF_STATIC_ARRAY(offsets) )
        handle_internal_error( "edge_intersected" );

    edge_point->edge_intersected = edge;
    edge_point->coord[0] = x1;
    edge_point->coord[1] = y1;
    edge_point->coord[2] = z1;
}

BICAPI  void  translate_from_edge_index(
    int       edge_index,
    int       offset[] )
{
    offset[0] = offsets[edge_index][0];
    offset[1] = offsets[edge_index][1];
    offset[2] = offsets[edge_index][2];
}

static  void  create_case(
    int          x,
    int          y,
    int          z,
    Case_types   case_flags[2][2][2],
    case_struct  *case_info )
{
    int          i, j, k, tx, ty, tz, x1, y1, z1, x2, y2, z2;
    int          poly, n_indices, ind, vertex;
    voxel_point_type tmp;
    int          translation_indices[2][2][2][3];
    int          edge_indices[MAX_INDICES_PER_VOXEL][2][3];
    Case_types   transformed_case[2][2][2];
    int          sizes[MAX_POLYGONS_PER_VOXEL];
    BOOLEAN      left_handed;

    for_less( i, 0, 2 )
    for_less( j, 0, 2 )
    for_less( k, 0, 2 )
    {
        tx = i ^ x;
        ty = j ^ y;
        tz = k ^ z;
        translation_indices[i][j][k][0] = tx;
        translation_indices[i][j][k][1] = ty;
        translation_indices[i][j][k][2] = tz;
        transformed_case[i][j][k] = case_flags[tx][ty][tz];
    }

    case_info->n_polygons = lookup_case( transformed_case,
                                         sizes, edge_indices );
 
    if( case_info->n_polygons == 0 )
        return;

    ALLOC( case_info->sizes, case_info->n_polygons );
    n_indices = 0;
    for_less( poly, 0, case_info->n_polygons )
    {
        case_info->sizes[poly] = sizes[poly];
        n_indices += sizes[poly];
    }

    ALLOC( case_info->indices, n_indices );

    for_less( ind, 0, n_indices )
    {
        tx = edge_indices[ind][0][X];
        ty = edge_indices[ind][0][Y];
        tz = edge_indices[ind][0][Z];
        x1 = translation_indices[tx][ty][tz][X];
        y1 = translation_indices[tx][ty][tz][Y];
        z1 = translation_indices[tx][ty][tz][Z];

        tx = edge_indices[ind][1][X];
        ty = edge_indices[ind][1][Y];
        tz = edge_indices[ind][1][Z];
        x2 = translation_indices[tx][ty][tz][X];
        y2 = translation_indices[tx][ty][tz][Y];
        z2 = translation_indices[tx][ty][tz][Z];

        translate_to_edge_index( x1, y1, z1, x2, y2, z2,
                                 &case_info->indices[ind] );
    }

    left_handed = ((x + y + z) % 2) == 1;

    if( left_handed )
    {
        ind = 0;
        for_less( poly, 0, case_info->n_polygons )
        {
            for_less( vertex, 0, sizes[poly]/2 )
            {
                tmp = case_info->indices[ind+vertex];
                case_info->indices[ind+vertex] =
                             case_info->indices[ind+sizes[poly]-1-vertex];
                case_info->indices[ind+sizes[poly]-1-vertex] = tmp;
            }
            ind += sizes[poly];
        }
    }
}

static  void  delete_case(
    case_struct  *case_info );

BICAPI  int  get_tetra_isosurface_polygons(
    int               x,
    int               y,
    int               z,
    Real              corners[2][2][2],
    Real              isovalue,
    int               *sizes[],
    voxel_point_type  *points[] )
{
    int          c0, c1, c2, c3, c4, c5, c6, c7, xc, yc, zc;
    case_struct  *voxel_case;

    check_initialized();

    c0 = (corners[0][0][0] <= isovalue);
    c1 = (corners[0][0][1] <= isovalue);
    c2 = (corners[0][1][0] <= isovalue);
    c3 = (corners[0][1][1] <= isovalue);
    c4 = (corners[1][0][0] <= isovalue);
    c5 = (corners[1][0][1] <= isovalue);
    c6 = (corners[1][1][0] <= isovalue);
    c7 = (corners[1][1][1] <= isovalue);

    xc = x & 1;
    yc = y & 1;
    zc = z & 1;

    voxel_case = &cases[xc][yc][zc][c0][c1][c2][c3][c4][c5][c6][c7];

#ifdef DEBUG
{
    Case_types   case_flags[2][2][2];

    case_flags[0][0][0] = (Case_types) c0;
    case_flags[0][0][1] = (Case_types) c1;
    case_flags[0][1][0] = (Case_types) c2;
    case_flags[0][1][1] = (Case_types) c3;
    case_flags[1][0][0] = (Case_types) c4;
    case_flags[1][0][1] = (Case_types) c5;
    case_flags[1][1][0] = (Case_types) c6;
    case_flags[1][1][1] = (Case_types) c7;
    delete_case( voxel_case );
    create_case( xc, yc, zc, case_flags, voxel_case );
}
#endif

    *sizes = voxel_case->sizes;
    *points = voxel_case->indices;

    return( voxel_case->n_polygons );
}

static  void  delete_case(
    case_struct  *case_info )
{
    if( case_info->n_polygons > 0 )
    {
        FREE( case_info->sizes );
        FREE( case_info->indices );
    }
}

BICAPI  void  delete_tetra_marching_cubes_table( void )
{
    int          x, y, z;
    Case_types   c0, c1, c2, c3, c4, c5, c6, c7;

    if( initialized )
    {
        for_less( x, 0, 2 )
        for_less( y, 0, 2 )
        for_less( z, 0, 2 )
        for_enum( c0, MAX_CASES, Case_types )
        for_enum( c1, MAX_CASES, Case_types )
        for_enum( c2, MAX_CASES, Case_types )
        for_enum( c3, MAX_CASES, Case_types )
        for_enum( c4, MAX_CASES, Case_types )
        for_enum( c5, MAX_CASES, Case_types )
        for_enum( c6, MAX_CASES, Case_types )
        for_enum( c7, MAX_CASES, Case_types )
        {
            delete_case( &cases[x][y][z][c0][c1][c2][c3][c4][c5][c6][c7] );
        }
    }
}
