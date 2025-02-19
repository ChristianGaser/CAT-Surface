#include "bicpl_internal.h"

#include  "minimize_lsq_include.c"

BICAPI  Real   minimize_lsq_float(
    int              n_parameters,
    Real             constant_term,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[],
    Real             max_step_size,
    int              n_iters,
    Real             node_values[] )
{
    return( private_minimize_lsq( n_parameters, constant_term,
                                  linear_terms, square_terms,
                                  n_cross_terms, cross_parms, cross_terms,
                                  max_step_size, n_iters, node_values ) );
}

BICAPI  void  initialize_lsq_terms_float(
    int              n_parameters,
    Real             *constant_term,
    float            *linear_terms[],
    float            *square_terms[],
    int              *n_cross_terms[],
    int              **cross_parms[],
    float            **cross_terms[] )
{
    int   parm;

    ALLOC( *linear_terms, n_parameters );
    ALLOC( *square_terms, n_parameters );
    ALLOC( *n_cross_terms, n_parameters );
    ALLOC( *cross_parms, n_parameters );
    ALLOC( *cross_terms, n_parameters );

    *constant_term = 0.0;

    for_less( parm, 0, n_parameters )
    {
        (*linear_terms)[parm] = 0.0f;
        (*square_terms)[parm] = 0.0f;
        (*n_cross_terms)[parm] = 0;
    }
}

BICAPI  void  reset_lsq_terms_float(
    int              n_parameters,
    Real             *constant_term,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[] )
{
    int   parm, n;

    *constant_term = 0.0;

    for_less( parm, 0, n_parameters )
    {
        linear_terms[parm] = 0.0f;
        square_terms[parm] = 0.0f;
        for_less( n, 0, n_cross_terms[parm] )
            cross_terms[parm][n] = 0.0f;
    }
}

BICAPI  void  add_to_lsq_terms_float(
    int              n_parameters,
    Real             *constant_term,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[],
    int              n_in_list,
    int              list[],
    Real             weights[],
    Real             constant,
    int              alloc_increment )
{
    int   p, q, p1, p2, t;

    *constant_term += constant * constant;

    for_less( p, 0, n_in_list )
    {
        linear_terms[list[p]] += (float) (2.0 * weights[p] * constant);
        square_terms[list[p]] += (float) (weights[p] * weights[p]);
        for_less( q, p+1, n_in_list )
        {
            p1 = MIN( list[p], list[q] );
            p2 = MAX( list[p], list[q] );
            for_less( t, 0, n_cross_terms[p1] )
            {
                if( cross_parms[p1][t] == p2 )
                    break;
            }

            if( t >= n_cross_terms[p1] )
            {
                t = n_cross_terms[p1];
                SET_ARRAY_SIZE( cross_terms[p1], t, t+1, alloc_increment );
                SET_ARRAY_SIZE( cross_parms[p1], t, t+1, alloc_increment );
                cross_parms[p1][t] = p2;
                cross_terms[p1][t] = 0.0f;
                ++n_cross_terms[p1];
            }

            cross_terms[p1][t] += (float) (2.0 * weights[p] * weights[q]);
        }
    }
}

BICAPI  void  realloc_lsq_terms_float(
    int              n_parameters,
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[] )
{
    int   p;

    for_less( p, 0, n_parameters )
    {
        if( n_cross_terms[p] > 0 )
        {
            REALLOC( cross_terms[p], n_cross_terms[p] );
            REALLOC( cross_parms[p], n_cross_terms[p] );
        }
    }
}

BICAPI  void  delete_lsq_terms_float(
    int              n_parameters,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[] )
{
    int  p;

    for_less( p, 0, n_parameters )
    {
        if( n_cross_terms[p] > 0 )
        {
            FREE( cross_parms[p] );
            FREE( cross_terms[p] );
        }
    }

    FREE( linear_terms );
    FREE( square_terms );
    FREE( n_cross_terms );
    FREE( cross_terms );
    FREE( cross_parms );
}

BICAPI  void  create_lsq_hypersurface_float(
    STRING           filename,
    int              parm1,
    int              parm2,
    int              x_size,
    int              y_size,
    Real             x_min,
    Real             x_max,
    Real             y_min,
    Real             y_max,
    Real             scale,
    int              n_parameters,
    Real             constant,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[],
    Real             parameters[] )
{
    Real             save1, save2, p1, p2, val;
    object_struct    *object;
    quadmesh_struct  *quadmesh;
    Point            point;
    Vector           normal;
    int              x, y;

    object = create_object( QUADMESH );
    quadmesh = get_quadmesh_ptr( object );

    initialize_quadmesh( quadmesh, WHITE, NULL, x_size, y_size );

    save1 = parameters[parm1];
    save2 = parameters[parm2];

    for_less( x, 0, x_size )
    for_less( y, 0, x_size )
    {
        p1 = INTERPOLATE( (Real) x / (Real) (x_size-1), x_min, x_max );
        p2 = INTERPOLATE( (Real) y / (Real) (y_size-1), y_min, y_max );
        parameters[parm1] = p1;
        parameters[parm2] = p2;
        val = evaluate_fit( n_parameters, constant, linear_terms, square_terms,
                            n_cross_terms, cross_parms, cross_terms,
                            parameters );
        val *= scale;

        fill_Point( point, p1, p2, val );
        fill_Point( normal, 0.0, 0.0, 1.0 );
        set_quadmesh_point( quadmesh, x_size - 1 - x, y, &point, &normal );
    }

    parameters[parm1] = save1;
    parameters[parm2] = save2;

    compute_quadmesh_normals( quadmesh );

    (void) output_graphics_file( filename, BINARY_FORMAT, 1, &object );

    delete_object( object );
}
