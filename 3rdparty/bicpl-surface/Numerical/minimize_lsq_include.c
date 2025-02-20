
#include "bicpl_internal.h"
#define   DEFAULT_RATIO      1.0
#define Real float
#define LSQ_TYPE Real

static  Real  evaluate_fit(
    int              n_parameters,
    Real             constant_term,
    LSQ_TYPE         *linear_terms,
    LSQ_TYPE         *square_terms,
    int              n_cross_terms[],
    int              *cross_parms[],
    LSQ_TYPE         *cross_terms[],
    Real             parm_values[] )
{
    int   parm, n;
    Real  fit, parm_val;

    fit = constant_term;

    for_less( parm, 0, n_parameters )
    {
        parm_val = parm_values[parm];
        fit += parm_val * ((Real) linear_terms[parm] +
                            parm_val * (Real) square_terms[parm]);

        for_less( n, 0, n_cross_terms[parm] )
            fit += parm_val * parm_values[cross_parms[parm][n]] *
                   (Real) cross_terms[parm][n];
    }

    return( fit );
}

static  void  evaluate_fit_derivative(
    int              n_parameters,
    Real             constant_term,
    LSQ_TYPE         *linear_terms,
    LSQ_TYPE         *square_terms,
    int              n_cross_terms[],
    int              *cross_parms[],
    LSQ_TYPE         *cross_terms[],
    Real             parm_values[],
    Real             derivatives[] )
{
    int        parm, n, neigh_parm, n_terms, *cp;
    Real       term, parm_val, deriv_p;
    LSQ_TYPE   *ct;

    for_less( parm, 0, n_parameters )
        derivatives[parm] = 0.0;

    for_less( parm, 0, n_parameters )
    {
        parm_val = parm_values[parm];

        deriv_p = (Real) linear_terms[parm] +
                         2.0 * parm_val * (Real) square_terms[parm];

        ct = cross_terms[parm];
        cp = cross_parms[parm];
        n_terms = n_cross_terms[parm];
        for_less( n, 0, n_terms )
        {
            neigh_parm = *cp++;
            term = (Real) (*ct++);
            deriv_p += term * parm_values[neigh_parm];
            derivatives[neigh_parm] += term * parm_val;
        }
        derivatives[parm] += deriv_p;
    }
}

static  void  evaluate_fit_along_line(
    int              n_parameters,
    Real             constant_term,
    LSQ_TYPE         *linear_terms,
    LSQ_TYPE         *square_terms,
    int              n_cross_terms[],
    int              *cross_parms[],
    LSQ_TYPE         *cross_terms[],
    Real             parm_values[],
    Real             negative_gradient[],
    Real             line_coefs[],
    Real             *a_ptr,
    Real             *b_ptr )
{
    int        parm, n, n_terms, *cp;
    Real       a, b, line_coef;
    Real       a_inc;
    LSQ_TYPE   *ct;

    a = 0.0;
    b = 0.0;

    for_less( parm, 0, n_parameters )
    {
        line_coef = line_coefs[parm];
        a_inc = line_coef * (Real) square_terms[parm];

        ct = cross_terms[parm];
        cp = cross_parms[parm];
        n_terms = n_cross_terms[parm];
        for_less( n, 0, n_terms )
            a_inc += (Real) (*ct++) * line_coefs[*cp++];

        a += line_coef * a_inc;
        b += line_coef * (-negative_gradient[parm]);
    }

    *a_ptr = a;
    *b_ptr = b;
}

static  void  minimize_along_line(
    int              n_parameters,
    Real             constant_term,
    LSQ_TYPE         *linear_terms,
    LSQ_TYPE         *square_terms,
    int              n_cross_terms[],
    int              *cross_parms[],
    LSQ_TYPE         *cross_terms[],
    Real             max_step_size,
    Real             parm_values[],
    Real             negative_gradient[],
    Real             line_coefs[] )
{
    int     parm;
    Real    a, b, t, step_size;
    static  Real     ratio;
    static  BOOLEAN  first = TRUE;

    if( first )
    {
        first = FALSE;
        if( getenv( "LSQ_STEP_RATIO" ) == 0 ||
            sscanf( getenv( "LSQ_STEP_RATIO" ), "%g", &ratio ) != 1 )
            ratio = DEFAULT_RATIO;
    }

    evaluate_fit_along_line( n_parameters, constant_term, linear_terms,
                             square_terms, n_cross_terms, cross_parms,
                             cross_terms, parm_values, negative_gradient,
                             line_coefs, &a, &b );

    if( a == 0.0 )
        return;

    t = ratio * -b / (2.0 * a);

    if( max_step_size >= 0.0 )
    {
        step_size = 0.0;
        for_less( parm, 0, n_parameters )
            step_size += t * t * line_coefs[parm] * line_coefs[parm];

        step_size = sqrt( step_size );
        if( step_size > max_step_size )
            t *= max_step_size / step_size;
    }

    for_less( parm, 0, n_parameters )
        parm_values[parm] += t * line_coefs[parm];
}

static  Real   private_minimize_lsq(
    int              n_parameters,
    Real             constant_term,
    LSQ_TYPE         *linear_terms,
    LSQ_TYPE         *square_terms,
    int              n_cross_terms[],
    int              *cross_parms[],
    LSQ_TYPE         *cross_terms[],
    Real             max_step_size,
    int              n_iters,
    Real             parm_values[] )
{
    Real              fit, len, gg, dgg, gam;
    int               iter, p;
    int               update_rate;
    Real              *unit_dir, *g, *h, *xi;
    Real              last_update_time, current_time;

    ALLOC( g, n_parameters );
    ALLOC( h, n_parameters );
    ALLOC( xi, n_parameters );
    ALLOC( unit_dir, n_parameters );

    fit = evaluate_fit( n_parameters, constant_term, linear_terms,
                        square_terms, n_cross_terms, cross_parms, cross_terms,
                        parm_values );

    print( "Initial  %g\n", fit );
    (void) flush_file( stdout );

    evaluate_fit_derivative( n_parameters, constant_term, linear_terms,
                             square_terms, n_cross_terms,
                             cross_parms, cross_terms,
                             parm_values, xi );

    for_less( p, 0, n_parameters )
    {
        g[p] = -xi[p];
        h[p] = g[p];
        xi[p] = g[p];
    }

    update_rate = 1;
    last_update_time = current_cpu_seconds();

    for_less( iter, 0, n_iters )
    {
        len = 0.0;
        for_less( p, 0, n_parameters )
            len += xi[p] * xi[p];

        if( len != 0.0 )
        {
            len = sqrt( len );
            for_less( p, 0, n_parameters )
                unit_dir[p] = xi[p] / len;

            minimize_along_line( n_parameters, constant_term, linear_terms,
                                 square_terms, n_cross_terms, cross_parms,
                                 cross_terms,
                                 max_step_size, parm_values, g, unit_dir );
        }

        if( ((iter+1) % update_rate) == 0 || iter == n_iters - 1 )
        {
            fit =  evaluate_fit( n_parameters, constant_term, linear_terms,
                                 square_terms, n_cross_terms, cross_parms,
                                 cross_terms, parm_values );

            print( "%d: %g\n", iter+1, fit );
            (void) flush_file( stdout );
            current_time = current_cpu_seconds();
            if( current_time - last_update_time < 1.0 )
                update_rate *= 10;
            last_update_time = current_time;
        }

        evaluate_fit_derivative( n_parameters, constant_term, linear_terms,
                                 square_terms, n_cross_terms,
                                 cross_parms, cross_terms,
                                parm_values, xi );

        gg = 0.0;
        dgg = 0.0;
        for_less( p, 0, n_parameters )
        {
            gg += g[p] * g[p];
            dgg += (xi[p] + g[p]) * xi[p];
/*
            dgg += xi[p] * xi[p];
*/
        }

        if( gg == 0.0 )
            break;

        gam = dgg / gg;

        for_less( p, 0, n_parameters )
        {
            g[p] = -xi[p];
            h[p] = g[p] + gam * h[p];
            xi[p] = h[p];
        }
    }

    FREE( g );
    FREE( h );
    FREE( xi );
    FREE( unit_dir );

    return( fit );
}
