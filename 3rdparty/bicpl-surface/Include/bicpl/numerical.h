#ifndef  _NUMERICAL_H
#define  _NUMERICAL_H

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


#include  <volume_io.h>
#include  <bicpl/objects.h>
#include  <bicpl/amoeba.h>
#include  <bicpl/histogram.h>
#include  <bicpl/minimization.h>
#include  <bicpl/statistics.h>

typedef  struct
{
    int   n_parameters;
    Real  **second_derivs;
    Real  *constants;
} linear_least_squares;

typedef struct
{
    int   degrees_freedom;
    int   n_points;
    Real  max_dist;
    Real  *cumulative_probs;

} t_stat_struct;

#include  <bicpl/numeric_prototypes.h>

#endif
