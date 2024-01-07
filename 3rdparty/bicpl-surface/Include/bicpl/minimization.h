#ifndef  DEF_MINIMIZATION
#define  DEF_MINIMIZATION

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


typedef struct
{
    void            *evaluation_ptr;
    double          (*evaluate_fit_function)(void);
    double          **p;
    double          *y;
    double          *psum;
    double          *ptry;
} downhill_simplex_struct;

typedef  enum  { DOWNHILL_SIMPLEX, ONE_PARAMETER_MINIMIZATION }
                 Minimization_methods;

#endif
