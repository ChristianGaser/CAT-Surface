#ifdef  GLOBALS_LOOKUP_NAME
#undef  DEF_GLOBALS
#endif

#ifndef  DEF_GLOBALS
#define  DEF_GLOBALS

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

/*  ----------------- globals are being defined --------------------- */

#ifdef  GLOBALS_LOOKUP_NAME

#define  START_GLOBALS
#define  END_GLOBALS

#define  DEF_GLOBAL( variable, type, initial ) \
    type  variable = initial;

#define  DEF_GLOBAL3( variable, type, initial1, initial2, initial3 ) \
    type  variable = { initial1, initial2, initial3 };

#define  DEF_GLOBAL5( variable, type, initial1, initial2, initial3, initial4, initial5 ) \
    type  variable = { initial1, initial2, initial3, initial4, initial5 };

#define  DEF_GLOBAL_COLOUR( variable ) \
    Colour   variable;

#include  <global_variables.h>

#undef  START_GLOBALS
#undef  END_GLOBALS
#undef  DEF_GLOBAL
#undef  DEF_GLOBAL3
#undef  DEF_GLOBAL5
#undef  DEF_GLOBAL_COLOUR

/*  ----- database of correspondence between globals and string names ----- */

#include  <bicpl/global_lookup.h>

#define  START_GLOBALS   static  global_struct  GLOBALS_LOOKUP_NAME[] = {
#define  END_GLOBALS                                        };

#define  DEF_GLOBAL( variable, type, initial ) \
       { \
           (void *) &variable, CREATE_STRING(variable), GLUE(type,_type), FALSE\
       },

#define  DEF_GLOBAL3( variable, type, initial1, initial2, initial3 ) \
         DEF_GLOBAL( variable, type, initial1 )

#define  DEF_GLOBAL5( variable, type, i1, i2, i3, i4, i5 ) \
         DEF_GLOBAL( variable, type, i1 )

#define  DEF_GLOBAL_COLOUR( variable ) \
         DEF_GLOBAL( variable, Colour, 0 )

#include  <global_variables.h>

#undef  START_GLOBALS
#undef  END_GLOBALS
#undef  DEF_GLOBAL
#undef  DEF_GLOBAL3
#undef  DEF_GLOBAL5
#undef  DEF_GLOBAL_COLOUR

#else   /* GLOBALS_LOOKUP_NAME */

/*  ------ globals are being externed from application code ------------- */

#define  START_GLOBALS
#define  END_GLOBALS

#define  DEF_GLOBAL( variable, type, initial ) \
         extern  type  variable;

#define  DEF_GLOBAL3( variable, type, i1, i2, i3 ) \
         DEF_GLOBAL( variable, type, i1 )

#define  DEF_GLOBAL5( variable, type, i1, i2, i3, i4, i5 ) \
         DEF_GLOBAL( variable, type, i1 )

#define  DEF_GLOBAL_COLOUR( variable ) \
         extern  Colour  variable;

#include  <global_variables.h>

#undef  START_GLOBALS
#undef  END_GLOBALS
#undef  DEF_GLOBAL
#undef  DEF_GLOBAL3
#undef  DEF_GLOBAL5
#undef  DEF_GLOBAL_COLOUR

#endif   /* GLOBALS_LOOKUP_NAME */

#undef GLOBALS_LOOKUP_NAME

#endif
