#ifndef  DEF_QUEUE
#define  DEF_QUEUE

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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : queue.h
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Macros for maintaining queues of any element type.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */


#include  <volume_io/arrays.h>

#define  SHIFT_QUEUE_AT      100
#define  SHIFT_QUEUE_RATIO     1

/* ----------------------------- MNI Header -----------------------------------
@NAME       : QUEUE_STRUCT
@INPUT      : type
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Macro to define a queue of elements of type 'type'.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  QUEUE_STRUCT( type )                                                 \
         struct                                                               \
         {                                                                    \
             int      head, tail;                                             \
             int      n_queue_alloced;                                        \
             type     *queue;                                                 \
         }

/* ----------------------------- MNI Header -----------------------------------
@NAME       : INITIALIZE_QUEUE
@INPUT      : q
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Macro to initializes queue to empty.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  INITIALIZE_QUEUE( q )                                                \
         {                                                                    \
             (q).head = 0;                                                    \
             (q).tail = 0;                                                    \
             (q).n_queue_alloced = 0;                                         \
             (q).queue = 0;                                                   \
         }

/* ----------------------------- MNI Header -----------------------------------
@NAME       : INSERT_IN_QUEUE
@INPUT      : q
            : entry              - element to insert
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Macro to insert an entry in the queue.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  INSERT_IN_QUEUE( q, entry )                                          \
         {                                                                    \
             ADD_ELEMENT_TO_ARRAY_WITH_SIZE( (q).queue, (q).n_queue_alloced,  \
                          (q).head, entry, DEFAULT_CHUNK_SIZE )               \
         }

/* ----------------------------- MNI Header -----------------------------------
@NAME       : NUMBER_IN_QUEUE
@INPUT      : q
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Returns the number of entries in the queue.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  NUMBER_IN_QUEUE( q )  ((q).head - (q).tail)

/* ----------------------------- MNI Header -----------------------------------
@NAME       : IS_QUEUE_EMPTY
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns TRUE if the queue is empty
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  IS_QUEUE_EMPTY( q )  ((q).tail >= (q).head)

/* ----------------------------- MNI Header -----------------------------------
@NAME       : REMOVE_FROM_QUEUE
@INPUT      : q
@OUTPUT     : entry
@RETURNS    : 
@DESCRIPTION: Removes the first entry from the queue, and stores it in 'entry'.
            : Rather than shift the queue down every time, it only shifts
            : every once in a while, causing a little wasted space at the
            : beginning of the array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  REMOVE_FROM_QUEUE( q, entry )                                        \
         {                                                                    \
             if( !IS_QUEUE_EMPTY(q) )                                         \
             {                                                                \
                 entry = (q).queue[(q).tail];                                 \
                 ++(q).tail;                                                  \
                                                                              \
                 if( (q).tail == (q).head )                                   \
                 {                                                            \
                     (q).head = 0;                                            \
                     (q).tail = 0;                                            \
                 }                                                            \
                 else if( (q).tail > SHIFT_QUEUE_AT &&                        \
                          (q).tail/((q).head - (q).tail) > SHIFT_QUEUE_RATIO )\
                 {                                                            \
                     int  _i_;                                                \
                     for_less( _i_, (q).tail, (q).head )                      \
                     {                                                        \
                         (q).queue[_i_-(q).tail] = (q).queue[_i_];            \
                     }                                                        \
                     (q).head -= (q).tail;                                    \
                     (q).tail = 0;                                            \
                 }                                                            \
             }                                                                \
         }

/* ----------------------------- MNI Header -----------------------------------
@NAME       : DELETE_QUEUE
@INPUT      : q
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the queue.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#define  DELETE_QUEUE( q )                                                    \
         {                                                                    \
             if( (q).n_queue_alloced > 0 )                                    \
             {                                                                \
                 FREE( (q).queue );                                           \
             }                                                                \
         }


#endif
