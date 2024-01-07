#ifndef  DEF_BINTREE
#define  DEF_BINTREE

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


#include  <volume_io/geometry.h>

/**
 * A bintree is a binary space partition.
 * The code is in Data_structures/bintree.c.
 **/

typedef  unsigned char   node_info_type;

#define  LEFT_CHILD               0
#define  RIGHT_CHILD              1

#define  LEAF_SIGNAL              3
#define  SUBDIVISION_AXIS_BITS    3

#define  CHILDREN_BITS           12
#define  LEFT_CHILD_EXISTS       (1 << 2)
#define  RIGHT_CHILD_EXISTS      (1 << 3)

#define  NODE_INFO_OBJECTS_SHIFT  2
#define  MAX_NODE_INFO_OBJECTS   (255 >> NODE_INFO_OBJECTS_SHIFT)

typedef  struct  bintree_node_struct
{
    node_info_type   node_info;   /* --- lower two bits are split axis */
                                  /* --- 0, 1, 2, mean split on x, y, z, 3
                                         means leaf */

                                  /* --- next 2 lower bits are
                                         left and right child existence or
                                         remaining bits are number objects in 
                                         leaf */

    float            split_position;

    union
    {
    struct  bintree_node_struct    *children[1];
                          /* --- if not leaf then this is one or two
                              pointers to bintree_node_struct */
    int     object_list[1]; /* --- if leaf then this is number of objects
                                 followed by object indices list */
    }
    data;

} bintree_node_struct;

typedef  struct
{
    float    limits[N_DIMENSIONS][2];
} range_struct;

typedef  struct
{
    range_struct         range;
    int                  n_nodes;
    bintree_node_struct  *root;
} bintree_struct;

typedef  bintree_struct  *bintree_struct_ptr;

#endif
