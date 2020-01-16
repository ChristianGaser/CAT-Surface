/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   Create a surface around the white matter. This surface must be
   topologically equivalent to a sphere and will be used as the 
   starting surface for the white surface. For this to work, the
   white matter must be simply connected (all white voxels must 
   be connected together).

   Author: Claude Lepage, May 2011
*/


#include "CAT_NiftiIO.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"

#include <volume_io.h>
#include <time_stamp.h> 

#define CENTER 0x00000000

#define NODE0  0x00000001
#define NODE1  0x00000002
#define NODE2  0x00000004
#define NODE3  0x00000008
#define NODE4  0x00000010
#define NODE5  0x00000020
#define NODE6  0x00000040
#define NODE7  0x00000080

#define FACE0  0x00000100
#define FACE1  0x00000200
#define FACE2  0x00000400
#define FACE3  0x00000800
#define FACE4  0x00001000
#define FACE5  0x00002000

#define EDGE0  0x00010000
#define EDGE1  0x00020000
#define EDGE2  0x00040000
#define EDGE3  0x00080000
#define EDGE4  0x00100000
#define EDGE5  0x00200000
#define EDGE6  0x00400000
#define EDGE7  0x00800000
#define EDGE8  0x01000000
#define EDGE9  0x02000000
#define EDGE10 0x04000000
#define EDGE11 0x08000000

int NodeMask = NODE0 | NODE1 | NODE2 | NODE3 | NODE4 | NODE5 | NODE6 | NODE7;
int EdgeMask = EDGE0 | EDGE1 | EDGE2 | EDGE3 | EDGE4 | EDGE5 |
               EDGE6 | EDGE7 | EDGE8 | EDGE9 | EDGE10 | EDGE11;
int FaceMask = FACE0 | FACE1 | FACE2 | FACE3 | FACE4 | FACE5;
  

int neighbour[] = { NODE0, EDGE0, NODE1, EDGE1, FACE0, EDGE2, NODE4, EDGE3, NODE5,
                    EDGE4, FACE1, EDGE5, FACE2, CENTER, FACE3, EDGE6, FACE4, EDGE7,
                    NODE3, EDGE8, NODE2, EDGE9, FACE5, EDGE10, NODE7, EDGE11, NODE6 };


// Check if a voxel is safe to be removed.
// Edges:  EDGE0 = NODE0 + NODE1
//         EDGE1 = NODE0 + NODE4
//         EDGE2 = NODE1 + NODE5
//         EDGE3 = NODE4 + NODE5
//         EDGE4 = NODE0 + NODE3
//         EDGE5 = NODE1 + NODE2
//         EDGE6 = NODE4 + NODE7
//         EDGE7 = NODE5 + NODE6
//         EDGE8 = NODE2 + NODE3
//         EDGE9 = NODE3 + NODE7
//         EDGE10 = NODE2 + NODE6
//         EDGE11 = NODE6 + NODE7
// Faces:  FACE0 = EDGE0 + EDGE1 + EDGE2 + EDGE3  (opp FACE5)
//         FACE1 = EDGE0 + EDGE4 + EDGE5 + EDGE8  (opp FACE4)
//         FACE2 = EDGE1 + EDGE4 + EDGE6 + EDGE9  (opp FACE3)
//         FACE3 = EDGE2 + EDGE5 + EDGE7 + EDGE10  (opp FACE2)
//         FACE4 = EDGE3 + EDGE6 + EDGE7 + EDGE11  (opp FACE1)
//         FACE5 = EDGE8 + EDGE9 + EDGE10 + EDGE11  (opp FACE0)
//

int check_voxel( int flag ) {

  int  count = 0;
  int  node_mask = 0;
  int  edge_mask = 0;

  if( flag & FaceMask ) {
    if( flag & FACE0 ) {
      node_mask |= ( NODE0 | NODE1 | NODE4 | NODE5 );
      edge_mask |= ( EDGE0 | EDGE1 | EDGE2 | EDGE3 );
      count++;
    }
    if( flag & FACE1 ) {
      node_mask |= ( NODE0 | NODE1 | NODE2 | NODE3 );
      edge_mask |= ( EDGE0 | EDGE4 | EDGE5 | EDGE8 );
      count++;
    }
    if( flag & FACE2 ) {
      node_mask |= ( NODE0 | NODE3 | NODE4 | NODE7 );
      edge_mask |= ( EDGE1 | EDGE4 | EDGE6 | EDGE9 );
      count++;
    }
    if( flag & FACE3 ) {
      node_mask |= ( NODE1 | NODE2 | NODE5 | NODE6 );
      edge_mask |= ( EDGE2 | EDGE5 | EDGE7 | EDGE10 );
      count++;
    }
    if( flag & FACE4 ) {
      node_mask |= ( NODE4 | NODE5 | NODE6 | NODE7 );
      edge_mask |= ( EDGE3 | EDGE6 | EDGE7 | EDGE11 );
      count++;
    }
    if( flag & FACE5 ) {
      node_mask |= ( NODE2 | NODE3 | NODE6 | NODE7 );
      edge_mask |= ( EDGE8 | EDGE9 | EDGE10 | EDGE11 );
      count++;
    }

    // entities not touching found faces, otherwise a cut may
    // be made in the surface.
    node_mask = NodeMask ^ node_mask;
    edge_mask = EdgeMask ^ edge_mask;
    if( flag & node_mask ) return 0;
    if( flag & edge_mask ) return 0;

    // Special case for opposite faces: do not make a hole nor
    // break up a thin structure.
    if( count == 2 ) {
      if( flag & FACE0 && flag & FACE5 ) return 0;
      if( flag & FACE1 && flag & FACE4 ) return 0;
      if( flag & FACE2 && flag & FACE3 ) return 0;
    }
    if( count == 4 ) {
      if( !( flag & FACE0 ) && !( flag & FACE5 ) ) return 0;
      if( !( flag & FACE1 ) && !( flag & FACE4 ) ) return 0;
      if( !( flag & FACE2 ) && !( flag & FACE3 ) ) return 0;
    }
    return 1;

  } else {
    return 0;  // it's possible to have no face, don't delete, wait till later.
  }
  
  return 1;
}

//
//    val = 0 : voxel is outside sphere
//    val = 1 : voxel is inside sphere, but outside white matter surface
//    val > 1 : voxel is white matter and inside white matter surface
//
int make_sphere( int sizes[MAX_DIMENSIONS], short * val ) {

    int     i, j, k, ii, jj, iii;
    int     changed;
    int     di, dj, dk;

    int     n_voxels = sizes[0] * sizes[1] * sizes[2];

    int stencil = 1;
    int level_count = 0;

    // Create connectivity levels around the white surface until
    // all of the inside of the sphere (val>=1) has been labelled.

    short *levels  = (short *) malloc(sizeof(short) * n_voxels);
    if( !levels ) {
      printf( "Error allocating memory for levels.\n" );
      exit( 1 );
    }

    short *levels2  = (short *) malloc(sizeof(short) * n_voxels);
    if( !levels2 ) {
      printf( "Error allocating memory for levels2.\n" );
      exit( 1 );
    }

    int *level_ptr  = (int *) malloc(sizeof(int) * n_voxels);
    if( !level_ptr ) {
      printf( "Error allocating memory for level_ptr.\n" );
      exit( 1 );
    }
    
    int*level_ptr_tmp  = (int *) malloc(sizeof(int) * n_voxels/4);
    if( !level_ptr_tmp ) {
      printf( "Error allocating memory for level_ptr_tmp.\n" );
      exit( 1 );
    }

    level_count = 0;
    for( ii = 0; ii < n_voxels; ii++ ) {
      // ensure that border layer of voxels is outside
      i = ii / ( sizes[1] * sizes[2] );
      j = ( ii - i * sizes[1] * sizes[2] ) / sizes[2];
      k = ii - ( i * sizes[1] + j ) * sizes[2];
      if( i == 0 || i == sizes[0]-1 ||
          j == 0 || j == sizes[1]-1 ||
          k == 0 || k == sizes[2]-1 ) {
        val[ii] = 0;
      }
      if( val[ii] > 1 ) {
        levels[ii] = 0;    // white
        level_ptr[level_count] = ii;
        level_count++;
      } else {
        levels[ii] = -1;   // outside white (and outside sphere)
      }
    }

    int nlevel = 0;
    do {
      changed = 0;
      for( iii = 0; iii < level_count; iii++ ) {
        ii = level_ptr[iii];
        i = ii / ( sizes[1] * sizes[2] );
        j = ( ii - i * sizes[1] * sizes[2] ) / sizes[2];
        k = ii - ( i * sizes[1] + j ) * sizes[2];
        for( di = -1; di <= 1; di++ ) {
          for( dj = -1; dj <= 1; dj++ ) {
            for( dk = -1; dk <= 1; dk++ ) {
              if( ABS(di) + ABS(dj) + ABS(dk) <= stencil ) {
                jj = ( (i+di) * sizes[1] + (j+dj) ) * sizes[2] + k+dk;
                if( levels[jj] == -1 && val[jj] == 1 ) {
                  levels[jj] = nlevel+1;
                  level_ptr_tmp[changed] = jj;
                  changed++;
                }
              }
            }
          }
        }
      }
      nlevel++;
      level_count = changed;
      for( iii = 0; iii < level_count; iii++ ) {
        level_ptr[iii] = level_ptr_tmp[iii];
      }

      if( changed > 0 ) printf( "Outer level %d with %d voxels\n", 
                                nlevel, changed );
    } while( changed > 0 );

    // At this point, the outer surface of the labelled voxels
    // should be topologically equivalent to a sphere, after
    // propagation from white to outside.
    // Now propagate nlevel from outside to inside.

    level_count = 0;
    for( ii = 0; ii < n_voxels; ii++ ) {
      if( val[ii] == 0 ) {
        // These are all the voxels outside the sphere, so some
        // of them will touch the volume border. So in the loop
        // below, check for borders to avoid seg fault.
        levels2[ii] = 0;    // outside sphere
        level_ptr[level_count] = ii;
        level_count++;
      } else {
        levels2[ii] = -1;   // inside sphere (and inside white)
      }
    }

    stencil = 1;
    short nlevel2 = 0;
    do {
      changed = 0;
      for( iii = 0; iii < level_count; iii++ ) {
        ii = level_ptr[iii];
        i = ii / ( sizes[1] * sizes[2] );
        j = ( ii - i * sizes[1] * sizes[2] ) / sizes[2];
        k = ii - ( i * sizes[1] + j ) * sizes[2];
        for( di = -1; di <= 1; di++ ) {
          for( dj = -1; dj <= 1; dj++ ) {
            for( dk = -1; dk <= 1; dk++ ) {
              if( ABS(di) + ABS(dj) + ABS(dk) <= stencil ) {
                if( i+di >= 0 && i+di < sizes[0] &&
                    j+dj >= 0 && j+dj < sizes[1] &&
                    k+dk >= 0 && k+dk < sizes[2] ) {
                  jj = ( (i+di) * sizes[1] + (j+dj) ) * sizes[2] + k+dk;
                  if( levels2[jj] == -1 && val[jj] == 1 ) {
                    levels2[jj] = nlevel2+1;
                    level_ptr_tmp[changed] = jj;
                    changed++;
                  }
                }
              }
            }
          }
        }
      }

      nlevel2++;
      level_count = changed;
      for( iii = 0; iii < level_count; iii++ ) {
        level_ptr[iii] = level_ptr_tmp[iii];
      }

      if( changed > 0 ) printf( "Inner level %d with %d voxels\n", 
                                nlevel2, changed );
    } while( changed > 0 );

    // Blend levels and levels2 to make a better iteration pattern
    // for removal of the voxels.

    short W1 = 3;
    short W2 = 1;
    for( ii = 0; ii < n_voxels; ii++ ) {
      levels[ii] = W1 * levels[ii] + W2 * (  nlevel2 - levels2[ii] + 1 );
    }
    nlevel = W1 * nlevel + W2 * nlevel2;

    // Do a bit of blurring on the levels. This helps to prevent bridges.

    stencil = 1;
    for( iii = 0; iii < 10; iii++ ) {
      printf( "blurring iteration %d\n", iii );
      for( ii = 0; ii < n_voxels; ii++ ) {
        levels2[ii] = levels[ii];
        if( val[ii] == 1 ) {  // ignore borders
          i = ii / ( sizes[1] * sizes[2] );
          j = ( ii - i * sizes[1] * sizes[2] ) / sizes[2];
          k = ii - ( i * sizes[1] + j ) * sizes[2];
          int count = 0;
          int sum = 0;
          for( di = -1; di <= 1; di++ ) {
            for( dj = -1; dj <= 1; dj++ ) {
              for( dk = -1; dk <= 1; dk++ ) {
                if( ABS(di) + ABS(dj) + ABS(dk) <= stencil ) {
                  jj = ( (i+di) * sizes[1] + (j+dj) ) * sizes[2] + k+dk;
                  if( val[jj] == 1 ) {
                    sum += levels[jj];
                    count++;
                  }
                }
              }
	    }
	  }
          if( count ) {
            levels2[ii] = (short)rint( (float)sum / count );
          }   
        }
      }
      for( ii = 0; ii < n_voxels; ii++ ) {
        levels[ii] = levels2[ii];
      }
    }
    free(levels2);

    // Initialize voxels we can eliminate on this iteration. Proceed
    // by levels.

    int level_changed;
    stencil = 3;

    level_count = 0;
    for( ii = 0; ii < n_voxels; ii++ ) {
      if( levels[ii] == nlevel && val[ii] == 1 ) {
        levels2[level_count] = ii;
        level_count++;
      }
    }

    do {
      // Find current level around outer boundary.

      level_changed = 0;

      printf( "Processing level %d (%d):", nlevel, level_count );

      do {
        changed = 0;
      
        for( iii = 0; iii < level_count; iii++ ) {

          ii = level_ptr[iii];

          if( val[ii] == 1 ) {
            int flag = 0;
            i = ii / ( sizes[1] * sizes[2] );
            j = ( ii - i * sizes[1] * sizes[2] ) / sizes[2];
            k = ii - ( i * sizes[1] + j ) * sizes[2];
  
            short idx = 0;
  
            for( di = -1; di <= 1; di++ ) {
              for( dj = -1; dj <= 1; dj++ ) {
                for( dk = -1; dk <= 1; dk++ ) {
                  if( ABS(di) + ABS(dj) + ABS(dk) <= stencil ) {
                    jj = ( (i+di) * sizes[1] + (j+dj) ) * sizes[2] + k+dk;
                    if( val[jj] == 0 ) {
                      flag |= neighbour[idx];
                    }
                  }
                  idx++;
                }
              }
            }
            if( check_voxel( flag ) ) {
              val[ii] = 0;
              changed++;
            }
          }
        }
        if( changed > 0 ) printf( " %d", changed );
        level_changed += changed;
      } while( changed );

      // All unassigned voxels at this level will go to the next level.
      int unassigned = 0;
      for( iii = 0; iii < level_count; iii++ ) {
        ii = level_ptr[iii];
        if( val[ii] == 1 ) {
          unassigned++;
          levels[ii] = nlevel-1;
        }
      }
      printf( " (%d)\n", unassigned );
      nlevel--;
      level_count = 0;
      for( ii = 0; ii < n_voxels; ii++ ) {
        if( levels[ii] == nlevel && val[ii] == 1 ) {
          level_ptr[level_count] = ii;
          level_count++;
        }
      }

    // } while( level_changed || nlevel >= W1 + W2 );
    } while( level_changed || nlevel >= 1 );

    changed = 0;
    for( ii = 0; ii < n_voxels; ii++ ) {
      if( val[ii] == 1 ) changed++;
    }
    printf( "%d non-white voxels inside surface\n", changed );

    free(levels);
    free(level_ptr);
    free(level_ptr_tmp);

    return( OK );
}


// Save surface from voxels.

void save_voxel_surface( int sizes[3], Real dh[3], short * val,
                         Volume volume, char * obj_file ) {

    int   i, j, k, ii, jj;

    int n_nodes = (sizes[0]+1)*(sizes[1]+1)*(sizes[2]+1);
    int *nodes  = (int *) malloc(sizeof(int) * n_nodes);
    if( !nodes ) {
      printf( "Error allocating memory for nodes.\n" );
      exit( 1 );
    }
    for( i = 0; i < n_nodes; i++ ) {
      nodes[i] = -1;
    }

    FILE * fp = fopen( obj_file, "w" );

    int count = 0;
    int num_bndy = 0;
    int free_voxel = 0;

    short * old_im1, * old_i, * old_ip1;

    for( i = 0;  i < sizes[0];  ++i ) {
      // define pointers to the active planes.
      if( i > 0 ) {
        old_im1 = &(val[(i-1)*sizes[1]*sizes[2]]);
      }
      old_i = &(val[i*sizes[1]*sizes[2]]);
      if( i+1 < sizes[0] ) {
        old_ip1 = &(val[(i+1)*sizes[1]*sizes[2]]);
      }
      // loop on (y,z) in the x plane.
      for( j = 0;  j < sizes[1];  ++j ) {
        for( k = 0;  k < sizes[2];  ++k ) {

          if( old_i[j*sizes[2]+k] == 0 ) continue;

          // Count how many boundaries for this voxel.
          int  nn = 0;
          if( old_ip1[j*sizes[2]+k] == 0 ) nn++;
          if( old_im1[j*sizes[2]+k] == 0 ) nn++;
          if( old_i[(j+1)*sizes[2]+k] == 0 ) nn++;
          if( old_i[(j-1)*sizes[2]+k] == 0 ) nn++;
          if( old_i[j*sizes[2]+k+1] == 0 ) nn++;
          if( old_i[j*sizes[2]+k-1] == 0 ) nn++;
          if( nn == 6 ) continue;   // isolated voxel (should not occur now)
          if( nn == 0 ) continue;   // all surrounded by color - no boundary

          // Find the 8 nodes of the voxel.

          int  n0, n1, n2, n3, n4, n5, n6, n7;

          n0 = i*(sizes[1]+1)*(sizes[2]+1) + j*(sizes[2]+1) + k;
          n1 = (i+1)*(sizes[1]+1)*(sizes[2]+1) + j*(sizes[2]+1) + k;
          n2 = (i+1)*(sizes[1]+1)*(sizes[2]+1) + (j+1)*(sizes[2]+1) + k;
          n3 = i*(sizes[1]+1)*(sizes[2]+1) + (j+1)*(sizes[2]+1) + k;
          n4 = n0+1;
          n5 = n1+1;
          n6 = n2+1;
          n7 = n3+1;

          // Check the 6 face neighbours.

          if( old_ip1[j*sizes[2]+k] == 0 ) {
            // face (n1,n2,n6,n5):
            if( nodes[n1] == -1 ) {
              nodes[n1] = count;
              count++;
            }
            if( nodes[n2] == -1 ) {
              nodes[n2] = count;
              count++;
            }
            if( nodes[n6] == -1 ) {
              nodes[n6] = count;
              count++;
            }
            if( nodes[n5] == -1 ) {
              nodes[n5] = count;
              count++;
            }
            fprintf( fp, "%d %d %d %d %d %d\n", nodes[n2], nodes[n1], nodes[n6],
                                                nodes[n6], nodes[n1], nodes[n5] );
            num_bndy++;
          }

          if( old_im1[j*sizes[2]+k] == 0 ) {
            // face (n0,n4,n7,n3):
            if( nodes[n0] == -1 ) {
              nodes[n0] = count;
              count++;
            }
            if( nodes[n4] == -1 ) {
              nodes[n4] = count;
              count++;
            }
            if( nodes[n7] == -1 ) {
              nodes[n7] = count;
              count++;
            }
            if( nodes[n3] == -1 ) {
              nodes[n3] = count;
              count++;
            }
            fprintf( fp, "%d %d %d %d %d %d\n", nodes[n4], nodes[n0], nodes[n3],
                                                nodes[n7], nodes[n4], nodes[n3] );
            num_bndy++;
          }

          if( old_i[j*sizes[2]+k+1] == 0 ) {
            // face (n4,n5,n6,n7):
            if( nodes[n4] == -1 ) {
              nodes[n4] = count;
              count++;
            }
            if( nodes[n5] == -1 ) {
              nodes[n5] = count;
              count++;
            }
            if( nodes[n6] == -1 ) {
              nodes[n6] = count;
              count++;
            }
            if( nodes[n7] == -1 ) {
              nodes[n7] = count;
              count++;
            }
            fprintf( fp, "%d %d %d %d %d %d\n", nodes[n5], nodes[n4], nodes[n6],
                                                nodes[n6], nodes[n4], nodes[n7] );
            num_bndy++;
          }

          if( old_i[j*sizes[2]+k-1] == 0 ) {
            // face (n0,n3,n2,n1):
            if( nodes[n0] == -1 ) {
              nodes[n0] = count;
              count++;
            }
            if( nodes[n3] == -1 ) {
              nodes[n3] = count;
              count++;
            }
            if( nodes[n2] == -1 ) {
              nodes[n2] = count;
              count++;
            }
            if( nodes[n1] == -1 ) {
              nodes[n1] = count;
              count++;
            }
            fprintf( fp, "%d %d %d %d %d %d\n", nodes[n3], nodes[n0], nodes[n1],
                                                nodes[n3], nodes[n1], nodes[n2] );
            num_bndy++;
          }

          if( old_i[(j+1)*sizes[2]+k] == 0 ) {
            // face (n3,n7,n6,n2):
            if( nodes[n3] == -1 ) {
              nodes[n3] = count;
              count++;
            }
            if( nodes[n2] == -1 ) {
              nodes[n2] = count;
              count++;
            }
            if( nodes[n6] == -1 ) {
              nodes[n6] = count;
              count++;
            }
            if( nodes[n7] == -1 ) {
              nodes[n7] = count;
              count++;
            }
            fprintf( fp, "%d %d %d %d %d %d\n", nodes[n7], nodes[n3], nodes[n6],
                                                nodes[n6], nodes[n3], nodes[n2] );
            num_bndy++;
          }

          if( old_i[(j-1)*sizes[2]+k] == 0 ) {
            // face (n0,n1,n5,n4):
            if( nodes[n0] == -1 ) {
              nodes[n0] = count;
              count++;
            }
            if( nodes[n1] == -1 ) {
              nodes[n1] = count;
              count++;
            }
            if( nodes[n5] == -1 ) {
              nodes[n5] = count;
              count++;
            }
            if( nodes[n4] == -1 ) {
              nodes[n4] = count;
              count++;
            }
            fprintf( fp, "%d %d %d %d %d %d\n", nodes[n1], nodes[n0], nodes[n5],
                                                nodes[n5], nodes[n0], nodes[n4] );
            num_bndy++;
          }
        }
      }
    }
    fclose( fp );

    num_bndy *= 2;
    printf( "Found %d boundary faces\n", num_bndy );
    printf( "Found %d boundary nodes\n", count );

    // Renumbering of the boundary nodes and coordinates.

    int *index  = (int *) malloc(sizeof(int) * count);
    if( !index ) {
      printf( "Error allocating memory for index.\n" );
      exit( 1 );
    }

    for( i = 0; i < n_nodes; i++ ) {
      if( nodes[i] >= 0 ) {
        index[nodes[i]] = i;
      }
    }

    free(nodes);

    // Read the connectivity in exact memory.

    fp = fopen( obj_file, "r" );
    int *connec  = (int *) malloc(sizeof(int) * 3*num_bndy);
    if( !connec ) {
      printf( "Error allocating memory for connec.\n" );
      exit( 1 );
    }
    for( ii = 0; ii < num_bndy/2; ii++ ) {
      fscanf( fp, "%d %d %d %d %d %d\n", &connec[6*ii], &connec[6*ii+1],
              &connec[6*ii+2], &connec[6*ii+3], &connec[6*ii+4], &connec[6*ii+5] );
    }
    fclose( fp );
    remove( obj_file );

    // Remove points with 3 triangle neighbours, to simplify connectivity.
    int *n_ngh  = (int *) malloc(sizeof(int) * count);
    if( !n_ngh ) {
      printf( "Error allocating memory for n_ngh.\n" );
      exit( 1 );
    }
    int *ngh  = (int *) malloc(sizeof(int) * 3*count);
    if( !ngh ) {
      printf( "Error allocating memory for ngh.\n" );
      exit( 1 );
    }
    for( i = 0; i < count; i++ ) {
      n_ngh[i] = 0;
    }
    for( ii = 0; ii < num_bndy; ii++ ) {
      int n0 = connec[3*ii];
      int n1 = connec[3*ii+1];
      int n2 = connec[3*ii+2];
      if( n_ngh[n0] < 3 ) ngh[3*n0+n_ngh[n0]] = ii;
      if( n_ngh[n1] < 3 ) ngh[3*n1+n_ngh[n1]] = ii;
      if( n_ngh[n2] < 3 ) ngh[3*n2+n_ngh[n2]] = ii;
      n_ngh[n0]++;
      n_ngh[n1]++;
      n_ngh[n2]++;
    }

    int n0, n1, n2, idx[3];

    for( i = 0; i < count; i++ ) {
      if( n_ngh[i] == 3 ) {
        idx[0] = ngh[3*i];
        idx[1] = ngh[3*i+1];
        idx[2] = ngh[3*i+2];
        index[i] = -1;
        n0 = connec[3*idx[0]];
        n1 = connec[3*idx[0]+1];
        n2 = connec[3*idx[0]+2];
        for( j = 0; j < 3; j++ ) {
          ii = connec[3*idx[1]+j];
          connec[3*idx[1]+j] = -1;
          if( ii != n0 && ii != n1 && ii != n2 ) {
            if( n0 == i ) connec[3*idx[0]] = ii;
            if( n1 == i ) connec[3*idx[0]+1] = ii;
            if( n2 == i ) connec[3*idx[0]+2] = ii;
          }
          ii = connec[3*idx[2]+j];
          connec[3*idx[2]+j] = -1;
          if( ii != n0 && ii != n1 && ii != n2 ) {
            if( n0 == i ) connec[3*idx[0]] = ii;
            if( n1 == i ) connec[3*idx[0]+1] = ii;
            if( n2 == i ) connec[3*idx[0]+2] = ii;
          }
        }
      }
    }
    free(ngh);
    free(n_ngh);

    // Renumber the nodes after removal of some triangles.

    int * renum  = (int *) malloc(sizeof(int) * count);
    if( !renum ) {
      printf( "Error allocating memory for renum.\n" );
      exit( 1 );
    }

    jj = 0; 
    for( i = 0; i < count; i++ ) {
      renum[i] = -1;
      if( index[i] != -1 ) {
        renum[i] = jj;
        index[jj] = index[i];
        jj++;
      }
    }
    count = jj;
    printf( "Now have %d nodes\n", count );

    // Concatenate connectivity after removal of some triangles.

    jj = 0; 
    for( ii = 0; ii < num_bndy; ii++ ) {
      if( connec[3*ii] != -1 ) {
        connec[3*jj] = renum[connec[3*ii]];
        connec[3*jj+1] = renum[connec[3*ii+1]];
        connec[3*jj+2] = renum[connec[3*ii+2]];
        jj++;
      }
    }
    num_bndy = jj;
    printf( "Now have %d triangles\n", num_bndy );

    free(renum);

    fp = fopen( obj_file, "w" );
    fprintf( fp, "P 0.3 0.3 0.4 10 1 %d\n", count );

    Real       voxel[3], voxel_000[3], voxel_100[3], voxel_010[3], voxel_001[3];

    convert_world_to_voxel( volume, 0.0, 0.0, 0.0, voxel_000 );
    convert_world_to_voxel( volume, 1.0, 0.0, 0.0, voxel_100 );
    convert_world_to_voxel( volume, 0.0, 1.0, 0.0, voxel_010 );
    convert_world_to_voxel( volume, 0.0, 0.0, 1.0, voxel_001 );

    printf( "origin = %g %g %g\n", voxel_000[0], voxel_000[1], voxel_000[2] );
    printf( "x-axis = %g %g %g\n", voxel_100[0], voxel_100[1], voxel_100[2] );
    printf( "y-axis = %g %g %g\n", voxel_010[0], voxel_010[1], voxel_010[2] );
    printf( "z-axis = %g %g %g\n", voxel_001[0], voxel_001[1], voxel_001[2] );

    // orientation of axes for (x,y,z)
    int XX, YY, ZZ;
    if( voxel_100[0] != voxel_000[0] ) XX = 0;
    if( voxel_100[1] != voxel_000[1] ) XX = 1;
    if( voxel_100[2] != voxel_000[2] ) XX = 2;
    if( voxel_010[0] != voxel_000[0] ) YY = 0;
    if( voxel_010[1] != voxel_000[1] ) YY = 1;
    if( voxel_010[2] != voxel_000[2] ) YY = 2;
    if( voxel_001[0] != voxel_000[0] ) ZZ = 0;
    if( voxel_001[1] != voxel_000[1] ) ZZ = 1;
    if( voxel_001[2] != voxel_000[2] ) ZZ = 2;

    // offset from center to corner of voxel
    // coord = (i - origin - 1/2) * dh
    //       = (i - (origin + 1/2) * dh
    // So add dh/2 to origin (not subtract)
    voxel_000[0] += 0.5;
    voxel_000[1] += 0.5;
    voxel_000[2] += 0.5;

    // Determine the coordinates.
    float * coords  = (float *) malloc(sizeof(float) * 3 * count);
    if( !coords ) {
      printf( "Error allocating memory for coords.\n" );
      exit( 1 );
    }
    for( ii = 0; ii < count; ii++ ) {
      i = index[ii] / ( ( sizes[1]+1 ) * ( sizes[2]+1 ) );
      j = ( index[ii] - i * ( sizes[1]+1 ) * ( sizes[2]+1 ) ) / ( sizes[2]+1 );
      k = index[ii] - i * ( sizes[1]+1 ) * ( sizes[2]+1 ) - j * ( sizes[2]+1 );
      coords[3*ii+XX] = ( (float)i - voxel_000[0] ) * dh[0];
      coords[3*ii+YY] = ( (float)j - voxel_000[1] ) * dh[1];
      coords[3*ii+ZZ] = ( (float)k - voxel_000[2] ) * dh[2];
    }

    // Do a little bit of smoothing on the coordinates (simple averaging).

    float * new_coords  = (float *) malloc(sizeof(float) * 3 * count);
    if( !new_coords ) {
      printf( "Error allocating memory for new_coords.\n" );
      exit( 1 );
    }

    for( ii = 0; ii < count; ii++ ) {
      index[ii] = 0;
    }
    for( ii = 0; ii < num_bndy; ii++ ) {
      int n0 = connec[3*ii];
      int n1 = connec[3*ii+1];
      int n2 = connec[3*ii+2];
      index[n0]++;
      index[n1]++;
      index[n2]++;
    }

    float relax = 0.80;

    // for( int iter = 0; iter < 50; iter++ ) {
    for( int iter = 0; iter < 20; iter++ ) {
      if( (iter%10) == 0 ) printf( "Smoothing iteration %d...\n", iter );
      for( ii = 0; ii < 3*count; ii++ ) {
        new_coords[ii] = 0.0;
      }
      for( ii = 0; ii < num_bndy; ii++ ) {
        int n0 = connec[3*ii];
        int n1 = connec[3*ii+1];
        int n2 = connec[3*ii+2];
        for( j = 0; j < 3; j++ ) {
          float xc = coords[3*n0+j] + coords[3*n1+j] + coords[3*n2+j];
          new_coords[3*n0+j] += xc;
          new_coords[3*n1+j] += xc;
          new_coords[3*n2+j] += xc;
        }
      }
      for( ii = 0; ii < 3*count; ii++ ) {
        new_coords[ii] = 0.5 * ( new_coords[ii] / (float)(index[ii/3]) - coords[ii] );
        coords[ii] = relax * coords[ii] + ( 1.0 - relax ) * new_coords[ii];
      }
    }
    free(new_coords);
    free(index);

    // Write out the coordinates.
    for( ii = 0; ii < count; ii++ ) {
      fprintf( fp, "%g %g %g\n", coords[3*ii+0], coords[3*ii+1],
               coords[3*ii+2] );
    }

    // The normals.
    float * normals  = (float *) malloc(sizeof(float) * 3 * count);
    if( !normals ) {
      printf( "Error allocating memory for normals.\n" );
      exit( 1 );
    }
    for( ii = 0; ii < 3*count; ii++ ) {
      normals[ii] = 0.0;
    }
    for( ii = 0; ii < num_bndy; ii++ ) {
      int n0 = connec[3*ii];
      int n1 = connec[3*ii+1];
      int n2 = connec[3*ii+2];

      float dx1 = coords[3*n1] - coords[3*n0];
      float dy1 = coords[3*n1+1] - coords[3*n0+1];
      float dz1 = coords[3*n1+2] - coords[3*n0+2];

      float dx2 = coords[3*n2] - coords[3*n0];
      float dy2 = coords[3*n2+1] - coords[3*n0+1];
      float dz2 = coords[3*n2+2] - coords[3*n0+2];

      float nx = dy1 * dz2 - dy2 * dz1;
      float ny = dz1 * dx2 - dz2 * dx1;
      float nz = dx1 * dy2 - dx2 * dy1;
      float mag = sqrt( nx*nx + ny*ny + nz*nz );
      if( mag <= 1.0e-20 ) mag = 1.0e10;
      nx /= mag;
      ny /= mag;
      nz /= mag;

      normals[3*n0] += nx;
      normals[3*n0+1] += ny;
      normals[3*n0+2] += nz;

      normals[3*n1] += nx;
      normals[3*n1+1] += ny;
      normals[3*n1+2] += nz;

      normals[3*n2] += nx;
      normals[3*n2+1] += ny;
      normals[3*n2+2] += nz;
    }
    free(coords);

    fprintf( fp, "\n" );
    for( ii = 0; ii < count; ii++ ) {
      float mag = sqrt( normals[3*ii] * normals[3*ii] +
                        normals[3*ii+1] * normals[3*ii+1] +
                        normals[3*ii+2] * normals[3*ii+2] );
      if( mag <= 1.0e-20 ) mag = 1.0e10;
      fprintf( fp, "%g %g %g\n", normals[3*ii]/mag, normals[3*ii+1]/mag,
               normals[3*ii+2]/mag );
    }
    free(normals);

    // The connectivity - part 1.
    fprintf( fp, "\n" );
    fprintf( fp, "%d\n", num_bndy );
    fprintf( fp, "0 1 1 1 1\n\n" );

    for( ii = 1; ii <= num_bndy; ii++ ) {
      fprintf( fp, "%d ", 3*ii );
      if( ii%8 == 0 ) fprintf( fp, "\n" );
    }

    // The connectivity - part 2.

    for( ii = 0; ii < 3*num_bndy; ii++ ) {
      if( ii%8 == 0 ) fprintf( fp, "\n" );
      fprintf( fp, "%d ", connec[ii] );
    }
    fclose( fp );
    free(connec);
}


void save_voxel_surface2( int sizes[3], Real dh[3], short * val,
                         Volume volume, char * obj_file ) {

    int   i, j, k, ii, jj;
    object_struct      *object;
    polygons_struct		 *polygons;
    Surfprop           spr;

/*    object  = create_object( POLYGONS );
    polygons = get_polygons_ptr( object );

    Surfprop_a(spr) = 0.3f;
    Surfprop_d(spr) = 0.6f;
    Surfprop_s(spr) = 0.6f;
    Surfprop_se(spr) = 30.0f;
    Surfprop_t(spr) = 1.0f;
    initialize_polygons( polygons, WHITE, &spr );
*/
    int n_nodes = (sizes[0]+1)*(sizes[1]+1)*(sizes[2]+1);
    int *nodes  = (int *) malloc(sizeof(int) * n_nodes);
    if( !nodes ) {
      printf( "Error allocating memory for nodes.\n" );
      exit( 1 );
    }
    for( i = 0; i < n_nodes; i++ ) {
      nodes[i] = -1;
    }

    int count = 0;
    int num_bndy = 0;
    int free_voxel = 0;

    short * old_im1, * old_i, * old_ip1;

    for( i = 0;  i < sizes[0];  ++i ) {
      // define pointers to the active planes.
      if( i > 0 ) {
        old_im1 = &(val[(i-1)*sizes[1]*sizes[2]]);
      }
      old_i = &(val[i*sizes[1]*sizes[2]]);
      if( i+1 < sizes[0] ) {
        old_ip1 = &(val[(i+1)*sizes[1]*sizes[2]]);
      }
      // loop on (y,z) in the x plane.
      for( j = 0;  j < sizes[1];  ++j ) {
        for( k = 0;  k < sizes[2];  ++k ) {

          if( old_i[j*sizes[2]+k] == 0 ) continue;

          // Count how many boundaries for this voxel.
          int  nn = 0;
          if( old_ip1[j*sizes[2]+k] == 0 ) nn++;
          if( old_im1[j*sizes[2]+k] == 0 ) nn++;
          if( old_i[(j+1)*sizes[2]+k] == 0 ) nn++;
          if( old_i[(j-1)*sizes[2]+k] == 0 ) nn++;
          if( old_i[j*sizes[2]+k+1] == 0 ) nn++;
          if( old_i[j*sizes[2]+k-1] == 0 ) nn++;
          if( nn == 6 ) continue;   // isolated voxel (should not occur now)
          if( nn == 0 ) continue;   // all surrounded by color - no boundary

          // Find the 8 nodes of the voxel.

          int  n0, n1, n2, n3, n4, n5, n6, n7;

          n0 = i*(sizes[1]+1)*(sizes[2]+1) + j*(sizes[2]+1) + k;
          n1 = (i+1)*(sizes[1]+1)*(sizes[2]+1) + j*(sizes[2]+1) + k;
          n2 = (i+1)*(sizes[1]+1)*(sizes[2]+1) + (j+1)*(sizes[2]+1) + k;
          n3 = i*(sizes[1]+1)*(sizes[2]+1) + (j+1)*(sizes[2]+1) + k;
          n4 = n0+1;
          n5 = n1+1;
          n6 = n2+1;
          n7 = n3+1;

          // Check the 6 face neighbours.

          if( old_ip1[j*sizes[2]+k] == 0 ) {
            // face (n1,n2,n6,n5):
            if( nodes[n1] == -1 ) {
              nodes[n1] = count;
              count++;
            }
            if( nodes[n2] == -1 ) {
              nodes[n2] = count;
              count++;
            }
            if( nodes[n6] == -1 ) {
              nodes[n6] = count;
              count++;
            }
            if( nodes[n5] == -1 ) {
              nodes[n5] = count;
              count++;
            }
            num_bndy++;
          }

          if( old_im1[j*sizes[2]+k] == 0 ) {
            // face (n0,n4,n7,n3):
            if( nodes[n0] == -1 ) {
              nodes[n0] = count;
              count++;
            }
            if( nodes[n4] == -1 ) {
              nodes[n4] = count;
              count++;
            }
            if( nodes[n7] == -1 ) {
              nodes[n7] = count;
              count++;
            }
            if( nodes[n3] == -1 ) {
              nodes[n3] = count;
              count++;
            }
            num_bndy++;
          }

          if( old_i[j*sizes[2]+k+1] == 0 ) {
            // face (n4,n5,n6,n7):
            if( nodes[n4] == -1 ) {
              nodes[n4] = count;
              count++;
            }
            if( nodes[n5] == -1 ) {
              nodes[n5] = count;
              count++;
            }
            if( nodes[n6] == -1 ) {
              nodes[n6] = count;
              count++;
            }
            if( nodes[n7] == -1 ) {
              nodes[n7] = count;
              count++;
            }
            num_bndy++;
          }

          if( old_i[j*sizes[2]+k-1] == 0 ) {
            // face (n0,n3,n2,n1):
            if( nodes[n0] == -1 ) {
              nodes[n0] = count;
              count++;
            }
            if( nodes[n3] == -1 ) {
              nodes[n3] = count;
              count++;
            }
            if( nodes[n2] == -1 ) {
              nodes[n2] = count;
              count++;
            }
            if( nodes[n1] == -1 ) {
              nodes[n1] = count;
              count++;
            }
            num_bndy++;
          }

          if( old_i[(j+1)*sizes[2]+k] == 0 ) {
            // face (n3,n7,n6,n2):
            if( nodes[n3] == -1 ) {
              nodes[n3] = count;
              count++;
            }
            if( nodes[n2] == -1 ) {
              nodes[n2] = count;
              count++;
            }
            if( nodes[n6] == -1 ) {
              nodes[n6] = count;
              count++;
            }
            if( nodes[n7] == -1 ) {
              nodes[n7] = count;
              count++;
            }
            num_bndy++;
          }

          if( old_i[(j-1)*sizes[2]+k] == 0 ) {
            // face (n0,n1,n5,n4):
            if( nodes[n0] == -1 ) {
              nodes[n0] = count;
              count++;
            }
            if( nodes[n1] == -1 ) {
              nodes[n1] = count;
              count++;
            }
            if( nodes[n5] == -1 ) {
              nodes[n5] = count;
              count++;
            }
            if( nodes[n4] == -1 ) {
              nodes[n4] = count;
              count++;
            }
            num_bndy++;
          }
        }
      }
    }

    num_bndy *= 2;
    printf( "Found %d boundary faces\n", num_bndy );
    printf( "Found %d boundary nodes\n", count );

    // Renumbering of the boundary nodes and coordinates.

    int *index  = (int *) malloc(sizeof(int) * count);
    if( !index ) {
      printf( "Error allocating memory for index.\n" );
      exit( 1 );
    }

    for( i = 0; i < n_nodes; i++ ) {
      if( nodes[i] >= 0 ) {
        index[nodes[i]] = i;
      }
    }

    free(nodes);

    int *connec  = (int *) malloc(sizeof(int) * 3*num_bndy);
    if( !connec ) {
      printf( "Error allocating memory for connec.\n" );
      exit( 1 );
    }

    ii = 0;
    for( i = 0;  i < sizes[0];  ++i ) {
      // define pointers to the active planes.
      if( i > 0 ) {
        old_im1 = &(val[(i-1)*sizes[1]*sizes[2]]);
      }
      old_i = &(val[i*sizes[1]*sizes[2]]);
      if( i+1 < sizes[0] ) {
        old_ip1 = &(val[(i+1)*sizes[1]*sizes[2]]);
      }
      // loop on (y,z) in the x plane.
      for( j = 0;  j < sizes[1];  ++j ) {
        for( k = 0;  k < sizes[2];  ++k ) {

          if( old_i[j*sizes[2]+k] == 0 ) continue;

          // Count how many boundaries for this voxel.
          int  nn = 0;
          if( old_ip1[j*sizes[2]+k] == 0 ) nn++;
          if( old_im1[j*sizes[2]+k] == 0 ) nn++;
          if( old_i[(j+1)*sizes[2]+k] == 0 ) nn++;
          if( old_i[(j-1)*sizes[2]+k] == 0 ) nn++;
          if( old_i[j*sizes[2]+k+1] == 0 ) nn++;
          if( old_i[j*sizes[2]+k-1] == 0 ) nn++;
          if( nn == 6 ) continue;   // isolated voxel (should not occur now)
          if( nn == 0 ) continue;   // all surrounded by color - no boundary

          // Find the 8 nodes of the voxel.

          int  n0, n1, n2, n3, n4, n5, n6, n7;

          n0 = i*(sizes[1]+1)*(sizes[2]+1) + j*(sizes[2]+1) + k;
          n1 = (i+1)*(sizes[1]+1)*(sizes[2]+1) + j*(sizes[2]+1) + k;
          n2 = (i+1)*(sizes[1]+1)*(sizes[2]+1) + (j+1)*(sizes[2]+1) + k;
          n3 = i*(sizes[1]+1)*(sizes[2]+1) + (j+1)*(sizes[2]+1) + k;
          n4 = n0+1;
          n5 = n1+1;
          n6 = n2+1;
          n7 = n3+1;

          // Check the 6 face neighbours.

          if( old_ip1[j*sizes[2]+k] == 0 ) {
            connec[6*ii]   = nodes[n2];
            connec[6*ii+1] = nodes[n1];
            connec[6*ii+2] = nodes[n6];
            connec[6*ii+3] = nodes[n6];
            connec[6*ii+4] = nodes[n1];
            connec[6*ii+5] = nodes[n5];

            ii++;
          }

          if( old_im1[j*sizes[2]+k] == 0 ) {
            connec[6*ii]   = nodes[n4];
            connec[6*ii+1] = nodes[n0];
            connec[6*ii+2] = nodes[n3];
            connec[6*ii+3] = nodes[n7];
            connec[6*ii+4] = nodes[n4];
            connec[6*ii+5] = nodes[n3];

            ii++;
          }

          if( old_i[j*sizes[2]+k+1] == 0 ) {
            connec[6*ii]   = nodes[n5];
            connec[6*ii+1] = nodes[n4];
            connec[6*ii+2] = nodes[n6];
            connec[6*ii+3] = nodes[n6];
            connec[6*ii+4] = nodes[n4];
            connec[6*ii+5] = nodes[n7];
            ii++;
          }

          if( old_i[j*sizes[2]+k-1] == 0 ) {
            connec[6*ii]   = nodes[n3];
            connec[6*ii+1] = nodes[n0];
            connec[6*ii+2] = nodes[n1];
            connec[6*ii+3] = nodes[n3];
            connec[6*ii+4] = nodes[n1];
            connec[6*ii+5] = nodes[n2];

            ii++;
          }

          if( old_i[(j+1)*sizes[2]+k] == 0 ) {
            connec[6*ii]   = nodes[n7];
            connec[6*ii+1] = nodes[n3];
            connec[6*ii+2] = nodes[n6];
            connec[6*ii+3] = nodes[n6];
            connec[6*ii+4] = nodes[n3];
            connec[6*ii+5] = nodes[n2];

            ii++;
          }

          if( old_i[(j-1)*sizes[2]+k] == 0 ) {
            connec[6*ii]   = nodes[n1];
            connec[6*ii+1] = nodes[n0];
            connec[6*ii+2] = nodes[n5];
            connec[6*ii+3] = nodes[n5];
            connec[6*ii+4] = nodes[n0];
            connec[6*ii+5] = nodes[n4];

            ii++;
          }
        }
      }
    }

    // Remove points with 3 triangle neighbours, to simplify connectivity.
    int *n_ngh  = (int *) malloc(sizeof(int) * count);
    if( !n_ngh ) {
      printf( "Error allocating memory for n_ngh.\n" );
      exit( 1 );
    }
    int *ngh  = (int *) malloc(sizeof(int) * 3*count);
    if( !ngh ) {
      printf( "Error allocating memory for ngh.\n" );
      exit( 1 );
    }
    for( i = 0; i < count; i++ ) {
      n_ngh[i] = 0;
    }
    for( ii = 0; ii < num_bndy; ii++ ) {
      int n0 = connec[3*ii];
      int n1 = connec[3*ii+1];
      int n2 = connec[3*ii+2];
      if( n_ngh[n0] < 3 ) ngh[3*n0+n_ngh[n0]] = ii;
      if( n_ngh[n1] < 3 ) ngh[3*n1+n_ngh[n1]] = ii;
      if( n_ngh[n2] < 3 ) ngh[3*n2+n_ngh[n2]] = ii;
      n_ngh[n0]++;
      n_ngh[n1]++;
      n_ngh[n2]++;
    }

    int n0, n1, n2, idx[3];

    for( i = 0; i < count; i++ ) {
      if( n_ngh[i] == 3 ) {
        idx[0] = ngh[3*i];
        idx[1] = ngh[3*i+1];
        idx[2] = ngh[3*i+2];
        index[i] = -1;
        n0 = connec[3*idx[0]];
        n1 = connec[3*idx[0]+1];
        n2 = connec[3*idx[0]+2];
        for( j = 0; j < 3; j++ ) {
          ii = connec[3*idx[1]+j];
          connec[3*idx[1]+j] = -1;
          if( ii != n0 && ii != n1 && ii != n2 ) {
            if( n0 == i ) connec[3*idx[0]] = ii;
            if( n1 == i ) connec[3*idx[0]+1] = ii;
            if( n2 == i ) connec[3*idx[0]+2] = ii;
          }
          ii = connec[3*idx[2]+j];
          connec[3*idx[2]+j] = -1;
          if( ii != n0 && ii != n1 && ii != n2 ) {
            if( n0 == i ) connec[3*idx[0]] = ii;
            if( n1 == i ) connec[3*idx[0]+1] = ii;
            if( n2 == i ) connec[3*idx[0]+2] = ii;
          }
        }
      }
    }
    free(ngh);
    free(n_ngh);

    // Renumber the nodes after removal of some triangles.

    int * renum  = (int *) malloc(sizeof(int) * count);
    if( !renum ) {
      printf( "Error allocating memory for renum.\n" );
      exit( 1 );
    }

    jj = 0; 
    for( i = 0; i < count; i++ ) {
      renum[i] = -1;
      if( index[i] != -1 ) {
        renum[i] = jj;
        index[jj] = index[i];
        jj++;
      }
    }
    count = jj;
    printf( "Now have %d nodes\n", count );

    // Concatenate connectivity after removal of some triangles.

    jj = 0; 
    for( ii = 0; ii < num_bndy; ii++ ) {
      if( connec[3*ii] != -1 ) {
        connec[3*jj] = renum[connec[3*ii]];
        connec[3*jj+1] = renum[connec[3*ii+1]];
        connec[3*jj+2] = renum[connec[3*ii+2]];
        jj++;
      }
    }
    num_bndy = jj;
    printf( "Now have %d triangles\n", num_bndy );

    free(renum);

    FILE * fp = fopen( obj_file, "w" );
    fprintf( fp, "P 0.3 0.3 0.4 10 1 %d\n", count );

    Real       voxel[3], voxel_000[3], voxel_100[3], voxel_010[3], voxel_001[3];

    convert_world_to_voxel( volume, 0.0, 0.0, 0.0, voxel_000 );
    convert_world_to_voxel( volume, 1.0, 0.0, 0.0, voxel_100 );
    convert_world_to_voxel( volume, 0.0, 1.0, 0.0, voxel_010 );
    convert_world_to_voxel( volume, 0.0, 0.0, 1.0, voxel_001 );

    printf( "origin = %g %g %g\n", voxel_000[0], voxel_000[1], voxel_000[2] );
    printf( "x-axis = %g %g %g\n", voxel_100[0], voxel_100[1], voxel_100[2] );
    printf( "y-axis = %g %g %g\n", voxel_010[0], voxel_010[1], voxel_010[2] );
    printf( "z-axis = %g %g %g\n", voxel_001[0], voxel_001[1], voxel_001[2] );

    // orientation of axes for (x,y,z)
    int XX, YY, ZZ;
    if( voxel_100[0] != voxel_000[0] ) XX = 0;
    if( voxel_100[1] != voxel_000[1] ) XX = 1;
    if( voxel_100[2] != voxel_000[2] ) XX = 2;
    if( voxel_010[0] != voxel_000[0] ) YY = 0;
    if( voxel_010[1] != voxel_000[1] ) YY = 1;
    if( voxel_010[2] != voxel_000[2] ) YY = 2;
    if( voxel_001[0] != voxel_000[0] ) ZZ = 0;
    if( voxel_001[1] != voxel_000[1] ) ZZ = 1;
    if( voxel_001[2] != voxel_000[2] ) ZZ = 2;

    // offset from center to corner of voxel
    // coord = (i - origin - 1/2) * dh
    //       = (i - (origin + 1/2) * dh
    // So add dh/2 to origin (not subtract)
    voxel_000[0] += 0.5;
    voxel_000[1] += 0.5;
    voxel_000[2] += 0.5;

    // Determine the coordinates.
    float * coords  = (float *) malloc(sizeof(float) * 3 * count);
    if( !coords ) {
      printf( "Error allocating memory for coords.\n" );
      exit( 1 );
    }
    for( ii = 0; ii < count; ii++ ) {
      i = index[ii] / ( ( sizes[1]+1 ) * ( sizes[2]+1 ) );
      j = ( index[ii] - i * ( sizes[1]+1 ) * ( sizes[2]+1 ) ) / ( sizes[2]+1 );
      k = index[ii] - i * ( sizes[1]+1 ) * ( sizes[2]+1 ) - j * ( sizes[2]+1 );
      coords[3*ii+XX] = ( (float)i - voxel_000[0] ) * dh[0];
      coords[3*ii+YY] = ( (float)j - voxel_000[1] ) * dh[1];
      coords[3*ii+ZZ] = ( (float)k - voxel_000[2] ) * dh[2];
    }

    // Do a little bit of smoothing on the coordinates (simple averaging).

    float * new_coords  = (float *) malloc(sizeof(float) * 3 * count);
    if( !new_coords ) {
      printf( "Error allocating memory for new_coords.\n" );
      exit( 1 );
    }

    for( ii = 0; ii < count; ii++ ) {
      index[ii] = 0;
    }
    for( ii = 0; ii < num_bndy; ii++ ) {
      int n0 = connec[3*ii];
      int n1 = connec[3*ii+1];
      int n2 = connec[3*ii+2];
      index[n0]++;
      index[n1]++;
      index[n2]++;
    }

    float relax = 0.80;

    // for( int iter = 0; iter < 50; iter++ ) {
    for( int iter = 0; iter < 20; iter++ ) {
      if( (iter%10) == 0 ) printf( "Smoothing iteration %d...\n", iter );
      for( ii = 0; ii < 3*count; ii++ ) {
        new_coords[ii] = 0.0;
      }
      for( ii = 0; ii < num_bndy; ii++ ) {
        int n0 = connec[3*ii];
        int n1 = connec[3*ii+1];
        int n2 = connec[3*ii+2];
        for( j = 0; j < 3; j++ ) {
          float xc = coords[3*n0+j] + coords[3*n1+j] + coords[3*n2+j];
          new_coords[3*n0+j] += xc;
          new_coords[3*n1+j] += xc;
          new_coords[3*n2+j] += xc;
        }
      }
      for( ii = 0; ii < 3*count; ii++ ) {
        new_coords[ii] = 0.5 * ( new_coords[ii] / (float)(index[ii/3]) - coords[ii] );
        coords[ii] = relax * coords[ii] + ( 1.0 - relax ) * new_coords[ii];
      }
    }
    free(new_coords);
    free(index);

    // Write out the coordinates.
    for( ii = 0; ii < count; ii++ ) {
      fprintf( fp, "%g %g %g\n", coords[3*ii+0], coords[3*ii+1],
               coords[3*ii+2] );
    }

    // The normals.
    float * normals  = (float *) malloc(sizeof(float) * 3 * count);
    if( !normals ) {
      printf( "Error allocating memory for normals.\n" );
      exit( 1 );
    }
    for( ii = 0; ii < 3*count; ii++ ) {
      normals[ii] = 0.0;
    }
    for( ii = 0; ii < num_bndy; ii++ ) {
      int n0 = connec[3*ii];
      int n1 = connec[3*ii+1];
      int n2 = connec[3*ii+2];

      float dx1 = coords[3*n1] - coords[3*n0];
      float dy1 = coords[3*n1+1] - coords[3*n0+1];
      float dz1 = coords[3*n1+2] - coords[3*n0+2];

      float dx2 = coords[3*n2] - coords[3*n0];
      float dy2 = coords[3*n2+1] - coords[3*n0+1];
      float dz2 = coords[3*n2+2] - coords[3*n0+2];

      float nx = dy1 * dz2 - dy2 * dz1;
      float ny = dz1 * dx2 - dz2 * dx1;
      float nz = dx1 * dy2 - dx2 * dy1;
      float mag = sqrt( nx*nx + ny*ny + nz*nz );
      if( mag <= 1.0e-20 ) mag = 1.0e10;
      nx /= mag;
      ny /= mag;
      nz /= mag;

      normals[3*n0] += nx;
      normals[3*n0+1] += ny;
      normals[3*n0+2] += nz;

      normals[3*n1] += nx;
      normals[3*n1+1] += ny;
      normals[3*n1+2] += nz;

      normals[3*n2] += nx;
      normals[3*n2+1] += ny;
      normals[3*n2+2] += nz;
    }
    free(coords);

    fprintf( fp, "\n" );
    for( ii = 0; ii < count; ii++ ) {
      float mag = sqrt( normals[3*ii] * normals[3*ii] +
                        normals[3*ii+1] * normals[3*ii+1] +
                        normals[3*ii+2] * normals[3*ii+2] );
      if( mag <= 1.0e-20 ) mag = 1.0e10;
      fprintf( fp, "%g %g %g\n", normals[3*ii]/mag, normals[3*ii+1]/mag,
               normals[3*ii+2]/mag );
    }
    free(normals);

    // The connectivity - part 1.
    fprintf( fp, "\n" );
    fprintf( fp, "%d\n", num_bndy );
    fprintf( fp, "0 1 1 1 1\n\n" );

    for( ii = 1; ii <= num_bndy; ii++ ) {
      fprintf( fp, "%d ", 3*ii );
      if( ii%8 == 0 ) fprintf( fp, "\n" );
    }

    // The connectivity - part 2.

    for( ii = 0; ii < 3*num_bndy; ii++ ) {
      if( ii%8 == 0 ) fprintf( fp, "\n" );
      fprintf( fp, "%d ", connec[ii] );
    }
    fclose( fp );

//    (void) output_graphics_any_format( obj_file, ASCII_FORMAT, 1, &object, NULL);

    free(connec);
}

int surface( int sizes[MAX_DIMENSIONS], Real dh[MAX_DIMENSIONS],
             Volume mask, char * obj_file, double min_threshold ) {

    int    i, j, k, ii;
    Real   fval;

    int count;

    short * val = NULL;
    int * nodes = NULL;

    // Copy the volume data in temporary arrays, for speed.
    // This is faster than calling get_volume_real_value
    // many times.

    val  = (short *) malloc(sizeof(short) * sizes[0]*sizes[1]*sizes[2]);
    if( !val ) {
      printf( "Error allocating memory for val.\n" );
      exit( 1 );
    }

    count = 0;
    for( int i = 0;  i < sizes[0];  ++i ) {
      for( int j = 0;  j < sizes[1];  ++j ) {
        for( int k = 0;  k < sizes[2];  ++k ) {
          fval = get_volume_real_value( mask, i, j, k, 0, 0 );
          if (fval > min_threshold)
            val[count] = 255;
          else val[count] = 0;
          count++;
        }
      }
    }

    // Process the image.

    int ret = OK;
    ret = make_sphere( sizes, val );

    // Save surface

    save_voxel_surface( sizes, dh, val, mask, obj_file );

    free(val);

    return( ret );
}



private  void  usage(
    STRING   executable )
{
			STRING  usage_str = "\n\
	Usage: CAT_MarchingCubesSphere input.nii output_surface_file [threshold]\n\
	\n\
			 Copyright Alan C. Evans\n\
			 Professor of Neurology\n\
			 McGill University\n\n";

    print_error( usage_str, executable );
}

int  main( int argc, char* argv[] ) {
	
    STRING               input_filename, output_filename;
	  double               min_threshold;
  
    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_real_argument( 0.5, &min_threshold );

    // Read the white matter volume. 
    Volume in_volume;
    if ( input_volume_all( input_filename, 3, NULL, 
                       MI_ORIGINAL_TYPE, 0, 0, 0,
                       TRUE, &in_volume, NULL ) != OK ) {
      printf("Error: cannot read white matter volume file %s\n",input_filename);
      return 1;
    }

    if ( get_volume_n_dimensions( in_volume ) != 3 ) {
      printf("Error: volume in %s does not have three dimensions\n",input_filename);
      return 1;
    }

    int sizes[MAX_DIMENSIONS];
    get_volume_sizes( in_volume, sizes );

    Real       dh[3];
    get_volume_separations( in_volume, dh );

    int ret = surface( sizes, dh, in_volume, output_filename, min_threshold );

    delete_volume( in_volume );

    return ( ret != OK );

}

