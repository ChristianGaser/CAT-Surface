/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2012: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is dual-licensed as follows:                                 *
*                                                                           *
* (1) You may use TMesh as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use TMesh as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "heap.h"

namespace T_MESH
{

abstractHeap::abstractHeap(int size)
{
 heap = new void *[size+1];
 numels = 0;
 maxels = size;
 positions = NULL;
}

abstractHeap::~abstractHeap()
{
 delete(heap);
}

int abstractHeap::upheap(int k)
{
 if (k < 2) return k;
 
 void *t = heap[k];
 int fk = (k%2)?((k-1)/2):(k/2);
 void *f = heap[fk];

 if (compare(t, f) <= 0)
 {
  heap[k] = f;
  heap[fk] = t;
  if (positions != NULL)
  {
   positions[(intptr_t)f] = k;
   positions[(intptr_t)t] = fk;
  }
  return upheap(fk);
 }
 return k;
}

int abstractHeap::downheap(int k)
{
 int j;
 
 void *t = heap[k];
 int fk = (numels%2)?((numels-1)/2):(numels/2);
 if (k > fk) return k;
 
 j = k+k;
 if (j < numels && compare(heap[j], heap[j+1]) >= 0) j++;
 void *f = heap[j];
 if (compare(t, f) >= 0)
 {
  heap[k] = f;
  heap[j] = t;
  if (positions != NULL)
  {
   positions[(intptr_t)f] = k;
   positions[(intptr_t)t] = j;
  }
  return downheap(j);
 }

 return k;
}

int abstractHeap::insert(void *t)
{
 if (numels == maxels) return -1;

 heap[++numels] = t;
 if (positions != NULL) positions[(intptr_t)t] = numels;
 return upheap(numels);
}

void *abstractHeap::removeHead()
{
 void *t = heap[1];
 if (positions != NULL) positions[(intptr_t)t] = 0;
 heap[1] = heap[numels--];
 if (numels)
 {
  if (positions != NULL) positions[(intptr_t)heap[1]] = 1;
  downheap(1);
 }

 return t;
}

} //namespace T_MESH
