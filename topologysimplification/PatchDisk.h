#pragma once
#include "SurfaceObj.h"

class PatchDisk
{
	void _Init();
	void _Alloc(int which_patch);
public:
	Surface disk;
	Loop ring,init_ring;
	int *vtrans ;
	int *ftrans ;

	PatchDisk();
	PatchDisk(int which_patch);
	PatchDisk(const string s):disk(s),init_ring(10){
		vtrans = new int[disk.nvertices];
		ftrans = new int[disk.nfaces];
		_Init();
	};
	~PatchDisk(void);

	void Init();
	void Create(int which_patch);
};


//create + _Alloc + constructeurs
//main
// surface.cpp : knitpatch et cutloop (pdisk)
