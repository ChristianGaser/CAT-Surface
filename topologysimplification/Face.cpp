#include "Face.h"

Face::Face(void)
{
	marked=0;

#if (!FACE_IS_TRIANGLE)
	// by default, a face is a triangle
	nv = 3;
	v = new int[nv];
	f = new int[nv];
#endif
}

Face::~Face(void)
{
#if (!FACE_IS_TRIANGLE)
	if(v) delete [] v;
	if(f) delete [] f;
#endif
}

const Face &Face::operator=(const Face &face)
{
#if (!FACE_IS_TRIANGLE)
	if(nv < face.nv){
		delete [] v;
		delete [] f;
		v = new int[face.nv];
		f = new int[face.nv];
	}
	nv=face.nv;
	for(int n = 0 ; n < nv ; n++){
		v[n] = face.v[n];
		f[n] = face.f[n];
	}
#else
	for(int n = 0 ; n < 3 ; n++){
		v[n] = face.v[n];
		f[n] = face.f[n];
	}
#endif
	marked=face.marked;	
	nx=face.nx; ny=face.ny; nz=face.nz; 
	x=face.x; y=face.y; z=face.z; 

	return face;
}
