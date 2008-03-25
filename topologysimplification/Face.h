#pragma once

#define FACE_IS_TRIANGLE 1

class Face
{
public:

#if (!FACE_IS_TRIANGLE)
	int nv; // number of vertices 
	int *v; // list of vertices
	int *f; // list of neighboring faces
#else
	int v[3]; // list of vertices
	int f[3]; // list of neighboring faces
#endif

	int marked;	//for computational purposes

	double nx,ny,nz; // the face normal
	double x,y,z; // the coordinates of the face

	//constructor/destructor
	Face(void);
	~Face(void);
	const Face& operator=(const Face &face);
};
