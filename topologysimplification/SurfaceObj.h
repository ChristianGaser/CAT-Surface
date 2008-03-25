#pragma once

#include "Vertex.h"
#include "Face.h"
#include "Loop.h"

extern "C"
{
	#include  <bicpl.h>
}

#include <string>
using namespace std;

#define UNKNOWN_TYPE_OF_SURFACE 0 
#define CLOSED_SURFACE 1 
#define OPEN_SURFACE 2 
#define PATCH 3

class PatchDisk;

class Surface
{
private:
	int _InitSurfaceConnectivity(void);
	bool _InitFaceConnectivity(void);
	int _InitFaceCoordinates(void);
	int _FindFace(int n1, int n2, int fn) const;
	int _Allocate(int nv, int nf);
	void _OverAlloc(int nextrav,int nextraf);
	double _FaceDistance(int fdst, int fsrc);

public:
	// number of vertices, edges, and faces, and Euler characteristic
	int nvertices, nedges, nfaces, euler;
	// list of vertices
	int maxvertices;
	Vertex *vertices;
	// list of faces
	int maxfaces;
	Face *faces;
	//type of surface
	int type_of_surface;


	//constructor/destructor
	Surface(void);
	Surface(int nv, int nf);
	Surface(const string s);
	~Surface(void);

	Surface *Clone() const;
	void Expand(int nextrav,int nextraf);
	void Center();
	void scale(double scaling_factor);
	int OpenFile(const string s,int verbose=0);
	int ObjToSurface(object_struct* object);
	object_struct *SurfaceToObj();
	int WriteFile(const string s, int verbose = -1) const;
	int WriteCurvature(const string s, int verbose) const;
	int GetDefectLabels(const string s);
	int OpenCurvatureFile(const string s);
	bool IsSurfaceValid(int verbose = 0);
	void PrintDefectInfo(int ndefect=-1);
	int InitSurface(void);
	int GetEuler(int &nv, int &ne, int &nf, int mark = -1);
	int GetEuler(int mark = -1);
	int GetEuler(const int *list_of_faces, int nfs);
	Surface *ExtractPatch(int mark,int nextravertices=0, int nextrafaces=0);
	void SetMarks(int mark);
	void SetMarks(const int *v, int nv, int mark);
	void Smooth(int niters);
	void Smooth(int niters, const int *v,int nv);
	void SmoothMarked(int niters,int mark);
	void ExpandMarks(int niters,int mark);
	
	///////////////////////////////////////////////////////////////
	//
	//       For The Topology Correction
	//
	///////////////////////////////////////////////////////////////
	
	//when the surface is a patch extracted from another surface
	int *vtrans_to,*ftrans_to,*vtrans_from,*ftrans_from;
	Surface *surface_source;
	PatchDisk *disk;

	double GetLoopLength(Loop &loop);
	void OpenLoop(Loop &loop);
	void CutLoop(Loop &loop);
	bool LoopValid(Loop &loop);
	void KnitPatch(Loop &loop , PatchDisk *pdisk);
	void IncreaseEuler(int nattempts,int maxinitface = -1);
	void CorrectTopology();
	int CutPatch(int seed=-1, int maxinitface = -1, int nattempts = 10);


};


