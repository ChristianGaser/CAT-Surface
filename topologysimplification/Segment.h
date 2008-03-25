#pragma once

#define NUMBER_OF_POINTS 10 
#define INCREASE_NUMBER_OF_POINTS 1.2

class Segment
{
private:
	int npoints,maxpoints;
	int* points;

	int nvertices,nfaces,nedges,euler;
	int marked;

	void _ReallocSegment(int new_maxpoints=-1);

public:
	// constructor / destructor
	Segment();
	~Segment();

	int size() const;
	void clear();
	void AddPoint(int pt);
	void AddSegment(Segment *s);
	void Transfer(Segment &b);
	int GetEuler();
	const int* GetPointList(void) const;
	int GetMark() const;
	void SetMark(int m);
};
