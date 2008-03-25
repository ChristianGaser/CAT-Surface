#include "Segment.h"
#include <iostream>
#include <cstring>

Segment::Segment(){
	npoints=0;
	marked=0; 
	maxpoints=NUMBER_OF_POINTS;
	points = new int[maxpoints];
}

Segment::~Segment(){
	if(points) delete [] points;
}

int Segment::GetMark() const{ 
	return marked;
}
void Segment::SetMark(int m){
	marked=m;
}

void Segment::_ReallocSegment(int new_maxpoints){
	if(new_maxpoints < 0 || new_maxpoints < maxpoints) 
		new_maxpoints = int(maxpoints*INCREASE_NUMBER_OF_POINTS+1);
	
	int *new_points = new int[new_maxpoints];
	for(int i = 0 ; i < maxpoints ; i++)
		new_points[i]=points[i];

	maxpoints=new_maxpoints;
	delete [] points;
	points=new_points;
}

void Segment::AddPoint(int pt){
	if(npoints == maxpoints)
		_ReallocSegment();
	points[npoints++]=pt;
}

int Segment::size() const{
	return npoints;
}
void Segment::clear(){
	npoints=0;
}

void Segment::AddSegment(Segment *s){
	int new_maxpoints = npoints+s->npoints;
	if(new_maxpoints > maxpoints) _ReallocSegment(new_maxpoints);
	
	//copy the segment s 
	memcpy(&points[npoints],s->points,s->npoints*sizeof(int));
	npoints += s->npoints;

	s->clear();
}

const int* Segment::GetPointList(void) const
{
	return points;
}

void Segment::Transfer(Segment &b){
	npoints=b.npoints;
	maxpoints=b.maxpoints;
	delete [] points;
	points=b.points;
	b.maxpoints=0;
	b.npoints=0;
	b.points=NULL;
}

