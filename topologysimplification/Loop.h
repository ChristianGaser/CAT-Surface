#pragma once

class Loop
{
private:
	void _ReAlloc(int maxpts=-1);
public:
	int npoints,maxpoints;
	int *points;
public:
	//constructor/destructor
	Loop(void);
	Loop(int maxpts);
	~Loop(void);

	void Alloc(int maxpts);
	void AddPoint(int pt);
	void Print() const;
	int End();
	void Pop();
	int Replace(int pt, int new_pt);
	int operator[](int n){
		//Assert((n >= 0 ) && (n < npoints));	
		return points[n];
	}
	const Loop & operator=(const Loop& loop);
};
