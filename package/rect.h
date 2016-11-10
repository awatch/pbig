#ifndef __RECT_H_
#define __RECT_H_

typedef unsigned int CType;

typedef struct CoordInfo {
	CType xmin, xmax; 
	CType ymin, ymax;
	CType zmin, zmax;
} CoordInfo;


typedef struct RectInfo {
	CoordInfo coordinate;
    unsigned int rid;   
} RectInfo;

#endif
