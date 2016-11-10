#ifndef __RECT_H_
#define __RECT_H_

typedef unsigned int CType;

typedef struct CoordInfo {
	CType xmin, xmax; 
	CType ymin, ymax;
#if DIMS == 3
	CType zmin, zmax;
#endif
} CoordInfo;


typedef struct RectInfo {
	CoordInfo coordinate;
    unsigned int rid;   
} RectInfo;

#endif
