#ifndef __OBJECT_H
#define __OBJECT_H

using namespace std;
#include <set>

struct object {

	//unsigned int id;
    //unsigned int x;
    //unsigned int y;
    //unsigned int type;
    //int sizeIndex;
#if DIMS == 2
    //unsigned int rectSizes[2];
#elif DIMS == 3
    //unsigned int rectSizes[3];
#endif
    //unsigned int move;
    //unsigned int dir;
	unsigned int xmin;
	unsigned int xmax;
	unsigned int ymin;
	unsigned int ymax;
#if DIMS == 3
    //unsigned int z;
	unsigned int zmin;
	unsigned int zmax;
#endif
	//set<unsigned int> overlap;
};

unsigned int pair_intersect(struct object *, struct object *);
void showObjInfo(struct object *, int);
void swapObj(struct object*, struct object *);




#endif
