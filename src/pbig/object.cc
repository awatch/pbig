using namespace std;
#include <iostream>
#include "object.h"


#if 0
unsigned int pair_intersect(struct object *i, struct object *j){


	if ( i->id == j->id )
		return 0;

#if DIMS == 2
	if ( i->xmin <= j->xmax && 
		 j->xmin <= i->xmax && 
		 i->ymin <= j->ymax && 
		 j->ymin <= i->ymax ) 
#elif DIMS == 3
	if ( i->xmin <= j->xmax && 
		 j->xmin <= i->xmax && 
		 i->ymin <= j->ymax && 
		 j->ymin <= i->ymax && 
		 i->zmin <= j->zmax && 
		 j->zmin <= i->zmax ) 
#endif

	{
		//allOverlapSets[(long long int)i->id*NumObjs+(long long int)j->id] = 1;
		//i->overlap.insert(j->id);	
		return j->id;
	}
	else
		return 0;

}
#endif

#if 0
void swapObj(struct object *o1, struct object *o2) {

	struct object t = *o1;
	unsigned int rid1, rid2;

	rid1 = o1->id;
	rid2 = o2->id;

	*o1 = *o2;
	*o2 = t;

	o1->id = rid1;
	o2->id = rid2;
}

void showObjInfo(struct object *objects, int i)
{

	printf("Object %d pos(%6d,%6d) rect [%6d %6d][%6d %6d] mov_dist %d regSizes [%d %d] dir %d\n", i, 
			                                                           objects[i].x, 
																	   objects[i].y, 
																	   objects[i].xmin,
																	   objects[i].xmax,
																	   objects[i].ymin,
																	   objects[i].ymax,
																	   objects[i].zmin,
																	   objects[i].zmax,
																	   objects[i].move, 
																	   objects[i].regSizes[0], 
																	   objects[i].regSizes[1], 
																	   objects[i].regSizes[2], 
																	   objects[i].dir);
}
#endif
