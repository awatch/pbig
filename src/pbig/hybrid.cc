using namespace std;
#include <iostream>
#include "hybrid.h"

#if USE_PBIG || USE_PBIG_LIST
#include "pbig.h"
#endif



void Init_Hybrid(unsigned int numObjs, unsigned int range, unsigned int cs, unsigned int device, unsigned int hasProfSteps)
{


#if USE_PBIG
	int device_init(unsigned int);
	int deviceID = device_init( device);
	pbig_init( numObjs, range, cs, deviceID, hasProfSteps);
#endif
}

void Intersect_Hybrid(unsigned int num, RectInfo *input, struct report *r)
{



#if USE_PBIG 
	pbig_check(input, r);
#endif
}

void Finalize_Hybrid() {


#if USE_PBIG 
	pbig_finalize();
#endif
}

