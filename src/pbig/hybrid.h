#ifndef __HYBRID_H_
#define __HYBRID_H_

using namespace std;
#include <vector>
#include <set>

#include <time.h>
#include <unistd.h>
#include <sys/time.h>
extern struct timeval begin, end;
extern float intersect_time;


extern unsigned int HAVE_STEPS;

struct report;
void Intersect_Hybrid(unsigned int, struct RectInfo *, struct report*);
void Init_Hybrid(unsigned int, unsigned int, unsigned int, unsigned int,unsigned int); 
void Finalize_Hybrid();


#if USE_SEQ
	typedef vector<unsigned int> ObjectSet;
	//typedef unsigned int* ObjectSet;

	extern ObjectSet *cellObjSets;
	void p_hybrid_init(unsigned int);
	void *t_run(void *); 
#endif




#endif
