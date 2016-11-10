#ifndef __PBIG_H_
#define __PBIG_H_

using namespace std;
#include <vector>
#include <set>

struct report;
void pbig_init(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int);
void pbig_check(struct RectInfo*, struct report*);
void pbig_finalize();


#endif

