#ifndef __REPORT_H
#define __REPORT_H

#include "rect.h"

struct report {

public:
	report();
	void operator()(unsigned int a, unsigned int b){ sum++; };
	void found(const RectInfo *a, const RectInfo *b);
	void print();
	void reset() { sum = 0; }
private:
	int sum;
};

#endif 
