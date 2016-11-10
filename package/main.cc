#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include "pbig.h"
#include "report.h"
#include "rect.h"


int main()
{
	int numRects = 1000000;
	int maxrange = 10000;
	int deviceID = 0;
	struct timeval begin, end;

	struct RectInfo *input = (struct RectInfo*) malloc(sizeof(struct RectInfo)*numRects);
	printf("GPU Initialization ...\n");
	pbig_init( numRects, maxrange, 0, deviceID, 0);


	for(int j=0;j<10;j++) {

		// Prepare input
		for(int i=0;i<numRects;i++) {
			input[i].rid = i;
			input[i].coordinate.xmin = rand()%10000;
			input[i].coordinate.xmax = input[i].coordinate.xmin + 100;
			input[i].coordinate.ymin = rand()%10000;
			input[i].coordinate.ymax = input[i].coordinate.ymin + 100;
			input[i].coordinate.zmin = rand()%10000;
			input[i].coordinate.zmax = input[i].coordinate.zmin + 100;
		}

		report r;
		printf("Intersection Checking ...\n");
		gettimeofday(&begin, NULL);
			pbig_check(input, &r);
		gettimeofday(&end, NULL);

		// Show number of intersecting pairs and elapsed time
		r.print();
		printf("%ld us\n", 1000000*(end.tv_sec-begin.tv_sec)+end.tv_usec-begin.tv_usec);
	}



	delete input;
	pbig_finalize();
}
