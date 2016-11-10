//#define LOCALCOMM

using namespace std;

#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <pthread.h>
#include <math.h>
#include "object.h"
#include "rect.h"
#include "report.h"

#define STDRATIO 0.2



#if USE_PBIG
#include "hybrid.h"
#include "pbig.h"
#endif







//#define OVERLAP_BITMAP 1
#ifdef OVERLAP_BITMAP
#include <bitset>
unsigned int *allOverlapSets;
#else
#endif


RectInfo *input;


struct timeval begin, end;
struct timeval beginCell, endCell;
struct timeval cutBegin, cutEnd;
float total_time = 0;
float total_timeCore = 0;
float intersect_time = 0;
uint numObjs=0;
uint cellSize=0;



int NThreads;

//long long total_matchs = 0;
//long long real_matchs = 0;
//long long repeat_matchs = 0;
//long long unnece_matchs = 0;








/* the function array is related to enum func_e in Endecode.h */



//struct object *objects;
unsigned int maxRange=0;
int TEST_STEP=-1;
int device=-1;
int InputRectSize = 0;
int PossRegCount = 0;
char *SimCase;
//unsigned int RegsUsed = 0;

unsigned int hasProfSteps = 0;


void ParseArg(int argc, char *argv[]) {

    int c;

    while (1) {
        //int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[] = {
            {"cs", 1, 0, 'c'},
            {"device", 1, 0, 'd'},
            {"regs", 1, 0, 'f'},
            {"ds", 1, 0, 'g'},
            {"numobj", 1, 0, 'n'},
            {"step", 1, 0, 'p'},
            {"numthreads", 1, 0, 'r'},
            {"rectSizes", 1, 0, 'S'},
            {"regtype", 1, 0, 't'},
            {"psteps", 1, 0, 'x'},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "c:d:f:g:n:p:r:S:t:x:",
                long_options, &option_index);


        if (c == -1)
            break;

        switch (c) {
            case 0:
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'c':
                printf ("Grid Cell Size with value '%s'\n", optarg);
#if USE_PBIG
                cellSize = atoi(optarg);
#endif
                break;

            case 'd':
                printf ("deivce with value '%s'\n", optarg);
                device = atoi(optarg);
                break;

            case 'f':
                //printf ("registers with value '%s'\n", optarg);
                //RegsUsed = atoi(optarg);
                break;

			case 'g':
                printf ("Dimension Size with value '%s'\n", optarg);
                maxRange = atoi(optarg);
                break;

            case 'n':
                printf ("numObjs with value '%s'\n", optarg);
                numObjs = atoi(optarg);
                break;

            case 'p':
                printf ("Steps with value '%s'\n", optarg);
                TEST_STEP = atoi(optarg);
                break;

            case 'r':
                printf ("Threads with value '%s'\n", optarg);
                NThreads = atoi(optarg);
                break;

            case 'S':
                {
                    char *tok;
                    printf ("InputRectSize with value '%s'\n", optarg);
					InputRectSize = atoi(optarg);
#if 0
                    tok = strtok( optarg, "_");
                    while ( tok != NULL ) { 
                        PossRectSizes[PossRegCount] = atoi(tok);
                        printf(" PossRectSizes[%d]=%d", PossRegCount, PossRectSizes[PossRegCount++]);
                        tok = strtok( NULL, "_");
                    }
                    printf("\n");
#endif
                }
                break;

            case 't':
                printf ("SimCase with value '%s'\n", optarg);
                SimCase = optarg;
                break;

            case 'x':
                printf ("Profiling Steps '%s'\n", optarg);
                hasProfSteps = atoi(optarg);
                break;

            case '?':
                break;

            default:
                printf ("?? getopt returned character code 0%o ??\n", c);
        }
    }

    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
    }

}

double getUniformRnd(unsigned int max)
{
	//return max * (((double)rand()) / double(RAND_MAX));
	return rand() % (unsigned int) max;
}

double getNormalRnd(unsigned int mean, unsigned int std) 
{
	
	double u = rand() / (double) RAND_MAX;
	double v = rand() / (double) RAND_MAX;
	double x = sqrt( -2 * log(u)) * cos( 2 * M_PI * v ) * std + mean;

	return x;
}

#define Inc(X, Y, MAX) (X+Y) >= (MAX) ? (MAX-1) : (X)+(Y);
#define Dec(X, Y, MIN) (X-Y) <  (MIN) ? (MIN) : (X)-(Y);

void InitTestPerformance() {

    uint i;

    srand(100);
    

    // decide the sizes of regions of objects and the locations of objects
	printf("Object Size = %ld\n", sizeof(struct object));
	printf("numObjs=%d, Total Memory Size for Objects = %ld\n", numObjs, numObjs*sizeof(struct object));
    //objects = new struct object[numObjs];

#if 1
    for(i=0;i<numObjs;i++) {

		//objects[i].id = i;
		int posX, posY, posZ;
		int rectSizes[3];	

        // set the rectangle size
        if ( SimCase[0] == 'U' ) {
			rectSizes[0] = getUniformRnd(  InputRectSize  );
			rectSizes[1] = getUniformRnd(  InputRectSize  );
#if DIMS == 3
			rectSizes[2] = getUniformRnd(  InputRectSize  );
#endif
#if 0
            objects[i].sizeIndex = 0;
            objects[i].rectSizes[0] = rand() % PossRectSizes[objects[i].sizeIndex];
            objects[i].rectSizes[1] = rand() % PossRectSizes[objects[i].sizeIndex];
#if DIMS == 3
            objects[i].rectSizes[2] = rand() % PossRectSizes[objects[i].sizeIndex];
#endif
#endif
        }
		else if ( SimCase[0] == 'N' ) {
			rectSizes[0] = getNormalRnd(  InputRectSize , InputRectSize*STDRATIO );
			rectSizes[1] = getNormalRnd(  InputRectSize , InputRectSize*STDRATIO );
#if DIMS == 3
			rectSizes[2] = getNormalRnd(  InputRectSize , InputRectSize*STDRATIO );
#endif
		}
        else {
            printf("Error SimCase\n");
            exit(-1);
        }

        // set the location according to the size
        if ( SimCase[1] == 'U' ) {
			//objects[i].x = getUniformRnd( maxRange );	
			//objects[i].y = getUniformRnd( maxRange ); 
			posX = getUniformRnd( maxRange );	
			posY = getUniformRnd( maxRange ); 
#if DIMS ==  3
			//objects[i].z = getUniformRnd( maxRange ); 
			posZ = getUniformRnd( maxRange ); 
#endif
		}
		else if ( SimCase[1] == 'N' ) {
			//objects[i].x = getNormalRnd(  maxRange/2 , (maxRange/2)*STDRATIO );	
			//objects[i].y = getNormalRnd(  maxRange/2 , (maxRange/2)*STDRATIO ); 
			posX = getNormalRnd(  maxRange/2 , (maxRange/2)*STDRATIO );	
			posY = getNormalRnd(  maxRange/2 , (maxRange/2)*STDRATIO ); 
#if DIMS == 3 
			//objects[i].z = getNormalRnd(  maxRange/2 , (maxRange/2)*STDRATIO ); 
			posZ = getNormalRnd(  maxRange/2 , (maxRange/2)*STDRATIO ); 
#endif
		}
		else {
            printf("Error SimCase\n");
            exit(-1);
		}

#if 0
#if DIMS == 2
        objects[i].move = ( objects[i].rectSizes[0] + objects[i].rectSizes[1] ) / 2;
#elif DIMS == 3
        objects[i].move = ( objects[i].rectSizes[0] + objects[i].rectSizes[1] + objects[i].rectSizes[2] ) / 3;
#endif
        objects[i].dir  = 0;
        if ( i % 2 == 0 )
            objects[i].type = 1;
        else
            objects[i].type = 2;
#endif

		//objects[i].xmin = Dec( objects[i].x, objects[i].regSize, 0);
		//objects[i].xmin = posX; //objects[i].x;
		//objects[i].xmax = Inc( posX/*objects[i].x*/, rectSizes[0], maxRange);
		input[i].coordinate.xmin = posX; //objects[i].x;
		input[i].coordinate.xmax = Inc( posX/*objects[i].x*/, rectSizes[0], maxRange);

		//objects[i].ymin = Dec( objects[i].y, objects[i].regSize, 0);
		//objects[i].ymin = posY;//objects[i].y;
		//objects[i].ymax = Inc( posY/*objects[i].y*/, rectSizes[1], maxRange);
		input[i].coordinate.ymin = posY;//objects[i].y;
		input[i].coordinate.ymax = Inc( posY/*objects[i].y*/, rectSizes[1], maxRange);

#if DIMS ==  3
		//objects[i].zmin = posZ;//objects[i].z;
		//objects[i].zmax = Inc( posZ/*objects[i].z*/, rectSizes[2], maxRange);
		input[i].coordinate.zmin = posZ;//objects[i].z;
		input[i].coordinate.zmax = Inc( posZ/*objects[i].z*/, rectSizes[2], maxRange);
#endif
		input[i].rid = i;

		//showObjInfo( objects, i);
    }
#endif



}


#if 0
void objMove(int step)
{

    int i;


    for(i=0;i<numObjs;i++) {

#if 1 
        //objects[i].dir = rand() % 4 + 1;
        objects[i].dir = 0;

        if ( objects[i].dir == 0 )
            continue;
        if ( objects[i].dir == 1 ) // east
            objects[i].x = (objects[i].x+objects[i].move) >= maxRange ? maxRange-1 : objects[i].x+objects[i].move;
        if ( objects[i].dir == 2 ) // west
            objects[i].x = (int)(objects[i].x-objects[i].move) < 0 ? 0 : objects[i].x-objects[i].move;
        if ( objects[i].dir == 3 ) // south
            objects[i].y = (objects[i].y+objects[i].move) >= maxRange ? maxRange-1 : objects[i].y+objects[i].move;
        if ( objects[i].dir == 4 ) // north
            objects[i].y = (int)(objects[i].y-objects[i].move) < 0 ? 0 : objects[i].y-objects[i].move;
#if DIMS ==  3
        if ( objects[i].dir == 5 ) // up 
            objects[i].z = (objects[i].z+objects[i].move) >= maxRange ? maxRange-1 : objects[i].z+objects[i].move;
        if ( objects[i].dir == 6 ) // down 
            objects[i].z = (int)(objects[i].z-objects[i].move) < 0 ? 0 : objects[i].z-objects[i].move;
#endif

		//objects[i].xmin = Dec( objects[i].x, objects[i].regSize, 0);
		objects[i].xmin = objects[i].x;
		objects[i].xmax = Inc( objects[i].x, objects[i].rectSizes[0], maxRange);

		//objects[i].ymin = Dec( objects[i].y, objects[i].regSize, 0);
		objects[i].ymin = objects[i].y;
		objects[i].ymax = Inc( objects[i].y, objects[i].rectSizes[1], maxRange);

#if DIMS ==  3
		objects[i].zmin = objects[i].z;
		objects[i].zmax = Inc( objects[i].z, objects[i].rectSizes[2], maxRange);
#endif

		
#endif
    }

}
#endif


static void statistics(int step)
{

    cout << "\t=====================================" << endl;
    cout << "\t*Results report                      " << endl;
    cout << "\t=====================================" << endl;
    cout<<"\tmatch latency = "<< total_time << " us" << endl;
    cout<<"\tmatch latency Core = "<< intersect_time << " us" << endl;
    cout<<"\t" << "step" << step << "," << (total_time)  << endl;


    total_time = 0.0;
    total_timeCore = 0.0;
	intersect_time = 0.0;

}

void PrepareData() {


#if 0
    for(int i=0;i<numObjs;i++) {
#if DIMS == 2
		input[i].coordinate.xmin = objects[i].xmin;
		input[i].coordinate.xmax = objects[i].xmax;
		input[i].coordinate.ymin = objects[i].ymin;
		input[i].coordinate.ymax = objects[i].ymax;
#elif DIMS == 3
		input[i].coordinate.xmin = objects[i].xmin;
		input[i].coordinate.xmax = objects[i].xmax;
		input[i].coordinate.ymin = objects[i].ymin;
		input[i].coordinate.ymax = objects[i].ymax;
		input[i].coordinate.zmin = objects[i].zmin;
		input[i].coordinate.zmax = objects[i].zmax;
#endif
		input[i].rid  = i;//objects[i].id; 
    }
#endif
}

// check answer
#if 0 
void CheckMatchResult() {



}
#endif


void Intersect()
{

	//printf("Prepare Data\n");
	//PrepareData();

	printf("Running ...\n");


#if USE_PBIG 
    gettimeofday(&begin, NULL);
		report r;
		Intersect_Hybrid(numObjs, input, &r);
		r.print();
		r.reset();
    gettimeofday(&end, NULL);
    total_time += 1000000*(end.tv_sec-begin.tv_sec)+end.tv_usec-begin.tv_usec;
#endif


}


void TestMainLoop(int countDimension)
{


    /******************************************/
    /*         main loop                      */
    /******************************************/
    cout << __func__ << endl;

    cout << "Allocate memory for arrayDimension" << endl;
    int step;

    for(step=0;step<=TEST_STEP;step++) {

        // adjust object's position
        //objMove(step);

		Intersect();
#if 0
		CheckMatchResult();
#endif

        statistics(step);

    }

}


void TestPerformance()
{


    int countDimension = 2;
    
    InitTestPerformance();

    TestMainLoop( countDimension);


}


void InitAlgo()
{

	input = (RectInfo*) malloc(sizeof(RectInfo)*numObjs);


#if USE_PBIG 
	Init_Hybrid(numObjs, maxRange, cellSize, device, hasProfSteps);
#endif



#ifdef OVERLAP_BITMAP
    allOverlapSets = new unsigned int[numObjs*numObjs];
    for(int i=0;i<numObjs*numObjs;i++) 
        allOverlapSets[i] = 0;
#endif


}


int main(int argc, char **argv)
{


    ParseArg(argc,argv);

    InitAlgo();

    TestPerformance();


#if USE_PBIG 
	Finalize_Hybrid();
#endif

	
}





/*
 * vim: ts=4 sts=4 sw=4
 */
