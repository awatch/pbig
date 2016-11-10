#ifndef __PBIG_INTERNAL_H
#define __PBIG_INTERNAL_H

#ifdef RICE_COMPRESSION
typedef unsigned long long int OTYPE;
#define OTYPE_SHIFT 3
#else
typedef unsigned int OTYPE;
#define OTYPE_SHIFT 2
#endif

#define MAXBLOCKS 64000
#define MAXJOBS   5000000
//#define MAX_Grid_Recs (MAXALLOCRECTS*6)
//#define MATCHING_MB 240
#define MaxBatchRecs (480*NTPB)


#ifndef COMPRESSION
#define OVERLAP_LEN 8
#else
#define OVERLAP_LEN 4
#endif


#if PROFILING == 1 
#define B_PROFILING(n) \
	if ( hasProfSteps >= n ) { 
#define E_PROFILING(n) \
	}
#define ASSERT(n) \
{ \
	if ( threadIdx.x >= NTPB ) \
		log[2]+=n; \
}
#else
#define B_PROFILING(n) 
#define E_PROFILING(n) 
#define ASSERT(n) 
#endif

#define STime(timer) \
{                                                                                   \
	cutResetTimer(timer);                                                          \
	cutStartTimer(timer);                                                          \
}

#define ETime(timer) \
{                                                                                   \
	cutilSafeCall( cutilDeviceSynchronize() );                                      \
	cutStopTimer(timer);                                                           \
}

#define GTime(timer) cutGetTimerValue(timer)

#endif
