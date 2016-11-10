
#ifndef CUDA_COMM_H_
#define CUDA_COMM_H_

#include "log.h"
#include <cuda.h>



#define UNUSED 0x0fffff11

#define MAXLOG 3

#define CHECK_CUDA_ERR \
	{ \
	cudaError_t cerr = cudaGetLastError(); \
	if ( cerr != cudaSuccess ) { \
		printf("\n"); \
		printf("\t\t\t\t\t\t\t\t\t[CUDA] %d %s", __LINE__, cudaGetErrorString(cerr)); \
		for(int i=0;i<MAXLOG;i++) \
			printf(" %d", h_log[i]); \
		exit(-1);  cuda_finalize(); } \
	}

#define CUDA_COPY_LOG \
	{ \
	cudaMemcpy( h_log, d_log, MAXLOG*sizeof(unsigned int), cudaMemcpyDeviceToHost); \
		for(int i=0;i<1;i++) \
			LOG(" %d", h_log[i]); \
	}


#define CUDA_INIT_LOG \
	{ \
	for(int i=0;i<MAXLOG;i++) h_log[i]=0; \
	cudaMemcpy( d_log, h_log, MAXLOG*sizeof(unsigned int), cudaMemcpyHostToDevice); \
	}

#define CUDA_EVENT_START \
	cudaEventRecord(e_start, 0);

#define CUDA_EVENT_END \
	cudaEventRecord(e_stop, 0); \
	cudaEventSynchronize(e_stop); \
	cudaEventElapsedTime(&elapsedTime, e_start, e_stop);

#define MEM_STATUS(s) \
	do { \
		size_t free_mem, total_mem; \
		cudaMemGetInfo( &free_mem, &total_mem); \
		free_mem /= 1024*1024; \
		total_mem /= 1024*1024; \
		LOG(" + %s\n", s);  \
		LOG(" ++++ Total Device Memory free: %ld MB, used: %ld MB, total: %ld MB\n", free_mem, total_mem-free_mem, total_mem); \
	} while(0)


#define CHECK_AVAIL_MEM(size) \
	do { \
		size_t free_mem, total_mem; \
		cudaMemGetInfo( &free_mem, &total_mem); \
		if ( size >= free_mem ) { \
			printf("Required: %d, Free: %d\n", size, free_mem); \
			CUDA_ERR("out of device memory"); \
		} \
	} while(0)

#define IS_AVAIL_MEM(avail,size) \
	do { \
		size_t free_mem, total_mem; \
		cudaMemGetInfo( &free_mem, &total_mem); \
		if ( size >= free_mem ) { \
			*avail = 0; \
		} \
		else \
			*avail = 1; \
	} while(0)


#ifdef __CUDACC__
#define ALIGN(x) __align__(x)
#else
#define ALIGN(x) __attribute__ ((aligned (x)))
#endif


#endif

