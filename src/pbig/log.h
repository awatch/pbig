#ifndef __LOG_H
#define __LOG_H

#define CUDA_ERR(fmt, rest...) \
	{ \
			printf("ERROR:"); \
			printf(fmt, ## rest); \
			if (RELEASE != 1) \
			printf(" at file %s line %d\n", __FILE__,  __LINE__); \
			cuda_finalize(); \
			exit(-1); \
	}

#define ERR(fmt, rest...) \
	{ \
			printf("ERROR:"); \
			printf(fmt, ## rest); \
			exit(-1); \
	}

//#define RELEASE_VERSION  1
#if RELEASE == 1
#define LOG(fmt, rest...) 
#else
#define LOG(fmt, rest...) \
	printf( fmt, ## rest) 
#endif

#define DBG_MSG 0
#if DBG_MSG == 1
#define DBG(fmt, rest...) \
	printf( fmt, ## rest) 
#else
#define DBG(fmt, rest...) 
#endif


#endif
