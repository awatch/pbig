

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <thrust/version.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include <helper_cuda.h>
#include <helper_timer.h>

#include "rect.h"
#include "report.h"
#include "cuda_comm.h"
#include "pbig_internal.h"
#include "pbig.h"

struct report *_report;

void allocMemForInitAndRects();

void allocMemForGrid();
void deallocMemForGrid();

void allocMemForGridRecs();
void deallocMemForGridRecs();

void allocMemForOverlap();
void allocMemForJob();

#define USE_ORDER 0

#define HEADER_LEN 4

#if DIMS == 2
typedef int2 GridIdx;
#elif DIMS == 3
typedef int3 GridIdx;
#endif

int MAX_Grid_Recs;
int MAX_OVERLAP_DATA;

float time_copy_out = 0.0;
float time_decode = 0.0;
float time_check = 0.0;
float time_order = 0.0;
float total_cuda_time = 0.0;

typedef unsigned uint;

StopWatchInterface *hTimerCopyOut = NULL;
StopWatchInterface *hTimerCheck = NULL;
StopWatchInterface *hTimerOrder = NULL;

unsigned int _assignCs = 0;
unsigned int _NumObjs = 0;
unsigned int _cellSize = 0;
unsigned int _maxRange = 0;
unsigned int _hasProfSteps = 0;

//extern unsigned int RegsUsed;

float4 f4;

unsigned int cuDeviceID = 0;

unsigned int nRowCells = 0;
unsigned int nCells = 0;

unsigned long long int *h_sum;
unsigned int *h_log;
unsigned long long int *d_sum;
unsigned int *d_log;

unsigned long long int numOverlaps;
unsigned long long int numCompPairs;
unsigned int max_nlaps = 0;

// all region information
RectInfo *d_rectInfoList;
RectInfo *h_rectInfoList;

struct CoordInfo *d_orderRectCoordInfoList;


unsigned int *h_hasProfSteps;
unsigned int *d_hasProfSteps;


// overlapped information
unsigned int *h_overlapDataCnt;
unsigned int *d_overlapDataCnt;
OTYPE *h_dst_overlap;
OTYPE *d_dst_overlap;

unsigned int *d_Grid_NRs;
unsigned int *d_Grid_Merged;
unsigned int *h_Grid_NRs;
unsigned int *h_Grid_Merged;
unsigned int *d_Grid_StartI;
unsigned int *h_Grid_StartI;
unsigned int *d_Grid_Recs;
unsigned int *h_Grid_Recs;

cudaArray* cuMergedArray;
texture<unsigned int, cudaTextureType1D, cudaReadModeElementType> texMergedRef;


struct Job {

	unsigned int gid_s;
	unsigned int gid_e;
	unsigned int seg_s;
	unsigned int seg_e;
};

struct Job *h_jobs;
struct Job *d_jobs;

cudaEvent_t e_start, e_stop;
float elapsedTime;

void cuda_finalize() {

    sdkStopTimer(&hTimerCopyOut);
    sdkStopTimer(&hTimerCheck);
    sdkStopTimer(&hTimerOrder);

    cudaEventDestroy(e_start);
    cudaEventDestroy(e_stop);
	cudaThreadExit();
	LOG("Close CUDA system\n");
}

void pbig_finalize() {
	cuda_finalize();
}




#if 1
__global__ void d_clearData(unsigned int nCells, unsigned int *d_Grid_NRs, unsigned int *d_Grid_StartI,unsigned int *d_Grid_Merged, unsigned int *d_overlapDataCnt, unsigned int *log)
{

	int tid = blockIdx.x * blockDim.x + threadIdx.x;


	if ( tid < nCells ) {
		d_Grid_NRs[tid] = 0;
		d_Grid_StartI[tid] = 0;
		d_Grid_Merged[tid] = UNUSED;
		__threadfence();
	}

#if 1
    if ( !threadIdx.x && !blockIdx.x ) {
		*d_overlapDataCnt = 0;
		__threadfence();
	}
#endif

	
}
#endif


__global__ void d_precounting( unsigned int nIter, RectInfo *allRectInfoList, unsigned int *d_Grid_NRs, unsigned int *d_Grid_Merged, int _NumObjs, unsigned int _cellSize, unsigned int shift, unsigned int nCells, unsigned int nRowCells, unsigned int *log) {

	unsigned int rid;
	RectInfo rectInfo;;


    rid = (blockIdx.x+nIter*MAXBLOCKS) * blockDim.x + threadIdx.x;
	if ( rid < _NumObjs ) {
		
		rectInfo = allRectInfoList[rid];

#if DIMS == 2
		unsigned int f_l = ((unsigned int)rectInfo.coordinate.xmin) >> shift;
		unsigned int f_u = ((unsigned int)rectInfo.coordinate.xmax) >> shift;
		unsigned int s_l = ((unsigned int)rectInfo.coordinate.ymin) >> shift;
		unsigned int s_u = ((unsigned int)rectInfo.coordinate.ymax) >> shift;
		for(int j=s_l;j<=s_u;j++) {
			unsigned int product = nRowCells * j;
			for(int i=f_l;i<=f_u;i++) {
				unsigned int gid = product + i;
				if ( gid >= nCells )
					continue;
				atomicAdd( &d_Grid_NRs[gid], 1);
			}
		}
#elif DIMS == 3
		unsigned int f_l = ((unsigned int)rectInfo.coordinate.xmin) >> shift;
		unsigned int f_u = ((unsigned int)rectInfo.coordinate.xmax) >> shift;
		unsigned int s_l = ((unsigned int)rectInfo.coordinate.ymin) >> shift;
		unsigned int s_u = ((unsigned int)rectInfo.coordinate.ymax) >> shift;
		unsigned int t_l = ((unsigned int)rectInfo.coordinate.zmin) >> shift;
		unsigned int t_u = ((unsigned int)rectInfo.coordinate.zmax) >> shift;
		for(int k=t_l;k<=t_u;k++) {
			unsigned int product = nRowCells * nRowCells * k;
			for(int j=s_l;j<=s_u;j++) {
				unsigned int m_gid = product + nRowCells * j;
				for(int i=f_l;i<=f_u;i++) {
					unsigned int gid = m_gid + i;
					if ( gid >= nCells )
						continue;
					atomicAdd( &d_Grid_NRs[gid], 1);
				}
			}
		}
#endif

	}
}

__global__ void d_postcounting( unsigned int nIter, RectInfo *allRectInfoList, unsigned int *d_Grid_NRs, unsigned int *d_Grid_Merged, int _NumObjs, unsigned int _cellSize, unsigned int shift, unsigned int nCells, unsigned int nRowCells, unsigned int *log) {

	unsigned int rid;
	RectInfo rectInfo;;
	unsigned int nRowCells2 = nRowCells * nRowCells;


    rid = (blockIdx.x+nIter*MAXBLOCKS) * blockDim.x + threadIdx.x;
	if ( rid < _NumObjs ) {
		
		rectInfo = allRectInfoList[rid];

#if DIMS == 2
		unsigned int f_l = ((unsigned int)rectInfo.coordinate.xmin) >> shift;
		unsigned int f_u = ((unsigned int)rectInfo.coordinate.xmax) >> shift;
		unsigned int s_l = ((unsigned int)rectInfo.coordinate.ymin) >> shift;
		unsigned int s_u = ((unsigned int)rectInfo.coordinate.ymax) >> shift;

		int ymin = -1;
		for(int j=s_l;j<=s_u;j++) {
			unsigned int product = nRowCells * j;
			int xmin = -1;
			//int t_ymin = -1;
			for(int i=f_l;i<=f_u;i++) {
				unsigned int gid = product + i;
				if ( gid >= nCells ) {
					continue;
				}

				if ( d_Grid_Merged[gid] != UNUSED ) {
					uint mgid = d_Grid_Merged[gid]; 
					GridIdx gridIdx;
					gridIdx.x = mgid % nRowCells;
					gridIdx.y = mgid / nRowCells;
					//t_ymin = max( t_ymin, gridIdx.y);
					if ( gridIdx.x > xmin && gridIdx.y > ymin ) {
						gid = mgid; //d_Grid_Merged[gid];
					}
					else
						continue;
				}
				xmin = i;

				atomicAdd( &d_Grid_NRs[gid], 1);
			}
			//if ( t_ymin > ymin )
			//	ymin = t_ymin;
			ymin = j;
		}
#elif DIMS == 3
		unsigned int f_l = ((unsigned int)rectInfo.coordinate.xmin) >> shift;
		unsigned int f_u = ((unsigned int)rectInfo.coordinate.xmax) >> shift;
		unsigned int s_l = ((unsigned int)rectInfo.coordinate.ymin) >> shift;
		unsigned int s_u = ((unsigned int)rectInfo.coordinate.ymax) >> shift;
		unsigned int t_l = ((unsigned int)rectInfo.coordinate.zmin) >> shift;
		unsigned int t_u = ((unsigned int)rectInfo.coordinate.zmax) >> shift;

		int zmin = -1;
		for(int k=t_l;k<=t_u;k++) {
			unsigned int product = nRowCells * nRowCells * k;
			int ymin = -1;
			for(int j=s_l;j<=s_u;j++) {
				unsigned int m_gid = product + nRowCells * j;
				int xmin = -1;
				for(int i=f_l;i<=f_u;i++) {
					unsigned int gid = m_gid + i;
					if ( gid >= nCells ) {
						continue;
					}

					if ( d_Grid_Merged[gid] != UNUSED ) {
						uint mgid = d_Grid_Merged[gid]; 
						GridIdx gridIdx;
						gridIdx.x = mgid % nRowCells;
						gridIdx.y = (mgid % nRowCells2) / nRowCells;
						gridIdx.z = mgid / nRowCells2;
						if ( gridIdx.x > xmin && gridIdx.y > ymin && gridIdx.z > zmin ) {
							gid = mgid; //d_Grid_Merged[gid];
						}
						else
							continue;

					}
					xmin = i;
					atomicAdd( &d_Grid_NRs[gid], 1);
				}
				ymin = j;
			}
			zmin = k;
		}
#endif


	}
			
}


__global__ void d_mapping( unsigned int nIter, RectInfo *allRectInfoList,unsigned int *d_Grid_StartI, unsigned int *d_Grid_Merged, unsigned int *d_Grid_Recs,int _NumObjs, unsigned int _cellSize, unsigned int shift, unsigned int nCells, unsigned int nRowCells, unsigned int *log) {

	unsigned int tid, rid;
	RectInfo rectInfo;;
	unsigned int nRowCells2 = nRowCells * nRowCells;


    tid = (blockIdx.x+nIter*MAXBLOCKS) * blockDim.x + threadIdx.x;
	if ( tid < _NumObjs ) {
		
		rid = allRectInfoList[tid].rid;
		rectInfo = allRectInfoList[tid];

#if DIMS == 2
		unsigned int f_l = ((unsigned int)rectInfo.coordinate.xmin) >> shift;
		unsigned int f_u = ((unsigned int)rectInfo.coordinate.xmax) >> shift;
		unsigned int s_l = ((unsigned int)rectInfo.coordinate.ymin) >> shift;
		unsigned int s_u = ((unsigned int)rectInfo.coordinate.ymax) >> shift;

		int ymin = -1;
		for(int j=s_l;j<=s_u;j++) {
			unsigned int product = nRowCells * j;
			int xmin = -1;
			//int t_ymin = -1;
			for(int i=f_l;i<=f_u;i++) {
				unsigned int gid = product + i;
				if ( gid >= nCells )
					continue;


#if MERGED == 1
				if ( d_Grid_Merged[gid] != UNUSED ) {
					uint mgid = d_Grid_Merged[gid]; 
					GridIdx gridIdx;
					gridIdx.x = mgid % nRowCells;
					gridIdx.y = mgid / nRowCells;
					//t_ymin = max( gridIdx.y, t_ymin);
					if ( gridIdx.x > xmin && gridIdx.y > ymin ) {
						gid = mgid; //d_Grid_Merged[gid];
					}
					else
						continue;
				}
				xmin = i;
#endif
				unsigned int idx = atomicAdd( &d_Grid_StartI[gid], 1);


				d_Grid_Recs[idx] = rid;
			}
			//if ( t_ymin > ymin )
			//	ymin = t_ymin;
			ymin = j;
		}

#elif DIMS == 3
		unsigned int f_l = ((unsigned int)rectInfo.coordinate.xmin) >> shift;
		unsigned int f_u = ((unsigned int)rectInfo.coordinate.xmax) >> shift;
		unsigned int s_l = ((unsigned int)rectInfo.coordinate.ymin) >> shift;
		unsigned int s_u = ((unsigned int)rectInfo.coordinate.ymax) >> shift;
		unsigned int t_l = ((unsigned int)rectInfo.coordinate.zmin) >> shift;
		unsigned int t_u = ((unsigned int)rectInfo.coordinate.zmax) >> shift;

		int zmin = -1;
		for(int k=t_l;k<=t_u;k++) {
			unsigned int product = nRowCells2 * k;
			int ymin = -1;
			for(int j=s_l;j<=s_u;j++) {
				unsigned int m_gid = product + nRowCells * j;
#if MERGED == 1
				int xmin = -1;
#endif
				for(int i=f_l;i<=f_u;i++) {
					unsigned int gid = m_gid + i;
					if ( gid >= nCells )
						continue;

#if MERGED == 1
					if ( d_Grid_Merged[gid] != UNUSED ) {
						uint mgid = d_Grid_Merged[gid]; 
						GridIdx gridIdx;
						gridIdx.x = mgid % nRowCells;
						gridIdx.y = (mgid % nRowCells2) / nRowCells;
						gridIdx.z = mgid / nRowCells2;
						if ( gridIdx.z > zmin && gridIdx.y > ymin && gridIdx.x > xmin ) {
							gid = mgid;
						}
						else
							continue;

					}
					xmin = i;

#endif
					unsigned int idx = atomicAdd( &d_Grid_StartI[gid], 1);

					d_Grid_Recs[idx] = rid;
				}
				ymin = j; 
			}
			zmin = k;
		}
#endif
	}


}

__device__ uint warp_scan(uint *ptr, uint idx)
{

	const uint lane = idx & 31;

	if ( lane >= 1  ) ptr[idx] = ptr[idx- 1] + ptr[idx];
	if ( lane >= 2  ) ptr[idx] = ptr[idx- 2] + ptr[idx];
	if ( lane >= 4  ) ptr[idx] = ptr[idx- 4] + ptr[idx];
	if ( lane >= 8  ) ptr[idx] = ptr[idx- 8] + ptr[idx];
	if ( lane >= 16 ) ptr[idx] = ptr[idx-16] + ptr[idx];

	//return  ptr[idx];
	return (lane>0) ? ptr[idx-1] : 0;
}

__device__ void block_scan(uint *ptr, uint idx)
{

	const uint lane = idx & 31;
	const uint warpid = idx >> 5;

	//  1.
	uint val = warp_scan( ptr, idx);
	__syncthreads();

	// 2. 
	if ( lane == 31 ) 
		ptr[warpid] = ptr[idx];
	__syncthreads();

	// 3.
	if ( warpid == 0 )
		warp_scan( ptr, idx);
	__syncthreads();

	// 4. 
	if ( warpid > 0 )
		val = ptr[warpid-1] + val;
	__syncthreads();

	// 5. 
	ptr[idx] = val;
	__syncthreads();
}

		

__device__ void prefix_scan(unsigned int hasProfSteps, unsigned int gid, unsigned int cnt, OTYPE *overlaps, const unsigned int MAX_OVERLAP_DATA, OTYPE *d_dst_overlap, unsigned int *d_overlapDataCnt, unsigned int *log) {


	unsigned int q;
	__shared__ unsigned int g_cnt;
	__shared__ unsigned int sh_idx[NTPB];
	__shared__ unsigned int g_idx;
	unsigned int thid;
	unsigned int last = NTPB - 1;


	// 6 op
	thid = threadIdx.x;                 // 1
	
#if 1
	sh_idx[thid] = cnt;
	if ( thid == last ) {               // 2
		g_idx = 0;
	}
	__syncthreads();
	block_scan( sh_idx, thid);
	if ( thid == last )
		g_cnt = sh_idx[thid]+cnt;
	__syncthreads(); // need this since the following code will change the sh_idx[last]
	if ( g_cnt == 0 )
		return;

	if ( thid == last ) {
		g_cnt = sh_idx[thid] + cnt + 2; // unused field + gid field
		g_idx = atomicAdd( d_overlapDataCnt, g_cnt);
		// TODO
		if ( g_idx > MAX_OVERLAP_DATA ) {
			// out of memory
			// TODO
			g_idx = 0;
		}
	}
#endif




	// local barrier
	__syncthreads();
	// 5 op
	unsigned int local_idx = sh_idx[thid] + g_idx + 2; 

	ASSERT(local_idx)
	B_PROFILING(6)

		// 1 op
		if ( thid == last ) 
		{
			// 4 op + 2 gst + 2 gld
			d_dst_overlap[g_idx] = UNUSED;
			d_dst_overlap[g_idx+1] = gid;
		}

		// 2 op
		for(q=0;q<cnt;q++) {
			// cnt * 3 op + cnt * gld + cnt * gst
			d_dst_overlap[local_idx+q] = overlaps[q];
		}

	E_PROFILING(6)


}


// TODO: for 1D texture, len cannot be greater than 32768
void initTexTable(cudaArray *cuMergedArray, unsigned int* mergedTable, unsigned int len)
{
	cudaChannelFormatDesc channelDesc_dist = cudaCreateChannelDesc(  32, 0, 0, 0, cudaChannelFormatKindUnsigned);

	checkCudaErrors( cudaMallocArray( &cuMergedArray, &channelDesc_dist, len, 0) );

	cudaMemcpyToArray( cuMergedArray, 0, 0, mergedTable, len*sizeof(unsigned int), cudaMemcpyHostToDevice);

	texMergedRef.normalized = 0;
	texMergedRef.filterMode = cudaFilterModePoint;
	texMergedRef.addressMode[0] = cudaAddressModeWrap;

	cudaBindTextureToArray( texMergedRef, cuMergedArray,  channelDesc_dist);	
}

__device__ __inline__ unsigned int calcOverlaps( CoordInfo c1, CoordInfo c2 ) {

#if DIMS == 2
	return	 (c1.xmin <= c2.xmax) + 
			 (c2.xmin <= c1.xmax) + 
			 (c1.ymin <= c2.ymax) +
			 (c2.ymin <= c1.ymax) ;
#elif DIMS == 3 
	return	 (c1.xmin <= c2.xmax) + 
			 (c2.xmin <= c1.xmax) + 
			 (c1.ymin <= c2.ymax) + 
			 (c2.ymin <= c1.ymax) + 
			 (c1.zmin <= c2.zmax) + 
			 (c2.zmin <= c1.zmax) ;
#endif
}

__device__ __inline__ unsigned int calcRightToWrite(unsigned int *d_Grid_Merged, unsigned int gid, CoordInfo c1, CoordInfo c2, const unsigned int nRowCells, unsigned int _cellSize, unsigned int shift ) {

	unsigned int first_cid;
#if DIMS == 2
	unsigned int f_l = (unsigned int) max( c1.xmin, c2.xmin) >> shift;
	unsigned int s_l = (unsigned int) max( c1.ymin, c2.ymin) >> shift;
	first_cid = nRowCells * s_l + f_l;                 
#elif DIMS == 3
	unsigned int f_l = (unsigned int) max( c1.xmin, c2.xmin) >> shift;
	unsigned int s_l = (unsigned int) max( c1.ymin, c2.ymin) >> shift;
	unsigned int t_l = (unsigned int) max( c1.zmin, c2.zmin) >> shift;
	first_cid = nRowCells * nRowCells * t_l + nRowCells * s_l + f_l;                    // 2
#endif

#if MERGED == 1
	if ( d_Grid_Merged[first_cid] != UNUSED ) {
		first_cid = d_Grid_Merged[first_cid];
	}
	//if ( tex1D( texMergedRef, first_cid) != UNUSED ) {
	//	first_cid = tex1D( texMergedRef, first_cid);
	//}
#endif

	return (first_cid == gid);
}

__global__ void d_ordering(RectInfo *allRectInfoList, struct CoordInfo *orderRectCoordInfoList, unsigned int *d_Grid_Recs, unsigned int kstart, unsigned int kend, unsigned int *log) {
//__global__ void d_ordering(RectInfo *allRectInfoList, uint *xminList, uint *xmaxList, uint *yminList, uint *ymaxList, uint *zminList, uint *zmaxList, unsigned int *d_Grid_Recs, unsigned int kstart, unsigned int kend, unsigned int *log) {

	CoordInfo coord;
	__shared__ uint data[2048];

	uint tid = blockDim.x * blockIdx.x + threadIdx.x;
	if ( tid + kstart < kend ) {

		uint rid = d_Grid_Recs[kstart+tid];
		orderRectCoordInfoList[tid] = allRectInfoList[rid].coordinate;
		//coord = allRectInfoList[rid].coordinate;
		//xminList[tid] = coord.xmin;
		//xmaxList[tid] = coord.xmax;
		//yminList[tid] = coord.ymin;
		//ymaxList[tid] = coord.ymax;
		//zminList[tid] = coord.zmin;
		//zmaxList[tid] = coord.zmax;

	}
	
	if ( threadIdx.x > NTPB )  {
		data[tid] = coord.xmin;
		log[0] = data[tid];
	}
		
}


#if DIMS == 2
#define MATCHED_ANS 5
#elif DIMS == 3
#define MATCHED_ANS 7
#endif
#define Incr(cur,bound,incr) ( (((cur) + (incr)) >= (bound)) ? (bound) - (cur) : (incr) )

extern __shared__ unsigned int sh_rect[];
__global__ void d_matching( unsigned int *d_hasProfSteps, RectInfo *allRectInfoList, struct CoordInfo *orderRectCoordInfoList, unsigned int *d_Grid_NRs, unsigned int *d_Grid_StartI,unsigned int *d_Grid_Merged, unsigned int *d_Grid_Recs, unsigned int globalJobOffset, unsigned int kstart, struct Job *jobs, unsigned int _cellSize, unsigned int shift, const unsigned int nRowCells, const unsigned int MAX_OVERLAP_DATA, OTYPE *d_dst_overlap, unsigned int *d_overlapDataCnt, unsigned int *log) {

	CoordInfo rect1Coord, rect2Coord;

	unsigned int rid1, rid2;

	unsigned int area_type = 0;
	unsigned int area_x, area_y;
	unsigned int nRects_x, nRects_y;

	unsigned int hasProfSteps;
	unsigned int overcnt;
#ifdef COMPRESSION
	OTYPE loc_overlaps[NTPB];
	unsigned char *overlaps = (unsigned char *)loc_overlaps;
	unsigned short *overlaps2 = (unsigned short*)loc_overlaps;
#else
	unsigned int overlaps[NTPB*2];
#endif

	__shared__ struct Job job;


	__shared__ unsigned int bid;
	unsigned int numLocalRects, firstRectIdx;
	unsigned int encode_unit;
	unsigned int gid;

	CoordInfo *shCoords;
	unsigned int *shIds;

	// 2 op + 3 gld
	if ( threadIdx.x == 0 ) {
		bid = blockIdx.x + globalJobOffset;
		job = jobs[bid];
	}

	hasProfSteps = d_hasProfSteps[0];
	shCoords = (CoordInfo*) sh_rect;
#if DIMS == 2
	shIds = sh_rect + 4 * NTPB;
#elif DIMS == 3
	shIds = sh_rect + 6 * NTPB;
#endif

	// local barrier
	__syncthreads();

	unsigned int job_cell_start = job.seg_s;
	unsigned int job_cell_end   = job.seg_e;
	unsigned int start_gid  = job.gid_s;
	unsigned int end_gid    = job.gid_e;

	area_y = job_cell_start / NTPB;

	//ASSERT(job_cell_start+start_gid+end_gid+useless[threadIdx.x])
	ASSERT(job_cell_start+start_gid+end_gid)

	B_PROFILING(1)
	// 2 op
	//for(gid=start_gid;gid<=end_gid;gid++) { 
	gid = start_gid;

		// 5 op + 2 gld
		numLocalRects = d_Grid_NRs[gid];
		firstRectIdx = d_Grid_StartI[gid];
		nRects_y = Incr(job_cell_start,numLocalRects,NTPB);
#ifdef RICE_COMPRESSION
		encode_unit  = 1;
#else
		encode_unit  = ( numLocalRects > 256 ) + 1;
#endif


		// 2 op
		if ( job_cell_start + threadIdx.x < numLocalRects ) {
			// 7 op + 2 gld with random access
#if USE_ORDER == 1
			unsigned int assign_loc = firstRectIdx+job_cell_start+threadIdx.x;
			rid1 = d_Grid_Recs[assign_loc];
			uint pos = assign_loc - kstart;
			rect1Coord = orderRectCoordInfoList[pos];
			//rect1Coord.xmin = xminList[pos];
			//rect1Coord.xmax = xmaxList[pos];
			//rect1Coord.ymin = yminList[pos];
			//rect1Coord.ymax = ymaxList[pos];
			//rect1Coord.zmin = zminList[pos];
			//rect1Coord.zmax = zmaxList[pos];


			//if ( rect1Coord.xmin != allRectInfoList[rid1].coordinate.xmin )
			//	atomicAdd( &log[0], 1);
#else
			unsigned int assign_loc = firstRectIdx+job_cell_start+threadIdx.x;
			rid1 = d_Grid_Recs[assign_loc];
			rect1Coord = allRectInfoList[rid1].coordinate;
#endif

		}

		ASSERT(rid1+rect1Coord.xmin+rect1Coord.xmax+rect1Coord.ymin+rect1Coord.ymax)
		B_PROFILING(2)
		// 2 op
		for(unsigned int i=0;i<numLocalRects;i+=NTPB) {
#if BALANCE == 1
			area_x = i / NTPB;
			if ( area_x == area_y )
				area_type = 1;
			else if ( (area_y > area_x && ((area_x+area_y)&1)==1) || (area_y < area_x && ((area_x+area_y)&1)==0) )
				area_type = 2;
			else
				continue;
			nRects_x = Incr(i,numLocalRects,NTPB);
#endif
			// 2 op
			if ( threadIdx.x + i < numLocalRects ) {
				// 7 op + 2 gld
				shIds[threadIdx.x] = d_Grid_Recs[firstRectIdx+threadIdx.x+i];
#if USE_ORDER == 1
				uint pos = firstRectIdx+threadIdx.x+i-kstart;
				shCoords[threadIdx.x] = orderRectCoordInfoList[pos];
				//shCoords[threadIdx.x].xmin = xminList[pos];
				//shCoords[threadIdx.x].xmax = xmaxList[pos];
				//shCoords[threadIdx.x].ymin = yminList[pos];
				//shCoords[threadIdx.x].ymax = ymaxList[pos];
				//shCoords[threadIdx.x].zmin = zminList[pos];
				//shCoords[threadIdx.x].zmax = zmaxList[pos];
				//if ( shCoords[threadIdx.x].xmin != allRectInfoList[shIds[threadIdx.x]].coordinate.xmin )
				//	atomicAdd( &log[0], 1);
#else
				shCoords[threadIdx.x] = allRectInfoList[shIds[threadIdx.x]].coordinate;
#endif
			}

			ASSERT(shIds[threadIdx.x]+shCoords[threadIdx.x].xmax)

			// local barrier
			__syncthreads();


			B_PROFILING(3)

			unsigned int seg_size = min( NTPB, 256);
			unsigned int seg_size_shift = __log2f(seg_size);
			for(unsigned int seg=0;seg<NTPB;seg+=seg_size) { 

#ifdef COMPRESSION
#ifdef RICE_COMPRESSION
				overcnt = 6;
#else
				if ( encode_unit == 2 )
					overcnt = 2;
				else
					overcnt = 4;
#endif
#else
				overcnt = 2;
#endif


#if BALANCE == 1
				int srange, nElm;
				if ( area_type == 1 ) {
					srange = threadIdx.x + 1 + seg;
					if ( threadIdx.x < (nRects_x>>1) )
						nElm = Incr(seg,(nRects_x>>1),seg_size);
					else
						nElm = Incr(seg,((nRects_x-1)>>1),seg_size);
				}
				else if ( area_type == 2 ) {
					srange = seg;
					nElm = Incr(seg,nRects_x,seg_size);
				}
#endif

				if ( job_cell_start + threadIdx.x >= numLocalRects )
					goto SYNC2;

				unsigned int k, proc_elm;
#if BALANCE == 1
				for(k=srange,proc_elm=0;proc_elm<nElm;k++,proc_elm++) {

					if ( k >= nRects_x )
						k = 0;
#else
				int srange = threadIdx.x > (i+seg) % NTPB ? threadIdx.x : (i+seg) % NTPB;
				int erange = i+seg_size+seg > numLocalRects ? numLocalRects % (NTPB) : seg_size+seg;
				for(k=srange;k<erange;k++) {
#endif

					// 4 op
					rid2 = shIds[k];
					rect2Coord = shCoords[k];

					unsigned int results, hasRight;
					results = calcOverlaps( rect1Coord, rect2Coord);
					hasRight = calcRightToWrite( d_Grid_Merged, gid, rect1Coord, rect2Coord, nRowCells, _cellSize, shift);
#if BALANCE == 0
					results -= ( k == threadIdx.x && rid1 <= rid2 );
#endif
					results += hasRight;                                    // 2
						
					ASSERT(results)
					// 1 op
					if ( results==MATCHED_ANS ) 
					{
						B_PROFILING(4)
#ifdef COMPRESSION
#ifdef RICE_COMPRESSION
						overlaps[overcnt++] = (i+k) & ( seg_size - 1 );
#else
						// 4 op + 1 gst + 1 gld
						if ( encode_unit == 2 ) {
							overlaps2[overcnt++] = i+k;
						}
						else 
							overlaps[overcnt++] = i+k;
#endif
#else
						// 4 op + 2 gst + 2 gld
						overlaps[overcnt++] = rid1;
						overlaps[overcnt++] = rid2;
#endif
						E_PROFILING(4)
					}
				}

SYNC2:
				// 2 ops + local barrier
				if ( i+NTPB < numLocalRects ) 
					__syncthreads();

				ASSERT(overlaps[overcnt-1]) 
				B_PROFILING(5)


#ifdef COMPRESSION
#ifdef RICE_COMPRESSION
				if ( overcnt == 6 ) 
					overcnt = 0;
				else {
					overlaps[0] = (job_cell_start+threadIdx.x) >> seg_size_shift;
					overlaps[1] = (i+seg) >> seg_size_shift;
					overlaps2[1] = threadIdx.x & ( seg_size - 1 );
					overlaps2[2] = overcnt-6;
				}
				prefix_scan( hasProfSteps, gid, (overcnt+sizeof(OTYPE)-1)>>OTYPE_SHIFT, loc_overlaps, MAX_OVERLAP_DATA, d_dst_overlap, d_overlapDataCnt, log);
#else
				// 5 op + 1 gst + 1 gld + prefix_scan
				if ( encode_unit == 2 ) {
					if ( overcnt == 2 ) 
						overcnt = 0;
					else {
						overlaps2[0] = job_cell_start+threadIdx.x;
						overlaps2[1] = overcnt-2;
					}
					prefix_scan( hasProfSteps, gid, (overcnt+1)>>1, loc_overlaps, MAX_OVERLAP_DATA, d_dst_overlap, d_overlapDataCnt, log);
				}
				else {
					if ( overcnt == 4 ) 
						overcnt = 0;
					else {
						overlaps[0] = job_cell_start+threadIdx.x;
						overlaps[1] = overcnt-4;
					}
					prefix_scan( hasProfSteps, gid, (overcnt+3)>>2, loc_overlaps, MAX_OVERLAP_DATA, d_dst_overlap, d_overlapDataCnt, log);
				}
#endif
#else
				// 2 op + 1 gst + 1 gld + prefix_scan
				if ( overcnt == 2 ) {
					overcnt = 0;
				}
				else {
					overlaps[0] = rid1;
					overlaps[1] = overcnt-2;
				}
				prefix_scan( hasProfSteps, gid, overcnt, overlaps, MAX_OVERLAP_DATA, d_dst_overlap, d_overlapDataCnt, log);
#endif
				overcnt = 0;

				E_PROFILING(5)
			}
			E_PROFILING(3)

		}
		E_PROFILING(2);
	//}
	E_PROFILING(1);

}

#if 0
void do_sth(unsigned int rid1, unsigned int rid2)
{
#if VERIFY == 1
	*h_sum += min( rid1, rid2);
	//*h_sum = *h_sum + 1;
#endif
}
#endif

#if 1
void printRec(unsigned int i) {

#if DIMS == 2
	LOG("rid %7d [%5d,%5d] [%5d,%5d]\n", i, 
										   h_rectInfoList[i].coordinate.xmin, 	
										   h_rectInfoList[i].coordinate.xmax, 	
										   h_rectInfoList[i].coordinate.ymin, 	
										   h_rectInfoList[i].coordinate.ymax);
#elif DIMS == 3
	LOG("rid %7d [%5d,%5d] [%5d,%5d]\n", i, 
										   h_rectInfoList[i].coordinate.xmin, 	
										   h_rectInfoList[i].coordinate.xmax, 	
										   h_rectInfoList[i].coordinate.ymin, 	
										   h_rectInfoList[i].coordinate.ymax,
										   h_rectInfoList[i].coordinate.zmin, 	
										   h_rectInfoList[i].coordinate.zmax);
#endif

}
#endif

#if 1
typedef pair<unsigned int, unsigned int> opair;
set<opair> overlapPairs; 

bool __check( unsigned int k, unsigned int q) {
	
	if ( k < q ) {
		if ( overlapPairs.find(make_pair(k,q)) != overlapPairs.end() )
			return false;
		overlapPairs.insert(make_pair(k,q));
		//LOG("Pair: %8d %8d\n", k, q);
	}
	else {
		if ( overlapPairs.find(make_pair(q,k)) != overlapPairs.end()  ) 
			return false;
		overlapPairs.insert(make_pair(q,k));
		//LOG("Pair: %8d %8d\n", q, k);
	}

	return true;
}

bool checkOverlaps( unsigned int k, unsigned int q) {

	bool isOverlap = false;
	
	//if ( k >= MAXALLOCRECTS || q >= MAXALLOCRECTS ) {
	//	CUDA_ERR("Error: in rid%d and rid%d\n", k, q);
	//}

#if DIMS == 2
	if ( 
		 (h_rectInfoList[k].coordinate.xmin <= h_rectInfoList[q].coordinate.xmax) && 
		 (h_rectInfoList[q].coordinate.xmin <= h_rectInfoList[k].coordinate.xmax) && 
		 (h_rectInfoList[k].coordinate.ymin <= h_rectInfoList[q].coordinate.ymax) && 
		 (h_rectInfoList[q].coordinate.ymin <= h_rectInfoList[k].coordinate.ymax) ) {
		isOverlap = __check(k,q);
	}
#elif DIMS == 3
	if (
		 (h_rectInfoList[k].coordinate.xmin <= h_rectInfoList[q].coordinate.xmax) && 
		 (h_rectInfoList[q].coordinate.xmin <= h_rectInfoList[k].coordinate.xmax) && 
		 (h_rectInfoList[k].coordinate.ymin <= h_rectInfoList[q].coordinate.ymax) && 
		 (h_rectInfoList[q].coordinate.ymin <= h_rectInfoList[k].coordinate.ymax) && 
		 (h_rectInfoList[k].coordinate.zmin <= h_rectInfoList[q].coordinate.zmax) && 
		 (h_rectInfoList[q].coordinate.zmin <= h_rectInfoList[k].coordinate.zmax) ) {
		isOverlap = __check(k,q);
	}
#endif

	if ( isOverlap == false ) {
		printRec(k);
		printRec(q);
		LOG("error in calculting overlaps\n");
	}

	return isOverlap;
}
#endif

//void _report(unsigned int a, unsigned int b)
//{
//	(*h_sum)++;
//}




unsigned int check_boundary(unsigned int dir, unsigned int gid)
{
#if DIMS == 2
	if ( dir == 0 ) {
		int col = gid % nRowCells;
		return col == 0 ? 0 : 1;
	}
	else if ( dir == 1 ) {
		int row = gid / nRowCells;
		return row == nRowCells ? 0 : 1;
	}	
	else
		return 0;
#elif DIMS == 3
	if ( dir == 0 ) {
		int col = gid % nRowCells;
		return col == 0 ? 0 : 1;
	}
	else if ( dir == 1 ) {
		int row = (gid % ( nRowCells * nRowCells )) / nRowCells;
		return row == 0 ? 0 : 1;
	}	
	else if ( dir == 2 ) {
		int high = gid / ( nRowCells * nRowCells );
		return high == nRowCells ? 0 : 1;
	}	
	else
		return 0;
#endif
}


unsigned int try_merge(unsigned int gid, GridIdx gid_idx, unsigned int dir, unsigned int *h_nr, unsigned int *v_nr, unsigned int *u_nr)
{

#define M_Cond NTPB

	unsigned int acc_nr = h_Grid_NRs[gid];
	unsigned int base_gid;

	DBG("\n=========before try_merge result: gid%d acc_nr %d\n", gid, acc_nr);
	
#if DIMS == 2 
	unsigned int step = 0;
	unsigned int merged_nr;
	// 1. check boundary for all merged grids along with horoizontal- or vertical-direction
	step = dir == 0 ? nRowCells : 1;
	merged_nr = dir == 0 ? *v_nr : *h_nr ;
	if ( dir == 0 )
		base_gid = gid + *h_nr; 
	else
		base_gid = gid + (*v_nr) * nRowCells;
	DBG("try_merge gid%d dir %d step %d merged_nr %d base_gid %d\n", gid, dir, step, merged_nr, base_gid);
	for(int i=0;i<merged_nr;i++) {
		if ( h_Grid_Merged[base_gid+i*step] != UNUSED )
			return 0;
	}

	// 2. check the accumulated number of rectangles < M_Cond
	for(int i=0;i<merged_nr;i++) {
		acc_nr += h_Grid_NRs[base_gid+i*step];
		if ( acc_nr > M_Cond )
			return 0;
	}

	// 3. if satisfy both conditions above, then do merge
	for(int i=0;i<merged_nr;i++) {
		h_Grid_Merged[base_gid+i*step] = gid;
		DBG("\t\tmerging gid%d and its nr %d\n", base_gid+i*step, h_Grid_NRs[base_gid+i*step]);
	}

	// 4. update the left-top gid
	h_Grid_Merged[gid] = gid;
	h_Grid_NRs[gid] = acc_nr;
	if ( dir == 0 )
		(*h_nr)++;
	else if ( dir == 1 )
		(*v_nr)++;


#elif DIMS == 3
	unsigned int step1, step2, merged_nr1, merged_nr2;
	// 1. check boundary for all merged grids along with horoizontal- or vertical-direction
	switch (dir) {
		case 0:
			step1 = nRowCells*nRowCells;
			step2 = nRowCells;
			merged_nr1 = *u_nr; 
			merged_nr2 = *v_nr; 
			break;
		case 1:
			step1 = nRowCells*nRowCells;
			step2 = 1;
			merged_nr1 = *u_nr; 
			merged_nr2 = *h_nr; 
			break;
		case 2:
			step1 = nRowCells;
			step2 = 1;
			merged_nr1 = *v_nr; 
			merged_nr2 = *h_nr; 
			break;
		default:
			exit(-1);
			break;
	}
	if ( dir == 0 )
		base_gid = gid + *h_nr; 
	else if ( dir == 1 )
		base_gid = gid + (*v_nr) * nRowCells;
	else if ( dir == 2 )
		base_gid = gid + (*u_nr) * nRowCells * nRowCells;

	DBG("try_merge gid%d dir %d step1 %d step2 %d merged_nr1 %d merged_nr2 %d base_gid %d\n", gid, dir, step1, step2, merged_nr1, merged_nr2, base_gid);
	for(int i=0;i<merged_nr1;i++) {
		int m_gid = base_gid + i * step1;
		for(int j=0;j<merged_nr2;j++) {
			if ( h_Grid_Merged[m_gid+j*step2] != UNUSED )
				return 0;
		}
	}

	// 2. check the accumulated number of rectangles < M_Cond
	for(int i=0;i<merged_nr1;i++) {
		int m_gid = base_gid + i * step1;
		for(int j=0;j<merged_nr2;j++) {
			acc_nr += h_Grid_NRs[m_gid+j*step2];
			if ( acc_nr > M_Cond )
				return 0;
		}
	}

	// 3. if satisfy both conditions above, then do merge
	for(int i=0;i<merged_nr1;i++) {
		int m_gid = base_gid + i * step1;
		for(int j=0;j<merged_nr2;j++) {
			h_Grid_Merged[m_gid+j*step2] = gid;
			//if ( m_gid+j*step2 >= nCells ) {
			//	CUDA_ERR("ERR: the merged cell %d is out of boundary\n", m_gid+j*step2);
			//}
		}
	}

	// 4. update the left-top gid
	h_Grid_Merged[gid] = gid;
	h_Grid_NRs[gid] = acc_nr;
	if ( dir == 0 )
		(*h_nr)++;
	else if ( dir == 1 )
		(*v_nr)++;
	else if ( dir == 2 )
		(*u_nr)++;
#endif

	return 1;
}

void h_merge_grids()
{
	unsigned int dir = 0;
	
	for(int i=0;i<nCells;i++) {

		unsigned int h_nr = 1; 
		unsigned int v_nr = 1; 
		unsigned int u_nr = 1; 

		// skip a merged grid
		if ( h_Grid_Merged[i] != UNUSED )
			continue;

		GridIdx gid_idx;
#if DIMS == 2
		gid_idx.x = i % nRowCells;
		gid_idx.y = i / nRowCells;
#elif DIMS == 3
		gid_idx.x = i % nRowCells;
		gid_idx.y = (i % (nRowCells*nRowCells)) / nRowCells;
		gid_idx.z = i / ( nRowCells * nRowCells );
#endif

		unsigned int progress;
		do {

			progress = 0;

			// horizontal
			dir = 0;	
			if ( gid_idx.x + h_nr < nRowCells )
				progress = try_merge( i, gid_idx, dir, &h_nr, &v_nr, &u_nr);

			// vertical
			dir = 1;
			if ( gid_idx.y + v_nr < nRowCells )
				progress += try_merge( i, gid_idx, dir, &h_nr, &v_nr, &u_nr);

#if DIMS == 3
			// up
			dir = 2;
			if ( gid_idx.z + u_nr < nRowCells )
				progress += try_merge( i, gid_idx, dir, &h_nr, &v_nr, &u_nr);
#endif
		} while ( progress != 0 );



#if 0
		// for debugging
		if ( h_Grid_Merged[i] != UNUSED && h_Grid_Merged[i] >= nCells ) {
			LOG("Error: in %s\n", __func__);
			cuda_finalize();
			exit(-1);
		}
#endif

	}

}

#if 1
struct avgX : public unary_function<RectInfo,RectInfo>
{
	__host__ __device__ float operator()(const RectInfo &x) const
	{
		return x.coordinate.xmax - x.coordinate.xmin;
	}
};

struct avgY : public unary_function<RectInfo,RectInfo>
{
	__host__ __device__ float operator()(const RectInfo &x) const
	{
		return x.coordinate.ymax - x.coordinate.ymin;
	}
};

#if DIMS == 3
struct avgZ : public unary_function<RectInfo,RectInfo>
{
	__host__ __device__ float operator()(const RectInfo &x) const
	{
		return x.coordinate.zmax - x.coordinate.zmin;
	}
};
#endif

#if DIMS == 2
#define NMR_BOUND 4
#elif DIMS == 3
#define NMR_BOUND 8
#endif
uint decideCellSize()
{
			
	float x,y,z;
	int nSamples;
	int nRecsPerCell = _NumObjs;
	int cs = 1;
	long long int nc; 
	long long int nrc;
	double nmr;

	while ( cs*2 < _maxRange ) {
		cs *= 2;
	}

	nSamples = 1000 > _NumObjs ? _NumObjs : 1000;

	LOG("Sampling ...\n");
	x = thrust::transform_reduce( h_rectInfoList, h_rectInfoList+nSamples, avgX(), 0, thrust::plus<float>()) / nSamples;  
	y = thrust::transform_reduce( h_rectInfoList, h_rectInfoList+nSamples, avgY(), 0, thrust::plus<float>()) / nSamples;  
#if DIMS == 2
	LOG("avgX = %lf, avgY = %lf\n", x, y);
#elif DIMS == 3
	z = thrust::transform_reduce( h_rectInfoList, h_rectInfoList+nSamples, avgZ(), 0, thrust::plus<float>()) / nSamples;  
	LOG("avgX = %lf, avgY = %lf, avgZ = %lf\n", x, y, z);
#endif


	while ( 1 ) {

		nrc = ( _maxRange + cs - 1 ) / cs;
		nc  = nrc * nrc;
		
		nmr = ( 1 + x / cs ) * ( 1 + y / cs );
#if DIMS == 3
		nc *= nrc;
		nmr *= ( 1 + z / cs );
#endif
		nRecsPerCell = nmr * _NumObjs / nc;
		LOG("nrc=%d, nc=%d, nmr=%lf, cs=%d, nRecsPerCell=%d\n", nrc, nc, nmr, cs, nRecsPerCell);
		if ( cs < x )
			break;
		if ( nRecsPerCell < 32 || nmr > NMR_BOUND )
			break;
		if ( cs == 1 ) 
			break;

		nRowCells = (_maxRange + cs -1 ) / cs;

		int isAvail;
#if DIMS == 3 
		IS_AVAIL_MEM( &isAvail, (nRowCells * nRowCells * nRowCells) * sizeof(unsigned int)*3);
#elif DIMS == 2 
		IS_AVAIL_MEM( &isAvail, (nRowCells * nRowCells) * sizeof(unsigned int)*3);
#endif
		if (  isAvail == 0 )
			break;

#if 0
#if DIMS == 3 
		if ( MAXGRIDS < nRowCells * nRowCells * nRowCells )
#elif DIMS == 2 
		if ( MAXGRIDS < nRowCells * nRowCells )
#endif
			break;
#endif

		cs /= 2;
	}

	return cs*2;

}
#endif


static void plan_phase(RectInfo *input)
{


    unsigned int numBlocks;
	total_cuda_time = 0;

	h_rectInfoList = input;	

	CUDA_EVENT_START
		if ( _assignCs == 0 )
			_cellSize = decideCellSize();
		else
			_cellSize = _assignCs;
		nRowCells = (_maxRange + _cellSize -1) / _cellSize;
#if DIMS == 2
		nCells = nRowCells * nRowCells;
#elif DIMS == 3
		nCells = nRowCells * nRowCells * nRowCells;
#endif
	CUDA_EVENT_END
	LOG("host[decide _cellSize]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;

	allocMemForGrid();

	//if ( MAXGRIDS < nCells ) {
	//	LOG("Error: Due to MAXGRIDS (%d) < nCells (%d)\n", MAXGRIDS, nCells);
	//	exit(0);
	//}

	LOG(" _hasProfSteps = %d _cellSize = %d nCells = %d nRowCells %d\n", _hasProfSteps, _cellSize, nCells, nRowCells);
	LOG("%s %d objs NTPB %d NTPB %d\n", __func__, _NumObjs, NTPB, NTPB);

	*h_sum = 0;
	for(int k=0;k<MAXLOG;k++)
		h_log[k] = 0;

	/* copy input data into device mem */
	CUDA_EVENT_START
		cudaMemcpy(  d_rectInfoList,  h_rectInfoList, sizeof(RectInfo)*_NumObjs, cudaMemcpyHostToDevice);	
		cudaMemcpy( d_sum, h_sum, sizeof(unsigned long long int), cudaMemcpyHostToDevice);	
		cudaMemcpy( d_log, h_log, MAXLOG*sizeof(unsigned int), cudaMemcpyHostToDevice);	
	CUDA_EVENT_END
	LOG("kernel[copy in]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;


	for(int i=0;i<nCells;i++) {
		h_Grid_Merged[i] = UNUSED;
	}


	int totalThreads = nCells;	
	numBlocks = ( totalThreads + 512 -1 ) / 512 ;
	LOG("nCells=threads = %d, numBlock = %lld\n", totalThreads, numBlocks);
	CUDA_EVENT_START
		d_clearData<<<numBlocks,512>>>(nCells, d_Grid_NRs, d_Grid_StartI, d_Grid_Merged, d_overlapDataCnt, d_log);
	CUDA_EVENT_END
	CUDA_COPY_LOG
	CHECK_CUDA_ERR
	LOG("kernel[d_clearData]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;

	overlapPairs.clear();
}


static void map_phase()
{
    unsigned long long int numBlocks;
	const unsigned int blockSizeForMapping = 512;

	/*  device: d_precounting */
    numBlocks = ( _NumObjs + blockSizeForMapping -1 ) / blockSizeForMapping;
    LOG("_NumObjs = %d, numBlock = %lld\n", _NumObjs, numBlocks);
	CUDA_EVENT_START
		for(int nIter=0;nIter<(numBlocks+MAXBLOCKS-1)/MAXBLOCKS;nIter++) {
			int nb = MAXBLOCKS*(nIter+1) < numBlocks ? MAXBLOCKS : numBlocks%MAXBLOCKS ; 
			d_precounting<<<nb,blockSizeForMapping>>>(   
														 nIter,
														 d_rectInfoList,
														 d_Grid_NRs,
														 d_Grid_Merged,
														 _NumObjs,
														 _cellSize,
														 (unsigned int) log2((float)_cellSize),
														 nCells,
														 nRowCells,
														 d_log);
		}
		cudaMemcpy( h_Grid_NRs, d_Grid_NRs, sizeof(unsigned int)*nCells, cudaMemcpyDeviceToHost);	
	CUDA_EVENT_END
	CHECK_CUDA_ERR
	CUDA_COPY_LOG
	LOG("device[d_precounting]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;

#if MERGED == 1
	/*  host: h_merge_grids */
	CUDA_EVENT_START
		h_merge_grids();
		cudaMemcpy( d_Grid_Merged, h_Grid_Merged, sizeof(unsigned int)*nCells, cudaMemcpyHostToDevice);	
		cudaMemset( d_Grid_NRs, 0, sizeof(unsigned int)*nCells);
		//initTexTable( cuMergedArray, h_Grid_Merged, nCells);
	CUDA_EVENT_END
	LOG("host[h_merge_grids]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;

	/*  device: d_postcounting */
    numBlocks = ( _NumObjs + blockSizeForMapping - 1 ) / blockSizeForMapping;
    LOG("_NumObjs = %d, numBlock = %lld\n", _NumObjs, numBlocks);
	CUDA_EVENT_START
		for(int nIter=0;nIter<(numBlocks+MAXBLOCKS-1)/MAXBLOCKS;nIter++) {
			int nb = MAXBLOCKS*(nIter+1) < numBlocks ? MAXBLOCKS : numBlocks%MAXBLOCKS ; 
			d_postcounting<<<nb,blockSizeForMapping>>>(   
														 nIter,
														 d_rectInfoList,
														 d_Grid_NRs,
														 d_Grid_Merged,
														 _NumObjs,
														 _cellSize,
														 (unsigned int) log2((float)_cellSize),
														 nCells,
														 nRowCells,
														 d_log);
		}
		cudaMemcpy( h_Grid_NRs, d_Grid_NRs, sizeof(unsigned int)*nCells, cudaMemcpyDeviceToHost);	
	CUDA_EVENT_END
	CHECK_CUDA_ERR
	CUDA_COPY_LOG
	LOG("device[d_postcounting]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;
#endif


	/*  deivce: d_exclusive_scan */
	CUDA_EVENT_START
		thrust::exclusive_scan( h_Grid_NRs, h_Grid_NRs + nCells, h_Grid_StartI);
		cudaMemcpy( d_Grid_StartI, h_Grid_StartI, sizeof(unsigned int)*nCells, cudaMemcpyHostToDevice);	
	CUDA_EVENT_END
	CHECK_CUDA_ERR
	//LOG("cpu time %ld\n", (cuda_end.tv_sec*1000000+cuda_end.tv_usec)-(cuda_start.tv_sec*1000000+cuda_start.tv_usec));
	LOG("device[d_exclusive_scan]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;

	LOG("Last item in d_Grid_Recs is at %d\n", h_Grid_StartI[nCells-1]+h_Grid_NRs[nCells-1]-1);
	MAX_Grid_Recs = h_Grid_StartI[nCells-1]+h_Grid_NRs[nCells-1];
	CUDA_EVENT_START
		allocMemForGridRecs();
	CUDA_EVENT_END
	LOG("device[d_alloc_GridRecs]: elapsed time %lf\n", elapsedTime*1000);
	if ( h_Grid_StartI[nCells-1]+h_Grid_NRs[nCells-1]-1 >= MAX_Grid_Recs ) {
		LOG("Error: the number of rectangles in grids >= MAX_Grid_Recs(%ld)\n", MAX_Grid_Recs); 
	}

	/*  device: d_mapping */
    numBlocks = ( _NumObjs + blockSizeForMapping -1 ) / blockSizeForMapping;
    LOG("_NumObjs = %d, numBlock = %lld\n", _NumObjs, numBlocks);
	CUDA_EVENT_START
		for(int nIter=0;nIter<(numBlocks+MAXBLOCKS-1)/MAXBLOCKS;nIter++) {
			int nb = MAXBLOCKS*(nIter+1) < numBlocks ? MAXBLOCKS : numBlocks%MAXBLOCKS ; 
			d_mapping<<<nb,blockSizeForMapping>>>(   
														 nIter,
														 d_rectInfoList,
														 d_Grid_StartI,
														 d_Grid_Merged,
														 d_Grid_Recs,
														 _NumObjs,
														 _cellSize,
														 (unsigned int) log2((float)_cellSize),
														 nCells,
														 nRowCells,
														 d_log);
		}
		cudaMemcpy( h_Grid_Recs, d_Grid_Recs, sizeof(unsigned int)*(h_Grid_StartI[nCells-1]+h_Grid_NRs[nCells-1]), cudaMemcpyDeviceToHost);	
		cudaMemcpy( d_Grid_StartI, h_Grid_StartI, sizeof(unsigned int)*nCells, cudaMemcpyHostToDevice);	
	CUDA_EVENT_END
	CHECK_CUDA_ERR
	CUDA_COPY_LOG
	LOG("device[d_mapping]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;


	unsigned int avgRecsCell = 0;
	unsigned int maxRecs = 0;
	unsigned int maxRecsCell = 0;
	for(int i=0;i<nCells;i++) {
		avgRecsCell += h_Grid_NRs[i];
		if ( h_Grid_NRs[i] > maxRecs ) {
			maxRecs = h_Grid_NRs[i];
			maxRecsCell = i;
		}
	}
	unsigned int numMappedRecs = 0;
	for(int i=0;i<nCells;i++) {
		numMappedRecs += h_Grid_NRs[i];
	}

#if 0
	LOG("\n\n\n");
	for(int i=0;i<nCells;i++) {

		int nregs = h_Grid_NRs[i];
		int pos = h_Grid_StartI[i];
		for(int j=0;j<nregs;j++) {
			if ( h_Grid_Recs[pos+j] == 4078 || h_Grid_Recs[pos+j] == 19072 ) {
				LOG("rid%d at cell %d\n", h_Grid_Recs[pos+j], i);
			}
		}
	}
	LOG("\n\n\n");
	unsigned int rec_histgram[100] = { 0 };
	for(int i=0;i<100;i++)
		rec_histgram[i] = 0;
	for(int i=0;i<nCells;i++) {
		rec_histgram[h_Grid_NRs[i]*100/maxRecs]++;
	}
	for(int i=0;i<100;i++) {
		LOG("bucket%d = %d\n", i, rec_histgram[i]); 
	}
#endif
	LOG("Sampling AvgRecs Per Cell %.2f, maxRecs %d on cell %d, numMappedRecs %d\n", (float)avgRecsCell/(float)nCells, maxRecs, maxRecsCell, numMappedRecs);
}

unsigned int numBatches = 0;
struct Batch {
	uint job_s, job_e;
};
struct Batch jobBatches[10000];

static unsigned int sched_phase()
{


	numCompPairs = 0;
	numBatches = 0;

	uint job_i = 0;
	uint accRecs = 0;
	uint job_s = 0;

#if 1
	CUDA_EVENT_START
	for(int i=0;i<nCells;i++) {

	
		if ( h_Grid_NRs[i] <= 1 ) { 
			accRecs += h_Grid_NRs[i];
			continue;
		}

		if ( h_Grid_NRs[i] >= MaxBatchRecs )
			CUDA_ERR("h_Grid_NRs[i]=%d > MaxBatchRecs=%d\n", h_Grid_NRs[i], MaxBatchRecs);

		if ( accRecs + h_Grid_NRs[i] <= MaxBatchRecs ) {
			// pack into the current batch
			accRecs += h_Grid_NRs[i];
		}
		else {
			// new batch
			jobBatches[numBatches].job_s = job_s; 
			jobBatches[numBatches].job_e = job_i; 
			accRecs = h_Grid_NRs[i];
			numBatches++;
			job_s = job_i;
		}
		
		numCompPairs += ((unsigned long long int)h_Grid_NRs[i]*h_Grid_NRs[i]/2);

		for(int j=0;j<(h_Grid_NRs[i]+NTPB-1)/NTPB;j++) {
			h_jobs[job_i].gid_s = i;
			h_jobs[job_i].gid_e = i;
			h_jobs[job_i].seg_s = j * NTPB;
			h_jobs[job_i].seg_e = (j + 1) * NTPB >= h_Grid_NRs[i] ? h_Grid_NRs[i] : (j+1)*NTPB ;
			job_i++;
		}

	}
	DBG("last batch %d %d %d\n", numBatches, job_s, job_i);
	jobBatches[numBatches].job_s = job_s; 
	jobBatches[numBatches].job_e = job_i; 
	numBatches++;

	for(int i=0;i<numBatches;i++) {
		DBG("batch %5d job(%6d,%6d) recidx(%6d,%6d)\n",
					   i,
				       jobBatches[i].job_s,
					   jobBatches[i].job_e,
					   h_Grid_StartI[h_jobs[jobBatches[i].job_s].gid_s],
					   h_Grid_StartI[h_jobs[jobBatches[i].job_e-1].gid_s]+
					      h_Grid_NRs[h_jobs[jobBatches[i].job_e-1].gid_s]
					   );

		if (h_Grid_StartI[h_jobs[jobBatches[i].job_e-1].gid_s] +
			h_Grid_NRs[h_jobs[jobBatches[i].job_e-1].gid_s]   -
			h_Grid_StartI[h_jobs[jobBatches[i].job_s].gid_s] > MaxBatchRecs )
			CUDA_ERR("batch size > %d\n", MaxBatchRecs); 
	}
	DBG("Number of Jobs = %d\n", job_i);
#endif

	//h_block_jobs_e[job_i] = nCells-1;
	//h_block_jobs_offset[job_i] = 0;
	//job_i++;

	if ( MAXJOBS < job_i ) {
		LOG("Error: Due to MAXJOBS (%d) < job_i (%d)\n", MAXJOBS , job_i);
		exit(0);
	}

	//LOG("cpu[block_scheduling]: elapsed time %ld\n", (cuda_end.tv_sec*1000000+cuda_end.tv_usec)-(cuda_start.tv_sec*1000000+cuda_start.tv_usec));
	//total_cuda_time += (cuda_end.tv_sec*1000000+cuda_end.tv_usec)-(cuda_start.tv_sec*1000000+cuda_start.tv_usec);
	cudaMemcpy( d_jobs, h_jobs, sizeof(struct Job)*job_i, cudaMemcpyHostToDevice);

	CUDA_EVENT_END
	CHECK_CUDA_ERR
	LOG("host[block_scheduling]: elapsed time %lf\n", elapsedTime*1000);
	total_cuda_time += elapsedTime*1000;
	LOG("Total %d jobs AvgCells per job %.2f\n", job_i, (float)(nCells)/job_i);
	for(int i=(job_i>=10 ? job_i-10 : 0);i<job_i;i++) {
		DBG("scheduling blocks %d %d %d %d\n", i, h_jobs[i].gid_s, h_jobs[i].gid_e, h_jobs[i].seg_s);
	}


	return job_i;

}



void decode(unsigned int);

static void match_phase(unsigned int job_i)
{

    unsigned long long int numBlocks;

/*  kernel running: matching */
#if 1 
	*h_sum = 0;
	cudaMemcpy( d_sum, h_sum, sizeof(unsigned long long int), cudaMemcpyHostToDevice);	

	unsigned int sharedmem_size = sizeof(RectInfo)*NTPB;

#if 0
	if ( NStreams != NOuts ) {
		LOG("Error: NStreams != NOuts \n");
		cuda_finalize();
		exit(0);
	}
#endif

	cudaStream_t streams[NStreams+USE_ORDER*2];
	cudaEvent_t exec_events[NStreams+USE_ORDER*2];
	cudaEvent_t data_events[NStreams+USE_ORDER*2];
	for(int s=0;s<NStreams+USE_ORDER*2;s++) {
		cudaStreamCreate( &(streams[s]));
		cudaEventCreate( &(exec_events[s]));
		cudaEventCreate( &(data_events[s]));
		//cudaEventCreateWithFlags( &(exec_events[s]), cudaEventDisableTiming);
		//cudaEventCreateWithFlags( &(data_events[s]), cudaEventDisableTiming);
		CHECK_CUDA_ERR
	}
	
	uint orderDepth = 0, streamDepth = 0;
	unsigned int cur = 0, prev = 1;
	unsigned long long int output_size = 0;
	int extraloops;
	unsigned int prefetch_output = 0;

#if USE_ORDER == 1
	orderDepth = 1;
#endif
#if NStreams == 2
	streamDepth = 1;
#endif
	extraloops = orderDepth*1 + streamDepth*2;

	if ( NStreams == 1 ) {
		prev = 0;
	}

	for(int i=0;i<NStreams;i++) 
		h_overlapDataCnt[i] = 0;
	cudaMemcpy( d_overlapDataCnt, h_overlapDataCnt, sizeof(unsigned int)*NStreams, cudaMemcpyHostToDevice);
	time_copy_out = 0.0;
	time_decode = 0.0;
	time_check = 0.0;
	time_order = 0.0;
	numOverlaps = 0;

	*h_hasProfSteps = _hasProfSteps;
	cudaMemcpy( d_hasProfSteps, h_hasProfSteps, sizeof(unsigned int), cudaMemcpyHostToDevice);
	CHECK_CUDA_ERR

	LOG("numBatches = %d\n", numBatches);
	uint orderLoopS		= 1 + orderDepth*0,					orderLoopE	= numBatches + orderDepth*0;
	uint checkLoopS		= 1 + orderDepth*1,					checkLoopE	= numBatches + orderDepth*1;
	uint  copyLoopS		= 1 + orderDepth*1+streamDepth*1,	copyLoopE	= numBatches + orderDepth*1+streamDepth*1;
	uint decodeLoopS	= 1 + orderDepth*1+streamDepth*2,	decodeLoopE = numBatches + orderDepth*1+streamDepth*2;
	CUDA_EVENT_START
	for(int loop=1;loop<=numBatches+extraloops;loop++) {

		DBG(" --- loop %d --- \n", loop);
#if USE_ORDER == 1
#if NStreams == 1
		STime(hTimerOrder);
#endif
		if ( loop >= orderLoopS && loop <= orderLoopE ) {
			uint kstart, kend;
			uint job_s, job_e;
			job_s = jobBatches[loop-1].job_s;
			job_e = jobBatches[loop-1].job_e;
			numBlocks = job_e - job_s;
			kstart = h_Grid_StartI[h_jobs[job_s].gid_s];
			kend   = h_Grid_StartI[h_jobs[job_e-1].gid_s]+h_Grid_NRs[h_jobs[job_e-1].gid_s];
			DBG("Ordering loops = %d InB=%d numBlocks=%d\n", loop, (loop-1)%2, numBlocks);
			DBG("kend(%d) - kstart(%d) = %d\n", kend, kstart, kend-kstart);
			//cudaEventSynchronize(exec_events[2+((loop-1)%2)]);
			d_ordering<<<numBlocks,NTPB,0,streams[orderDepth+streamDepth]>>>  
				(
									 d_rectInfoList,
									 d_orderRectCoordInfoList+MaxBatchRecs*((loop-1)%2),
									 d_Grid_Recs,
									 kstart,
									 kend,
									 d_log);
		}
#if NStreams == 1
		ETime(hTimerOrder);
		time_order += GTime(hTimerOrder) * 1000.0;
#endif
#endif

#if 1
#if NStreams == 1
		STime(hTimerCheck);
#endif
		//CHECK_CUDA_ER
		if ( loop >= checkLoopS && loop <= checkLoopE ) {
			uint kstart, kend;
			uint job_s, job_e;
			job_s = jobBatches[loop-orderDepth-1].job_s;
			job_e = jobBatches[loop-orderDepth-1].job_e;
			numBlocks = job_e - job_s;
			kstart = h_Grid_StartI[h_jobs[job_s].gid_s];
			kend   = h_Grid_StartI[h_jobs[job_e-1].gid_s] + h_Grid_NRs[h_jobs[job_e-1].gid_s];
			DBG("Matching loops = %d InB=%d, OutB=%d, numBlocks=%d\n", loop, cur, cur, numBlocks);
			DBG("kstart(%d), kend(%d)\n", kstart, kend);


				d_matching<<<numBlocks,NTPB,sharedmem_size,streams[cur]>>>  
					(
													 d_hasProfSteps,
													 d_rectInfoList,
													 d_orderRectCoordInfoList+MaxBatchRecs*((loop)%2),
													 d_Grid_NRs,
													 d_Grid_StartI,
													 d_Grid_Merged,
													 d_Grid_Recs,
													 job_s,
													 kstart,
													 d_jobs,
													 _cellSize,
													 log2((float)_cellSize),
													 nRowCells,
													 MAX_OVERLAP_DATA,
													 d_dst_overlap+cur*MAX_OVERLAP_DATA,
													 d_overlapDataCnt+cur,
													 d_log);
				//cudaEventRecord( exec_events[2+(loop%2)], 0);


		}
#endif
		//CUDA_COPY_LOG
		//CHECK_CUDA_ERR
#if NStreams == 1
		ETime(hTimerCheck);
		time_check += GTime(hTimerCheck) * 1000.0;
#endif


#if 1
#if NStreams == 2
		/* prefetch */
		if ( loop >= copyLoopS && loop <= copyLoopE && prefetch_output != 0 ) {
			DBG("Prefetching loops=%d, OutB=%d, prefetch_output=%d\n", loop, prev, prefetch_output);
			cudaMemcpyAsync( h_dst_overlap+prev*MAX_OVERLAP_DATA, 
							 d_dst_overlap+prev*MAX_OVERLAP_DATA, 
							 //((unsigned int*)d_dst_overlap4)+prev*MAX_OVERLAP_DATA, 
							 sizeof(OTYPE)*prefetch_output, cudaMemcpyDeviceToHost, streams[prev]);
		}
		//CHECK_CUDA_ERR

		/* decode */
		if ( loop >= decodeLoopS && loop <= decodeLoopE ) {
			cudaEventSynchronize(data_events[cur]);
			DBG("Decoding loops=%d, OutB=%d, decode data h_overlapDataCnt[%d]=%d\n", loop, cur, cur, h_overlapDataCnt[cur]);
			if ( h_overlapDataCnt[cur] != 0 ) {
				if ( h_overlapDataCnt[cur] >= MAX_OVERLAP_DATA ) {
					LOG("Error: loop %d decode data h_overlapDataCnt[%d]=%d\n", loop, cur, h_overlapDataCnt[cur]);
					exit(0);
				}
				decode(cur);
				h_overlapDataCnt[cur] = 0;
			}
		}
		//CHECK_CUDA_ERR
#endif

		/* get the output size */
		if ( loop >= copyLoopS && loop <= copyLoopE ) {
			cudaMemcpyAsync( h_overlapDataCnt+prev, d_overlapDataCnt+prev, 
							 sizeof(unsigned int), cudaMemcpyDeviceToHost, streams[prev]);

			cudaMemsetAsync( d_overlapDataCnt+prev, 0, sizeof(unsigned int), streams[prev]);
			cudaEventRecord( exec_events[prev], 0);
		}
		//CHECK_CUDA_ERR


		/* get the output */
		if ( loop >= copyLoopS && loop <= copyLoopE ) {
			cudaEventSynchronize(exec_events[prev]);
			DBG("Coping loops=%d, OutB=%d, h_overlapDataCnt[%d]=%d prefetch_output=%d\n", loop, prev, prev, h_overlapDataCnt[prev], prefetch_output);
			if ( h_overlapDataCnt[prev] >= MAX_OVERLAP_DATA ) {
				CUDA_ERR("h_overlapDataCnt[%d]=%d > MAX_OVERLAP_DATA(%d)\n", prev, h_overlapDataCnt[prev], MAX_OVERLAP_DATA);
			}
			output_size += h_overlapDataCnt[prev];
			if ( h_overlapDataCnt[prev] > prefetch_output ) {
				unsigned int rest = h_overlapDataCnt[prev] - prefetch_output ; 
#if NStreams == 1
				STime(hTimerCopyOut);
#endif
				cudaMemcpyAsync( h_dst_overlap+prev*MAX_OVERLAP_DATA+prefetch_output, 
								 d_dst_overlap+prev*MAX_OVERLAP_DATA+prefetch_output, 
								 //((unsigned int*)d_dst_overlap4)+prev*MAX_OVERLAP_DATA+prefetch_output, 
								 sizeof(OTYPE)*rest, cudaMemcpyDeviceToHost, streams[prev]);
#if NStreams == 1
				ETime(hTimerCopyOut);
				time_copy_out += GTime(hTimerCopyOut) * 1000.0;
#endif
			}
#if NStreams == 2
			prefetch_output = (float) h_overlapDataCnt[prev]*1.1;
#endif
			cudaEventRecord( data_events[prev], 0);
		}
		//CHECK_CUDA_ERR


#if NStreams == 1
		/* decode */
		if ( loop >= decodeLoopS && loop <= decodeLoopE ) {
			cudaEventSynchronize(data_events[0]);
			DBG("loop %d: decode data h_overlapDataCnt[%d]=%d\n", loop, cur, h_overlapDataCnt[cur]);
			if ( h_overlapDataCnt[cur] != 0 ) {
				if ( h_overlapDataCnt[cur] >= MAX_OVERLAP_DATA ) {
					CUDA_ERR("Error: loop %d decode data h_overlapDataCnt[%d]=%d > MAX_OVERLAP_DATA\n", loop, cur, h_overlapDataCnt[cur]);
				}
				decode(cur);
				h_overlapDataCnt[cur] = 0;
			}
		}
		//CHECK_CUDA_ERR
#endif

		cur  =  (cur+1) % NStreams;
		prev = (prev+1) % NStreams;
#endif
	}
	CUDA_EVENT_END
	CUDA_COPY_LOG
	CHECK_CUDA_ERR

	LOG("kernel[d_MDD]: elapsed time %lf\n", elapsedTime*1000);
	LOG("kernel[d_order]: elapsed time %lf\n", time_order);
	LOG("kernel[d_matching]: elapsed time %lf\n", time_check);
	LOG("kernel[copy out]: elapsed time %lf\n", time_copy_out);
	LOG("kernel[decode]: elapsed time %lf\n", time_decode);
	LOG("Output Size = %ld\n", output_size);
	LOG("Maximum number of overlaps for a box = %d\n", max_nlaps);
	LOG("Number of Intersections = %lld\n", numOverlaps);
	LOG("Number of Compared Pairs = %lld\n", numCompPairs);
	total_cuda_time += elapsedTime*1000;
	CHECK_CUDA_ERR

	for(int s=0;s<NStreams;s++) {
		cudaStreamDestroy( streams[s]);
	}
#endif



    DBG("\n\tCHECK Region ID sum=%lld\n", *h_sum);


}

void pbig_check(RectInfo *input, struct report *r) 
{

	_report = r;

	plan_phase( input);

	map_phase();

	unsigned int job_i = sched_phase();

	match_phase( job_i);

	LOG("\nhost+kernel: elapsed time %lf\n\n", total_cuda_time);

	total_cuda_time = 0.0;

	deallocMemForGrid();
	deallocMemForGridRecs();
}


void allocMemForInitAndRects()
{

	/* space for profiling */
	h_sum = (unsigned long long int*) malloc(sizeof(unsigned long long int));
	h_log = (unsigned int*) malloc(MAXLOG*sizeof(unsigned int));
	cudaHostAlloc( &h_hasProfSteps, sizeof(unsigned int), 0);
	// in device
	cudaMalloc( (void**) &d_sum, sizeof(unsigned long long int));
	cudaMalloc( (void**) &d_log, MAXLOG*sizeof(unsigned int));
	cudaMalloc( (void**) &d_hasProfSteps, sizeof(unsigned int));
	LOG("\t\t\t\t\t\t\t\t\t[CUDA] %d %s\n", __LINE__, cudaGetErrorString(cudaGetLastError()));


	/* space for input data */
	CHECK_AVAIL_MEM( _NumObjs*sizeof(RectInfo) + 2*MaxBatchRecs*sizeof(struct CoordInfo));
	cudaMalloc( (void**) &d_rectInfoList, _NumObjs*sizeof(RectInfo));
	cudaMalloc( (void**) &d_orderRectCoordInfoList, 2*MaxBatchRecs*sizeof(struct CoordInfo));
	MEM_STATUS("RectInfo");

}

void allocMemForGrid() {


	/* space for grid */
	//CHECK_AVAIL_MEM( MAXGRIDS * sizeof(unsigned int)*3 + MAX_Grid_Recs*sizeof(unsigned int));
	//cudaHostAlloc( &h_Grid_NRs, sizeof(unsigned int)*MAXGRIDS,0);
	//cudaHostAlloc( &h_Grid_StartI, sizeof(unsigned int)*MAXGRIDS,0);
	//cudaHostAlloc( &h_Grid_Merged, sizeof(unsigned int)*MAXGRIDS,0);
	//cudaMalloc( (void**) &d_Grid_NRs, MAXGRIDS*sizeof(unsigned int));
	//cudaMalloc( (void**) &d_Grid_StartI, MAXGRIDS*sizeof(unsigned int));
	//cudaMalloc( (void**) &d_Grid_Merged, MAXGRIDS*sizeof(unsigned int));
	//cudaHostAlloc( &h_Grid_Recs, sizeof(unsigned int)*MAX_Grid_Recs,0);
	//cudaMalloc( (void**) &d_Grid_Recs, MAX_Grid_Recs*sizeof(unsigned int));
	//MEM_STATUS("GridInfo");

	CHECK_AVAIL_MEM( nCells * sizeof(unsigned int)*3);
	cudaHostAlloc( &h_Grid_NRs, sizeof(unsigned int)*nCells,0);
	cudaHostAlloc( &h_Grid_StartI, sizeof(unsigned int)*nCells,0);
	cudaHostAlloc( &h_Grid_Merged, sizeof(unsigned int)*nCells,0);
	cudaMalloc( (void**) &d_Grid_NRs, nCells*sizeof(unsigned int));
	cudaMalloc( (void**) &d_Grid_StartI, nCells*sizeof(unsigned int));
	cudaMalloc( (void**) &d_Grid_Merged, nCells*sizeof(unsigned int));
	MEM_STATUS("GridInfo");
}


void deallocMemForGrid() {


	cudaFreeHost(h_Grid_NRs);
	cudaFreeHost(h_Grid_StartI);
	cudaFreeHost(h_Grid_Merged);
	cudaFree( d_Grid_NRs);
	cudaFree( d_Grid_StartI);
	cudaFree( d_Grid_Merged);
	MEM_STATUS("Dealloc GridInfo");
}


void allocMemForGridRecs() {

	CHECK_AVAIL_MEM( MAX_Grid_Recs*sizeof(unsigned int));
	//cudaHostAlloc( &h_Grid_Recs, sizeof(unsigned int)*MAX_Grid_Recs,0);
	h_Grid_Recs = (unsigned int*) malloc(sizeof(unsigned int)*MAX_Grid_Recs);
	cudaMalloc( (void**) &d_Grid_Recs, MAX_Grid_Recs*sizeof(unsigned int));
}

void deallocMemForGridRecs() {

	//cudaFreeHost( h_Grid_Recs);
	free(h_Grid_Recs);
	cudaFree( d_Grid_Recs);
	MEM_STATUS("Dealloc GridRecs");
}


void allocMemForOverlap() {


	/* space for overlap results */
	CHECK_AVAIL_MEM( MAX_OVERLAP_DATA*NStreams*sizeof(OTYPE));
	cudaHostAlloc( &h_overlapDataCnt, sizeof(unsigned int)*NStreams, 0);
	cudaMalloc( (void**) &d_overlapDataCnt, sizeof(unsigned int)*NStreams);

	LOG("\t\t\t\t\t\t\t\t\t[CUDA] %d %s\n", __LINE__, cudaGetErrorString(cudaGetLastError()));

	cudaHostAlloc( &h_dst_overlap, sizeof(OTYPE)*MAX_OVERLAP_DATA*NStreams, 0);
	cudaMalloc( (void**) &d_dst_overlap, MAX_OVERLAP_DATA*NStreams*sizeof(OTYPE));
	LOG("\t\t\t\t\t\t\t\t\t[CUDA] %d %s\n", __LINE__, cudaGetErrorString(cudaGetLastError()));
	MEM_STATUS("OverlapInfo");
}

void allocMemForJob() {
	/* space for job scheduling */
	CHECK_AVAIL_MEM( MAXJOBS*sizeof(struct Job));
	cudaHostAlloc( &h_jobs, sizeof(struct Job)*MAXJOBS, 0);
	cudaMalloc( (void**) &d_jobs, sizeof(struct Job)*MAXJOBS);
	LOG("\t\t\t\t\t\t\t\t\t[CUDA] %d %s\n", __LINE__, cudaGetErrorString(cudaGetLastError()));
	MEM_STATUS("JobInfo");
}

int device_init(unsigned int deviceSymbol)
{

	//checkCudaErrors(cudaSetDevice(0));
	//checkCudaErrors(cudaDeviceEnablePeerAccess(1, 0));
	//checkCudaErrors(cudaSetDevice(1));
	//checkCudaErrors(cudaDeviceEnablePeerAccess(0, 0));

	if ( deviceSymbol == 480 )
		cuDeviceID = 0;
	else if ( deviceSymbol == 2070 )
		cuDeviceID = 0;
	else {
		cuDeviceID = 0;
	}
	LOG("Use cuda deivce %d\n", cuDeviceID);

	return cuDeviceID;

#if 0
	if ( deviceSymbol == 480 || deviceSymbol == 2070 ) {

		if ( option == 0 ) {
			if ( NTPB >= 32 ) {
				LOG("Set L1 Cache 16 KB, shared memory 48 KB\n");
				cudaFuncSetCacheConfig( d_matching, cudaFuncCachePreferShared);
			}
			else {
				LOG("Set L1 Cache 48 KB, shared memory 16 KB\n");
				cudaFuncSetCacheConfig( d_matching, cudaFuncCachePreferL1);
			}
		}
		else if ( option == 16 ) {
				LOG("Set L1 Cache 16 KB, shared memory 48 KB\n");
				cudaFuncSetCacheConfig( d_matching, cudaFuncCachePreferShared);
		}
		else if ( option == 48 ) { 
				LOG("Set L1 Cache 48 KB, shared memory 16 KB\n");
				cudaFuncSetCacheConfig( d_matching, cudaFuncCachePreferL1);
		}
		else {

			LOG("Wrong options in setCacheSize \n");
			exit(-1);
		}
	}
	else {
		LOG("Default shared memory 16 KB\n");
	}
#endif
}

void pbig_init(unsigned int num, unsigned range, unsigned int cs, unsigned int deviceID, unsigned int hasProfSteps)
{

	_NumObjs = num;
	_maxRange = range;
	_hasProfSteps = hasProfSteps;
	_assignCs = cs;

#if RELEASE == 1
	if ( cs != 0 || hasProfSteps != 0 ) {
		LOG("The 3rd and 5th arguments in pbig_init() must be zero\n");
		_assignCs = 0;
		_hasProfSteps = 0;
	}
#endif

	MAX_Grid_Recs = _NumObjs * 4;
	//MAX_OVERLAP_DATA = _NumObjs * OVERLAP_LEN;
	//MAX_OVERLAP_DATA = _NumObjs;
	MAX_OVERLAP_DATA = 10000000*4;

	//LOG("cuda_init _NumObjs %d MAXALLOCRECTS %d MAXGRIDS %d MAX_OVERLAP_DATA %d\n", 
	//		_NumObjs, MAXALLOCRECTS, MAXGRIDS, MAX_OVERLAP_DATA);
	LOG("size:\n\tRectInfo \t%d bytes\n\n", 
			sizeof(RectInfo));


	// show optimization settings: USE_SM, MAP_ALL, NStreams, NStreams, Compression 
	LOG("Using shared memory\n");
#ifdef MAP_ALL 
	LOG("Using MAP_ALL\n");
#endif
#ifdef RICE_COMPRESSION
	LOG("Compression=%d Rice Compression=%d\n", COMPRESSION, RICE_COMPRESSION);
#endif
	LOG("NStreams=%d, NStreams=%d\n", NStreams, NStreams);
	
	cuDeviceID = deviceID;
	checkCudaErrors(cudaSetDevice(deviceID));
	cudaSetDeviceFlags( cudaDeviceBlockingSync); 

    cudaEventCreate(&e_start);
    cudaEventCreate(&e_stop);
	MEM_STATUS("Init");

    sdkCreateTimer(&hTimerCopyOut);
    sdkCreateTimer(&hTimerCheck);
    sdkCreateTimer(&hTimerOrder);

	//if ( _NumObjs > MAXALLOCRECTS ) {
	//	LOG("NumObjs(%ld) > MAXALLOCRECTS(%ld)\n", _NumObjs, MAXALLOCRECTS);
	//	exit(-1);
	//}

	// memory allocation
	allocMemForInitAndRects();
	allocMemForOverlap();
	allocMemForJob();

	*h_sum = 0;
	for(int k=0;k<MAXLOG;k++)
		h_log[k] = 0;




	CUDA_EVENT_START
		cudaMemset( d_sum, 0, sizeof(unsigned long long int));
		cudaMemset( d_log, 0, sizeof(unsigned int)*MAXLOG);
	CUDA_EVENT_END
	CHECK_CUDA_ERR
	LOG("kernel[CudaInit]: elapsed time %lf\n", elapsedTime*1000);

#if 1
	CHECK_CUDA_ERR
#endif


}


