#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "cuda_comm.h"
#include "pbig_internal.h"
#include "report.h"

extern int MAX_OVERLAP_DATA;
extern unsigned long long int numOverlaps;
extern struct RectInfo *h_rectInfoList;
extern unsigned int *h_overlapDataCnt;
extern OTYPE *h_dst_overlap;
extern unsigned int *h_Grid_StartI;
extern unsigned int *h_Grid_Recs;
extern unsigned int *h_Grid_NRs;
extern struct report *_report;
extern float time_decode;
extern unsigned int max_nlaps;

#define min(a,b) (a)>(b)?(b):(a)
#define max(a,b) (a)>(b)?(a):(b)

void decode(unsigned int cur) {


	struct timeval decode_start, decode_end;
	gettimeofday( &decode_start, NULL);
	DBG("h_overlapDataCnt[%d] = %d\n", cur, h_overlapDataCnt[cur]);
	OTYPE *p = h_dst_overlap+cur*MAX_OVERLAP_DATA;
	unsigned int seg_size = min(NTPB,256);
	unsigned int seg_size_shift = log2((float)seg_size);
	unsigned int l_rid, nlaps, base1,base2;
	unsigned int seg1, seg2;
#if 1
	for(int i=2;i<h_overlapDataCnt[cur];i+=2) {
		//unsigned int magic = p[i-2];
		unsigned int gid = p[i-1];
		

#ifdef RICE_COMPRESSION
#else
		unsigned int encode_unit;
		if ( h_Grid_NRs[gid] > 256 )
			encode_unit = 2;
		else
			encode_unit = 1;
#endif

		do {
#ifdef COMPRESSION
#ifdef RICE_COMPRESSION
			seg1 =  ((unsigned char*)p)[(i<<OTYPE_SHIFT)+0];
			seg2 =  ((unsigned char*)p)[(i<<OTYPE_SHIFT)+1];
			//l_rid = ((unsigned char*)p)[(i<<2)+2];
			//nlaps = ((unsigned char*)p)[(i<<2)+3];
			l_rid = ((unsigned short*)p)[(i<<(OTYPE_SHIFT-1))+1];
			nlaps = ((unsigned short*)p)[(i<<(OTYPE_SHIFT-1))+2];
			base1 = h_Grid_StartI[gid] + (seg1 << seg_size_shift);
			base2 = h_Grid_StartI[gid] + (seg2 << seg_size_shift);
#else
			if ( encode_unit == 2 ) {
				l_rid = ((unsigned short*)p)[(i<<1)+0];
				nlaps = ((unsigned short*)p)[(i<<1)+1];
			}
			else {
				l_rid = ((unsigned char*)p)[(i<<2)+0];
				nlaps = ((unsigned char*)p)[(i<<2)+1];
			}
			base1 = base2 = h_Grid_StartI[gid];
#endif
			unsigned int rid1 = h_Grid_Recs[base1+l_rid];
			//if ( gid == 0 ) 
			//	printRec( rid1);
			for(int k=0;k<nlaps;k++) {
				unsigned int t;
#ifdef RICE_COMPRESSION
				t = ((unsigned char *)p)[((i+1)<<OTYPE_SHIFT)+2+k];
				//t = ((unsigned char *)p)[((i+1)<<2)+k];
#else
				if ( encode_unit == 2 )
					t = ((unsigned short *)p)[((i+1)<<1)+k];
				else
					t = ((unsigned char *)p)[((i+1)<<2)+k];
#endif

				//if ( !checkOverlaps( rid1, h_Grid_Recs[base2+t]) ) {
				//	LOG("i=%d gid=%d numLocalRects=%d nlaps=%d seg1 %d seg2 %d offset1 %d offset2 %d rid1 %d rid2 %d\n", i, gid, h_Grid_NRs[gid], nlaps, seg1, seg2, l_rid, t, h_Grid_Recs[base1+l_rid], h_Grid_Recs[base2+t]);
				//	printRec( rid1);
				//	printRec( h_Grid_Recs[base2+t]);
				//	exit(0);
				//}
				//if ( gid == 0 ) {
				//	LOG("\t");
				//	printRec( h_Grid_Recs[base+t]);
				//}
				//(*_report).found( &(h_rectInfoList[rid1]), &(h_rectInfoList[h_Grid_Recs[base2+t]])); 
				(*_report)( rid1, h_Grid_Recs[base2+t]); 
			}
			//  i + header_size + numer of overlaps in 4bytes
#ifdef RICE_COMPRESSION 
			i = i + ((6 + nlaps + sizeof(OTYPE)-1) >> OTYPE_SHIFT);
			//i = i + 1 + ((nlaps + 3) >> 2);
#else
			if ( encode_unit == 2 ) 
				i = i + 1 + ((nlaps + 1) >> 1);
			else
				i = i + 1 + ((nlaps + 3) >> 2);
#endif
#else // non-compression
			unsigned int rid1 = p[i];
			unsigned int nlaps = p[i+1];
			for(int k=0;k<nlaps;k+=2) {
				//if ( !checkOverlaps( p[i+2+k], p[i+2+k+1] )) {
				//	LOG("i=%d gid=%d nlaps=%d rid1 %d rid2 %d\n", i, gid, nlaps, p[i+2+k], p[i+2+k+1]);
				//	exit(0);
				//}
				//(*_report).found( &(h_rectInfoList[p[i+2+k]]), &(h_rectInfoList[p[i+2+k+1]])); 
				(*_report)( p[i+2+k], p[i+2+k+1]); 
			}
			i = i + 2 + nlaps;
			nlaps /= 2;
#endif
			max_nlaps = max( max_nlaps, nlaps);

			numOverlaps += nlaps;
		} while ( p[i] != UNUSED && i < h_overlapDataCnt[cur]);


	}
#endif
	gettimeofday( &decode_end, NULL);
	time_decode += (decode_end.tv_sec*1000000+decode_end.tv_usec)-(decode_start.tv_sec*1000000+decode_start.tv_usec);
	//printf("kernel[decode]: elapsed time %lf\n", time_decode);

}
