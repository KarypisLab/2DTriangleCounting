/*!
  \file
  \brief The various triangle counting routines
 */

#include "tc.h"
#include <mpi.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#ifndef DYNAMIC_CHUNK
#define DYNAMIC_CHUNK 16
#endif

#ifndef PPT_DYNAMIC_CHUNK
#define PPT_DYNAMIC_CHUNK 500
#endif

#ifndef factor
#define factor 0.5
#endif

#ifndef BLK_SZ
#define BLK_SZ 1
#endif

#ifndef NPROCS
#define NPROCS 1024
#endif

#define adj_t int64_t
#define MPI_ADJ_INT MPI_LONG_LONG_INT

#define blk_t int
#define MPI_BLK_CNT_INT MPI_INT
//#define ADJNCY_ADM_MAX 1431655766.0
//#define XADJ_ADM_MAX 134217728
#define ADJNCY_ADM_MAX 1431655766.0
#define XADJ_ADM_MAX 1431655766.0

/* OSX timer includes */
 #ifdef __MACH__
   #include <mach/mach.h>
   #include <mach/mach_time.h>
 #endif

//double tct, pkt, unpkt, commt;
//double dcntt, ineedt1, ineedt2, ineedt3, ineedt4;
//double ult1, ult2, ult3, ult4, sortt;
//double end1, end2, end3;

 /**
 * @brief Return the number of seconds since an unspecified tim
 e (e.g., Unix
 *        epoch). This is accomplished with a high-resolution m
 onotonic timer,
 *        suitable for performance timing.
 *
 * @return The number of seconds.
 */
static inline double monotonic_seconds()
{
 #ifdef __MACH__
   /* OSX */
   static mach_timebase_info_data_t info;
   static double seconds_per_unit;
   if(seconds_per_unit == 0) {
     mach_timebase_info(&info);
     seconds_per_unit = (info.numer / info.denom) / 1e9;
   }
   return seconds_per_unit * mach_absolute_time();
 #else
   /* Linux systems */
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC, &ts);
   return ts.tv_sec + ts.tv_nsec * 1e-9;
 #endif
}


static void print_time(double const seconds) {
        printf("partial time: %0.04fs\n", seconds);
}

/*************************************************************************/
/*! Determine the iperm for the key order using counting sort.
 */
/*************************************************************************/
ssize_t *gk_sskvipermi(vtxs_t n, ssize_t *xadj)
{
  ssize_t i; 
  ssize_t *iperm, *perm;

  int32_t *counts; 
  gk_i32kv_t *cand;
  int range;

  /* TODO: check if we need to change this? */

  cand = gk_i32kvmalloc(n, "cand");
  for (i=0; i<n; i++) {
    //printf(" %zd ", (xadj[i+1]-xadj[i]));
    cand[i].key = (int32_t)(xadj[i+1]-xadj[i]);
    cand[i].val = i;
  }
  //printf("\n");

  for (range=0, i=0; i<n; i++) {
    if (cand[i].key > range)
      range = cand[i].key;
  }
  range++;

  counts = gk_i32smalloc(range+1, 0, "counts");
  for (i=0; i<n; i++)
    counts[cand[i].key]++;
  MAKECSR(i, range, counts);

  //iperm = gk_i32smalloc(n, 0, "iperm");
  iperm = (ssize_t *) malloc(n*sizeof(ssize_t));
  perm = (ssize_t *) malloc(n*sizeof(ssize_t));
  memset(iperm, 0, n*sizeof(ssize_t));
  for (i=0; i<n; i++)
    iperm[counts[cand[i].key]++] = i;

  for (i=0; i<n; i++) {
    perm[iperm[i]] = i;
  }

  gk_free((void **)&counts, &cand, &iperm, LTERM);

  return perm;
}

int mod(int a, int b) {
  int r = a%b;
  return r < 0 ? r + b : r;
}

int ceil1(double x) {
  if((int)x * (double)(1.0) == x)
    return x;
  else
    return (int)x + 1;
}

int64_t mpitc_kernelJIKBlockCyclic7(
  int shift,
  matrix_info_tn *la_matrix, 
  blob_info_t *rblob, 
  blob_info_t *cblob) 
{

  /* TODO: Include checks that will ensure that it doesn't 
    count triangles for sub blocks of the form i<j<k. */
  matrix_info_tn *ra_matrix, *rb_matrix;

  int mype, npes, rootp; 
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  int64_t ntriangles=0;
  int32_t ltriangles=0; /* Do we need to change this to int64_t too? */

  ssize_t la_vi, la_vj;
  ssize_t vi, vj, vk;

  double start;
  double ltime, gtime, atime, mtime;  

  adj_t l, nc, hmsize, *hmap;

  adj_t rb_vtx;
  ssize_t *la_xadj, *ra_xadj, *rb_xadj;
  adj_t *la_adjncy, *ra_adjncy, *rb_adjncy; 

  rootp = (int) sqrt(npes);

  /* unpkt */

  /* Unpack 1. */
  ra_matrix = mpitc_unpackBlobBlockCyclic(rblob);
  /* Unpack 2. */
  rb_matrix = mpitc_unpackBlobBlockCyclic(cblob);

  la_xadj = la_matrix->xadj; 
  ra_xadj = ra_matrix->xadj; 
  rb_xadj = rb_matrix->xadj; 

  la_adjncy = la_matrix->adjncy;
  ra_adjncy = ra_matrix->adjncy; 
  rb_adjncy = rb_matrix->adjncy; 

  //MPI_Barrier(MPI_COMM_WORLD);
  //start = MPI_Wtime();

  /* Figuring out the size of hmap. */
  for (hmsize=0, vi=0; vi<rb_matrix->g_nvtxs; vi++) {
    hmsize = gk_max(hmsize, (adj_t)(rb_xadj[vi+1]-rb_xadj[vi]));
  }

  for (l=1; hmsize>(1<<l); l++);
  hmsize = (1<<(l+4))-1;

  hmap = gk_i64smalloc(hmsize+1, -1, "hmap");

  /* For each communicated ra-chunk, find the number of 
      residing triangles. */
  int32_t g_off = BLK_SZ*rootp;
  adj_t ch_vtx = -1;

  /* tct */
  adj_t nnz_idx; 
  /* nvtxs corresponds to either unvtxs/lnvtxs. */
  /* This helps you navigate the nnz idxs (verify). */
  vtxs_t vtxs = rb_matrix->nvtxs;

  ssize_t lcnt; 
  for (la_vi=0; la_vi<rb_matrix->nvtxs; la_vi++) {

    nnz_idx = rb_matrix->nnz_idx[la_vi];

    if (la_xadj[nnz_idx+1]==la_xadj[nnz_idx])
      continue;

    for (lcnt=0, la_vj=la_xadj[nnz_idx]; la_vj<la_xadj[nnz_idx+1]; la_vj++) {
      lcnt += ra_xadj[la_adjncy[la_vj]+1]-ra_xadj[la_adjncy[la_vj]];
    }

    if (!lcnt)
      continue;

    if (rb_adjncy[rb_xadj[nnz_idx+1]-1]-rb_adjncy[rb_xadj[nnz_idx]]+1 > hmsize) { /*Possible collisions.*/

      /* Hash rb adjncy to perform lookups, i-th chunk in focus. */
      for (nc=0, vi=rb_xadj[nnz_idx+1]-1; vi>=rb_xadj[nnz_idx]; vi--) {

        vj = rb_adjncy[vi];
        for (l=vj&hmsize; hmap[l]!=-1; l=((l+1)&hmsize), nc++);

        hmap[l] = vj;
      }

      ch_vtx = rb_adjncy[rb_xadj[nnz_idx]];

      if (nc > 0) {
        for (ltriangles=0, la_vj=la_xadj[nnz_idx]; la_vj<la_xadj[nnz_idx+1]; la_vj++) {

          rb_vtx = la_adjncy[la_vj];

          for (vj=ra_xadj[rb_vtx+1]-1; vj>=ra_xadj[rb_vtx]; vj--) {

            vk = ra_adjncy[vj];

            if (vk<ch_vtx)
              break; 

            for (l=(vk&hmsize); hmap[l]!=vk && hmap[l]!=-1; l=((l+1)&hmsize)) {
            }

            if (hmap[l]==vk) {
              ltriangles++;
            }
          }
        }
        ntriangles += ltriangles;

        for (vi=rb_xadj[nnz_idx+1]-1; vi>=rb_xadj[nnz_idx]; vi--) {

          vj = rb_adjncy[vi];
          for (l=vj&hmsize; hmap[l]!=vj; l=((l+1)&hmsize)) {
          }
          hmap[l] = -1;
        }
      } else {
        for (ltriangles=0, la_vj=la_xadj[nnz_idx]; la_vj<la_xadj[nnz_idx+1]; la_vj++) {

          rb_vtx = la_adjncy[la_vj];

          /* Intersection -- 
             only intersect in case of elements being persisted in hmap. */
          for (vj=ra_xadj[rb_vtx+1]-1; vj>=ra_xadj[rb_vtx]; vj--) {

            vk = ra_adjncy[vj];

            if (vk<ch_vtx)
              break; 

            if (hmap[vk&hmsize]==vk) {
              ltriangles++;
            }
          }
        }

        ntriangles += ltriangles;

        for (vi=rb_xadj[nnz_idx+1]-1; vi>=rb_xadj[nnz_idx]; vi--) {

          vj = rb_adjncy[vi];
          hmap[vj&hmsize] = -1;
        }
      }
    } else { /*no collisions.*/

      /* Hash rb adjncy to perform lookups, i-th chunk in focus. */
      for (vi=rb_xadj[nnz_idx+1]-1; vi>=rb_xadj[nnz_idx]; vi--) {

        vj = rb_adjncy[vi];
        hmap[vj&hmsize] = vj;
      }

      ch_vtx = rb_adjncy[rb_xadj[nnz_idx]];

      for (ltriangles=0, la_vj=la_xadj[nnz_idx]; la_vj<la_xadj[nnz_idx+1]; la_vj++) {

        rb_vtx = la_adjncy[la_vj];

        /* Intersection -- 
           only intersect in case of elements being persisted in hmap. */

        for (vj=ra_xadj[rb_vtx+1]-1; vj>=ra_xadj[rb_vtx]; vj--) {

          vk = ra_adjncy[vj];

          if (vk<ch_vtx)
            break; 

          if (hmap[vk&hmsize]==vk) {
            ltriangles++;
          }
        }
      }

      ntriangles += ltriangles;

      for (vi=rb_xadj[nnz_idx+1]-1; vi>=rb_xadj[nnz_idx]; vi--) {

        vj = rb_adjncy[vi];
        hmap[vj&hmsize] = -1;
      }
    }
  }

  /*
  MPI_Reduce(&ltime, &gtime, 1,
      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ltime, &mtime, 1,
      MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  //printf("rank: %d ltime: %f. \n", mype, ltime);

  if (mype==0) {
    atime = gtime/ (double) npes;
    printf("rank: %d tct max: %f, avg: %f, imb: %0.4f. \n", mype, mtime, atime, (mtime/atime));
  }

  MPI_Barrier(MPI_COMM_WORLD); 
  */

  free(hmap);
  free(ra_xadj);
  free(rb_xadj);
  free(ra_adjncy);
  free(rb_adjncy);
  free(ra_matrix->nnz_idx);
  free(rb_matrix->nnz_idx);
  free(ra_matrix);
  free(rb_matrix);

  return ntriangles;
}

int intcmp(const void *x, const void *y)
{       
  return(*((int32_t *)x) - *((int32_t *)y));
}

void SetUpRecvInfo(SendInfoTypeGK *sinfo, RecvInfoTypeGK *rinfo) {

	int npes, mype;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  ssize_t i, j;

	MPI_Request send_req[NPROCS], recv_req[NPROCS];
	MPI_Status statuses[NPROCS];

	int32_t *recv_row = (int32_t *) malloc(npes*sizeof(int32_t));
	int32_t *send_row = (int32_t *) malloc(npes*sizeof(int32_t));
	for(i=0; i<npes; i++)
		send_row[i] = 0;

	for(i=0; i<sinfo->n_nbrs; i++) {
		send_row[sinfo->nbr_ids[i]] = (int32_t) sinfo->ptr[i+1] - sinfo->ptr[i];
	}

	/* I decide what I want to request - that's in the SEND row, but it's actually the recv row. 
   * I will get need to populate the RECV row, but, that's the send row.. */
	MPI_Alltoall(send_row, 1, MPI_INT, recv_row, 1, MPI_INT, MPI_COMM_WORLD);

	/* What will I receive? */
	/* Setup. */

  /* ineedt */
	ssize_t nrecv;
	rinfo->n_nbrs = 0;
	for(nrecv=i=0; i<npes; i++) {
		if(recv_row[i]) {
			nrecv+=recv_row[i];
			rinfo->n_nbrs++;
		}
	}

	rinfo->ptr = (ssize_t *) malloc((rinfo->n_nbrs+1) * sizeof(ssize_t));
	rinfo->nbr_ids = (int *) malloc(rinfo->n_nbrs * sizeof(int));
	rinfo->ind = (adj_t *) malloc(nrecv * sizeof(adj_t));

	ssize_t n_inds = 0;
	rinfo->n_nbrs = 0;
	rinfo->ptr[0] = 0;

	for(i=0; i<npes; i++) {
		if(recv_row[i]) {
			rinfo->nbr_ids[rinfo->n_nbrs++] = (int32_t) i;
			n_inds += recv_row[i];
			rinfo->ptr[rinfo->n_nbrs] = n_inds;
		}
	}

	int tag1; 
  tag1 = 2;
	/* All non-blocking calls - helps not deadlock the processing. */

  #ifdef all2all
  int *rcount = (int *) malloc(npes * sizeof(int));
  int *rdispls = (int *) malloc(npes * sizeof(int));
  int *scount = (int *) malloc(npes * sizeof(int));
  int *sdispls = (int *) malloc(npes * sizeof(int));

  //if (mype==0)
    //gk_startwctimer(vault->timer_3);
  j=0;
  for(i=0; i<npes; i++) {
    if(j<rinfo->n_nbrs && rinfo->nbr_ids[j] == i) {
      rcount[i] = rinfo->ptr[j+1] - rinfo->ptr[j];
      j++;
    } else {
      rcount[i] = 0;
    }

    if(i==0) 
      rdispls[i] = 0;
    else {
      rdispls[i] = rcount[i-1] + rdispls[i-1];
    }
  }

  j=0;
  for(i=0; i<npes; i++) {
    if(j<sinfo->n_nbrs && sinfo->nbr_ids[j] == i) {
      scount[i] = sinfo->ptr[j+1] - sinfo->ptr[j];
      j++;
    } else {
      scount[i] = 0;
    }

    if(i==0) 
      sdispls[i] = 0;
    else {
      sdispls[i] = scount[i-1] + sdispls[i-1];
    }
  }

  MPI_Alltoallv(sinfo->ind, scount, sdispls, MPI_INT, rinfo->ind, rcount, rdispls, MPI_INT, MPI_COMM_WORLD); 

  free(rcount);
  free(rdispls);
  free(scount);
  free(sdispls); 
  #endif  
  //if (mype==0)
    //gk_stopwctimer(vault->timer_3);

	/* Recv. */
	for(i=0; i<rinfo->n_nbrs; i++) {
		MPI_Irecv(
      rinfo->ind+rinfo->ptr[i], 
      rinfo->ptr[i+1]-rinfo->ptr[i], 
      MPI_ADJ_INT, rinfo->nbr_ids[i], tag1, 
      MPI_COMM_WORLD, recv_req+i);
	}

	/* Sends. */
	for(i=0; i<sinfo->n_nbrs; i++) {
		MPI_Isend(
      sinfo->ind+sinfo->ptr[i], 
      sinfo->ptr[i+1]-sinfo->ptr[i], 
      MPI_ADJ_INT, sinfo->nbr_ids[i], tag1, 
      MPI_COMM_WORLD, send_req+i);
	}

	MPI_Waitall(rinfo->n_nbrs, recv_req, statuses);
	MPI_Waitall(sinfo->n_nbrs, send_req, statuses);
  #ifdef isr
  #endif
	free(recv_row);
	free(send_row);
}

/*
n: Total unique nbrs in processor p_i
m: size of unique vtxids - (computed before hand) - xadj_list for local list. 
 */
void SetUpSendInfo(
    SendInfoTypeGK *sinfo, 
    ssize_t m, adj_t *nbr_list, 
    ssize_t *vtxdist,
    AccumPRDetails *req_data, 
    vtxs_t nvtxs) 
{ 

  int npes, mype; 
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  ssize_t offset, until;
	ssize_t i, j, k, l, nunique;

  offset = vtxdist[mype]; 
  until = vtxdist[mype+1];
  /* if (mype==npes-1)
    until = nvtxs;
  else
    until = partition*(mype+1); */

	sinfo->ind = (adj_t *) malloc(m * sizeof(adj_t));
	sinfo->nbr_ids = (int *) malloc(npes * sizeof(int));
	sinfo->ptr = (ssize_t *) malloc((npes+1) * sizeof(ssize_t));

	/*
	   Step 1: Determine which vertex ranks need to be received from other processors. 
     It's the unfortunately populated in the SEND row, but, it's the RECV row. 
	 */
  /* This is what I will recv - unique number of nvtxs, 
  * m (nunique) is traversed to find the vtx ids outside. */
	req_data->nvtxs = 0;
  //printf("mype: %d nbr list: ", mype);
	for(i=0; i<m; i++) {
    //printf(" %"PRId32" ", nbr_list[i]);
		if(nbr_list[i]<offset || nbr_list[i]>=until) {
			req_data->vtx_ids[req_data->nvtxs] = nbr_list[i];		
			sinfo->ind[req_data->nvtxs++] = nbr_list[i];		
		}
	}
  //printf("\n");

	/* j has the count of unique entries. */
	nunique = req_data->nvtxs;

  /* printf("mype: %d in setup sinfo, ind ", mype);
  for (i=0; i<nunique; i++) {
    printf(" %"PRId32" ", sinfo->ind[i]);
  }
  printf("\n"); */

	sinfo->n_nbrs = 0;
	sinfo->ptr[0] = 0;

	ssize_t supposed_until;
	supposed_until = 0;

	for(j=i=0; i<npes; i++) {
		l=j; //some random index to be happy.
		/* if(i==npes-1)
			supposed_until = nvtxs;
		else
			supposed_until = supposed_offset+supposed_chunk; */ 
    supposed_until = vtxdist[i+1];

		for(; j<nunique; j++) {
			if(sinfo->ind[j] >= supposed_until)
				break;
		}

		if(j>l) {
			sinfo->nbr_ids[sinfo->n_nbrs] = i; /* id of the ith nbr. */
      sinfo->n_nbrs++;
			sinfo->ptr[sinfo->n_nbrs] = j; /* i+1th entry gives the end index for nbr i. */
		}

		//supposed_offset = supposed_offset + supposed_chunk;
    //supposed_offset = vtxdist[i+1];
	}
}

/**
 * Requesting rank values for adjancent vertices that are not local. Also send to guys who are away from you. 
 */
void GetValues(
  ssize_t *local_vals, 
  ssize_t *remote_vals, 
  RecvInfoTypeGK *rinfo, 
  SendInfoTypeGK *sinfo) {

  int npes, mype;
  ssize_t i, j;

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  MPI_Request recv_reqs[NPROCS], send_reqs[NPROCS];
  MPI_Status statuses[NPROCS];

  /* 
     Rinfo values are what I want to get, what happens to sinfo? Needs to be sent as well. 
     Do the same isend/irecv - seems to be faster = but why? 
     */

  #ifdef all2allv
  int *rcount = (int *) malloc(npes * sizeof(int));
  int *rdispls = (int *) malloc((npes+1) * sizeof(int));
  int *scount = (int *) malloc(npes * sizeof(int));
  int *sdispls = (int *) malloc((npes+1) * sizeof(int));

  j=0;
  for(i=0; i<npes; i++) {
    if(j<rinfo->n_nbrs && rinfo->nbr_ids[j] == i) {
      rcount[i] = rinfo->ptr[j+1] - rinfo->ptr[j];
      j++;
    } else {
      rcount[i] = 0;
    }

    if(i==0) 
      rdispls[i] = 0;
    else {
      rdispls[i] = rcount[i-1] + rdispls[i-1];
    }
  }

  j=0;
  for(i=0; i<npes; i++) {
    if(j<sinfo->n_nbrs && sinfo->nbr_ids[j] == i) {
      scount[i] = sinfo->ptr[j+1] - sinfo->ptr[j];
      j++;
    } else {
      scount[i] = 0;
    }

    if(i==0) 
      sdispls[i] = 0;
    else {
      sdispls[i] = scount[i-1] + sdispls[i-1];
    }
  }

  MPI_Alltoallv(
    local_vals, rcount, 
    rdispls, MPI_LONG_LONG_INT, 
    remote_vals, scount, 
    sdispls, MPI_LONG_LONG_INT, 
    MPI_COMM_WORLD); 

  free(rcount);
  free(rdispls);
  free(scount);
  free(sdispls);   
  #endif

  /* Post the recvs. */
  for(i=0; i<sinfo->n_nbrs; i++) {
    MPI_Irecv(
      remote_vals+sinfo->ptr[i], 
      sinfo->ptr[i+1]-sinfo->ptr[i], 
      MPI_LONG_LONG_INT, 
      sinfo->nbr_ids[i], 0, 
      MPI_COMM_WORLD, send_reqs+i);
  } 

  /* Send the required values as well. */
  for(i=0; i<rinfo->n_nbrs; i++) {
    MPI_Isend(
      local_vals+rinfo->ptr[i], 
      rinfo->ptr[i+1]-rinfo->ptr[i], 
      MPI_LONG_LONG_INT, 
      rinfo->nbr_ids[i], 0, 
      MPI_COMM_WORLD, recv_reqs+i);
  }

  MPI_Waitall(sinfo->n_nbrs, send_reqs, statuses);
  MPI_Waitall(rinfo->n_nbrs, recv_reqs, statuses); 
  #ifdef isr
  #endif
}

matrix_info_ta *mpitc_ExtractUpperLowerAdjncy(
  gk_graph_t *graph, 
  ssize_t *perm)
{

  int mype, npes, rootp; 
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  ssize_t ej, ejend, nedges; 
  ssize_t hrange, lrange, drange, range;
  ssize_t x_offset, total, vi_val, vj_val, proc;
  ssize_t vi_off, vi_t, vj_t, nunique;

  ssize_t *xadj, *uxadj, *lxadj, *vtxdist;  
  ssize_t *all_uxadj, *all_lxadj, *all_ruxadj, *all_rlxadj; 
  ssize_t *urank_array, *rrank_array, *vi_offset, *remote_vals, *all_perm; 

  vtxs_t nvtxs; 
  int32_t vi, vj, vk, vl, rank;
 
  adj_t *all_uadjncy, *all_ladjncy, *all_ruadjncy, *all_rladjncy; 
  adj_t *adjncy, *uadjncy, *ladjncy, *dadjncy, *unadjncy;
  int32_t *counts; 
  chunk_t *iperm;

  ssize_t *nvtxs_count, *rnvtxs_count, *nvtxs_displs, *rnvtxs_displs; 

  //gk_i32kv_t *cand; 
  gk_i64kv_t *cand;
  matrix_info_ta *mat = (matrix_info_ta *) malloc(sizeof(*mat));

  /* Access to my partition ranks - have to send these to requesting processors. */
  nvtxs = graph->nvtxs;
  xadj = graph->xadj; 
  adjncy = graph->adjncy; 
  vtxdist = graph->vtxdist;

  ssize_t offset = (ssize_t) vtxdist[mype];
  ssize_t until = (ssize_t) vtxdist[mype+1];

  total = xadj[nvtxs];

  dadjncy = (adj_t *) malloc(total*sizeof(adj_t));
  
  /* sorting to find duplicates */
  #ifdef ancy_chngs
    cand = gk_i32kvmalloc(total, "cand");
  #endif
  cand = gk_i64kvmalloc(total, "cand");

  for (vi_t=0; vi_t<total; vi_t++) {
    cand[vi_t].key = adjncy[vi_t];
    cand[vi_t].val = (int64_t) vi_t;
  }

  #ifdef ancy_chngs
    gk_i32kvsorti(total, cand); 
  #endif
  gk_i64kvsorti(total, cand); 

  /* Ineed array here - get unique. */
  /* Has issues when #vtxs < #procs. */
  dadjncy[0] = (adj_t) cand[0].key; 
  for (nunique=vi_t=1; vi_t<total; vi_t++)
    if (cand[vi_t].key!=cand[vi_t-1].key)
      dadjncy[nunique++]=(adj_t)cand[vi_t].key;

  SendInfoTypeGK sinfo;
  RecvInfoTypeGK rinfo;
  AccumPRDetails sdata;

  sdata.vtx_ids = (adj_t *) malloc(nunique*sizeof(adj_t));
  
  /* Step 2: Setup send info - RECV info in reality. */
  SetUpSendInfo(&sinfo, nunique, dadjncy, vtxdist, &sdata, vtxdist[npes]);
  free(dadjncy);
  
  /* Step 3: Setup recv info - SEND info in reality. */
  SetUpRecvInfo(&sinfo, &rinfo);

  /* ineedt */
  remote_vals = (ssize_t *) malloc((sinfo.ptr[sinfo.n_nbrs]) * sizeof(ssize_t));

  sdata.vtx_vals = (ssize_t *) malloc(rinfo.ptr[rinfo.n_nbrs]*sizeof(ssize_t));

  for (vj=vi=0; vi<rinfo.ptr[rinfo.n_nbrs]; vi++)
    sdata.vtx_vals[vj++] = perm[rinfo.ind[vi]-offset];
  
  GetValues((&sdata)->vtx_vals, remote_vals, &rinfo, &sinfo);
  free(sdata.vtx_ids);
  free(sdata.vtx_vals);

  /* ineedt */
  for (vi_t=0, vj_t=0; vi_t<total; vi_t++) {
    vk = cand[vi_t].key;
    if (vk>=offset && vk<until)
      adjncy[cand[vi_t].val] = (adj_t) perm[(ssize_t) vk-offset];
    else {
      if (vk!=sinfo.ind[vj_t])
        vj_t++;
      adjncy[cand[vi_t].val] = (adj_t) remote_vals[vj_t];
    }
  }
  free(cand);
  free(remote_vals);

  ssize_t *lcounts;

  nvtxs_count   = (chunk_t *) malloc(npes*sizeof(chunk_t));
  rnvtxs_count  = (chunk_t *) malloc(npes*sizeof(chunk_t));
  nvtxs_displs  = (chunk_t *) malloc((npes+1)*sizeof(chunk_t));
  rnvtxs_displs = (chunk_t *) malloc((npes+1)*sizeof(chunk_t));
  lcounts       = (ssize_t *) malloc(npes*sizeof(ssize_t));

  memset(nvtxs_count, 0, npes*sizeof(ssize_t)); 
  memset(lcounts, 0, npes*sizeof(ssize_t));

  //mat->uvtxs_displs = rnvtxs_displs;

  /* SUB MATRIX INFO - 1 */
  rootp = sqrt((double) npes);

  /* ult */ 
  urank_array = (ssize_t *) malloc(nvtxs*rootp*sizeof(ssize_t));
  vi_offset = (ssize_t *) malloc(nvtxs*sizeof(ssize_t));

  /* Figuring out which set of rootp 
   * procs a vi will be assigned to - 
   * so I can send adj vertices to these procs. */ 
  for (vi=0; vi<nvtxs; vi++){
    vi_val = (perm[vi]/BLK_SZ)%rootp;
    vi_offset[vi] = (ssize_t) nvtxs_count[vi_val*rootp];
    nvtxs_count[vi_val*rootp]++;
    lcounts[vi_val*rootp]++;
  }

  for (vi=0; vi<npes; vi++) {
    nvtxs_count[vi] = nvtxs_count[(vi/rootp)*rootp];
    lcounts[vi] = lcounts[(vi/rootp)*rootp];
  }
  
  free(sinfo.ind);
  free(sinfo.nbr_ids);
  free(sinfo.ptr);
  free(rinfo.ind);
  free(rinfo.nbr_ids);
  free(rinfo.ptr);

  MPI_Alltoall(
    nvtxs_count, 1,
    MPI_CHUNK_INT,
    rnvtxs_count, 1,
    MPI_CHUNK_INT,
    MPI_COMM_WORLD);

  nvtxs_displs[0] = 0;
  for (vi=1; vi<=npes; vi++)
    nvtxs_displs[vi] = nvtxs_displs[vi-1]+nvtxs_count[vi-1]; 

  rnvtxs_displs[0] = 0;
  for (vi=1; vi<=npes; vi++)
    rnvtxs_displs[vi] = rnvtxs_displs[vi-1]+rnvtxs_count[vi-1]; 
  free(rnvtxs_count);

  /* ult */
  ssize_t p;
  for (vi=0; vi<nvtxs; vi++) {
    p = perm[vi];
    vi_val = (p/BLK_SZ)%rootp;
    for (vj=0; vj<rootp; vj++) {
      urank_array[nvtxs_displs[vi_val*rootp+vj]++] = p;
    }
  } 

  nvtxs_displs[0] = 0;
  for (vi=1; vi<=npes; vi++)
    nvtxs_displs[vi] = nvtxs_displs[vi-1]+nvtxs_count[vi-1];
  free(nvtxs_count);
  
  ssize_t *lcounts_cp;
  lcounts_cp = (ssize_t *) malloc(npes*sizeof(ssize_t));
  memcpy(lcounts_cp, lcounts, npes*sizeof(ssize_t));

  rrank_array = ancy_all2allv_xadj(lcounts, urank_array, MPI_COMM_WORLD);
  free(urank_array);

  /* ult */
  all_uxadj = (ssize_t *) malloc((nvtxs_displs[npes]+1)*sizeof(ssize_t));
  all_lxadj = (ssize_t *) malloc((nvtxs_displs[npes]+1)*sizeof(ssize_t));
  memset(all_uxadj, 0, ((nvtxs_displs[npes]+1)*sizeof(ssize_t)));
  memset(all_lxadj, 0, ((nvtxs_displs[npes]+1)*sizeof(ssize_t)));

  for (vi_t=0; vi_t<nvtxs; vi_t++) {
    p = perm[vi_t];
    vi_val = (p/BLK_SZ)%rootp;
    vi_off = vi_offset[vi_t];
    for (vj_t=xadj[vi_t]; vj_t<xadj[vi_t+1]; vj_t++) {
      rank = adjncy[vj_t];
      vj_val = (rank/BLK_SZ)%rootp; 
      proc = vi_val*rootp + vj_val; 

      if (p<rank) {
        all_uxadj[(ssize_t) nvtxs_displs[proc]+vi_off]++; 
      } else if (p>rank) {
        all_lxadj[(ssize_t) nvtxs_displs[proc]+vi_off]++;
      }
    }
  }

  //all_ruxadj = (ssize_t *) malloc((rnvtxs_displs[npes]+1)*sizeof(ssize_t));
  //all_rlxadj = (ssize_t *) malloc((rnvtxs_displs[npes]+1)*sizeof(ssize_t));

  memcpy(lcounts, lcounts_cp, npes*sizeof(ssize_t));

  all_ruxadj = ancy_all2allv_xadj(lcounts, all_uxadj, MPI_COMM_WORLD);
  all_rlxadj = ancy_all2allv_xadj(lcounts_cp, all_lxadj, MPI_COMM_WORLD);

  /* all_uxadj, all_lxadj makecsr */
  MAKECSR(vi, nvtxs_displs[npes], all_uxadj);
  MAKECSR(vi, nvtxs_displs[npes], all_lxadj);

  MAKECSR(vi, rnvtxs_displs[npes], all_ruxadj);
  MAKECSR(vi, rnvtxs_displs[npes], all_rlxadj);
 
  /* ult */
  all_uadjncy = (adj_t *) malloc(all_uxadj[nvtxs_displs[npes]]*sizeof(adj_t));
  all_ladjncy = (adj_t *) malloc(all_lxadj[nvtxs_displs[npes]]*sizeof(adj_t));

  ssize_t test_l, t_disp;
  t_disp = all_lxadj[nvtxs_displs[npes]];
  for (vi_t=0; vi_t<nvtxs; vi_t++) {
    p = perm[vi_t];
    vi_val = (p/BLK_SZ)%rootp;
    vi_off = vi_offset[vi_t];
    for (vj_t=xadj[vi_t]; vj_t<xadj[vi_t+1]; vj_t++) {
      rank = adjncy[vj_t];
      vj_val = (rank/BLK_SZ)%rootp; 
      proc = vi_val*rootp + vj_val; 

      if (p<rank) {
        test_l = all_uxadj[(ssize_t) nvtxs_displs[proc]+vi_off];
        all_uadjncy[test_l] = rank;
        all_uxadj[(ssize_t) nvtxs_displs[proc]+vi_off]++;
      } else if (p>rank) {
        test_l = all_lxadj[(ssize_t) nvtxs_displs[proc]+vi_off];
        all_ladjncy[test_l] = rank;
        all_lxadj[(ssize_t) nvtxs_displs[proc]+vi_off]++;
      }
    }
  }

  SHIFTCSR(vi, nvtxs_displs[npes], all_uxadj); 
  SHIFTCSR(vi, nvtxs_displs[npes], all_lxadj);

  ssize_t ru_sz = all_ruxadj[rnvtxs_displs[npes]];
  ssize_t rl_sz = all_rlxadj[rnvtxs_displs[npes]];

  //all_ruadjncy = (adj_t *) malloc(ru_sz*sizeof(adj_t));
  //all_rladjncy = (adj_t *) malloc(rl_sz*sizeof(adj_t));

  /* lcounts -- counts for each processor id. */
  for (vi=0; vi<npes; vi++) 
    lcounts[vi] = all_uxadj[nvtxs_displs[vi+1]]-all_uxadj[nvtxs_displs[vi]];  
  free(all_uxadj); 
  all_ruadjncy = ancy_all2allv(lcounts, all_uadjncy, MPI_COMM_WORLD);
  free(all_uadjncy);

  for (vi=0; vi<npes; vi++)
    lcounts[vi] = all_lxadj[nvtxs_displs[vi+1]]-all_lxadj[nvtxs_displs[vi]];
  free(all_lxadj);
  all_rladjncy = ancy_all2allv(lcounts, all_ladjncy, MPI_COMM_WORLD);
  free(all_ladjncy);

  free(lcounts);
  free(lcounts_cp);
  free(nvtxs_displs);

  /* arrays of interest - rrank_array, all_ruxadj, 
   * all_rlxadj, allruadjncy, allrladjncy, rnvtxs_displs. */
  
  /* sortt */
  /* sort by rank. */ 
  for (hrange=0, lrange=0, vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    if ((p=rrank_array[vi_t]) > hrange)
      hrange = p;
    else if (p < lrange)
      lrange = p;
  }
  drange = hrange-lrange;
  drange++; 

  counts = gk_i32smalloc(drange+1, 0, "counts");
  for (vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++)
    counts[rrank_array[vi_t]-lrange]++;
  MAKECSR(vi, drange, counts);

  //iperm = gk_i32smalloc(rnvtxs_displs[npes], 0, "iperm");
  iperm = (chunk_t *) malloc(rnvtxs_displs[npes]*sizeof(chunk_t));
  for (vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++)
    iperm[counts[rrank_array[vi_t]-lrange]++] = (chunk_t) vi_t;

  gk_free((void **)&counts, LTERM);
  free(rrank_array);

  mat->uxadj    = uxadj   = (ssize_t *) malloc((rnvtxs_displs[npes]+1)*sizeof(ssize_t));
  mat->uadjncy  = uadjncy = (adj_t *) malloc(all_ruxadj[rnvtxs_displs[npes]]*sizeof(adj_t));
  
  uxadj[0] = nedges = 0;
  //printf("mype: %d rnvtxs displs %"PRId32". \n", mype, rnvtxs_displs[npes]);
  for (vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    vj_t = iperm[vi_t];
    for (ej=all_ruxadj[vj_t], ejend=all_ruxadj[vj_t+1]; ej<ejend; ej++) {
      vk = all_ruadjncy[ej];/* keep only the upper part */
      uadjncy[nedges++] = vk;
    }
    uxadj[vi_t+1] = nedges;

    #ifdef ancy_chgs
    if (nedges-uxadj[vi_t] > 1)
      gk_i32sorti(nedges-uxadj[vi_t], uadjncy+uxadj[vi_t]);  /* sort adjncy list */
    #endif
    if (nedges-uxadj[vi_t] > 1)
      gk_i64sorti(nedges-uxadj[vi_t], uadjncy+uxadj[vi_t]);  /* sort adjncy list */
  }
  free(all_ruxadj);
  free(all_ruadjncy);
  
  mat->lxadj    = lxadj   = (ssize_t *) malloc((rnvtxs_displs[npes]+1)*sizeof(ssize_t));
  mat->ladjncy  = ladjncy = (adj_t *) malloc(all_rlxadj[rnvtxs_displs[npes]]*sizeof(adj_t));

  lxadj[0] = nedges = 0;
  for (vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    vj = iperm[vi_t];
    for (ej=all_rlxadj[vj], ejend=all_rlxadj[vj+1]; ej<ejend; ej++) {
      vk = all_rladjncy[ej];/* keep only the upper part */
      ladjncy[nedges++] = vk;
    }
    lxadj[vi_t+1] = nedges;

    #ifdef ancy_chgs
    if (nedges-lxadj[vi_t] > 1)
      gk_i32sorti(nedges-lxadj[vi_t], ladjncy+lxadj[vi_t]);  /* sort adjncy list */
    #endif
    if (nedges-lxadj[vi_t] > 1)
      gk_i64sorti(nedges-lxadj[vi_t], ladjncy+lxadj[vi_t]);  /* sort adjncy list */
    
  }
  free(all_rlxadj);
  free(all_rladjncy);

  mat->nvtxs = rnvtxs_displs[npes];

  ssize_t nnz_count; 
  adj_t *lnnz_idx, *unnz_idx; 

  for (nnz_count=0, vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    if (lxadj[vi_t+1]-lxadj[vi_t])
      nnz_count++;
  }

  lnnz_idx = (adj_t *) malloc(nnz_count*sizeof(adj_t));
  for (vj=0, vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    if (lxadj[vi_t+1]-lxadj[vi_t])
      lnnz_idx[vj++] = vi_t;
  }

  mat->lnnz_idx = lnnz_idx; 
  mat->lnvtxs = nnz_count; 
  
  for (nnz_count=0, vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    if (uxadj[vi_t+1]-uxadj[vi_t])
      nnz_count++;
  }

  unnz_idx = (adj_t *) malloc(nnz_count*sizeof(adj_t));
  for (vj=0, vi_t=0; vi_t<rnvtxs_displs[npes]; vi_t++) {
    if (uxadj[vi_t+1]-uxadj[vi_t])
      unnz_idx[vj++] = vi_t;
  }

  mat->unnz_idx = unnz_idx; 
  mat->unvtxs = nnz_count; 

  free(iperm);
  free(vi_offset);

  /*_ancy new*/
  free(rnvtxs_displs);

  return mat;
}

ssize_t *mpitc_DistrCountingSort(gk_graph_t *graph)
{
  int mype, npes;
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 
  
  int32_t nvtxs; 
  ssize_t *xadj; 
  adj_t *adjncy;
  
  ssize_t *deg, *perm; 
  ssize_t *lcount, *loffset, *gcount, *lgcount;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj; 

  int32_t vi, vj, vk; 
  ssize_t x_offset, range, g_range; 

  deg = (ssize_t *) malloc(nvtxs*sizeof(ssize_t));
  perm = (ssize_t *) malloc(nvtxs*sizeof(ssize_t));
  for (range=0, vi=0; vi<nvtxs; vi++) {
    deg[vi] = xadj[vi+1]-xadj[vi];
    if (range<deg[vi])
      range = deg[vi];
  }

  /*All reduce max value to find range.*/
  MPI_Allreduce(&range, &g_range, 1, MPI_LONG_LONG_INT, 
        MPI_MAX, MPI_COMM_WORLD);
  g_range++;

  loffset = (ssize_t *) malloc(nvtxs*sizeof(ssize_t));
  lcount = (ssize_t *) malloc(g_range*sizeof(ssize_t));
  gcount = (ssize_t *) malloc(g_range*sizeof(ssize_t));
  lgcount = (ssize_t *) malloc(g_range*sizeof(ssize_t));

  memset(lcount, 0, g_range*sizeof(ssize_t));
  memset(gcount, 0, g_range*sizeof(ssize_t));

  for (vi=0; vi<nvtxs; vi++) {
    loffset[vi] = lcount[deg[vi]];
    lcount[deg[vi]]++;
  }

  MPI_Scan(lcount, gcount, g_range, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  if(mype==npes-1)
    memcpy(lgcount, gcount, g_range*sizeof(ssize_t));

  MPI_Bcast(lgcount, g_range, MPI_LONG_LONG_INT, npes-1, MPI_COMM_WORLD);

  ssize_t sum=0; 
  gcount[0] -= lcount[0];
  for(vi=1; vi<g_range; vi++) {
    sum += lgcount[vi-1];
    gcount[vi] += sum-lcount[vi];
  }

  for (vi=0; vi<nvtxs; vi++) {
    perm[vi] = gcount[deg[vi]] + loffset[vi];
  }

  free(lgcount);
  free(loffset);
  free(lcount);
  free(gcount);
  free(deg);
  /* printf("mype: %d permutation array: ", mype);
  for (vi=0; vi<nvtxs; vi++)
    printf(" %zd ", perm[vi]);
  printf("\n"); */

  return perm; 
}

gk_graph_t *mpitc_cyclicGraphDistr(gk_graph_t *graph) {

  int npes, mype;
  
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
 
  gk_graph_t *t_graph; 
  int32_t vi, nvtxs, proc_id;
  ssize_t offset, coffset, until, vj_t;

  adj_t *adjncy, *cadjncy, *radjncy;
  ssize_t *gcount, *cdisp, *lcount, *rcount, *rcdisp; 
  ssize_t *xadj, *vtxdist, *cxadj, *rxadj;

  t_graph = malloc(sizeof(*t_graph));
  vtxdist = graph->vtxdist;
  offset  = vtxdist[mype];
  until   = vtxdist[mype+1];
  nvtxs   = graph->nvtxs; 

  xadj    = graph->xadj; 
  adjncy  = graph->adjncy;

  rcount  = (ssize_t *) malloc(npes*sizeof(ssize_t));
  gcount  = (ssize_t *) malloc((npes+1)*sizeof(ssize_t));
  cdisp   = (ssize_t *) malloc((npes+1)*sizeof(ssize_t));
  rcdisp  = (ssize_t *) malloc((npes+1)*sizeof(ssize_t));
  lcount  = (ssize_t *) malloc(npes*sizeof(ssize_t));
  memset(lcount, 0, npes*sizeof(ssize_t));

  for (proc_id=0, vi=offset; vi<until; vi++) {
    proc_id = vi%npes;
    lcount[proc_id]++;
  }
  free(graph->vtxdist);

  MPI_Scan(lcount, gcount, npes, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Alltoall(lcount, 1, MPI_LONG_LONG_INT, rcount, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
  memcpy(cdisp, lcount, npes*sizeof(ssize_t));
  MAKECSR(vi, npes, cdisp);
  memcpy(rcdisp, rcount, npes*sizeof(ssize_t));
  MAKECSR(vi, npes, rcdisp);

  MPI_Bcast(gcount, npes, MPI_LONG_LONG_INT, npes-1, MPI_COMM_WORLD);
  MAKECSR(vi, npes, gcount);

  t_graph->vtxdist  = gcount;
  t_graph->nvtxs    = gcount[mype+1]-gcount[mype];

  cxadj     = (ssize_t *) malloc((nvtxs+1)*sizeof(ssize_t));
  for (vi=0; vi<nvtxs; vi++) {
    proc_id = (vi+offset)%npes;
    cxadj[cdisp[proc_id]] = (xadj[vi+1]-xadj[vi]);
    cdisp[proc_id]++;
  }

  SHIFTCSR(vi, npes, cdisp);

  rxadj = ancy_all2allv_xadj(lcount, cxadj, MPI_COMM_WORLD);
  t_graph->xadj = rxadj; 

  MAKECSR(vi, nvtxs, cxadj);
  MAKECSR(vi, rcdisp[npes], rxadj);

  cadjncy = (adj_t *) malloc(xadj[nvtxs]*sizeof(adj_t));
  
  for (vi=0; vi<npes; vi++)
    lcount[vi] = cxadj[cdisp[vi+1]]-cxadj[cdisp[vi]];

  for (vi=0; vi<nvtxs; vi++) {
    proc_id = (vi+offset)%npes;
    coffset = cxadj[cdisp[proc_id]];
    for (vj_t=xadj[vi]; vj_t<xadj[vi+1]; vj_t++) {
      cadjncy[coffset++] = gcount[adjncy[vj_t]%npes]+(adjncy[vj_t]/npes);
    }
    cdisp[proc_id]++;
  }

  radjncy = ancy_all2allv(lcount, cadjncy, MPI_COMM_WORLD);
  t_graph->adjncy = radjncy;

  free(lcount);
  free(rcount);
  free(cdisp);
  free(rcdisp);
  free(cxadj);
  free(cadjncy);

  return t_graph;
}

matrix_info_tn *extractFromMat(matrix_info_ta *mat, int uflag) {

  matrix_info_tn *matrix_info; 
  matrix_info = malloc(sizeof(*matrix_info));

  if (uflag) {
    
    matrix_info->g_nvtxs = mat->nvtxs;
    matrix_info->xadj = mat->uxadj;
    matrix_info->nvtxs = mat->unvtxs;
    matrix_info->adjncy = mat->uadjncy;
    matrix_info->nnz_idx = mat->unnz_idx;
    //matrix_info->vi_chunk = mat->uvtxs_displs;
    //matrix_info->reverse_idx = mat->ureverse_idx;
  } else {
    
    matrix_info->g_nvtxs = mat->nvtxs;
    matrix_info->xadj = mat->lxadj;
    matrix_info->nvtxs = mat->lnvtxs;
    matrix_info->adjncy = mat->ladjncy;
    matrix_info->nnz_idx = mat->lnnz_idx;
    //matrix_info->vi_chunk = mat->uvtxs_displs;
    //matrix_info->reverse_idx = mat->lreverse_idx;
  }

  return matrix_info; 
}

/*************************************************************************/
/*! The hash-map-based triangle-counting routine that uses the JIK
  triangle enumeration scheme.
  This is the version used to get the first set of experiments and it
  does not use the tmap direct addressing array for the tail.
 */
/*************************************************************************/
int64_t mpitc_MapJIK2D(params_t *params, vault_t *vault)
{

  double ttimer=0.;

  int npes, mype;
  double start, end;

  ssize_t *perm, until, offset, *xadj, *vtxdist;
  ssize_t ei, eiend, ej, ejend;
 
  adj_t vi, vj, vk, vl;
  vtxs_t nvtxs;
  adj_t *adjncy; 
  int64_t ntriangles, g_ntriangles;

  gk_graph_t *graph=NULL;
  matrix_info_ta *mat=NULL; 
  blob_info_t *rblob=NULL, *cblob=NULL;

  /* MPI decl. */
  ntriangles = 0;
  g_ntriangles = 0;

  /* Obtain the number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  int rootp = (int) sqrt(npes);
  double *timers = (double *) malloc (rootp * sizeof(double)); 

  /* ttimer */
  MPI_Barrier(MPI_COMM_WORLD);
  if (mype==0)
    ttimer = MPI_Wtime();

  /* Adding timers for pp. */
  /* timer 1 */
  MPI_Barrier(MPI_COMM_WORLD);
  if (mype==0)
    gk_startwctimer(vault->timer_1);

  graph = vault->graph;

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;

  /* Cyclic distribution. */
  gk_graph_t *r_graph = NULL; 
  if (npes>1) {

    r_graph = mpitc_cyclicGraphDistr(graph);

    free(xadj);
    free(adjncy);

    graph->nvtxs    = r_graph->nvtxs; 
    graph->xadj     = r_graph->xadj;
    graph->adjncy   = r_graph->adjncy;
    graph->vtxdist  = r_graph->vtxdist; 
  }

  if (npes>1)
    perm = mpitc_DistrCountingSort(graph);
  else
    perm = gk_sskvipermi(nvtxs, graph->xadj);

  mat = mpitc_ExtractUpperLowerAdjncy(graph, perm);

  free(perm);
  free(graph->xadj);
  free(graph->adjncy);
  free(graph->vtxdist);
  free(graph);

  if (npes>1)
    free(r_graph);

  int *transpose = (int *) malloc(npes*sizeof(int));
  for (vi=0; vi<rootp; vi++) 
    for (vj=0; vj<rootp; vj++) 
      transpose[vi+vj*rootp] = vi*rootp+vj;

  matrix_info_tn *a_matrix_info, *b_matrix_info, *c_matrix_info;

  a_matrix_info = extractFromMat(mat, 1);

  /* Need to transpose b_matrix_info. */
  b_matrix_info = extractFromMat(mat, 1);
  rblob = transposeSubMatrix(b_matrix_info, MPI_COMM_WORLD, transpose[mype]);

  /* Need to also transpose c_matrix_info. */
  c_matrix_info = extractFromMat(mat, 0);
  cblob = transposeSubMatrix(c_matrix_info, MPI_COMM_WORLD, transpose[mype]);

  free(mat->lxadj);
  free(mat->ladjncy);
  free(mat->lnnz_idx);
  free(b_matrix_info);
  free(c_matrix_info);

  b_matrix_info = mpitc_unpackBlobBlockCyclic(rblob);
  c_matrix_info = mpitc_unpackBlobBlockCyclic(cblob);

  free(transpose);
  free(rblob->blob);
  free(cblob->blob);
  free(rblob);
  free(cblob);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mype==0)
    gk_stopwctimer(vault->timer_1);

  if (mype==0)
    gk_startwctimer(vault->timer_2);
  
  /* Cannon algorithm communication. - square-root p stages. */
  // Finding the square-root of the #procs over here.
  int row_rank, col_rank; 
  int row_color = mype/rootp;
  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, row_color, mype, &row_comm);
  MPI_Comm_rank(row_comm, &row_rank);

  int col_color = mype%rootp;
  MPI_Comm col_comm; 
  MPI_Comm_split(MPI_COMM_WORLD, col_color, mype, &col_comm);
  MPI_Comm_rank(col_comm, &col_rank);

  /* tct -- do we need this now? */
  //start = MPI_Wtime();
  int32_t g_off = BLK_SZ*rootp;
  for (vi=0; vi<c_matrix_info->xadj[c_matrix_info->g_nvtxs]; vi++)
    c_matrix_info->adjncy[vi] = 
      (c_matrix_info->adjncy[vi]/g_off)*BLK_SZ + 
      (c_matrix_info->adjncy[vi]%BLK_SZ); 

  rblob = 
    initialShiftBlockCyclic(
      a_matrix_info, row_comm, 
      row_color, row_rank);
  MPI_Barrier(row_comm);

  cblob = 
    initialShiftBlockCyclic(
      b_matrix_info, col_comm, 
      col_color, col_rank);
  MPI_Barrier(col_comm);

  free(a_matrix_info->nnz_idx);
  free(a_matrix_info->xadj);
  free(a_matrix_info->adjncy);
  free(a_matrix_info);

  free(b_matrix_info->nnz_idx);
  free(b_matrix_info->xadj);
  free(b_matrix_info->adjncy);
  free(b_matrix_info);

  ntriangles += 
    mpitc_kernelJIKBlockCyclic7(
      0, c_matrix_info, rblob, cblob);

  int shift;
  for (shift=1; shift<rootp; shift++) {

    rblob = 
      mpitc_upleftshift(row_comm, row_color, row_rank, rblob);
    cblob = 
      mpitc_upleftshift(col_comm, col_color, col_rank, cblob);

    ntriangles += mpitc_kernelJIKBlockCyclic7(shift, c_matrix_info, rblob, cblob);
  }

  MPI_Reduce(&ntriangles, &g_ntriangles, 1, MPI_LONG_LONG, MPI_SUM, 0,
           MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mype==0)
    gk_stopwctimer(vault->timer_2);
  
  free(mat);
  free(rblob->blob);
  free(cblob->blob);
  free(rblob);
  free(cblob);

  free(c_matrix_info->nnz_idx);
  free(c_matrix_info->xadj);
  free(c_matrix_info->adjncy);
  free(c_matrix_info);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mype==0)
    printf("Overall separate timer: %f. \n", (MPI_Wtime()-ttimer));

  return g_ntriangles;
}

matrix_info_tn *mpitc_unpackBlobBlockCyclic(blob_info_t *blob_info) 
{

  int mype, npes, rootp, vi, vj;  
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 

  unsigned char* rrblob;
  ssize_t rmembytes;
  rrblob = blob_info->blob;
  rmembytes = blob_info->membytes;

  //printf("mype: %d rmembytes: %"PRId64". \n", mype, rmembytes);

  /* Remote block details - we expect blocks to come from 
     every other row except the first rootp rows as they
     do not get transferred around. */
  ssize_t bytecnt = 0;
  ssize_t xadj_counter=0, nnz_counter=0, nnz=0;
  int32_t nvtxs, g_nvtxs; 

  /* Chunk size. */
  //memcpy(r_chunksizes, rrblob+bytecnt, (rootp+1)*sizeof(int32_t));
  memcpy(&nvtxs, rrblob+bytecnt, 1*sizeof(int32_t));
  bytecnt += sizeof(int32_t);

  memcpy(&g_nvtxs, rrblob+bytecnt, 1*sizeof(int32_t));
  bytecnt += sizeof(int32_t);

  xadj_counter = g_nvtxs;
  //printf("mype: %d xadj counter: %zd. \n", mype, xadj_counter);

  nnz_counter = *((ssize_t *) (rrblob+bytecnt+xadj_counter*sizeof(ssize_t)));
  //printf("mype: %d nnz counter: %zd. \n", mype, nnz_counter);

  //xadj_counter += 1;
  ssize_t *r_xadj = (ssize_t *) malloc((xadj_counter+1)*sizeof(ssize_t));
  memcpy(r_xadj, rrblob+bytecnt, (xadj_counter+1)*sizeof(ssize_t));
  bytecnt += (xadj_counter+1)*sizeof(ssize_t);

  adj_t *r_adjncy = (adj_t *) malloc(nnz_counter*sizeof(adj_t));
  memcpy(r_adjncy, rrblob+bytecnt, nnz_counter*sizeof(adj_t));
  bytecnt += nnz_counter*sizeof(adj_t);

  adj_t *r_nnz_idx = (adj_t *) malloc(nvtxs*sizeof(adj_t));
  memcpy(r_nnz_idx, rrblob+bytecnt, nvtxs*sizeof(adj_t));
  bytecnt += nvtxs*sizeof(adj_t);

  /* ssize_t *r_reverse_idx = (ssize_t *) malloc(g_nvtxs*sizeof(ssize_t));
  memcpy(r_reverse_idx, rrblob+bytecnt, g_nvtxs*sizeof(ssize_t));
  bytecnt += g_nvtxs*sizeof(ssize_t); */

  matrix_info_tn *r_matrix_info = (matrix_info_tn *) malloc(sizeof(*r_matrix_info)); 
  r_matrix_info->nvtxs = nvtxs; 
  r_matrix_info->g_nvtxs = g_nvtxs; 
  r_matrix_info->xadj = r_xadj;
  r_matrix_info->adjncy = r_adjncy;
  //r_matrix_info->reverse_idx = r_reverse_idx; 
  r_matrix_info->nnz_idx = r_nnz_idx;

  //printf("mype: %d REACHED HERE!. \n", mype);

  return r_matrix_info;
}

unsigned char *blockedBlobSendRecv(
    ssize_t lcount, ssize_t rcount, 
    int snd, int rcv, 
    unsigned char* larray, 
    MPI_Comm comm) {

  size_t i, j, k;
  ssize_t max_lcnt, c_lcnt;
  int size, mype;
  int global_size, global_mype; 

  blk_t *blk_cnts, *rblk_cnts;
  ssize_t *blk_displs, *rblk_displs;
  int *max_slcnts, *max_sldispls, *max_rlcnts, *max_rldispls;
  unsigned char *rarray;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_mype);

  max_lcnt = ceil1(lcount/(double) XADJ_ADM_MAX);

  /* This is for size procs times number of 
    alltoallv we will have to make. */
  blk_cnts = (blk_t *) malloc(max_lcnt*sizeof(blk_t));
  blk_displs = (ssize_t *) malloc((max_lcnt+1)*sizeof(ssize_t));
  memset(blk_cnts, 0, sizeof(blk_t)*max_lcnt);

  ssize_t accum;
  for (j=0; j<max_lcnt; j++) {
    if ((lcount-XADJ_ADM_MAX) > 0) {
      blk_cnts[j] = XADJ_ADM_MAX;
      lcount -= XADJ_ADM_MAX;
      accum += XADJ_ADM_MAX;
    } else { /* Should ideally be the last blk with entries. */
      blk_cnts[j] = lcount;
      lcount = 0;
    }
  }
  accum = 0;

  blk_displs[0]=0;
  for (i=1; i<=max_lcnt; i++) 
    blk_displs[i] = (ssize_t) blk_cnts[i-1]+blk_displs[i-1];
   
  MPI_Status status; 
  ssize_t max_rlcnt;
  /* First, tell the recvr that you will be sending 
    max_cnts so as to init an array of that size. */
  double start, ltime, mtime, gtime, atime;
 
  MPI_Send(&max_lcnt, 1, MPI_LONG_LONG_INT,
      snd, 124, comm);
  MPI_Recv(&max_rlcnt, 1, MPI_LONG_LONG_INT,
      rcv, 124, comm, &status);

  rblk_cnts = (blk_t *) malloc(max_rlcnt*sizeof(blk_t));

  MPI_Send(blk_cnts, max_lcnt, MPI_INT, 
      snd, 125, comm);
  MPI_Recv(rblk_cnts, max_rlcnt, MPI_INT,
      rcv, 125, comm, &status);

  rblk_displs = (ssize_t *) malloc((max_rlcnt+1)*sizeof(ssize_t));
  rblk_displs[0] = 0;
  for (i=1; i<=max_rlcnt; i++)
    rblk_displs[i] = (ssize_t) rblk_cnts[i-1]+rblk_displs[i-1];

  rarray = (unsigned char *) malloc(rblk_displs[max_rlcnt]*sizeof(int32_t));

  MPI_Request sreq[100], rreq[100];
  MPI_Status statuses[100];

  for (i=0; i<max_lcnt; i++) {
    MPI_Isend((void *)((int32_t *) larray+blk_displs[i]), blk_cnts[i], 
        MPI_INT, snd, (mype+i), comm, sreq+i);
  }

  for (i=0; i<max_rlcnt; i++) {
    MPI_Irecv((void *)((int32_t *) rarray+rblk_displs[i]), rblk_cnts[i], 
        MPI_INT, rcv, (rcv+i), comm, rreq+i);
  }
  
  MPI_Waitall(max_lcnt, sreq, statuses);
  MPI_Waitall(max_rlcnt, rreq, statuses);

  free(blk_cnts);
  free(rblk_cnts);
  free(blk_displs);
  free(rblk_displs);

  return rarray;
}

ssize_t mpitc_blobBytesBlockCyclic(matrix_info_tn *a_matrix_info) 
{
  int mype, npes, rootp, i;  
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 
  //rootp = (int) sqrt(npes);
  rootp = npes;

  ssize_t nbytes32, nbytes64, membytes;
  int32_t nvtxs; 
  ssize_t nnz;

  nbytes64 = 0;
  nbytes32 = 0;
  /* 1: Send over number of bytes to each processor. */
  /* to account for the xadj - 
     nvtxs, vi offset, vj offset per processor. */
  /* nnz to be sent. */
  //nbytes64 = a_matrix_info->vi_chunk[rootp]+1;/* nvtxs per processor. */
  nbytes64 = a_matrix_info->g_nvtxs+1; //+a_matrix_info->g_nvtxs; /* nvtxs per processor. */
  //nbytes32 = a_matrix_info->xadj[a_matrix_info->g_nvtxs];
  nbytes64 += a_matrix_info->xadj[a_matrix_info->g_nvtxs] + a_matrix_info->nvtxs; /* adjncy + nnz_idx */
  nbytes32 = 1+1; /* nnz to be sent. */
  membytes = nbytes32 + nbytes64*2;

  return membytes;
}

blob_info_t *mpitc_packBlobBlockCyclic(
matrix_info_tn *a_matrix_info) 
{
  int mype, npes;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 

  vtxs_t nvtxs, g_nvtxs;
  ssize_t nbytes32, nbytes64, membytes, bytecnt;
  ssize_t nnz; 

  ssize_t *xadj;
  adj_t *adjncy;
  unsigned char *blob; 

  /* New. */
  adj_t *nnz_idx;

  ssize_t n_ints; 
  
  /* #reverse idx global vtx, #nnz idx actual nnz present. */
  g_nvtxs   = a_matrix_info->g_nvtxs;
  nvtxs     = a_matrix_info->nvtxs; 
  xadj      = a_matrix_info->xadj;
  adjncy    = a_matrix_info->adjncy;
  nnz_idx   = a_matrix_info->nnz_idx;

  membytes = mpitc_blobBytesBlockCyclic(a_matrix_info);
  
  n_ints = membytes;
  //printf("mype: %d n ints: %zd %zd, %zd, %zd. \n", mype, membytes, sizeof(int32_t), sizeof(int64_t), sizeof(ssize_t));

  blob = (unsigned char*) malloc(n_ints*sizeof(int32_t));

  bytecnt = 0;

  /* x-axis dim in the chunk. */
  memcpy(blob+bytecnt, &nvtxs, 1*sizeof(vtxs_t));
  bytecnt+=1*sizeof(vtxs_t);

  memcpy(blob+bytecnt, &g_nvtxs, 1*sizeof(vtxs_t));
  bytecnt+=1*sizeof(vtxs_t);

  /* xadj for the chunk. */
  memcpy(blob+bytecnt, xadj, (g_nvtxs+1)*sizeof(ssize_t));
  bytecnt+=(g_nvtxs+1)*sizeof(ssize_t);

  /* nnz to be sent - for mem counting. */
  nnz = xadj[g_nvtxs];

  /* adjncy for the chunk - uses nnz. */
  memcpy(blob+bytecnt, adjncy, nnz*sizeof(adj_t));
  bytecnt+=nnz*sizeof(adj_t);

  //printf("mype: %d bytecnt: %zd membytes: %zd. \n", mype, bytecnt, membytes);
  memcpy(blob+bytecnt, nnz_idx, nvtxs*sizeof(adj_t));
  bytecnt+=nvtxs*sizeof(adj_t);

  //memcpy(blob+bytecnt, reverse_idx, g_nvtxs*sizeof(ssize_t));
  //bytecnt+=g_nvtxs*sizeof(ssize_t);

  blob_info_t *blob_info;
  blob_info = (blob_info_t *) malloc(sizeof(*blob_info));
  //blob_info->membytes = membytes;
  blob_info->membytes = n_ints;
  blob_info->blob = blob;

  return blob_info; 
}
 
blob_info_t *transposeSubMatrix(
    matrix_info_tn *a_matrix_info, MPI_Comm comm, 
    int r) 
{

  int mype, npes, rootp;  
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 
  rootp = (int) sqrt(npes);

  MPI_Status status;
  MPI_Request send_reqs[32], recv_reqs[32];
  MPI_Status statuses[32];

  vtxs_t nvtxs; 
  ssize_t nnz;
  ssize_t membytes, rmembytes;

  blob_info_t *rblob_info;
  //rblob_info = (blob_info_t *) malloc(sizeof(*rblob_info));

  if (mype != r) {

    unsigned char *blob, *rrblob, *rrblob1;
    blob_info_t *lblob_info;

    lblob_info = mpitc_packBlobBlockCyclic(a_matrix_info);

    blob = lblob_info->blob;
    membytes = lblob_info->membytes;

    MPI_Isend(&membytes, 1, MPI_LONG_LONG_INT, r, 123, MPI_COMM_WORLD, send_reqs);
    MPI_Irecv(&rmembytes, 1, MPI_LONG_LONG_INT, r, 123, MPI_COMM_WORLD, recv_reqs);
    MPI_Waitall(1, recv_reqs, statuses);
    MPI_Waitall(1, send_reqs, statuses);

    rrblob = blockedBlobSendRecv(membytes, rmembytes, r, r, blob, comm);

    #ifdef _ancy
    rrblob1 = (unsigned char *) malloc(rmembytes * sizeof(int32_t));

    /* Isend, Irecv */  
    MPI_Isend(
        (void *) blob, membytes, MPI_INT, r, 
        1, comm, send_reqs);
    MPI_Irecv(
        (void *) rrblob1, rmembytes, MPI_INT, r, 
        1, comm, recv_reqs); 
    MPI_Waitall(1, recv_reqs, statuses);
    MPI_Waitall(1, send_reqs, statuses); 
    #endif

    rblob_info = (blob_info_t *) malloc(sizeof(*rblob_info));
    rblob_info->blob = rrblob;
    rblob_info->membytes = rmembytes;

    free(blob);
    free(lblob_info);

  } else {

    /* Row / Col 0 will remain the same. */
    rblob_info = mpitc_packBlobBlockCyclic(a_matrix_info);
  }

  return rblob_info;
}

blob_info_t *initialShiftBlockCyclic0(
    blob_info_t *ablob, MPI_Comm comm, 
    int my_color, int my_rank) 
{

  int mype, npes, rootp;  
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 
  rootp = (int) sqrt(npes);

  MPI_Status status;
  MPI_Request send_reqs[32], recv_reqs[32];
  MPI_Status statuses[32];

  int s, r;
  vtxs_t nvtxs; 
  ssize_t nnz, membytes, rmembytes;

  blob_info_t *rblob_info = (blob_info_t *) malloc(sizeof(*rblob_info));

  if (my_color > 0) {

    unsigned char *blob, *rrblob;
    blob_info_t *lblob_info;

    /* pkt */
    //lblob_info = mpitc_packBlobBlockCyclic(a_matrix_info);

    //blob = lblob_info->blob;
    //membytes = lblob_info->membytes;

    blob      = ablob->blob;
    membytes  = ablob->membytes;

    /* commt */
    s = mod((my_rank-my_color), rootp);
    r = mod((my_rank+my_color), rootp);
    //printf("mype: %d, s %d, r %d. \n", mype, s, r);

    MPI_Isend(&membytes, 1, MPI_LONG_LONG_INT, s, 1, comm, send_reqs);
    MPI_Irecv(&rmembytes, 1, MPI_LONG_LONG_INT, r, 1, comm, recv_reqs);
    MPI_Waitall(1, recv_reqs, statuses);
    MPI_Waitall(1, send_reqs, statuses);

    rrblob = blockedBlobSendRecv(membytes, rmembytes, s, r, blob, comm);

    #ifdef _ancy
    //rrblob = (unsigned char *) malloc(rmembytes * sizeof(unsigned char));
    rrblob = (unsigned char *) malloc(rmembytes * sizeof(int32_t));

    /* Isend, Irecv */  
    MPI_Isend(
        (void *) blob, membytes, MPI_INT, s, 
        1, comm, send_reqs);
    MPI_Irecv(
        (void *) rrblob, rmembytes, MPI_INT, r, 
        1, comm, recv_reqs); 
    MPI_Waitall(1, recv_reqs, statuses);
    MPI_Waitall(1, send_reqs, statuses);
    #endif

    //free(blob);
    //free(lblob_info);
    //free(ablob->blob);
    //free(ablob);

    //rblob_info = (blob_info_t *) malloc(sizeof(*rblob_info));
    rblob_info->blob = rrblob;
    rblob_info->membytes = rmembytes;
  } else {
    
    /* Row / Col 0 will remain the same. */

    /* pkt */
    //rblob_info = mpitc_packBlobBlockCyclic(a_matrix_info);
    //rblob_info = ablob;
    rblob_info->membytes = ablob->membytes; 
    rblob_info->blob = gk_malloc((ablob->membytes)*sizeof(int32_t), "rblob");
    memcpy(rblob_info->blob, ablob->blob, (ablob->membytes)*sizeof(int32_t));
   
    /* commt */ 
  }

  return rblob_info;
}


blob_info_t *initialShiftBlockCyclic(
    matrix_info_tn *a_matrix_info, MPI_Comm comm, 
    int my_color, int my_rank) 
{

  int mype, npes, rootp;  
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes); 
  rootp = (int) sqrt(npes);

  MPI_Status status;
  MPI_Request send_reqs[32], recv_reqs[32];
  MPI_Status statuses[32];

  int s, r;
  vtxs_t nvtxs; 
  ssize_t nnz, membytes, rmembytes;

  blob_info_t *rblob_info;
  // = (blob_info_t *) malloc(sizeof(*rblob_info));

  if (my_color > 0) {

    unsigned char *blob, *rrblob;
    blob_info_t *lblob_info;

    /* pkt */
    lblob_info = mpitc_packBlobBlockCyclic(a_matrix_info);

    blob = lblob_info->blob;
    membytes = lblob_info->membytes;

    /* commt */
    s = mod((my_rank-my_color), rootp);
    r = mod((my_rank+my_color), rootp);
    //printf("mype: %d, s %d, r %d. \n", mype, s, r);

    MPI_Isend(&membytes, 1, MPI_LONG_LONG_INT, s, 1, comm, send_reqs);
    MPI_Irecv(&rmembytes, 1, MPI_LONG_LONG_INT, r, 1, comm, recv_reqs);
    MPI_Waitall(1, recv_reqs, statuses);
    MPI_Waitall(1, send_reqs, statuses);

    rrblob = blockedBlobSendRecv(membytes, rmembytes, s, r, blob, comm);

    #ifdef _ancy
    //rrblob = (unsigned char *) malloc(rmembytes * sizeof(unsigned char));
    rrblob = (unsigned char *) malloc(rmembytes * sizeof(int32_t));

    /* Isend, Irecv */  
    MPI_Isend(
        (void *) blob, membytes, MPI_INT, s, 
        1, comm, send_reqs);
    MPI_Irecv(
        (void *) rrblob, rmembytes, MPI_INT, r, 
        1, comm, recv_reqs); 
    MPI_Waitall(1, recv_reqs, statuses);
    MPI_Waitall(1, send_reqs, statuses);
    #endif

    free(blob);
    free(lblob_info);

    rblob_info = (blob_info_t *) malloc(sizeof(*rblob_info));
    rblob_info->blob = rrblob;
    rblob_info->membytes = rmembytes;
  } else {
    
    /* Row / Col 0 will remain the same. */

    /* pkt */
    rblob_info = mpitc_packBlobBlockCyclic(a_matrix_info);
   
    /* commt */ 
  }

  return rblob_info;
}

blob_info_t *mpitc_upleftshift(
  MPI_Comm comm, int my_color, 
  int my_rank, blob_info_t *sblob_info) 
{
  int tag1, tag2;
  
  int s, r;
  int npes, mype, rootp;

  blob_info_t *rblob_info;
  ssize_t smembytes, rmembytes;
  unsigned char *sblob, *rblob;
  double start; 
  double ltime, mtime, gtime, atime;

  MPI_Status status, statuses[32];
  MPI_Request recv_requests[1], send_requests[1];

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  rootp = (int) sqrt(npes);

  tag2 = 1, tag1 = 2;

  smembytes = sblob_info->membytes;
  sblob = sblob_info->blob;

  /* Send ra_matrix one left. */
  s = mod(my_rank-1, rootp);
  r = mod(my_rank+1, rootp);
  
  MPI_Send(&smembytes, 1, MPI_LONG_LONG_INT, s, tag1, comm);
  MPI_Recv(&rmembytes, 1, MPI_LONG_LONG_INT, r, tag1, comm, &status);

  rblob = blockedBlobSendRecv(smembytes, rmembytes, s, r, sblob, comm);

  #ifdef _ancy
  MPI_Isend(sblob, smembytes, MPI_INT, s, tag2, comm, send_requests);
  MPI_Irecv(rblob, rmembytes, MPI_INT, r, tag2, comm, recv_requests);
  MPI_Waitall(1, recv_requests, statuses);
  MPI_Waitall(1, send_requests, statuses);
  #endif

  free(sblob);
  free(sblob_info);

  rblob_info = (blob_info_t *) malloc(sizeof(*rblob_info));
  rblob_info->membytes = rmembytes;
  rblob_info->blob = rblob;
  return rblob_info;
}


