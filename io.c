/*!
\file
\brief This file contains routines for reading in the data set
*/

#include "tc.h"

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>

#include <mpi.h>
#include "make_graph.h"

#include "splittable_mrg.h"

#define blk_t int
#define adj_t int64_t
#define MPI_BLK_CNT_INT MPI_INT
#define MPI_ADJ_INT MPI_LONG_LONG_INT
#define ADJNCY_ADM_MAX 1431655766.0
#define XADJ_ADM_MAX 1431655766.0

typedef struct {
  int64_t key;
  int64_t val;
} kvi64_t;

int intcmp1(const void *x, const void *y)
{
  int64_t test; 
  if (((kvi64_t *)x)->key == ((kvi64_t *)y)->key) {
    test = (((kvi64_t *)x)->val - ((kvi64_t *)y)->val);
    if (test < 0)
      return -1; 
    else if (test > 0)
      return 1; 
    
    return 0; 
  }
  else {
    test = (((kvi64_t *)x)->key - ((kvi64_t *)y)->key);
    if (test < 0)
      return -1;
    else if (test > 0)
      return 1;

    return 0;   
  }
}

int _ceil(double x) {
  if((int)x * (double)(1.0) == x)
    return x;
  else
    return (int)x + 1;
}

int64_t _mod(int64_t a, int64_t b) {
  int64_t r = a%b;
  return r<0 ? r+b : r;
}

ssize_t *ancy_all2allv_xadj(ssize_t *lcounts, ssize_t *ladjncy, MPI_Comm comm) {

  double start, ltime, gtime, mtime; 

  int size, rank;
  ssize_t *radjncy; 

  size_t i, j, k;
  ssize_t max_lcnt, c_lcnt; 
  blk_t *blk_cnts, *rblk_cnts;
  ssize_t *blk_displs, *rblk_displs; 
  int *max_slcnts, *max_sldispls, *max_rlcnts, *max_rldispls; 

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (size==1) {
    radjncy = (ssize_t *) malloc((lcounts[0]+1)*sizeof(ssize_t)); 
    memcpy(radjncy, ladjncy, lcounts[0]*sizeof(ssize_t)); 
    return radjncy;  
  }

  for (max_lcnt=-1, i=0; i<size; i++) {
    c_lcnt = _ceil(lcounts[i]/(double) XADJ_ADM_MAX);
    if (max_lcnt < c_lcnt)
      max_lcnt = c_lcnt;
  }

  /* This is for size procs times number of 
     alltoallv we will have to make. */
  blk_cnts = (blk_t *) malloc(size*max_lcnt*sizeof(blk_t));
  blk_displs = (ssize_t *) malloc((size*max_lcnt+1)*sizeof(ssize_t));
  memset(blk_cnts, 0, sizeof(blk_t)*size*max_lcnt);

  ssize_t accum;
  for (i=0; i<size; i++) {
    for (j=0; j<max_lcnt; j++) {
      if ((lcounts[i]-XADJ_ADM_MAX) > 0) {
        blk_cnts[i*max_lcnt+j] = XADJ_ADM_MAX; 
        lcounts[i] -= XADJ_ADM_MAX; 
        accum += XADJ_ADM_MAX; 
      } else { /* Should ideally be the last blk with entries. */
        blk_cnts[i*max_lcnt+j] = lcounts[i];
        lcounts[i] = 0;
      }
    }
    accum = 0;
  }

  blk_displs[0]=0;
  for (i=1; i<=size*max_lcnt; i++) {
    blk_displs[i] = (ssize_t) blk_cnts[i-1]+blk_displs[i-1];
  }

  /* I will need local buffers for: 
     1) counts 
     2) displs 
     3) rcounts 
     4) rdispls. */

  max_slcnts = (int *) malloc(size*sizeof(int));
  max_sldispls = (int *) malloc((size+1)*sizeof(int));
  for (i=0; i<size; i++)
    max_slcnts[i] = (int) max_lcnt; 

  max_rlcnts = (int *) malloc(size*sizeof(int));
  max_rldispls = (int *) malloc((size+1)*sizeof(int));

  /* First, do an all2allv of max_cnts so as to init an array of that size. */
  MPI_Alltoall(max_slcnts, 1, MPI_INT, 
      max_rlcnts, 1, MPI_INT, comm);

  max_rldispls[0] = 0;
  max_sldispls[0] = 0;
  for (i=1; i<=size; i++) {
    max_rldispls[i] = max_rlcnts[i-1] + max_rldispls[i-1];
    max_sldispls[i] = max_slcnts[i-1] + max_sldispls[i-1];
  }

  rblk_cnts = (blk_t *) malloc(max_rldispls[size]*sizeof(blk_t));

  MPI_Alltoallv(
      blk_cnts, max_slcnts, max_sldispls, 
      MPI_BLK_CNT_INT, rblk_cnts, max_rlcnts, 
      max_rldispls, MPI_BLK_CNT_INT, comm);

  rblk_displs = (ssize_t *) malloc((max_rldispls[size]+1)*sizeof(ssize_t)); 
  rblk_displs[0] = 0;
  for (i=1; i<=max_rldispls[size]; i++) {
    rblk_displs[i] = (ssize_t) rblk_cnts[i-1]+rblk_displs[i-1];
  }

  radjncy = (ssize_t *) malloc((rblk_displs[max_rldispls[size]]+1)*sizeof(ssize_t));

  MPI_Request sreq[1024], rreq[1024];
  MPI_Status statuses[1024];

  int dst; 
  for (i=0; i<size; i++) {
    dst = _mod((i+rank), size);
    for (k=0, j=max_sldispls[dst]; j<max_sldispls[dst+1]; j++, k++) {
      MPI_Isend((void *) (ladjncy+blk_displs[j]), blk_cnts[j],
          MPI_LONG_LONG_INT, dst, (rank+k), comm, sreq+j);
    }
  }

  for (i=0; i<size; i++) {
    dst = _mod((i+rank), size);
    for (k=0, j=max_rldispls[dst]; j<max_rldispls[dst+1]; j++, k++)  {
      MPI_Irecv((void *) (radjncy+rblk_displs[j]), rblk_cnts[j], 
          MPI_LONG_LONG_INT, dst, (dst+k), comm, rreq+j);
    }
  }

  MPI_Waitall(max_sldispls[size], sreq, statuses);
  MPI_Waitall(max_rldispls[size], rreq, statuses);

  free(blk_cnts);
  free(blk_displs);
  free(rblk_cnts);
  free(rblk_displs);
  free(max_slcnts);
  free(max_rlcnts);
  free(max_sldispls);
  free(max_rldispls);
 
  return radjncy; 
}


int64_t *ancy_all2allv(ssize_t *lcounts, int64_t *ladjncy, MPI_Comm comm) {

  double start, ltime, gtime, mtime; 

  size_t i, j, k;
  ssize_t max_lcnt, c_lcnt; 
  int size, rank;
  blk_t *blk_cnts, *rblk_cnts;
  ssize_t *blk_displs, *rblk_displs; 
  int *max_slcnts, *max_sldispls, *max_rlcnts, *max_rldispls; 
  int64_t *radjncy; 

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (size==1) {
    radjncy = (int64_t *) malloc(lcounts[0]*sizeof(int64_t));
    memcpy(radjncy, ladjncy, lcounts[0]*sizeof(int64_t));
    return radjncy; 
  }

  for (max_lcnt=-1, i=0; i<size; i++) {
    c_lcnt = _ceil(lcounts[i]/(double) ADJNCY_ADM_MAX);
    if (max_lcnt < c_lcnt)
      max_lcnt = c_lcnt;
  }

  /* This is for size procs times number of 
     alltoallv we will have to make. */
  blk_cnts = (blk_t *) malloc(size*max_lcnt*sizeof(blk_t));
  blk_displs = (ssize_t *) malloc((size*max_lcnt+1)*sizeof(ssize_t));
  memset(blk_cnts, 0, sizeof(blk_t)*size*max_lcnt);

  ssize_t accum;
  for (i=0; i<size; i++) {
    for (j=0; j<max_lcnt; j++) {
      if ((lcounts[i]-ADJNCY_ADM_MAX) > 0) {
        blk_cnts[i*max_lcnt+j] = ADJNCY_ADM_MAX; 
        lcounts[i] -= ADJNCY_ADM_MAX; 
        accum += ADJNCY_ADM_MAX; 
      } else { /* Should ideally be the last blk with entries. */
        blk_cnts[i*max_lcnt+j] = lcounts[i];
        lcounts[i] = 0;
      }
    }
    accum = 0;
  }

  blk_displs[0]=0;
  for (i=1; i<=size*max_lcnt; i++) {
    blk_displs[i] = (ssize_t) blk_cnts[i-1]+blk_displs[i-1];
  }

  max_slcnts = (int *) malloc(size*sizeof(int));
  max_sldispls = (int *) malloc((size+1)*sizeof(int));
  for (i=0; i<size; i++)
    max_slcnts[i] = (int) max_lcnt; 

  max_rlcnts = (int *) malloc(size*sizeof(int));
  max_rldispls = (int *) malloc((size+1)*sizeof(int));

  /* First, do an all2allv of max_cnts so as to init an array of that size. */
  MPI_Alltoall(max_slcnts, 1, MPI_INT, 
      max_rlcnts, 1, MPI_INT, comm);

  max_rldispls[0] = 0;
  max_sldispls[0] = 0;
  for (i=1; i<=size; i++) {
    max_rldispls[i] = max_rlcnts[i-1] + max_rldispls[i-1];
    max_sldispls[i] = max_slcnts[i-1] + max_sldispls[i-1];
  }

  rblk_cnts = (blk_t *) malloc(max_rldispls[size]*sizeof(blk_t));

  MPI_Alltoallv(
      blk_cnts, max_slcnts, max_sldispls, 
      MPI_BLK_CNT_INT, rblk_cnts, max_rlcnts, 
      max_rldispls, MPI_BLK_CNT_INT, comm);

  rblk_displs = (ssize_t *) malloc((max_rldispls[size]+1)*sizeof(ssize_t)); 
  rblk_displs[0] = 0;
  for (i=1; i<=max_rldispls[size]; i++) {
    rblk_displs[i] = (ssize_t) rblk_cnts[i-1]+rblk_displs[i-1];
  }

  radjncy = gk_malloc(rblk_displs[max_rldispls[size]]*sizeof(int64_t), "radjncy within ancy_all2allv");

  MPI_Request sreq[1024], rreq[1024];
  MPI_Status statuses[1024];

  int dst;
  for (i=0; i<size; i++) {
    dst = _mod((i+rank), size);
    for (k=0, j=max_sldispls[dst]; j<max_sldispls[dst+1]; j++, k++) {
      MPI_Isend((void *) (ladjncy+blk_displs[j]), blk_cnts[j],
          MPI_LONG_LONG_INT, dst, (rank+k), comm, sreq+j);
    }
  }

  for (i=0; i<size; i++) {
    dst = _mod((i+rank), size);
    for (k=0, j=max_rldispls[dst]; j<max_rldispls[dst+1]; j++, k++)  {
      MPI_Irecv((void *) (radjncy+rblk_displs[j]), rblk_cnts[j], 
          MPI_LONG_LONG_INT, dst, (dst+k), comm, rreq+j);
    }
  }

  MPI_Waitall(max_sldispls[size], sreq, statuses);
  MPI_Waitall(max_rldispls[size], rreq, statuses);

  free(rblk_cnts);
  free(rblk_displs);
  free(blk_cnts);
  free(blk_displs);
  free(max_sldispls);
  free(max_rldispls);
  free(max_slcnts);
  free(max_rlcnts);

  return radjncy; 
}

/*************************************************************************
* This function generates an RMAT graph using the 
* grapp libraries. 
*
* For data structure tracking - 
*   BO - big ones. 
**************************************************************************/
void ParallelReadGraphRMAT(gk_graph_t *graph, int scale, MPI_Comm comm)
{
  int log_numverts;
  int size, rank;

  int *vcounts, *vdispls, *rvcounts, *rvdispls; 

  ssize_t my_edges;
  unsigned long global_edges;
  double start, stop;

  size_t i, j, k, trker;
  ssize_t *lcounts, *lvcounts; 

  int64_t nedges, t_key;
  int64_t key, val, p_key, p_val;
  ssize_t nvtxs, chunk, p_edges, adj_entries; 

  ssize_t *lxadj, *rxadj, *rxadj1, *vtxdist;
  int64_t *ladjncy, *radjncy;

  log_numverts = scale; /* In base 2 */

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) fprintf(stderr, "Graph size is %" PRId64 " vertices and %" PRId64 " edges\n", INT64_C(1) << log_numverts, INT64_C(16) << log_numverts);

  /* Start of graph generation timing */
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  packed_edge* result;

  /* BO - result. */
  make_graph(log_numverts, INT64_C(16) << log_numverts, 1, 2, &nedges, &result);
  /* End of graph generation timing */

  nvtxs = (INT64_C(1) << log_numverts);
  chunk = nvtxs/(ssize_t) size;
  my_edges = nedges;

  vtxdist = gk_malloc((size+1)*sizeof(ssize_t), "vtxdist allocation. \n");
  
  if (vtxdist == NULL) 
   printf("rank: %d vtxdist malloc didn't go through. \n", rank);
  
  for (i=0; i<size; i++) {
    vtxdist[i] = chunk*i;
  }
  vtxdist[i] = nvtxs;

  /* Creating an undirected version of the 
    edges here, j will tell me the total number 
    of edges. */

  lxadj     = (ssize_t *) malloc((nvtxs+1)*sizeof(ssize_t));
  memset(lxadj,   0, (nvtxs+1)*sizeof(ssize_t));

  p_key = -1; 
  p_val = -1;
  
  for (j=0, i=0; i<my_edges; i++) {
    key = (int64_t) get_v0_from_edge(&result[i]);
    val = (int64_t) get_v1_from_edge(&result[i]);

    if (key==val)
      continue;

    if (key==p_key && val==p_val)
      continue;

    lxadj[key]++;
    lxadj[val]++;
    j+=2;

    p_key = key;
    p_val = val;
  }

  MAKECSR(i, nvtxs, lxadj);

  p_key = -1;
  p_val = -1;

  ladjncy   = (int64_t *) malloc(j*sizeof(int64_t));
  for (i=0; i<my_edges; i++) {
    key = (int64_t) get_v0_from_edge(&result[i]);
    val = (int64_t) get_v1_from_edge(&result[i]);

    if (key==val)
      continue;

    if (key==p_key && val==p_val)
      continue;

    ladjncy[lxadj[key]++]=val;
    ladjncy[lxadj[val]++]=key;

    p_key = key;
    p_val = val;
  }
  free(result);

  SHIFTCSR(i, nvtxs, lxadj);

  for (i=0; i<nvtxs; i++)
    gk_i64sorti(lxadj[i+1]-lxadj[i], ladjncy+lxadj[i]);
  
  ssize_t *refined_lxadj; 
  refined_lxadj = gk_malloc((nvtxs+1)*sizeof(ssize_t), "refined lxadj");
  memset(refined_lxadj, 0, (nvtxs+1)*sizeof(ssize_t));

  /* Deduplicates, prepping for all2allv. */
  lcounts   = (ssize_t *) malloc(size*sizeof(ssize_t));
  memset(lcounts, 0, size*sizeof(ssize_t));

  p_key = -1; 
  p_val = -1; 

  for (p_edges=0, i=0; i<nvtxs; i++) {
    for (j=lxadj[i]; j<lxadj[i+1]; j++) {

      if (i==ladjncy[j])
        continue; 

      if (i==p_key && ladjncy[j]==p_val)
        continue; 

      t_key = i/chunk;
      if (t_key > size-1)
        t_key = size-1;

      lcounts[t_key]++;
      refined_lxadj[i]++;
      ladjncy[p_edges++] = ladjncy[j];
    
      p_key = i;
      p_val = ladjncy[j];
    }
  }
  free(lxadj);

  lvcounts = (ssize_t *) malloc(size*sizeof(ssize_t));
  for (i=0; i<size; i++) {
    lvcounts[i] = (ssize_t) vtxdist[i+1]-vtxdist[i];
  }

  rxadj = ancy_all2allv_xadj(lvcounts, refined_lxadj, MPI_COMM_WORLD);
  free(refined_lxadj);

  ssize_t exp_rvtxs; 
  exp_rvtxs = size*(vtxdist[rank+1]-vtxdist[rank]);  

  MAKECSR(i, exp_rvtxs, rxadj);
  adj_entries = rxadj[exp_rvtxs];

  /* Handling ladjncy now with new all2allv. */
  radjncy = ancy_all2allv(lcounts, ladjncy, MPI_COMM_WORLD); 
  free(ladjncy); 

  free(lcounts);
  free(lvcounts);

  ssize_t mychunk, idx;
  mychunk = vtxdist[rank+1]-vtxdist[rank];
  ssize_t *wqxadj;
  int64_t *wqadjncy; 
  
  wqxadj    = gk_malloc((mychunk+1)*sizeof(ssize_t), "without qsort");
  memset(wqxadj, 0, (mychunk+1)*sizeof(ssize_t));
  for (i=0; i<size; i++) {
    for (j=0; j<mychunk; j++) {
      idx = i*mychunk+j;
      wqxadj[j] += rxadj[idx+1]-rxadj[idx];
    }
  }
  MAKECSR(i, mychunk, wqxadj);
  
  wqadjncy  = gk_malloc(adj_entries*sizeof(int64_t), "adj entries");
  for (i=0; i<exp_rvtxs; i++) {
    memcpy(
      wqadjncy+wqxadj[i%mychunk], 
      radjncy+rxadj[i], 
      (rxadj[i+1]-rxadj[i])*sizeof(int64_t));
    wqxadj[i%mychunk] += rxadj[i+1]-rxadj[i]; 
  }
  free(rxadj);
  free(radjncy);

  SHIFTCSR(i, mychunk, wqxadj);

  for (i=0; i<mychunk; i++) {
    gk_i64sorti(wqxadj[i+1]-wqxadj[i], wqadjncy+wqxadj[i]);
  }

  ssize_t *lxadj1;
  lxadj1 = gk_malloc((mychunk+1)*sizeof(ssize_t), "lxadj1");
  memset(lxadj1, 0, (mychunk+1)*sizeof(ssize_t));

  p_key=-1;
  p_val=-1;

  for (p_edges=0, i=0; i<mychunk; i++) {
    for (j=wqxadj[i]; j<wqxadj[i+1]; j++) {

      key = (int64_t) i;
      val = (int64_t) wqadjncy[j];

      if (key+vtxdist[rank]==val)
        continue;

      if (key==p_key && val==p_val)
        continue;

      lxadj1[i]++;
      wqadjncy[p_edges] = val;
      p_edges++;

      p_key = key;
      p_val = val;
    }
  }
  free(wqxadj);

  MAKECSR(i, mychunk, lxadj1);

  ssize_t g_edges, l_edges, m_edges;
  g_edges = 0;
  m_edges = 0;
  l_edges = (ssize_t) lxadj1[mychunk];
  
  MPI_Reduce(&l_edges, &g_edges, 1, 
        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&l_edges, &m_edges, 1, 
        MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
  MPI_Barrier(MPI_COMM_WORLD);
  stop = MPI_Wtime();
  
  if (rank == 0) {
    printf("myrank: %d #(pruned edges): %zd. \n", rank, g_edges);
    printf("myrank: %d #(max #edge): %zd, #(avg edges): %zd. \n", rank, m_edges, (g_edges/size));
    printf("%lu edge%s generated in %fs (%f Medges/s on %d processor%s)\n", g_edges, (g_edges == 1 ? "" : "s"), (stop - start), g_edges / (stop - start) * 1.e-6, size, (size == 1 ? "" : "s"));
  }

  graph->xadj     = lxadj1;
  graph->nvtxs    = mychunk;
  graph->vtxdist  = vtxdist; 
  graph->adjncy   = wqadjncy; 
}

/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
void ParallelReadGraphMETIS(gk_graph_t *graph, char *filename, MPI_Comm comm)
{
  ssize_t i, k, l, pe;
  int npes, mype, ier;
  int32_t gnvtxs, nvtxs, your_nvtxs, maxnvtxs = -1, edge;
  ssize_t your_nedges, gnedges, maxnedges = -1;
  ssize_t *vtxdist, *xadj, *your_xadj;
  //int32_t *adjncy, *your_adjncy;
  adj_t *adjncy, *your_adjncy;
  MPI_Status stat;
  char *line = NULL, *oldstr, *newstr;
  FILE *fpin = NULL;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = (ssize_t *) malloc((npes+1)*sizeof(ssize_t));
  memset(vtxdist, 0, (npes+1)*sizeof(ssize_t)); 

  if (mype == npes-1) {
    ier = 0;
    fpin = fopen(filename, "r");

    if (fpin == NULL) {
      printf("COULD NOT OPEN FILE '%s' FOR SOME REASON!\n", filename);
      ier++;
    }

    MPI_Bcast(&ier, 1, MPI_INT, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    line = gk_cmalloc(MAXLINE+1, "line");

    while (fgets(line, MAXLINE, fpin) && line[0] == '%');

    sscanf(line, "%"PRId32" %zd", &gnvtxs, &gnedges);
    gnedges *=2;

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=gnvtxs; i<npes; i++) {
      l = k/(npes-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }

    MPI_Bcast((void *)vtxdist, npes+1, MPI_LONG_LONG_INT, npes-1, comm);
  }
  else {
    MPI_Bcast(&ier, 1, MPI_INT, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    MPI_Bcast((void *)vtxdist, npes+1, MPI_LONG_LONG_INT, npes-1, comm);
  }

  graph->vtxdist = vtxdist;

  nvtxs = graph->nvtxs = (int32_t) vtxdist[mype+1]-vtxdist[mype];
  xadj  = graph->xadj  = (ssize_t *) malloc((graph->nvtxs+1)*sizeof(ssize_t));

  /*******************************************/
  /* Go through first time and generate xadj */
  /*******************************************/
  if (mype == npes-1) {
    maxnvtxs = (int32_t) vtxdist[1];
    for (i=1; i<npes; i++) 
      maxnvtxs = (maxnvtxs < vtxdist[i+1]-vtxdist[i] ? vtxdist[i+1]-vtxdist[i] : maxnvtxs);

    your_xadj = (ssize_t *) malloc((maxnvtxs+1)*sizeof(ssize_t));

    maxnedges = 0;
    for (pe=0; pe<npes; pe++) {
      your_nvtxs = (int32_t) vtxdist[pe+1]-vtxdist[pe];

      for (i=0; i<your_nvtxs; i++) {
        your_nedges = 0;

        while (fgets(line, MAXLINE, fpin) && line[0] == '%'); /* skip lines with '#' */
        oldstr = line;
        newstr = NULL;

        for (;;) {
          edge  = strtoll(oldstr, &newstr, 10) -1;
          oldstr = newstr;

          if (edge < 0)
            break;

          your_nedges++;
        }
        your_xadj[i] = your_nedges;
      }

      MAKECSR(i, your_nvtxs, your_xadj);
      maxnedges = (maxnedges < your_xadj[your_nvtxs] ? your_xadj[your_nvtxs] : maxnedges);

      if (pe < npes-1) {
        MPI_Send((void *)your_xadj, your_nvtxs+1, MPI_LONG_LONG_INT, pe, 0, comm);
      }
      else {
        for (i=0; i<your_nvtxs+1; i++)
          xadj[i] = your_xadj[i];
      }
    }
    fclose(fpin);
    gk_free((void **)&your_xadj, LTERM);
  }
  else {
    MPI_Recv((void *)xadj, nvtxs+1, MPI_LONG_LONG_INT, npes-1, 0, comm, &stat);
  }

  adjncy = graph->adjncy = (adj_t *) malloc(xadj[nvtxs]*sizeof(adj_t));

  /***********************************************/
  /* Now go through again and record adjncy data */
  /***********************************************/
  if (mype == npes-1) {
    ier = 0;
    fpin = fopen(filename, "r");

    if (fpin == NULL){
      printf("COULD NOT OPEN FILE '%s' FOR SOME REASON!\n", filename);
      ier++;
    }

    MPI_Bcast(&ier, 1, MPI_INT, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    /* get first line again */
    while (fgets(line, MAXLINE, fpin) && line[0] == '%');

    your_adjncy = (adj_t *) malloc(maxnedges*sizeof(adj_t));

    for (pe=0; pe<npes; pe++) {
      your_nvtxs  = vtxdist[pe+1]-vtxdist[pe];
      your_nedges = 0;

      for (i=0; i<your_nvtxs; i++) {
        while (fgets(line, MAXLINE, fpin) && line[0] == '%');
        oldstr = line;
        newstr = NULL;

        for (;;) {
          edge   = strtoll(oldstr, &newstr, 10) -1;
          oldstr = newstr;

          if (edge < 0)
            break;

          your_adjncy[your_nedges] = edge;
          your_nedges++;
        }
      }
      if (pe < npes-1) {
        MPI_Send((void *)your_adjncy, your_nedges, MPI_ADJ_INT, pe, 0, comm);
      }
      else {
        for (i=0; i<your_nedges; i++)
          adjncy[i] = your_adjncy[i];
      }
    }
    fclose(fpin);
    gk_free((void **)&your_adjncy, &line, LTERM);
  }
  else {
    MPI_Bcast(&ier, 1, MPI_INT, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    MPI_Recv((void *)adjncy, xadj[nvtxs], MPI_ADJ_INT, npes-1, 0, comm, &stat);
  }
}



/**************************************************************************/
/*! Reads the input data */
/**************************************************************************/
vault_t *loadData(params_t *params)
{

  int mype, npes; 
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  vault_t *vault;
  vault = (vault_t *)gk_malloc(sizeof(vault_t), "loadData: vault");
  memset(vault, 0, sizeof(vault_t));

  gk_graph_t *graph = malloc(sizeof(*graph));

  switch (params->iftype) {
    case IFTYPE_METIS:
      /* read the graph */
      printf("Reading graph %s...\n", params->infile);
      GKASSERT(gk_fexists(params->infile)); 
      ParallelReadGraphMETIS(graph, params->infile, MPI_COMM_WORLD);
      break;

    case IFTYPE_RMAT:
      if(params->infile)
        printf("Ignoring given input file, RMAT graph with scale %d being generated..\n", params->scale);
      
      ParallelReadGraphRMAT(graph, params->scale, MPI_COMM_WORLD);
      break;

    default:
      errexit("Unknown iftype of %d\n", params->iftype);
  }

  vault->graph = graph;

  return vault;
}

