/*!
\file
\brief Data structures used in the program
*/

#ifndef _STRUCT_TC_H_
#define _STRUCT_TC_H_

#define adj_t int64_t
#define vtxs_t int32_t 
#define chunk_t ssize_t
#define MPI_CHUNK_INT MPI_LONG_LONG_INT

/*************************************************************************
* stores the current xadj[] and tail location information that is 
* initially used by mapjikv2
**************************************************************************/
typedef struct {
  ssize_t iend;      /* end index of the adjacency list */
  uint32_t len;      /* the length of the "active" adjancency list */
  uint32_t tlen;     /* the length of the "tail" part of the list */
} ste_indices_t;


/*************************************************************************
* the data vault 
**************************************************************************/
typedef struct {
  gk_graph_t *graph;      /* the graph */

  /* timers */
  double timer_global;
  double timer_1;
  double timer_2;
  double timer_3;
  double timer_4;
  double timer_5;
  double timer_6;
} vault_t; 


/*************************************************************************
* run parameters 
**************************************************************************/
typedef struct {
  int tctype;           /* The algorithm to use for triangle counting */
  int otype;            /* The type of graph ordering */
  int iftype;           /* The format of the input file */

  int seed;             /* Seed for the random number generator */
  int scale;   					/* Scale for generating the RMAT graphs */
  int dbglvl;           /* Debuging information */

  char *infile;         /* The file storing the input data */

  int argc;             /* For mpi runs. */ 
  char **argv;
} params_t;

typedef struct {
  int n_nbrs; /* Number of nbrs I will recv data from. */
  int *nbr_ids; /* Id of the nbr I will recv data from. */
  ssize_t *ptr; /* csr for each nbr. */
  adj_t *ind; /* actual vertex ids it needs to recv values from each processor during transfer. */
} SendInfoTypeGK;

typedef struct {
  int n_nbrs; /* Number of nbrs I will recv data from. */
  int *nbr_ids; /* Id of the nbr I will recv data from. */
  ssize_t *ptr; /* csr for each nbr. */
  adj_t *ind; /* actual vertex ids it needs to recv values from each processor during transfer. */
} RecvInfoTypeGK;

typedef struct {
  vtxs_t nvtxs;
  adj_t *vtx_ids;
  ssize_t *vtx_vals;
} AccumPRDetails;


typedef struct {
  vtxs_t nvtxs;
  int32_t *vtx_ids;
  ssize_t *vtx_adjncy_counts;
} RemoteNbrDetails;


/*************************************************************************
* sub matrix info 
**************************************************************************/
typedef struct matrix_info {
  vtxs_t nvtxs;

  int local_i; 
  int local_j;
  int global_i; 
  int global_j;

  chunk_t vi_offset; 
  chunk_t vj_offset;
  chunk_t vi_chunk;
  chunk_t vj_chunk;
  ssize_t *xadj;
  adj_t *adjncy;
} matrix_info_t;



/*************************************************************************
* sub matrix info 
**************************************************************************/
typedef struct smatrix_info_n {
  vtxs_t nvtxs;

  int32_t rc_blks;
  chunk_t *vi_chunk;
  chunk_t *vj_chunk;
  ssize_t *xadj;
  adj_t *adjncy;
} smatrix_info_tn;



/*************************************************************************
* sub matrix info 
**************************************************************************/
typedef struct matrix_info_n {
  vtxs_t nvtxs;
  vtxs_t g_nvtxs; 

  chunk_t *vi_offset; 
  chunk_t *vj_offset;
  chunk_t *vi_chunk;
  chunk_t *vj_chunk;

  ssize_t *xadj;
  adj_t *adjncy;
  ssize_t *reverse_idx;
  adj_t *nnz_idx;
} matrix_info_tn;


/*************************************************************************
* sub matrix info 
**************************************************************************/
typedef struct matrix_info_a {
  vtxs_t nvtxs;

  chunk_t *uvtxs_displs;

  int32_t unvtxs;
  ssize_t *uxadj;
  ssize_t *ureverse_idx;
  adj_t *uadjncy;
  adj_t *unnz_idx; 

  int32_t lnvtxs;
  ssize_t *lxadj;
  ssize_t *lreverse_idx;
  adj_t *ladjncy;
  adj_t *lnnz_idx;

} matrix_info_ta;



/*************************************************************************
* proc info 
**************************************************************************/
typedef struct proc_info {
  int i; 
  int j; 
  int m_i; 
  int m_j; 
  matrix_info_t *la_matrix; 
  matrix_info_t *lb_matrix; 
  matrix_info_t *ra_matrix; 
  matrix_info_t *rb_matrix;
} proc_info_t; 



/*************************************************************************
* proc info 
**************************************************************************/
typedef struct proc_info_n {
  int i; 
  int j; 
  int m_i; 
  int m_j; 
  matrix_info_tn *la_matrix; 
  matrix_info_tn *lb_matrix; 
  matrix_info_tn *ra_matrix; 
  matrix_info_tn *rb_matrix;
} proc_info_tn; 

/**
* blob + #bytes
*/
typedef struct blob_info {
  ssize_t membytes;
  unsigned char* blob;
} blob_info_t;

#endif 
