/*!
\file
\brief This file contains function prototypes
*/

#ifndef _PROTO_H_
#define _PROTO_H_

#include <mpi.h>

/* io.c */
vault_t *loadData(params_t *params);
void ParallelReadGraphRMAT(gk_graph_t *graph, int scale, MPI_Comm comm);
void ParallelReadGraphMETIS(gk_graph_t *graph, char *filename, MPI_Comm comm);

/* cmdline.c */
params_t *getcmdline_params(int argc, char *argv[]);

/* mpi_tc.c */
int32_t *gk_i32kvipermi(int32_t n, gk_i32kv_t *cand);
ssize_t *gk_sskvipermi(int32_t n, ssize_t *xadj);
matrix_info_tn *mpitc_prep2DGraph(gk_graph_t *graph);
matrix_info_tn *mpitc_prep2DGraphTranspose(gk_graph_t *graph);
int64_t mpitc_MapJIK2D(params_t *params, vault_t *vault);
int64_t mpitc_kernel(int shift, matrix_info_tn *la_matrix, blob_info_t *rblob, blob_info_t *cblob);
blob_info_t *initialShiftBlockCyclic(matrix_info_tn *a_matrix_info, MPI_Comm comm, int my_color, int my_rank);
blob_info_t *mpitc_upleftshift(MPI_Comm comm, int my_color, int my_rank, blob_info_t *sblob_info);
int64_t mpitc_kernelJIK(int shift, matrix_info_tn *la_matrix, blob_info_t *rblob, blob_info_t *cblob);
int64_t mpitc_blobBytesBlockCyclic(matrix_info_tn *a_matrix_info);
blob_info_t *mpitc_packBlobBlockCyclic(matrix_info_tn *a_matrix_info);
matrix_info_tn *mpitc_unpackBlobBlockCyclic(blob_info_t *blob_info);
blob_info_t *transposeSubMatrix(matrix_info_tn *a_matrix_info, MPI_Comm comm, int r);
ssize_t *ancy_all2allv_xadj(ssize_t *lcounts, ssize_t *ladjncy, MPI_Comm comm);
int64_t *ancy_all2allv(ssize_t *lcounts, int64_t *ladjncy, MPI_Comm comm);
#endif
