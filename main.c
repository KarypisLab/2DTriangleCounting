/*!
\file
\brief The entry point of the triangle counting code
*/

#include "tc.h"

/*************************************************************************
* The entry point 
**************************************************************************/
int main(int argc, char *argv[])
{
  int64_t ntriangles=0;
  params_t *params;
  vault_t *vault;

  params = getcmdline_params(argc, argv);

  /* Initialize MPI */
  MPI_Init(&params->argc, &params->argv);

  int mype, npes;
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  vault = loadData(params);

  setvbuf(stdout, NULL, _IOLBF, 0);

  if (mype == 0) {
    printf("\n-----------------\n");
    printf("infile: %s\n", params->infile);
    printf("#nvtxs: %d\n", vault->graph->nvtxs);
    printf("tctype: %s, otype: %s\n", tctypenames[params->tctype], otypenames[params->otype]);
    printf("\n");
  }

  /* perform various non-timed initializations */
  srand(params->seed);

  /* count triangles... */
  switch (params->tctype) {

    case TCTYPE_MAPJIK2D: 
      MPI_Barrier(MPI_COMM_WORLD);
      if (mype==0)
        gk_startwctimer(vault->timer_global);
      ntriangles = mpitc_MapJIK2D(params, vault);
      MPI_Barrier(MPI_COMM_WORLD);
      if (mype==0)
        gk_stopwctimer(vault->timer_global);
      break;

    default:
      errexit("Unknown tctype of %d\n", params->tctype);
  }


  if (mype == 0) {

    printf("#triangles: %"PRId64"; rate: %.4lf MT/sec\n", ntriangles, 
        ((double)ntriangles)/((double)1e6*gk_getwctimer(vault->timer_2))); 
    printf("wall clock time:\t%.3lfs\n", gk_getwctimer(vault->timer_global));
    if (gk_getwctimer(vault->timer_1) > 0)
      printf("pre-processing:\t\t%.3lfs\n", gk_getwctimer(vault->timer_1));
    if (gk_getwctimer(vault->timer_2) > 0)
      printf("triangle counting:\t%.3lfs\n", gk_getwctimer(vault->timer_2));
    if (gk_getwctimer(vault->timer_3) > 0)
      printf("   timer_3: %.3lfs\n", gk_getwctimer(vault->timer_3));
    if (gk_getwctimer(vault->timer_4) > 0)
      printf("   timer_4: %.3lfs\n", gk_getwctimer(vault->timer_4));
    if (gk_getwctimer(vault->timer_5) > 0)
      printf("   timer_5: %.3lfs\n", gk_getwctimer(vault->timer_5));
    if (gk_getwctimer(vault->timer_6) > 0)
      printf("   timer_6: %.3lfs\n", gk_getwctimer(vault->timer_6));
    printf("\n-----------------\n");

    printf("RUNRESULT: %s %s %d %"PRId64" %.4lf %.3lf %.3lf %.3lf\n",
        params->infile, 
        tctypenames[params->tctype], 
        params->scale, 
        ntriangles, 
        ((double)ntriangles)/((double)1e6*gk_getwctimer(vault->timer_2)),
        gk_getwctimer(vault->timer_global), 
        gk_getwctimer(vault->timer_1), 
        gk_getwctimer(vault->timer_2)
        );
  }

  MPI_Finalize();
}

