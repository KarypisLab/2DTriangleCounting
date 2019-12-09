/*!
\file  
\brief Parsing of command-line arguments
*/

#include "tc.h"

/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
  {"tctype",            1,      0,      CMD_TCTYPE},
  {"otype",             1,      0,      CMD_OTYPE},
  {"iftype",            1,      0,      CMD_IFTYPE},

  {"scale",             1,      0,      CMD_SCALE},
  {"dbglvl",            1,      0,      CMD_DBGLVL},
  {"help",              0,      0,      CMD_HELP},
  {0,                   0,      0,      0}
};

/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
static gk_StringMap_t tctype_options[] = {
  {"mapjik2d",     TCTYPE_MAPJIK2D},
  {NULL,                 0}
};

/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
static gk_StringMap_t otype_options[] = {
  {"none",      OTYPE_NONE},
  {"incd",      OTYPE_INCD},
  {"decd",      OTYPE_DECD},
  {NULL,        0}
};

static gk_StringMap_t iftype_options[] = {
  {"rmat",         IFTYPE_RMAT},
  {"metis",       IFTYPE_METIS},
  {NULL,                 0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: mpirun -np <#ranks> ./build/Linux-x86_64/mpitc [options] infile",
" ",
" ",        
" Options",
"  -tctype=text",
"     Specifies the type of triangle counting algorithm to use.",
"     Possible values are:",
"        mapjik2d      hash-map jik enumeration 2D [default]",
" ",
"  -iftype=text",
"     Specifies the format of the input file. ",
"     Possible values are:",
"        rmat    Generate an RMAT graph for a given scale",
"        metis   METIS format [default]",
" ",
"  -scale=int",
"     Specifies the scale for generating the RMAT graphs.",
"     Default value is 10.",
" ",
"  -dbglvl=int",
"     Specifies the level of debugging information to be displayed.",
"     Default value is 0.",
" ",
"  -help",
"     Prints this message.",
""
};

 

/*************************************************************************/
/*! This is the entry point of the command-line argument parser */
/*************************************************************************/
params_t *getcmdline_params(int argc, char *argv[])
{
  gk_idx_t i, j, k;
  int type=0;
  int c, option_index;
  params_t *params;

  params = (params_t *)gk_malloc(sizeof(params_t), "cmdline_parse: params");
  memset(params, 0, sizeof(params_t)); 

  /* print the command line */
  for (i=0; i<argc; i++)
    printf("%s ", argv[i]);
  printf("\n");

  /* initialize the params data structure */
  params->infile  = NULL;

  params->tctype   = TCTYPE_MAPJIK2D;
  params->otype    = OTYPE_INCD;
  params->iftype   = IFTYPE_METIS;
  params->dbglvl   = 0;
	params->scale 	 = 10;
  params->argc = argc;
  params->argv = argv;


  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case CMD_TCTYPE:
        if (gk_optarg) {
          if ((params->tctype = gk_GetStringID(tctype_options, gk_optarg)) == -1) 
            errexit("Invalid tctype of %s.\n", gk_optarg);
        }
        break; 

      case CMD_OTYPE:
        if (gk_optarg) {
          if ((params->otype = gk_GetStringID(otype_options, gk_optarg)) == -1) 
            errexit("Invalid otype of %s.\n", gk_optarg);
        }
        break; 

      case CMD_IFTYPE:
        if (gk_optarg) {
          if ((params->iftype = gk_GetStringID(iftype_options, gk_optarg)) == -1) 
            errexit("Invalid iftype of %s.\n", gk_optarg);
        }
        break; 

      case CMD_SCALE:
  			if (gk_optarg) {
    			if ((params->scale = atoi(gk_optarg)) < 1)
            errexit("The scale must be greater than 0.\n");
  			}
  			break;

      case CMD_DBGLVL:
        if (gk_optarg) {
          params->dbglvl = atoi(gk_optarg);
          if (params->dbglvl < 0) 
            errexit("The -dbglvl must be non-negative.\n");
        }
        break;

      case CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(EXIT_SUCCESS);
        break;

      default:
        printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  /* Get the operation to be performed */
  if (argc-gk_optind != 1 && params->iftype != 1) {
    printf("Missing required parameters.\n  Use %s -help for a summary of the options.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  params->infile = gk_strdup(argv[gk_optind++]);

  return params;
}


