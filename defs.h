/*!
\file
\brief This file contains various constant definitions
*/

#ifndef __DEF_H__
#define __DEF_H__

#define MAXLINE         64*1024*1024
#define MAX_STRLEN      1024*128


/* command-line options */
#define CMD_TCTYPE              1
#define CMD_OTYPE               2
#define CMD_IFTYPE              3
#define CMD_DBGLVL              100
#define CMD_HELP                200
#define CMD_SCALE               50

/* tctypes */
#define TCTYPE_MAPJIK2D         1

/* otypes */
#define OTYPE_NONE              1
#define OTYPE_INCD              2
#define OTYPE_DECD              3

/* iftype */
#define IFTYPE_RMAT              1
#define IFTYPE_METIS             2

/* The text labels for the different tctypes */
static char tctypenames[][20] = 
                {"", "mapjik2d", ""}; 


/* The text labels for the different otypes */
static char otypenames[][10] = 
                {"", "none", "incd", "decd", 
                 ""}; 

/* The text labels for the different iftypes */
static char iftypenames[][10] = 
                {"", "rmat", "metis", ""}; 

#endif

