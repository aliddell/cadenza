/*
 * Blue Harvest (working title)
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * blueharvest.h: header file for Blue Harvest
 */

#ifndef _BLUE_HARVEST_H
#define _BLUE_HARVEST_H

#include <getopt.h>

/***************************************
 * data structures from alphaCertified *
 ***************************************/
#include "alphaCertified.h"

/******************************
 * constants for Blue Harvest *
 *****************************/
#ifndef BH_USE_RATIONAL
#define BH_USE_RATIONAL 0
#endif
#ifndef BH_USE_FLOAT
#define BH_USE_FLOAT 1
#endif

#ifndef BH_PROGRAM_NAME
#define BH_PROGRAM_NAME "Blue Harvest"
#endif
#ifndef BH_AUTHORS
#define BH_AUTHORS "Jonathan D. Hauenstein, Alan C. Liddell, Jr., and Ian T. Haywood"
#endif
#ifndef BH_VERSION
#define BH_VERSION "0.0.1"
#endif
#ifndef BH_BUILD_DATE
#define BH_BUILD_DATE "Aug 01, 2013"
#endif
#ifndef BH_MAX_DATECHAR
#define BH_MAX_DATECHAR 25
#endif

#ifndef BH_EXIT_BADFILE
#define BH_EXIT_BADFILE 1
#endif

/**************************************
 * global variables for blueharvest.c *
 **************************************/
int verbose_flag, help_flag, arithmetic_type;
char *pointsfile, *sysfile;

/*******************************************
 * function declarations for blueharvest.c *
 *******************************************/
void getargs(int argc, char *argv[]);

/**********************************
 * function declarations for io.c *
 **********************************/
void print_error(char *msg);
void prog_info();
void usage();
void display_config();
polynomial parse_polynomial(char *input);

#endif
