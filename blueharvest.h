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

#include <errno.h>
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

/* verbosity flags */
#ifndef BH_TACITURN
#define BH_TACITURN 0
#endif
#ifndef BH_CHATTY
#define BH_CHATTY 1
#endif
#ifndef BH_VERBOSE
#define BH_VERBOSE 2
#endif
#ifndef BH_LOQUACIOUS
#define BH_LOQUACIOUS 3
#endif

/* meta constants */
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
#define BH_BUILD_DATE "Jul 26, 2013"
#endif

/* string width and other memory-related constants */
#ifndef BH_MAX_DATECHAR
#define BH_MAX_DATECHAR 25
#endif
#ifndef BH_TERMWIDTH
#define BH_TERMWIDTH 80
#endif

/* exit codes */
#ifndef BH_EXIT_SUCCESS
#define BH_EXIT_SUCCESS 0
#endif
/* can't open file for reading or writing */
#ifndef BH_EXIT_BADFILE
#define BH_EXIT_BADFILE 1
#endif
/* unexpected EOF in file read */
#ifndef BH_EXIT_BADREAD
#define BH_EXIT_BADREAD 2
#endif
/* general parse error */
#ifndef BH_EXIT_BADPARSE
#define BH_EXIT_BADPARSE 3
#endif

/**************************************
 * global variables for blueharvest.c *
 **************************************/
int verbosity, help_flag, default_precision, arithmetic_type;
char *pointsfile, *sysfile;

/*********************
 * function pointers *
 *********************/
polynomial_system (*read_system_file)(char *filename);
void (*free_system)(void *ps);

/*******************************************
 * function declarations for blueharvest.c *
 *******************************************/
void free_system_rational(void *ps);
void free_system_float(void *ps);
void getargs(int argc, char *argv[]);
void set_function_pointers();

/**********************************
 * function declarations for io.c *
 **********************************/
void print_error(char *msg);
void prog_info();
void usage();
void display_config();

/*******************************************
 * function declarations for io_rational.c *
 *******************************************/
polynomial_system read_system_file_rational(char *filename);
polynomial parse_polynomial_rational(FILE *sysfile, char *filename, int num_vars);
void parse_coeff_rational(char *str_coeff_real, char *str_coeff_imag, rational_complex_number c);

/****************************************
 * function declarations for io_float.c *
 ****************************************/
polynomial_system read_system_file_float(char *filename);
polynomial parse_polynomial_float(FILE *sysfile, char *filename, int num_vars);
void parse_coeff_float(char *str_coeff_real, char *str_coeff_imag, complex_number c);

#endif
