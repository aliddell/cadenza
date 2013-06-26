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
#ifndef BH_LACONIC
#define BH_LACONIC 0
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
#ifndef BH_MAX_FILENAME
#define BH_MAX_FILENAME 52
#endif
#ifndef BH_TERMWIDTH
#define BH_TERMWIDTH 80
#endif

/* exit codes */
#ifndef BH_EXIT_SUCCESS
#define BH_EXIT_SUCCESS 0
#endif
#ifndef BH_EXIT_BADFILE
#define BH_EXIT_BADFILE 1 /* can't open file for reading or writing */
#endif
#ifndef BH_EXIT_BADREAD
#define BH_EXIT_BADREAD 2 /* unexpected EOF in file read */
#endif
#ifndef BH_EXIT_BADPARSE
#define BH_EXIT_BADPARSE 3 /* general parse error */
#endif
#ifndef BH_EXIT_BADDEF
#define BH_EXIT_BADDEF 4 /* system is not square, not enough points */
#endif
#ifndef BH_EXIT_MEMORY
#define BH_EXIT_MEMORY 5 /* out of memory */
#endif

/**************************************
 * global variables for blueharvest.c *
 **************************************/
int verbosity, help_flag, default_precision, arithmetic_type;
char *pointsfile, *sysfile;

/*********************
 * function pointers *
 *********************/
void (*read_system_file)(char *filename, polynomial_system *system, void *v); /* reads in a polynomial system file, sets data */
int (*read_points_file)(char *filename, void **t, void **w, int num_var); /* reads in a points file, sets data, returns number of points */
void (*free_system)(void *system, void *v); /* frees dynamically-allocated memory for polynomial system */
void (*free_vector)(void *w, void *t, int num_points); /* frees dynamically-allocated memory for array of points vectors */
void (*test_pairwise)(polynomial_system *system, configurations *config, void *v, void *t, void *w, int num_points); /* applies f(x) + tv for all t */

/*******************************************
 * function declarations for blueharvest.c *
 *******************************************/
void free_system_rational(void *system, void *v); /* see free_system(void *) */
void free_system_float(void *system, void *v); /* see free_system(void *) */
void free_vector_rational(void *w, void *t, int num_points); /* see free_vector(void *, int) */
void free_vector_float(void *w, void *t, int num_points); /* see free_vector(void *, int) */
void getargs(int argc, char *argv[]); /* get command-line arguments, set flags and filenames */
void set_function_pointers(); /* set function pointers depending on arithmetic */

/**********************************
 * function declarations for io.c *
 **********************************/
void print_error(char *msg); /* prints msg to stderr */
void prog_info(); /* displays information about program, authors, libraries, &c on stderr */
void usage(); /* displays a helpful message about invocation on stderr */
void display_config(); /* displays the configuration (arithmetic type, precision, &c) on stderr */

/*******************************************
 * function declarations for io_rational.c *
 *******************************************/
void read_system_file_rational(char *filename, polynomial_system *system, void *v); /* see read_system_file(char *, void *) */
polynomial parse_polynomial_rational(FILE *sysfile, char *filename, int num_var); /* parses exponents of one polynomial, calls parse_complex */
void parse_complex_rational(char *str_real, char *str_imag, rational_complex_number c); /* parses coefficients of one polynomial */
int read_points_file_rational(char *filename, void **t, void **w, int num_var); /* see read_points_file(char *, void **, int) */
void print_points_rational(rational_complex_vector *points, int num_var);
void print_system_rational(polynomial_system *system);

/****************************************
 * function declarations for io_float.c *
 ****************************************/
void read_system_file_float(char *filename, polynomial_system *system, void *v); /* see read_system_file(char *, void *) */
polynomial parse_polynomial_float(FILE *sysfile, char *filename, int num_var); /* parses exponents of one polynomial, calls parse_complex */
void parse_complex_float(char *str_real, char *str_imag, complex_number c); /* parses coefficients of one polynomial */
int read_points_file_float(char *filename, void **t, void **w, int num_var); /* see read_points_file(char *, void **, int) */
void print_points_float(complex_vector *points, int num_var);
void print_system_float(polynomial_system *system);

/************************************************
 * function declarations for certify_rational.c *
 ************************************************/
void apply_tv_rational(polynomial_system *base, polynomial_system *F, mpq_t t, rational_complex_vector v); /* add (t_i)(v_i) to F */
void test_pairwise_rational(polynomial_system *system, configurations *config, void *v, void *t, void *w, int num_points); /* see test_pairwise(polynomial_system*, ...) */

/*********************************************
 * function declarations for certify_float.c *
 *********************************************/
void apply_tv_float(polynomial_system *base, polynomial_system *F, mpf_t t, complex_vector v); /* add (t_i)(v_i) to F */
void test_pairwise_float(polynomial_system *system, configurations *config, void *v, void *t, void *w, int num_points); /* see test_pairwise(polynomial_system*, ...) */

#endif
