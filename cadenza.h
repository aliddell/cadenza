/*
 * Cadenza
 *
 * Jonathan Hauenstein <jdhauens@ncsu.edu>
 * Alan Liddell <acliddel@ncsu.edu>
 * Ian Haywood <ithaywoo@ncsu.edu>
 *
 * cadenza.h: header file for Cadenza
 */

#ifndef _BLUE_HARVEST_H
#define _BLUE_HARVEST_H

/*******
 * MPI *
 *******/
#include "mpi.h"

#include <errno.h>
#include <getopt.h>
#include <sys/ioctl.h>

/*****************************************************
 * data structures and functions from alphaCertified *
 *****************************************************/
#include "alphaCertified.h"

/******************************
 * constants for Cadenza *
 *****************************/
#define BH_UNSET -1
#define BH_USE_RATIONAL 0
#define BH_USE_FLOAT 1
#define BH_SUB_TOLERANCE 100
#define BH_NEWT_TOLERANCE 20
#define BH_DISCONTINUOUS -1
#define BH_CONTUNKNOWN 0
#define BH_CONTINUOUS 1
#define BH_PREC_MIN 32
#define BH_DESCENDING 1
#define BH_ASCENDING -1

/* verbosity flags */
#define BH_LACONIC 0
#define BH_CHATTY 1
#define BH_VERBOSE 2
#define BH_LOQUACIOUS 3

/* meta constants */
#define BH_PROGRAM_NAME "Cadenza"
#define BH_AUTHORS "Jonathan D. Hauenstein, Ian T. Haywood, and Alan C. Liddell, Jr."
#define BH_VERSION "1.1.0"
#define BH_BUILD_DATE "Apr 29 2014"

/* filenames */
#define BH_FSUMMARY  "summary.out"
#define BH_FCONT    "continuous.out"
#define BH_FDISCONT "discontinuous.out"
#define BH_FUNCERTIFIED "uncertified.out"
#define BH_FPTSOUT  "points.out"

/* string width and other memory-related constants */
#define BH_MAX_DATECHAR 25
#define BH_MAX_FILENAME 1000
#define BH_MAX_STRING   1000
#define BH_TERMWIDTH    80

/* exit codes */
#define BH_EXIT_SUCCESS 0
#define BH_EXIT_BADFILE 1 /* can't open file for reading or writing */
#define BH_EXIT_BADREAD 2 /* unexpected EOF in file read */
#define BH_EXIT_BADPARSE 3 /* general parse error */
#define BH_EXIT_BADDEF 4 /* system is not square, not enough points, bad t-value, &c*/
#define BH_EXIT_MEMORY 5 /* out of memory */
#define BH_EXIT_OTHER 7 /* something else */

/**************************************
 * global variables for cadenza.c *
 **************************************/
int verbosity, help_flag, ver_flag, default_precision, arithmetic_type, newton_tolerance, subd_tolerance, termwidth, sort_order, sigdig;
char pointsfile[BH_MAX_FILENAME], sysfile[BH_MAX_FILENAME], configfile[BH_MAX_FILENAME], *error_string;

/*********************
 * function pointers *
 *********************/
void (*read_system_file)(polynomial_system *system, void *v); /* reads in a polynomial system file, sets data */
int (*read_points_file)(void **t, void **x, int num_var); /* reads in a points file, sets data, returns number of points */
void (*fprint_input)(FILE *outfile, polynomial_system *system, void *v, void *t, void *x, int num_points); /* prints input files to stdout for debugging */
void (*test_system)(polynomial_system *system, void *v, void *t, void *x, int num_points, void **t_final, void **x_final, void **sing, int *tested, int *succeeded, int *failed, int *num_sing); /* certify H(x, t) */
void (*fprint_solutions)(void *t, void *x, int num_points); /* print solutions to H(x, t) in a program-readable format */

/*******************************************
 * function declarations for cadenza.c *
 *******************************************/
void getargs(int argc, char *argv[]); /* get command-line arguments, set flags and filenames */
int set_termwidth(); /* set the terminal width */
void set_function_pointers(); /* set function pointers depending on arithmetic */
void free_v(void *v); /* free [rational_]complex_vector v */
void free_t(void *t, int num_points); /* free mp[qf]_t *t */
void free_x(void *x, int num_points); /* free [rational_]complex_vector *x */

/**********************************
 * function declarations for io.c *
 **********************************/
void print_error(char *msg, FILE *outfile); /* prints msg to outfile */
void prog_info(); /* displays information about program, authors, libraries, &c on stderr */
void usage(); /* displays a helpful message about invocation on stderr */
void read_config_file(); /* reads in a configuration file */
void display_config(); /* displays the configuration (arithmetic type, precision, &c) on stderr */
polynomial parse_polynomial(FILE *sysfh, int num_var);
void print_system(FILE *outfile, polynomial_system *system);
void initialize_output_files(polynomial_system *system, void *v, void *t, void *x, int num_points); /* creates empty output files in the working directory */
void summarize(int tested, int succeeded, int failed, int singularities); /* print a summary to stdout */

/****************************************
 * function declarations for parallel.c *
 ****************************************/
int send_complex_number(complex_number c, int to); /* send a complex_number to process `to' */
int recv_complex_number(complex_number c, int from); /* receive a complex_number from process `from' */
int send_complex_vector(complex_vector v, int to);
int recv_complex_vector(complex_number v, int from);
int send_rational_complex_number(rational_complex_number c, int to);
int recv_rational_complex_number(rational_complex_number c, int from);
int send_polynomial(polynomial *p, int to);
int recv_polynomial(polynomial *p, int from);
int send_polynomial_system(polynomial_system *F, int to);
int recv_polynomial_system(polynomial_system *F, int from);

/*******************************************
 * function declarations for io_rational.c *
 *******************************************/
int compare_mpq(const void *a, const void *b);
void sort_points_rational(mpq_t *t, rational_complex_vector *x, int num_points);
void read_system_file_rational(polynomial_system *system, void *v);
int read_points_file_rational(void **t, void **x, int num_var);
void fprint_input_rational(FILE *outfile, polynomial_system *system, void *v, void *t, void *x, int num_points);
void print_points_rational(FILE *outfile, rational_complex_vector points);
void fprint_continuous_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right);
void fprint_discontinuous_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right);
void fprint_uncertain_rational(mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t beta_sqr_left, mpq_t beta_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right);
void fprint_solutions_rational(void *t, void *x, int num_points);

/****************************************
 * function declarations for io_float.c *
 ****************************************/
int compare_mpf(const void *a, const void *b);
void sort_points_float(mpf_t *t, complex_vector *x, int num_points);
void read_system_file_float(polynomial_system *system, void *v); /* see read_system_file(char *, void *) */
int read_points_file_float(void **t, void **x, int num_var); /* see read_points_file(char *, void **, int) */
void fprint_input_float(FILE *outfile, polynomial_system *system, void *v, void *t, void *x, int num_points);
void print_points_float(FILE *outfile, complex_vector points);
void fprint_continuous_float(mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right);
void fprint_discontinuous_float(mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right);
void fprint_uncertain_float(mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t beta_left, mpf_t beta_right, mpf_t gamma_left, mpf_t gamma_right);
void fprint_solutions_float(void *t, void *x, int num_points);

/************************************************
 * function declarations for certify_rational.c *
 ************************************************/
int lte_sqr_rational(mpq_t X_sqr, mpq_t Y_sqr, mpq_t Z_sqr);
void compute_norm_sqr_Jv_rational(polynomial_system *F, rational_complex_vector v, rational_complex_vector x, mpq_t *norm_sqr_Jv);
int is_continuous_rational(rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right);
int is_not_continuous_rational(rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right);
int test_continuity_rational(rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpq_t alpha_sqr_left, mpq_t alpha_sqr_right, mpq_t gamma_sqr_left, mpq_t gamma_sqr_right);
void subdivide_segment_rational(polynomial_system *base, rational_complex_vector v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, mpq_t *t_mid, rational_complex_vector *x_mid, int num_var);
void apply_tv_rational(polynomial_system *base, polynomial_system *F, mpq_t t, rational_complex_vector v);
int compute_abg_sqr_rational(rational_complex_vector points, polynomial_system *F, mpq_t *alpha, mpq_t *beta, mpq_t *gamma);
void test_pairwise_rational(polynomial_system *system, rational_complex_vector *v, mpq_t t_left, mpq_t t_right, rational_complex_vector x_left, rational_complex_vector x_right, int num_var, int iter, mpq_t **t_final, rational_complex_vector **x_final, rational_complex_vector **sing, int *tested, int *succeeded, int *failed, int *num_sing, int check_left);
void test_system_rational(polynomial_system *system, void *v, void *t, void *x, int num_points, void **t_final, void **x_final, void **sing, int *tested, int *succeeded, int *failed, int *num_sing); /* see test_system(polynomial_system*, ...) */

/*********************************************
 * function declarations for certify_float.c *
 *********************************************/
void compute_norm_Jv_float(polynomial_system *F, complex_vector v, complex_vector x, mpf_t *norm_Jv);
int is_continuous_float(complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t gamma_left, mpf_t gamma_right);
int is_not_continuous_float(complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t gamma_left, mpf_t gamma_right);
int test_continuity_float(complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, polynomial_system F_left, polynomial_system F_right, mpf_t alpha_left, mpf_t alpha_right, mpf_t gamma_left, mpf_t gamma_right);
void subdivide_segment_float(polynomial_system *base, complex_vector v, mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, mpf_t *t_mid, complex_vector *x_mid, int num_var);
void apply_tv_float(polynomial_system *base, polynomial_system *F, mpf_t t, complex_vector v);
int compute_abg_float(complex_vector points, polynomial_system *F, mpf_t *alpha, mpf_t *beta, mpf_t *gamma);
void test_pairwise_float(polynomial_system *system, complex_vector *v, mpf_t t_left, mpf_t t_right, complex_vector x_left, complex_vector x_right, int num_var, int iter, mpf_t **t_final, complex_vector **x_final, complex_vector **sing, int *tested, int *succeeded, int *failed, int *num_sing, int check_left);
void test_system_float(polynomial_system *system, void *v, void *t, void *x, int num_points, void **t_final, void **x_final, void **sing, int *tested, int *succeeded, int *failed, int *num_sing); /* see test_system(polynomial_system*, ...) */

#define mpq_set_min(_setme, _prima, _secunda) { if (mpq_cmp(_prima, _secunda) <= 0) { mpq_set(_setme, _prima); } \
    else { mpq_set(_setme, _secunda); }}
#define mpq_set_max(_setme, _prima, _secunda) { if (mpq_cmp(_prima, _secunda) >= 0) { mpq_set(_setme, _prima); } \
    else { mpq_set(_setme, _secunda); }}
#define mpf_set_min(_setme, _prima, _secunda) { if (mpf_cmp(_prima, _secunda) <= 0) { mpf_set(_setme, _prima); } \
    else { mpf_set(_setme, _secunda); }}
#define mpf_set_max(_setme, _prima, _secunda) { if (mpf_cmp(_prima, _secunda) >= 0) { mpf_set(_setme, _prima); } \
    else { mpf_set(_setme, _secunda); }}

#define subtract_rational_vector(_diff, _minuend, _subtraend) { int _i; int _size = _minuend->size; for (_i = 0; _i < _size; _i++) \
    { subtract_rational(_diff->coord[_i], _minuend->coord[_i], _subtraend->coord[_i]); }}
#define subtract_vector(_diff, _minuend, _subtraend) { int _i; int _size = _minuend->size; for (_i = 0; _i < _size; _i++) \
    { subtract(_diff->coord[_i], _minuend->coord[_i], _subtraend->coord[_i]); }}

#endif
