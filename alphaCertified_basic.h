/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   alphaCertified_basic.h: basic header file for alphaCertified
*/

#ifndef _ALPHACERTIFIED_BASIC_H
#define _ALPHACERTIFIED_BASIC_H

// header files to include
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <mpfr.h>

// version and error information
#define VERSION_STRING "1.2.0"
#define DATE_STRING "August 15, 2011"
#define ERROR_MESSAGE "Cadenza must now exit due to this error.\n"

// length of random rational number
#define RATIONALDIGITLENGTH 10

// exit error codes
#define ERROR_OTHER 1
#define ERROR_WRITE_PRIVILEGE 2   // write privilege problems - e.g. not able to create and write to files
#define ERROR_FILE_NOT_EXIST 3    // file exist problems - e.g. expected files do not exist
#define ERROR_INVALID_SIZE 4      // size problems - e.g. expected size not the same as the current size
#define ERROR_MEMORY_ALLOCATION 5 // memory problems - e.g. unable to allocate memory
#define ERROR_CONFIGURATION 6     // configuration problems - e.g. function calls called with wrong input values
#define ERROR_INPUT_SYSTEM 7      // input errors
#define ERROR_INPUT_SYNTAX 8      // syntax errors for input system
#define ERROR_LOOP_FAILURE 9      // loop failed to exit properly - e.g. fail-safe for 'infinite loops'

// error code for LU decomposition/solving
#define ERROR_LU_DECOMP 1

// code for LU failure but exact solution
#define EXACT_SOLUTION_LU_ERROR -1

#define exponentiate_mpq(_zexp,_aexp,_dexp) { \
  if (_dexp == 0) { mpq_set_ui(_zexp, 1, 1); } \
  else if (_dexp == 1) { mpq_set(_zexp,_aexp); } \
  else if (_dexp == 2) { mpq_mul(_zexp,_aexp,_aexp); } \
  else if (_dexp > 2) { /* use binary decomposition of _dexp to compute (_aexp)^_dexp */ \
    int _i,_size; mpq_t _cexp; mpq_init(_cexp); mpq_set(_cexp,_aexp); mpq_set_ui(_zexp, 1, 1); \
    mpz_t _pow_int; mpz_init_set_ui(_pow_int, _dexp); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
    for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { mpq_mul(_zexp,_zexp,_cexp); } mpq_mul(_cexp,_cexp,_cexp); } \
    mpq_clear(_cexp); mpz_clear(_pow_int); free(_base2); } \
  else { printf("59\n"); printf("\nERROR: Invalid power (%d).\n",_dexp); errExit(ERROR_CONFIGURATION); }}

// declaration for errExit
void errExit(int errorCode);

// declaration for sqrt
int sqrt_upper(mpq_t sqrt_x, mpq_t x, int digits);

#endif

