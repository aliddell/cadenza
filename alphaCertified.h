/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   alphaCertified.h: header file for alphaCertified
*/

#ifndef _ALPHACERTIFIED_H
#define _ALPHACERTIFIED_H

// header files to include
#include "alphaCertified_basic.h"
#include <mpf2mpfr.h>

// basic complex number, vector, and matrix
typedef struct
{
  mpf_t re, im;
} _complex_number;

typedef struct
{
  mpq_t re, im;
} _rational_complex_number;

typedef _complex_number complex_number[1];
typedef _rational_complex_number rational_complex_number[1];

typedef struct 
{
  complex_number *coord;
  int alloc_size; // allocated size
  int curr_prec;  // current precision
  int size;       // size of the point
} _complex_vector;

typedef struct 
{
  rational_complex_number *coord;
  int alloc_size; // allocated size
  int size;       // size of the point
} _rational_complex_vector;

typedef _complex_vector complex_vector[1];
typedef _rational_complex_vector rational_complex_vector[1];

typedef struct 
{
  complex_number **entry;
  int alloc_rows; // allocated number of rows
  int alloc_cols; // allocated number of cols
  int curr_prec;  // current precision
  int rows;       // number of rows for the matrix
  int cols;       // number of cols for the matrix
} _complex_matrix;

typedef struct 
{
  rational_complex_number **entry;
  int alloc_rows; // allocated number of rows
  int alloc_cols; // allocated number of cols
  int rows;       // number of rows for the matrix
  int cols;       // number of cols for the matrix
} _rational_complex_matrix;

typedef _complex_matrix complex_matrix[1];
typedef _rational_complex_matrix rational_complex_matrix[1];

// structure for handling points
typedef struct
{
  complex_vector origX; // original point
  complex_vector x;     // current point
  complex_vector Nx;    // point after 1 newton iteration
  mpf_t norm_x;         // ||x||_2
  int isApproxSoln;     // 1 if approx solution, 0 if unknown
  int isActive;         // 1 if active, <= 0 if corresponds to same solution as -isActive
  int isReal;           // 1 if corresponds with a real solution
  complex_number origAlpha;
  complex_number origBeta;
  complex_number origGamma;
  complex_number alpha;
  complex_number beta;
  complex_number gamma;
} point_struct;

typedef struct
{
  rational_complex_vector origX; // original point
  rational_complex_vector x;     // current point
  rational_complex_vector Nx;    // point after 1 newton iteration
  mpq_t norm_sqr_x;              // ||x||_2^2
  int isApproxSoln;              // 1 if approx solution, 0 if unknown
  int isActive;                  // 1 if active, <= 0 if corresponds to same solution as -isActive
  int isReal;                    // 1 if corresponds with a real solution
  rational_complex_number origAlpha_sqr;
  rational_complex_number origBeta_sqr;
  rational_complex_number origGamma_sqr;
  rational_complex_number alpha_sqr;
  rational_complex_number beta_sqr;
  rational_complex_number gamma_sqr;
} rational_point_struct;

// polynomial
typedef struct
{
  int numVariables;                // number of variables
  int numTerms;                    // number of terms
  int degree;                      // degree of the polynomial
  int isReal;                      // determine if polynomial has real coefficients
  mpq_t norm_sqr;                  // norm^2 of the polynomial
  rational_complex_number *coeff;  // coefficients
  int **exponents;                 // monomial exponents
} polynomial;

// exponential function of the for y_i - g(beta_i * x_{j_i}) where g is one of exp, sin, or cos
typedef struct
{
  int xIndex;                    // variable index for x
  int yIndex;                    // variable index for y
  char expFunction;              // type of exponential function: 'X', 'C', or 'S'
  int isHyperbolic;              // determine if we are using sin & cos or sinh & cosh
  int isReal;                    // determine if beta is real
  rational_complex_number beta;  // constant in the exponential function: G(beta * X_i)
} exponential;

// polynomial system
typedef struct
{
  int numVariables;          // number of variables
  int numPolynomials;        // number of polynomials
  int maximumDegree;         // maximum degree
  int isReal;                // determine if N_f is a real map
  mpq_t norm_sqr;            // norm^2 of the polynomial system
  polynomial *polynomials;   // polynomials in the system
  int numExponentials;       // number of exponential functions
  exponential *exponentials; // exponentials in the system
} polynomial_system;

// configurations
typedef struct
{
  int arithmeticType;      // either rational (0) or floating point (1)
  int startingPrecision;   // starting precision, if using floating point arithmetic
  int algorithm;           // certify solutions (0) & distinct certify (1) & real distinct certify (2)
  int refineDigits;        // final number of digits to report certified solutions
  int numRandomSystems;    // number of random systems to use when dealing with overdetermined systems (>= 2)
  int randomDigits;        // number of digits identical to declare an overdetermined system soln a soln
  unsigned int randomSeed; // random seed
  int newtonOnly;          // whether to run Newton iterations only
  int newtonIts;           // number of Newton iterations to perform
  int realityCheck;        // -1 - assume real, 0 - only coefficient test, 1 - coeff & conj test, 2 - coeff, conj, and reality check
  int realityTest;         // 0 - local test, 1 - global test
} configurations;

// initialize number, vector, and matrix
#define initialize_number(_num) { mpf_init((_num)->re); mpf_init((_num)->im); }
#define initialize_number2(_num,_prec) { mpf_init2((_num)->re,_prec); mpf_init2((_num)->im,_prec); }

#define initialize_rational_number(_num) { mpq_init((_num)->re); mpq_init((_num)->im); }

#define initialize_vector2(_vec,_size,_prec) { int _i; (_vec)->coord = (complex_number *)malloc((_size) * sizeof(complex_number)); \
  for (_i = 0; _i < _size; _i++) initialize_number2((_vec)->coord[_i],_prec); (_vec)->curr_prec = _prec; (_vec)->alloc_size = (_vec)->size = _size; }
#define initialize_vector(_vec,_size) { int _p = mpf_get_default_prec(); initialize_vector2(_vec,_size,_p); }

#define initialize_rational_vector(_vec,_size) { int _i; (_vec)->coord = (rational_complex_number *)malloc((_size) * sizeof(rational_complex_number)); \
  for (_i = 0; _i < _size; _i++) initialize_rational_number((_vec)->coord[_i]); (_vec)->alloc_size = (_vec)->size = _size; }

#define initialize_matrix2(_mat,_rows,_cols,_prec) { int _i,_j; (_mat)->entry = (complex_number **)malloc((_rows) * sizeof(complex_number *)); \
  for (_i = 0; _i < _rows; _i++) { (_mat)->entry[_i] = (complex_number *)malloc((_cols) * sizeof(complex_number)); \
    for (_j = 0; _j < _cols; _j++) initialize_number2((_mat)->entry[_i][_j], _prec); } (_mat)->curr_prec = _prec; \
  (_mat)->alloc_rows = (_mat)->rows = _rows; (_mat)->alloc_cols = (_mat)->cols = _cols; }
#define initialize_matrix(_mat,_rows,_cols) { int _p = mpf_get_default_prec(); initialize_matrix2(_mat,_rows,_cols,_p); }

#define initialize_rational_matrix(_mat,_rows,_cols) { int _i,_j; (_mat)->entry = (rational_complex_number **)malloc((_rows) * sizeof(rational_complex_number)); \
  for (_i = 0; _i < _rows; _i++) { (_mat)->entry[_i] = (rational_complex_number *)malloc((_cols) * sizeof(rational_complex_number)); \
    for (_j = 0; _j < _cols; _j++) initialize_rational_number((_mat)->entry[_i][_j]); } \
  (_mat)->alloc_rows = (_mat)->rows = _rows; (_mat)->alloc_cols = (_mat)->cols = _cols; }

// change size of vector and matrix
#define change_size_vector(_vec,_new_size) { if ((_vec)->alloc_size != _new_size) { int _i; \
  for (_i = (_vec)->alloc_size - 1; _i >= 0; _i--) clear_number((_vec)->coord[_i]); \
  (_vec)->coord = (complex_number *)errRealloc((_vec)->coord, _new_size * sizeof(complex_number)); \
  for (_i = 0; _i < _new_size; _i++) initialize_number((_vec)->coord[_i]); (_vec)->size = (_vec)->alloc_size = _new_size; }}

#define change_size_rational_vector(_vec,_new_size) { if ((_vec)->alloc_size != _new_size) { int _i; \
  for (_i = (_vec)->alloc_size - 1; _i >= 0; _i--) clear_rational_number((_vec)->coord[_i]); \
  (_vec)->coord = (rational_complex_number *)errRealloc((_vec)->coord, _new_size * sizeof(rational_complex_number)); \
   for (_i = 0; _i < _new_size; _i++) initialize_rational_number((_vec)->coord[_i]); (_vec)->size = (_vec)->alloc_size = _new_size; }}

#define change_size_matrix(_mat,_new_rows,_new_cols) { \
  if ((_mat)->alloc_rows < _new_rows) { /* increase the rows first and then change the size of the columns */ \
    increase_rows_matrix(_mat,_new_rows); if ((_mat)->alloc_cols < _new_cols) { increase_cols_matrix(_mat,_new_cols); }} \
  else { /* change the size of the columns first and then change the rows */ \
    if ((_mat)->alloc_cols < _new_cols) { increase_cols_matrix(_mat,_new_cols); } \
    else if ((_mat)->alloc_cols > _new_cols) { decrease_cols_matrix(_mat,_new_cols); } \
    if ((_mat)->alloc_rows > _new_rows) { decrease_rows_matrix(_mat,_new_rows); }}}

#define change_size_rational_matrix(_mat,_new_rows,_new_cols) { \
  if ((_mat)->alloc_rows < _new_rows) { /* increase the rows first and then change the size of the columns */ \
    increase_rows_rational_matrix(_mat,_new_rows); if ((_mat)->alloc_cols < _new_cols) { increase_cols_rational_matrix(_mat,_new_cols); }} \
  else { /* change the size of the columns first and then change the rows */ \
    if ((_mat)->alloc_cols < _new_cols) { increase_cols_rational_matrix(_mat,_new_cols); } \
    else if ((_mat)->alloc_cols > _new_cols) { decrease_cols_rational_matrix(_mat,_new_cols); } \
    if ((_mat)->alloc_rows > _new_rows) { decrease_rows_rational_matrix(_mat,_new_rows); }}}

// increase rows & columns of a matrix
#define increase_rows_matrix(_mat,_new_rows) { if ((_mat)->alloc_rows < _new_rows) { int _i,_j; \
  (_mat)->entry = (complex_number **)errRealloc((_mat)->entry, (_new_rows) * sizeof(complex_number *)); \
  for (_i = (_mat)->alloc_rows; _i < _new_rows; _i++) { (_mat)->entry[_i] = (complex_number *)malloc((_mat)->alloc_cols * sizeof(complex_number)); \
    for (_j = 0; _j < (_mat)->alloc_cols; _j++) initialize_number2((_mat)->entry[_i][_j],(_mat)->curr_prec); } (_mat)->alloc_rows = (_mat)->rows = _new_rows; }}

#define increase_rows_rational_matrix(_mat,_new_rows) { if ((_mat)->alloc_rows < _new_rows) { int _i,_j; \
  (_mat)->entry = (rational_complex_number **)errRealloc((_mat)->entry, (_new_rows) * sizeof(rational_complex_number *)); \
  for (_i = (_mat)->alloc_rows; _i < _new_rows; _i++) { (_mat)->entry[_i] = (rational_complex_number *)malloc((_mat)->alloc_cols * sizeof(rational_complex_number)); \
    for (_j = 0; _j < (_mat)->alloc_cols; _j++) initialize_rational_number((_mat)->entry[_i][_j]); } (_mat)->alloc_rows = (_mat)->rows = _new_rows; }}

#define increase_cols_matrix(_mat,_new_cols) { if ((_mat)->alloc_cols < _new_cols) { int _i,_j; \
  for (_i = 0; _i < (_mat)->alloc_rows; _i++) { (_mat)->entry[_i] = (complex_number *)errRealloc((_mat)->entry[_i], (_new_cols) * sizeof(complex_number)); \
    for (_j = (_mat)->alloc_cols; _j  < _new_cols; _j++) initialize_number2((_mat)->entry[_i][_j], (_mat)->curr_prec); } \
  (_mat)->alloc_cols = (_mat)->cols = _new_cols; }}

#define increase_cols_rational_matrix(_mat,_new_cols) { if ((_mat)->alloc_cols < _new_cols) { int _i,_j; \
  for (_i = 0; _i < (_mat)->alloc_rows; _i++) { (_mat)->entry[_i] = (rational_complex_number *)errRealloc((_mat)->entry[_i], (_new_cols) * sizeof(rational_complex_number)); \
    for (_j = (_mat)->alloc_cols; _j  < _new_cols; _j++) initialize_rational_number((_mat)->entry[_i][_j]); } (_mat)->alloc_cols = (_mat)->cols = _new_cols; }}

// decrease rows & columns of a matrix
#define decrease_rows_matrix(_mat,_new_rows) { if ((_mat)->alloc_rows > _new_rows) { int _i,_j; /* clear extra rows */ \
  for (_i = _new_rows; _i < (_mat)->alloc_rows; _i++) if ((_mat)->entry[_i] != NULL) { \
    for (_j = 0; _j < (_mat)->alloc_cols; _j++) clear_number((_mat)->entry[_i][_j]); free((_mat)->entry[_i]); } \
  (_mat)->entry = (complex_number **)errRealloc((_mat)->entry, (_new_rows) * sizeof(complex_number)); (_mat)->alloc_rows = (_mat)->rows = _new_rows; }}

#define decrease_rows_rational_matrix(_mat,_new_rows) { if ((_mat)->alloc_rows > _new_rows) { int _i,_j; /* clear extra rows */ \
  for (_i = _new_rows; _i < (_mat)->alloc_rows; _i++) if ((_mat)->entry[_i] != NULL) { \
    for (_j = 0; _j < (_mat)->alloc_cols; _j++) clear_rational_number((_mat)->entry[_i][_j]); free((_mat)->entry[_i]); } \
  (_mat)->entry = (rational_complex_number **)errRealloc((_mat)->entry, (_new_rows) * sizeof(rational_complex_number)); (_mat)->alloc_rows = (_mat)->rows = _new_rows; }}

#define decrease_cols_matrix(_mat,_new_cols) { if ((_mat)->alloc_cols > _new_cols) { int _i,_j; /* clear extra cols */ \
  for (_i = 0; _i < (_mat)->alloc_rows; _i++) { for (_j = _new_cols; _j < (_mat)->alloc_cols; _j++) clear_number((_mat)->entry[_i][_j]); \
    (_mat)->entry[_i] = (complex_number *)errRealloc((_mat)->entry[_i], (_new_cols) * sizeof(complex_number)); } (_mat)->alloc_cols = (_mat)->cols = _new_cols; }}

#define decrease_cols_rational_matrix(_mat,_new_cols) { if ((_mat)->alloc_cols > _new_cols) { int _i,_j; /* clear extra cols */ \
  for (_i = 0; _i < (_mat)->alloc_rows; _i++) { for (_j = _new_cols; _j < (_mat)->alloc_cols; _j++) clear_rational_number((_mat)->entry[_i][_j]); \
    (_mat)->entry[_i] = (rational_complex_number *)errRealloc((_mat)->entry[_i], (_new_cols) * sizeof(rational_complex_number)); } (_mat)->alloc_cols = (_mat)->cols = _new_cols; }}

// set the precision for number, vector, and matrix
#define setPrec_number(_num,_prec) { mpf_set_prec((_num)->re, _prec); mpf_set_prec((_num)->im, _prec); }

#define setPrec_vector(_vec,_prec) { if ((_vec)->curr_prec != _prec) { int _i; \
  for (_i = 0; _i < (_vec)->alloc_size; _i++) setPrec_number((_vec)->coord[_i],_prec); (_vec)->curr_prec = _prec; }}

#define setPrec_matrix(_mat,_prec) { if ((_mat)->curr_prec != _prec) { int _i,_j; \
  for (_i = 0; _i < (_mat)->alloc_rows; _i++) for (_j = 0; _j < (_mat)->alloc_cols; _j++) \
    setPrec_number((_mat)->entry[_i][_j],_prec); (_mat)->curr_prec = _prec; }}

// clear number, vector, and matrix
#define clear_number(_num) { mpf_clear((_num)->re); mpf_clear((_num)->im); }
#define clear_rational_number(_num) { mpq_clear((_num)->re); mpq_clear((_num)->im); }

#define clear_vector(_vec) { int _i; for (_i = (_vec)->alloc_size - 1; _i >= 0; _i--) clear_number((_vec)->coord[_i]); free((_vec)->coord); \
  (_vec)->coord = NULL; (_vec)->alloc_size = (_vec)->size = 0; }

#define clear_rational_vector(_vec) { int _i; for (_i = (_vec)->alloc_size - 1; _i >= 0; _i--) clear_rational_number((_vec)->coord[_i]); free((_vec)->coord); \
  (_vec)->coord = NULL; (_vec)->alloc_size = (_vec)->size = 0; }

#define clear_matrix(_mat) { int _i,_j; for (_i = (_mat)->alloc_rows - 1; _i >= 0; _i--) { for (_j = (_mat)->alloc_cols - 1; _j >= 0; _j--) \
  clear_number((_mat)->entry[_i][_j]); free((_mat)->entry[_i]); (_mat)->entry[_i] = NULL; } free((_mat)->entry); (_mat)->entry = NULL; \
  (_mat)->alloc_rows = (_mat)->alloc_cols = (_mat)->rows = (_mat)->cols = 0; }

#define clear_rational_matrix(_mat) { int _i,_j; for (_i = (_mat)->alloc_rows - 1; _i >= 0; _i--) { for (_j = (_mat)->alloc_cols - 1; _j >= 0; _j--) \
  clear_rational_number((_mat)->entry[_i][_j]); free((_mat)->entry[_i]); (_mat)->entry[_i] = NULL; } free((_mat)->entry); (_mat)->entry = NULL; \
  (_mat)->alloc_rows = (_mat)->alloc_cols = (_mat)->rows = (_mat)->cols = 0; }

// copy number, vector & matrix
#define set_number(_a,_b) { mpf_set((_a)->re, (_b)->re); mpf_set((_a)->im, (_b)->im); }
#define set_rational_number(_a,_b) { mpq_set((_a)->re, (_b)->re); mpq_set((_a)->im, (_b)->im); }

#define copy_vector(_a,_b) { int _i,_size = (_b)->size; change_size_vector(_a,_size); for (_i = 0; _i < _size; _i++) set_number((_a)->coord[_i], (_b)->coord[_i]); }
#define copy_rational_vector(_a,_b) { int _i,_size = (_b)->size; change_size_rational_vector(_a,_size); for (_i = 0; _i < _size; _i++) \
  set_rational_number((_a)->coord[_i], (_b)->coord[_i]); }

#define copy_matrix(_a,_b) { int _i,_j,_rows = (_b)->rows,_cols = (_b)->cols; change_size_matrix(_a,_rows,_cols); \
  for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) set_number((_a)->entry[_i][_j], (_b)->entry[_i][_j]); }
#define copy_rational_matrix(_a,_b) { int _i,_j,_rows = (_b)->rows,_cols = (_b)->cols; change_size_rational_matrix(_a,_rows,_cols); \
  for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) set_rational_number((_a)->entry[_i][_j], (_b)->entry[_i][_j]); }

// set number, vector, and matrix to 0
#define set_zero_number(_num) { mpf_set_ui((_num)->re, 0); mpf_set_ui((_num)->im, 0); }
#define set_zero_rational_number(_num) { mpq_set_ui((_num)->re, 0, 1); mpq_set_ui((_num)->im, 0, 1); }

#define set_zero_vector(_vec) { int _i; for (_i = 0; _i < (_vec)->size; _i++) set_zero_number((_vec)->coord[_i]); }
#define set_zero_rational_vector(_vec) { int _i; for (_i = 0; _i < (_vec)->size; _i++) set_zero_rational_number((_vec)->coord[_i]); }

#define set_zero_matrix(_mat) { int _i,_j; for (_i = 0; _i < (_mat)->rows; _i++) for (_j = 0; _j < (_mat)->cols; _j++) \
  set_zero_number((_mat)->entry[_i][_j]); }

#define set_zero_rational_matrix(_mat) { int _i,_j; for (_i = 0; _i < (_mat)->rows; _i++) for (_j = 0; _j < (_mat)->cols; _j++) \
  set_zero_rational_number((_mat)->entry[_i][_j]); }

// set number and vector to 1
#define set_one_number(_num) { mpf_set_ui((_num)->re, 1); mpf_set_ui((_num)->im, 0); }
#define set_one_rational_number(_num) { mpq_set_ui((_num)->re, 1, 1); mpq_set_ui((_num)->im, 0, 1); }

#define set_one_vector(_vec) { int _i; for (_i = 0; _i < (_vec)->size; _i++) set_one_number((_vec)->coord[_i]); }
#define set_one_rational_vector(_vec) { int _i; for (_i = 0; _i < (_vec)->size; _i++) set_one_rational_number((_vec)->coord[_i]); }

// set matrix to identity
#define identity_matrix(_mat) { int _i,_j; for (_i = 0; _i < (_mat)->rows; _i++) for (_j = 0; _j < (_mat)->cols; _j++) \
  if (_i == _j) { set_one_number((_mat)->entry[_i][_j]); } else { set_zero_number((_mat)->entry[_i][_j]); }}

#define identity_rational_matrix(_mat) { int _i,_j; for (_i = 0; _i < (_mat)->rows; _i++) for (_j = 0; _j < (_mat)->cols; _j++) \
  if (_i == _j) { set_one_rational_number((_mat)->entry[_i][_j]); } else { set_zero_rational_number((_mat)->entry[_i][_j]); }}

// convert a rational to floating point 
#define convert_rational_number(_c, _c_rat) { mpf_set_q((_c)->re, (_c_rat)->re); mpf_set_q((_c)->im, (_c_rat)->im); }
#define convert_rational_vector(_v, _v_rat) { int _i, _size = (_v_rat)->size; change_size_vector(_v,_size); \
  for (_i = 0; _i < _size; _i++) convert_rational_number((_v)->coord[_i], (_v_rat)->coord[_i]); }

// negate a complex number
#define negate(_z,_a) { mpf_neg((_z)->re, (_a)->re); mpf_neg((_z)->im, (_a)->im); }
#define negate_rational(_z,_a) { mpq_neg((_z)->re, (_a)->re); mpq_neg((_z)->im, (_a)->im); }

// conjugate a complex number
#define conjugate(_z,_a) { mpf_set((_z)->re, (_a)->re); mpf_neg((_z)->im, (_a)->im); }
#define conjugate_rational(_z,_a) { mpq_set((_z)->re, (_a)->re); mpq_neg((_z)->im, (_a)->im); }

// conjugate a complex vector
#define conjugate_vector(_z,_a) { int _i,_size = (_a)->size; change_size_vector(_z,_size); \
  for (_i = 0; _i < _size; _i++) conjugate((_z)->coord[_i], (_a)->coord[_i]); }
#define conjugate_rational_vector(_z,_a) { int _i,_size = (_a)->size; change_size_rational_vector(_z,_size); \
  for (_i = 0; _i < _size; _i++) conjugate_rational((_z)->coord[_i], (_a)->coord[_i]); } 

// add two complex numbers
#define add(_z,_a,_b) { mpf_add((_z)->re, (_a)->re, (_b)->re); mpf_add((_z)->im, (_a)->im, (_b)->im); }
#define add_rational(_z,_a,_b) { mpq_add((_z)->re, (_a)->re, (_b)->re); mpq_add((_z)->im, (_a)->im, (_b)->im); }

// subtract two complex numbers
#define subtract(_z,_a,_b) { mpf_sub((_z)->re, (_a)->re, (_b)->re); mpf_sub((_z)->im, (_a)->im, (_b)->im); }
#define subtract_rational(_z,_a,_b) { mpq_sub((_z)->re, (_a)->re, (_b)->re); mpq_sub((_z)->im, (_a)->im, (_b)->im); }

// multiply two complex numbers
#define multiply(_z,_a,_b) { complex_number _c; initialize_number(_c); \
  mpf_mul((_c)->re, (_a)->re, (_b)->re); mpf_mul((_c)->im, (_a)->im, (_b)->im); mpf_sub((_c)->re, (_c)->re, (_c)->im); \
  mpf_mul((_c)->im, (_a)->re, (_b)->im); mpf_mul((_z)->im, (_a)->im, (_b)->re); mpf_add((_z)->im, (_z)->im, (_c)->im); \
  mpf_set((_z)->re, (_c)->re); clear_number(_c); }

#define multiply_rational(_z,_a,_b) { rational_complex_number _c; initialize_rational_number(_c); \
  mpq_mul((_c)->re, (_a)->re, (_b)->re); mpq_mul((_c)->im, (_a)->im, (_b)->im); mpq_sub((_c)->re, (_c)->re, (_c)->im); \
  mpq_mul((_c)->im, (_a)->re, (_b)->im); mpq_mul((_z)->im, (_a)->im, (_b)->re); mpq_add((_z)->im, (_z)->im, (_c)->im); \
  mpq_set((_z)->re, (_c)->re); clear_rational_number(_c); }

// multiply a complex number by a positive integer
#define multiply_positive_int(_z,_a,_b) { mpf_mul_ui((_z)->re, (_a)->re, _b); mpf_mul_ui((_z)->im, (_a)->im, _b); }
#define multiply_rational_positive_int(_z,_a,_b) { mpq_t _c; mpq_init(_c); mpq_set_ui(_c, _b, 1); \
  mpq_mul((_z)->re, (_a)->re, _c); mpq_mul((_z)->im, (_a)->im, _c); mpq_clear(_c); }

// multiply a complex number by a rational number
#define multiply_number(_z,_a,_b) { mpf_mul((_z)->re, (_a)->re, _b); mpf_mul((_z)->im, (_a)->im, _b); }
#define multiply_rational_number(_z,_a,_b) { mpq_mul((_z)->re, (_a)->re, _b); mpq_mul((_z)->im, (_a)->im, _b); }

// divide a complex number by another complex number
#define divide(_z,_a,_b) { complex_number _c; mpf_t _d; initialize_number(_c); mpf_init(_d); \
  mpf_mul((_c)->re, (_b)->re, (_b)->re); mpf_mul((_c)->im, (_b)->im, (_b)->im); mpf_add((_c)->re, (_c)->re, (_c)->im); mpf_ui_div((_d), 1, (_c)->re); \
  mpf_mul((_c)->re, (_a)->re, (_b)->re); mpf_mul((_c)->im, (_a)->im, (_b)->im); mpf_add((_c)->re, (_c)->re, (_c)->im); mpf_mul((_c)->re, (_c)->re, (_d)); \
  mpf_mul((_c)->im, (_a)->im, (_b)->re); mpf_mul((_z)->im, (_a)->re, (_b)->im); mpf_sub((_z)->im, (_c)->im, (_z)->im); mpf_mul((_z)->im, (_z)->im, (_d)); \
  mpf_set((_z)->re, (_c)->re); clear_number(_c); mpf_clear(_d); }

#define divide_rational(_z,_a,_b) { rational_complex_number _c; mpq_t _d; initialize_rational_number(_c); mpq_init(_d); \
  mpq_mul((_c)->re, (_b)->re, (_b)->re); mpq_mul((_c)->im, (_b)->im, (_b)->im); mpq_add((_c)->re, (_c)->re, (_c)->im); mpq_inv((_d), (_c)->re); \
  mpq_mul((_c)->re, (_a)->re, (_b)->re); mpq_mul((_c)->im, (_a)->im, (_b)->im); mpq_add((_c)->re, (_c)->re, (_c)->im); mpq_mul((_c)->re, (_c)->re, (_d)); \
  mpq_mul((_c)->im, (_a)->im, (_b)->re); mpq_mul((_z)->im, (_a)->re, (_b)->im); mpq_sub((_z)->im, (_c)->im, (_z)->im); mpq_mul((_z)->im, (_z)->im, (_d)); \
  mpq_set((_z)->re, (_c)->re); clear_rational_number(_c); mpq_clear(_d); }

// reciprocate a complex number
#define recip_number(_a,_b) { complex_number _c; initialize_number(_c); \
  mpf_mul((_c)->re, (_b)->re, (_b)->re); mpf_mul((_c)->im, (_b)->im, (_b)->im); mpf_add((_c)->re, (_c)->re, (_c)->im); \
  mpf_ui_div((_c)->re, 1, (_c)->re); mpf_mul((_a)->re, (_b)->re, (_c)->re); mpf_neg((_c)->re, (_c)->re); \
  mpf_mul((_a)->im, (_b)->im, (_c)->re); clear_number(_c); }

#define recip_rational_number(_a,_b) { rational_complex_number _c; initialize_rational_number(_c); \
  mpq_mul((_c)->re, (_b)->re, (_b)->re); mpq_mul((_c)->im, (_b)->im, (_b)->im); mpq_add((_c)->re, (_c)->re, (_c)->im); \
  mpq_inv((_c)->re, (_c)->re); mpq_mul((_a)->re, (_b)->re, (_c)->re); mpq_neg((_c)->re, (_c)->re); \
  mpq_mul((_a)->im, (_b)->im, (_c)->re); clear_rational_number(_c); }

// multiply two complex numbers and add to number
#define sum_multiply(_z,_a,_b) { complex_number _c; mpf_t _d; initialize_number(_c); mpf_init(_d); \
  mpf_mul((_c)->re, (_a)->re, (_b)->re); mpf_mul((_c)->im, (_a)->im, (_b)->im); mpf_sub((_c)->re, (_c)->re, (_c)->im); \
  mpf_mul((_c)->im, (_a)->re, (_b)->im); mpf_mul((_d), (_a)->im, (_b)->re); mpf_add((_c)->im, (_c)->im, (_d)); \
  add(_z,_z,_c); clear_number(_c); mpf_clear(_d); }

#define sum_multiply_rational(_z,_a,_b) { rational_complex_number _c; mpq_t _d; initialize_rational_number(_c); mpq_init(_d); \
  mpq_mul((_c)->re, (_a)->re, (_b)->re); mpq_mul((_c)->im, (_a)->im, (_b)->im); mpq_sub((_c)->re, (_c)->re, (_c)->im); \
  mpq_mul((_c)->im, (_a)->re, (_b)->im); mpq_mul((_d), (_a)->im, (_b)->re); mpq_add((_c)->im, (_c)->im, (_d)); \
  add_rational(_z,_z,_c); clear_rational_number(_c); mpq_clear(_d); }

// multiply two complex numbers and subtract from number
#define subtract_multiply(_z,_a,_b) { complex_number _c; mpf_t _d; initialize_number(_c); mpf_init(_d); \
  mpf_mul((_c)->re, (_a)->re, (_b)->re); mpf_mul((_c)->im, (_a)->im, (_b)->im); mpf_sub((_c)->re, (_c)->re, (_c)->im); \
  mpf_mul((_c)->im, (_a)->re, (_b)->im); mpf_mul((_d), (_a)->im, (_b)->re); mpf_add((_c)->im, (_c)->im, (_d)); \
  subtract(_z,_z,_c); clear_number(_c); mpf_clear(_d); }

#define subtract_multiply_rational(_z,_a,_b) { rational_complex_number _c; mpq_t _d; initialize_rational_number(_c); mpq_init(_d); \
  mpq_mul((_c)->re, (_a)->re, (_b)->re); mpq_mul((_c)->im, (_a)->im, (_b)->im); mpq_sub((_c)->re, (_c)->re, (_c)->im); \
  mpq_mul((_c)->im, (_a)->re, (_b)->im); mpq_mul((_d), (_a)->im, (_b)->re); mpq_add((_c)->im, (_c)->im, (_d)); \
  subtract_rational(_z,_z,_c); clear_rational_number(_c); mpq_clear(_d); }

// Householder update: _r = _a - 2 * _b * conj(_c)
#define householder_rational_multiplication(_r, _a, _b, _c) { mpq_t _s,_t,_u; mpq_init(_s); mpq_init(_t); mpq_init(_u); \
  mpq_mul(_s, (_b)->re, (_c)->re); mpq_mul(_t, (_b)->im, (_c)->im); mpq_add(_s, _s, _t); mpq_add(_s, _s, _s); \
  mpq_mul(_t, (_b)->im, (_c)->re); mpq_mul(_u, (_b)->re, (_c)->im); mpq_sub(_t, _t, _u); mpq_add(_t, _t, _t); \
  mpq_sub((_r)->re, (_a)->re, _s); mpq_sub((_r)->im, (_a)->im, _t); mpq_clear(_s); mpq_clear(_t); mpq_clear(_u); }

// square a number
#define square(_z,_a) { complex_number _b; initialize_number(_b); \
  mpf_mul((_b)->re, (_a)->re, (_a)->re); mpf_mul((_b)->im, (_a)->im, (_a)->im); mpf_sub((_b)->re, (_b)->re, (_b)->im); \
  mpf_mul((_b)->im, (_a)->re, (_a)->im); mpf_add((_z)->im, (_b)->im, (_b)->im); mpf_set((_z)->re, (_b)->re); \
  clear_number(_b); }

#define square_rational(_z,_a) { rational_complex_number _b; initialize_rational_number(_b); \
  mpq_mul((_b)->re, (_a)->re, (_a)->re); mpq_mul((_b)->im, (_a)->im, (_a)->im); mpq_sub((_b)->re, (_b)->re, (_b)->im); \
  mpq_mul((_b)->im, (_a)->re, (_a)->im); mpq_add((_z)->im, (_b)->im, (_b)->im); mpq_set((_z)->re, (_b)->re); \
  clear_rational_number(_b); }

// exponential of a number
#define exp_number(_z,_a) { complex_number _b; initialize_number(_b); \
  mpfr_exp((_b)->im, (_a)->re, __gmp_default_rounding_mode); mpfr_cos((_b)->re, (_a)->im, __gmp_default_rounding_mode); \
  mpf_mul((_z)->re, (_b)->im, (_b)->re); mpfr_sin((_z)->im, (_a)->im, __gmp_default_rounding_mode);  \
  mpf_mul((_z)->im, (_b)->im, (_z)->im); clear_number(_b); }

// cosine of a number
#define cos_number(_z,_a) { complex_number _b,_c; initialize_number(_b); initialize_number(_c); \
  mpfr_sin_cos((_b)->re, (_b)->im, (_a)->re, __gmp_default_rounding_mode); \
  mpfr_sinh_cosh((_c)->re, (_c)->im, (_a)->im, __gmp_default_rounding_mode); \
  mpf_mul((_z)->re, (_c)->im, (_b)->im); mpf_mul((_z)->im, (_c)->re, (_b)->re); mpf_neg((_z)->im, (_z)->im); \
  clear_number(_b); clear_number(_c); }

// sine of a number
#define sin_number(_z,_a) { complex_number _b,_c; initialize_number(_b); initialize_number(_c); \
  mpfr_sin_cos((_b)->re, (_b)->im, (_a)->re, __gmp_default_rounding_mode); \
  mpfr_sinh_cosh((_c)->re, (_c)->im, (_a)->im, __gmp_default_rounding_mode); \
  mpf_mul((_z)->re, (_c)->im, (_b)->re); mpf_mul((_z)->im, (_c)->re, (_b)->im); \
  clear_number(_b); clear_number(_c); }

// hyperbolic cosine of a number
#define cosh_number(_z,_a) { complex_number _b,_c; initialize_number(_b); initialize_number(_c); \
  mpfr_sinh_cosh((_b)->re, (_b)->im, (_a)->re, __gmp_default_rounding_mode); \
  mpfr_sin_cos((_c)->re, (_c)->im, (_a)->im, __gmp_default_rounding_mode); \
  mpf_mul((_z)->re, (_b)->im, (_c)->im); mpf_mul((_z)->im, (_b)->re, (_c)->re); \
  clear_number(_b); clear_number(_c); }

// hyperbolic sine of a number
#define sinh_number(_z,_a) { complex_number _b,_c; initialize_number(_b); initialize_number(_c); \
  mpfr_sinh_cosh((_b)->re, (_b)->im, (_a)->re, __gmp_default_rounding_mode); \
  mpfr_sin_cos((_c)->re, (_c)->im, (_a)->im, __gmp_default_rounding_mode); \
  mpf_mul((_z)->re, (_c)->im, (_b)->re); mpf_mul((_z)->im, (_c)->re, (_b)->im); \
  clear_number(_b); clear_number(_c); }

// exponentiate using binary
#define exponentiate(_zexp,_aexp,_dexp) { \
  if (_dexp == 0) { set_one_number(_zexp); } \
  else if (_dexp == 1) { set_number(_zexp,_aexp); } \
  else if (_dexp == 2) { square(_zexp,_aexp); } \
  else if (_dexp > 2) { /* use binary decomposition of _dexp to compute (_aexp)^_dexp */ \
    int _i,_size; complex_number _cexp; initialize_number(_cexp); set_number(_cexp,_aexp); set_one_number(_zexp); \
    mpz_t _pow_int; mpz_init_set_ui(_pow_int, _dexp); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
    for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { multiply(_zexp,_zexp,_cexp); } square(_cexp,_cexp); } \
    clear_number(_cexp); mpz_clear(_pow_int); free(_base2); } \
  else { printf("\nERROR: Invalid power (%d).\n",_dexp); errExit(ERROR_CONFIGURATION); }}

#define exponentiate_rational(_zexp,_aexp,_dexp) { \
  if (_dexp == 0) { set_one_rational_number(_zexp); } \
  else if (_dexp == 1) { set_rational_number(_zexp,_aexp); } \
  else if (_dexp == 2) { square_rational(_zexp,_aexp); } \
  else if (_dexp > 2) { /* use binary decomposition of _dexp to compute (_aexp)^_dexp */ \
    int _i,_size; rational_complex_number _cexp; initialize_rational_number(_cexp); set_rational_number(_cexp,_aexp); set_one_rational_number(_zexp); \
    mpz_t _pow_int; mpz_init_set_ui(_pow_int, _dexp); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
    for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { multiply_rational(_zexp,_zexp,_cexp); } square_rational(_cexp,_cexp); } \
    clear_rational_number(_cexp); mpz_clear(_pow_int); free(_base2); } \
  else { printf("\nERROR: Invalid power (%d).\n",_dexp); errExit(ERROR_CONFIGURATION); }}

// compute norm^2 & norm of a complex number & vector
#define norm_sqr_number(_n,_a) { mpf_t _n1; mpf_init(_n1); mpf_mul(_n, (_a)->re, (_a)->re); mpf_mul(_n1, (_a)->im, (_a)->im); mpf_add(_n, _n, _n1); mpf_clear(_n1); }
#define norm_sqr_rational_number(_n,_a) { mpq_t _n1; mpq_init(_n1); mpq_mul(_n, (_a)->re, (_a)->re); mpq_mul(_n1, (_a)->im, (_a)->im); mpq_add(_n, _n, _n1); mpq_clear(_n1); }

#define norm_number(_n,_a) { norm_sqr_number(_n,_a); mpf_sqrt(_n,_n); }

#define norm_sqr_vector(_n,_v) { int _i, _size = (_v)->size; mpf_t _n2; mpf_init(_n2); mpf_set_ui(_n, 0); \
  for (_i = 0; _i < _size; _i++) { norm_sqr_number(_n2,(_v)->coord[_i]); mpf_add(_n,_n,_n2); } mpf_clear(_n2); }

#define norm_sqr_rational_vector(_n,_v) { int _i, _size = (_v)->size; mpq_t _n2; mpq_init(_n2); mpq_set_ui(_n, 0, 1); \
  for (_i = 0; _i < _size; _i++) { norm_sqr_rational_number(_n2,(_v)->coord[_i]); mpq_add(_n,_n,_n2); } mpq_clear(_n2); }

#define norm_vector(_n,_v) { norm_sqr_vector(_n,_v); mpf_sqrt(_n,_n); }

// compute ||x||_1^2 & ||x||_1 for a vector
#define norm_one_sqr_vector(_n,_v) { int _i, _size = (_v)->size; mpf_t _n2; mpf_init(_n2); mpf_set_ui(_n, 1); \
  for (_i = 0; _i < _size; _i++) { norm_sqr_number(_n2,(_v)->coord[_i]); mpf_add(_n,_n,_n2); } mpf_clear(_n2); }

#define norm_one_sqr_rational_vector(_n,_v) { int _i, _size = (_v)->size; mpq_t _n2; mpq_init(_n2); mpq_set_ui(_n, 1, 1); \
  for (_i = 0; _i < _size; _i++) { norm_sqr_rational_number(_n2,(_v)->coord[_i]); mpq_add(_n,_n,_n2); } mpq_clear(_n2); }

#define norm_one_vector(_n,_v) { norm_one_sqr_vector(_n,_v); mpf_sqrt(_n,_n); }

// compute ||x - \pi(x)||_2^2 and ||x - \pi(x)||_2 for a vector
#define norm_sqr_imag_vector(_n,_v) { int _i, _size = (_v)->size; mpf_t _n1; mpf_init(_n1); mpf_set_ui(_n, 0); \
  for (_i = 0; _i < _size; _i++) { mpf_mul(_n1, (_v)->coord[_i]->im, (_v)->coord[_i]->im); mpf_add(_n, _n, _n1); } mpf_clear(_n1); }

#define norm_sqr_imag_rational_vector(_n,_v) { int _i, _size = (_v)->size; mpq_t _n1; mpq_init(_n1); mpq_set_ui(_n, 0, 1); \
  for (_i = 0; _i < _size; _i++) { mpq_mul(_n1, (_v)->coord[_i]->im, (_v)->coord[_i]->im); mpq_add(_n, _n, _n1); } mpq_clear(_n1); }

#define norm_imag_vector(_n,_v) { norm_sqr_imag_vector(_n,_v); mpf_sqrt(_n,_n); }

#endif

