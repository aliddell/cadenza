/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   norm.c: Compute the norm of a polynomial, polynomial system, and matrix
*/

#include "alphaCertified.h"

void norm_sqr_polynomial(mpq_t norm_sqr, polynomial *F)
/***************************************************************\
* USAGE: compute ||F||^2                                        *
\***************************************************************/
{
  int i, j, degree_term, numTerms = F->numTerms, numVars = F->numVariables, degree = F->degree;
  rational_complex_number abs_coeff_sqr;
  mpz_t degree_factorial, product, tempInt;

  // initialize memory
  initialize_rational_number(abs_coeff_sqr);
  mpz_init(degree_factorial);
  mpz_init(product);
  mpz_init(tempInt);

  // setup specific values
  mpq_set_ui(norm_sqr, 0, 1);
  mpz_fac_ui(degree_factorial, degree);

  for (i = 0; i < numTerms; i++)
  { // compute the normalizing constant
    degree_term = 0;
    mpz_set_ui(product, 1);
    for (j = 0; j < numVars; j++)
    {
      degree_term += F->exponents[i][j];
      if (F->exponents[i][j] > 1)
      {
        mpz_fac_ui(tempInt, F->exponents[i][j]);
        mpz_mul(product, product, tempInt);      
      }
    }
    mpz_fac_ui(tempInt, degree - degree_term);
    mpz_mul(product, product, tempInt);  

    // compute ||coeff||^2
    mpq_mul(abs_coeff_sqr->re, F->coeff[i]->re, F->coeff[i]->re);
    mpq_mul(abs_coeff_sqr->im, F->coeff[i]->im, F->coeff[i]->im);
    mpq_add(abs_coeff_sqr->re, abs_coeff_sqr->re, abs_coeff_sqr->im);

    // setup product/degree_factorial
    mpz_set(mpq_numref(abs_coeff_sqr->im), product);
    mpz_set(mpq_denref(abs_coeff_sqr->im), degree_factorial);
    mpq_canonicalize(abs_coeff_sqr->im);

    // multiply ||coeff||^2 * product / degree_factorial and then add on
    mpq_mul(abs_coeff_sqr->re, abs_coeff_sqr->re, abs_coeff_sqr->im);
    mpq_add(norm_sqr, norm_sqr, abs_coeff_sqr->re);
  }

  // clear memory
  clear_rational_number(abs_coeff_sqr);
  mpz_clear(degree_factorial);
  mpz_clear(product);
  mpz_clear(tempInt);

  return;
}

void norm_polynomial_system(complex_number norm, polynomial_system *F, int eval_prec)
/***************************************************************\
* USAGE: compute ||F|| in the given precision                   *
\***************************************************************/
{
  // set the default precision
  setPrec(eval_prec);

  // set precision on norm & initialize to 0
  setPrec_number(norm, eval_prec);
  set_zero_number(norm);

  // convert to floating point and compute sqr
  mpf_set_q(norm->re, F->norm_sqr);
  mpf_sqrt(norm->re, norm->re);  

  return;
}

void norm_frobenius(complex_number norm, complex_matrix A, int eval_prec)
/***************************************************************\
* USAGE: compute ||A||_F in the given precision                 *
\***************************************************************/
{
  int i, j, r = A->rows, c = A->cols;

  // set the default precision
  setPrec(eval_prec);

  // set precision on norm & initialize to 0
  setPrec_number(norm, eval_prec);
  set_zero_number(norm);
  
  for (i = 0; i < r; i++)
    for (j = 0; j < c; j++)
    {
      norm_sqr_number(norm->im, A->entry[i][j]);
      mpf_add(norm->re, norm->re, norm->im);
    }

  // compute sqrt of norm
  mpf_sqrt(norm->re, norm->re);
  mpf_set_ui(norm->im, 0);

  return;
}

void norm_sqr_frobenius_rational(rational_complex_number norm, rational_complex_matrix A)
/***************************************************************\
* USAGE: compute ||A||_F^2                                      *
\***************************************************************/
{
  int i, j, r = A->rows, c = A->cols;

  // initialize to 0
  set_zero_rational_number(norm);
  
  for (i = 0; i < r; i++)
    for (j = 0; j < c; j++)
    {
      norm_sqr_rational_number(norm->im, A->entry[i][j]);
      mpq_add(norm->re, norm->re, norm->im);
    }

  // set norm->im to 0
  mpq_set_ui(norm->im, 0, 1);

  return;
}


