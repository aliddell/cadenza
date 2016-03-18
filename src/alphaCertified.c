/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   alphaCertified.c: Main file for alphaCertified
*/

#include "alphaCertified.h"

int compute_beta(complex_number beta, polynomial_system *F, complex_vector x, int eval_prec)
/***************************************************************\
* USAGE: compute beta                                           *
\***************************************************************/
{
  int retVal = 0, *rowswaps = NULL;
  mpf_t pivot_tol, pivot_drop_tol;
  complex_vector Nx;
  complex_matrix LU;

  // set precision
  setPrec(eval_prec);

  // initialize
  mpf_init(pivot_tol); mpf_init(pivot_drop_tol);
  initialize_vector(Nx, 0);
  initialize_matrix(LU, 0, 0);

  // adjust precision on beta
  setPrec_number(beta, eval_prec);

  // determine pivot_tol & pivot_drop_tol
  determine_pivot_tolerances(pivot_tol, pivot_drop_tol, eval_prec);

  // compute a newton iteration and an LU decomposition of the Jacobian
  retVal = newton_iteration(Nx, beta, LU, &rowswaps, F, x, pivot_tol, pivot_drop_tol, eval_prec);

  if (retVal == 0)
  { // successful - compute beta
    mpf_sqrt(beta->re, beta->re);
    mpf_set_ui(beta->im, 0);
  }

  // clear memory
  free(rowswaps);
  mpf_clear(pivot_tol); mpf_clear(pivot_drop_tol);
  clear_vector(Nx);
  clear_matrix(LU);

  return retVal;
}

int compute_beta_sqr_rational(rational_complex_number beta_sqr, polynomial_system *F, rational_complex_vector x)
/***************************************************************\
* USAGE: compute beta^2                                         *
\***************************************************************/
{
  int retVal = 0, *rowswaps = NULL;
  rational_complex_vector Nx;
  rational_complex_matrix LU;

  // initialize
  initialize_rational_vector(Nx, 0);
  initialize_rational_matrix(LU, 0, 0);

  // compute a newton iteration and an LU decomposition of the Jacobian
  retVal = newton_iteration_rational(Nx, beta_sqr, LU, &rowswaps, F, x);

  // clear memory
  free(rowswaps);
  clear_rational_vector(Nx);
  clear_rational_matrix(LU);

  return retVal;
}

void compute_mu_from_LU(complex_number mu, mpf_t norm_x, complex_matrix LU, int *rowswaps, polynomial_system *F, complex_vector x, int eval_prec)
/***************************************************************\
* USAGE: essentially mu = max{1, ||JF^-1 * delta * ||F|| ||}    *
\***************************************************************/
{
  int i, rows = LU->rows, numPoly = F->numPolynomials;
  complex_number norm_F;
  complex_matrix A, B;

  // set precision
  setPrec(eval_prec);

  // initialize
  initialize_number(norm_F);
  initialize_matrix(A, rows, rows);
  initialize_matrix(B, rows, rows);

  // compute ||F||
  norm_polynomial_system(norm_F, F, eval_prec);

  // setup B to identity
  identity_matrix(B);

  // update the diagonal elements of B
  for (i = 0; i < numPoly; i++)
  {
    if (F->polynomials[i].degree > 1)
    { // multiply by sqrt(deg) * ||x||_1^{deg - 1}
      mpf_pow_ui(B->entry[i][i]->re, norm_x, F->polynomials[i].degree - 1);
      mpf_sqrt_ui(B->entry[i][i]->im, F->polynomials[i].degree);
      mpf_mul(B->entry[i][i]->re, B->entry[i][i]->re, B->entry[i][i]->im);
      mpf_set_ui(B->entry[i][i]->im, 0);
    }

    // multiply by ||F||
    multiply(B->entry[i][i], B->entry[i][i], norm_F);
  }

  // solve
  LUsolve(A, LU, rowswaps, B, eval_prec);

  // compute ||A||_F
  norm_frobenius(mu, A, eval_prec);

  // verify larger than 1
  if (mpf_cmp_ui(mu->re, 1) < 0)
  { // set to 1  
    mpf_set_ui(mu->re, 1);
  }

  // clear
  clear_number(norm_F);
  clear_matrix(A);
  clear_matrix(B);  

  return;
}

void exponentialBound(mpf_t bound, exponential *F, complex_vector vars)
/***************************************************************\
* USAGE: compute the bound associated with the exponential F    *
\***************************************************************/
{
  complex_number tempNum, tempNum2, beta, beta_sqr;

  // initialize
  initialize_number(tempNum);
  initialize_number(tempNum2);
  initialize_number(beta);
  initialize_number(beta_sqr);

  // convert beta to floating point
  convert_rational_number(beta, F->beta);

  // compute beta_sqr
  multiply(beta_sqr, beta, beta);

  if (F->expFunction == 'X')
  { // bound is max{|beta|, |beta^2 * exp(beta * x) / 2|}
    multiply(tempNum, beta, vars->coord[F->xIndex]);
    exp_number(tempNum, tempNum);
    multiply(tempNum, tempNum, beta_sqr);
    norm_number(bound, tempNum);
    mpf_div_ui(bound, bound, 2);

    // determine max
    norm_number(tempNum->re, beta);
    if (mpf_cmp(bound, tempNum->re) < 0)
      mpf_set(bound, tempNum->re);
  }
  else if (F->expFunction == 'C' || F->expFunction == 'S')
  { // bound is either max{|beta|,|beta^2 * sin(beta * x) / 2|,|beta^2 * cos(beta * x) / 2|}
    // OR max{|beta|,|beta^2 * sinh(beta * x) / 2|,|beta^2 * cosh(beta * x) / 2|}

    multiply(tempNum, beta, vars->coord[F->xIndex]);
    if (F->isHyperbolic)
    { // evaluate sinh & cosh
      cosh_number(tempNum2, tempNum);
      sinh_number(tempNum, tempNum);
    }
    else
    { // evaluate sin & cos
      cos_number(tempNum2, tempNum);
      sin_number(tempNum, tempNum);
    }
    multiply(tempNum, tempNum, beta_sqr);
    multiply(tempNum2, tempNum2, beta_sqr);

    // determine max for sin & cos or sinh & cosh
    norm_number(bound, tempNum);
    norm_number(tempNum->re, tempNum2);
    if (mpf_cmp(bound, tempNum->re) < 0)
      mpf_set(bound, tempNum->re);
    mpf_div_ui(bound, bound, 2);

    // determine max
    norm_number(tempNum->re, beta);
    if (mpf_cmp(bound, tempNum->re) < 0)
      mpf_set(bound, tempNum->re);
  }
  else
  { // error
    printf("ERROR: Invalid exponential function.\n\n");
    errExit(ERROR_INPUT_SYSTEM);
  }

  // clear
  clear_number(tempNum);
  clear_number(tempNum2);
  clear_number(beta);
  clear_number(beta_sqr);

  return;
}

void compute_gamma_from_LU(complex_number gamma, complex_matrix LU, int *rowswaps, polynomial_system *F, complex_vector x, int eval_prec)
/***************************************************************\
* USAGE: compute gamma                                          *
\***************************************************************/
{
  int i;
  mpf_t norm_x, polyPart, expPart;
  complex_number mu;

  // set precision
  setPrec(eval_prec);

  // initialize
  mpf_init(norm_x);
  mpf_init(polyPart);
  mpf_init(expPart);
  initialize_number(mu);

  // compute ||x||_1
  norm_one_vector(norm_x, x);

  // compute mu
  compute_mu_from_LU(mu, norm_x, LU, rowswaps, F, x, eval_prec);

  // compute the polynomial part: D^(3/2) / (2 * ||x||_1)
  mpf_set_ui(polyPart, F->maximumDegree);
  mpf_pow_ui(polyPart, polyPart, 3);
  mpf_sqrt(polyPart, polyPart);
  mpf_add(norm_x, norm_x, norm_x);
  mpf_div(polyPart, polyPart, norm_x);

  mpf_set_ui(expPart, 0);
  for (i = 0; i < F->numExponentials; i++)
  { // add on the bound for the ith exponential
    exponentialBound(norm_x, &F->exponentials[i], x);
    mpf_add(expPart, expPart, norm_x);
  }

  // compute gamma = mu * (polynomial bound + exponential bound)
  mpf_add(gamma->re, polyPart, expPart);
  mpf_mul(gamma->re, gamma->re, mu->re);
  mpf_set_ui(gamma->im, 0);

  // clear
  mpf_clear(norm_x);
  mpf_clear(polyPart);
  mpf_clear(expPart);
  clear_number(mu);

  return;
}

int compute_gamma(complex_number gamma, polynomial_system *F, complex_vector x, int eval_prec)
/***************************************************************\
* USAGE: compute gamma                                          *
\***************************************************************/
{
  int retVal = 0, *rowswaps = NULL;
  mpf_t pivot_tol, pivot_drop_tol;
  complex_vector f;
  complex_matrix J, LU;

  // set precision
  setPrec(eval_prec);

  // initialize
  mpf_init(pivot_tol); mpf_init(pivot_drop_tol);
  initialize_vector(f, 0);
  initialize_matrix(J, 0, 0);
  initialize_matrix(LU, 0, 0);

  // adjust precision on gamma
  setPrec_number(gamma, eval_prec);

  // evaluate the function and its Jacobian
  eval_polynomial_system(f, J, F, x, eval_prec);

  // determine pivot_tol & pivot_drop_tol
  determine_pivot_tolerances(pivot_tol, pivot_drop_tol, eval_prec);

  // compute an LU decomposition of the Jacobian
  retVal = LUdecomp(LU, &rowswaps, J, pivot_tol, pivot_drop_tol, eval_prec);

  if (retVal == 0)
  { // successful - compute gamma
    compute_gamma_from_LU(gamma, LU, rowswaps, F, x, eval_prec);
  }

  // clear memory
  free(rowswaps);
  mpf_clear(pivot_tol); mpf_clear(pivot_drop_tol);
  clear_vector(f);
  clear_matrix(J);
  clear_matrix(LU);

  return retVal;
}

void compute_mu_sqr_from_LU_rational(rational_complex_number mu_sqr, mpq_t norm_sqr_x, rational_complex_matrix LU, int *rowswaps, polynomial_system *F, rational_complex_vector x)
/***************************************************************\
* USAGE: compute mu^2 = max{1, ||f||^2 * ||Jv^-1 * delta||^2}   *
\***************************************************************/
{
  int i, j, deg, rows = LU->rows;
  mpq_t norm_sqr_F, norm_sqr_A, tempRat1, tempRat2;
  rational_complex_matrix A, B;

  // initialize
  mpq_init(norm_sqr_F);
  mpq_init(norm_sqr_A);
  mpq_init(tempRat1);
  mpq_init(tempRat2);
  initialize_rational_matrix(A, rows, rows);
  initialize_rational_matrix(B, rows, rows);

  // set ||F||^2
  mpq_set(norm_sqr_F, F->norm_sqr);  

  // setup B to identity
  identity_rational_matrix(B);

  // solve
  LUsolve_rational(A, LU, rowswaps, B);

  // compute the Frobenius norm over each column and multiply by the correct scaler
  mpq_set_ui(norm_sqr_A, 0, 1);  
  for (i = 0; i < rows; i++)
  { // compute norm of ith column
    mpq_set_ui(tempRat1, 0, 1);
    for (j = 0; j < rows; j++) 
    {
      norm_sqr_rational_number(tempRat2 , A->entry[j][i]);
      mpq_add(tempRat1, tempRat1, tempRat2);
    }
    // multiply by ith degree  
    mpq_set_ui(tempRat2, F->polynomials[i].degree, 1);
    mpq_mul(tempRat1, tempRat1, tempRat2);
    // compute ||x||_1^(2*(deg - 1))
    deg = F->polynomials[i].degree - 1;
    exponentiate_mpq(tempRat2, norm_sqr_x, deg);
    // multiply
    mpq_mul(tempRat1, tempRat1, tempRat2);

    // add on to norm_sqr_A
    mpq_add(norm_sqr_A, norm_sqr_A, tempRat1);
  }

  // compute ||f||^2 * ||A||_F^2
  set_zero_rational_number(mu_sqr);
  mpq_mul(mu_sqr->re, norm_sqr_F, norm_sqr_A);
 
  // verify larger than 1
  if (mpq_cmp_ui(mu_sqr->re, 1, 1) < 0)
  { // set to 1  
    mpq_set_ui(mu_sqr->re, 1, 1);
  }

  // clear
  mpq_clear(norm_sqr_F);
  mpq_clear(norm_sqr_A);  
  mpq_clear(tempRat1);
  mpq_clear(tempRat2);
  clear_rational_matrix(A);
  clear_rational_matrix(B);  

  return;
}

void compute_gamma_sqr_from_LU_rational(rational_complex_number gamma_sqr, rational_complex_matrix LU, int *rowswaps, polynomial_system *F, rational_complex_vector x)
/***************************************************************\
* USAGE: compute gamma^2                                        *
\***************************************************************/
{
  mpq_t norm_sqr_x, tempRat1, tempRat2;

  // initialize
  mpq_init(norm_sqr_x);
  mpq_init(tempRat1);
  mpq_init(tempRat2);

  // compute ||x||_1^2
  norm_one_sqr_rational_vector(norm_sqr_x, x);

  // compute mu^2
  compute_mu_sqr_from_LU_rational(gamma_sqr, norm_sqr_x, LU, rowswaps, F, x);

  // compute gamma^2 = mu^2 * D^3 / (4 * ||x||_1^2)
  mpq_div(gamma_sqr->re, gamma_sqr->re, norm_sqr_x);
  mpq_set_ui(tempRat1, 1, 4); 
  mpq_set_ui(tempRat2, F->maximumDegree, 1);
  mpq_mul(tempRat1, tempRat1, tempRat2);
  mpq_mul(tempRat1, tempRat1, tempRat2);
  mpq_mul(tempRat1, tempRat1, tempRat2);
  mpq_mul(gamma_sqr->re, gamma_sqr->re, tempRat1);
  mpq_set_ui(gamma_sqr->im, 0, 1);

  // clear
  mpq_clear(norm_sqr_x);
  mpq_clear(tempRat1);
  mpq_clear(tempRat2);

  return;
}

int compute_gamma_sqr(rational_complex_number gamma_sqr, polynomial_system *F, rational_complex_vector x)
/***************************************************************\
* USAGE: compute gamma^2                                        *
\***************************************************************/
{
  int retVal = 0, *rowswaps = NULL;
  rational_complex_vector f;
  rational_complex_matrix J, LU;

  // initialize
  initialize_rational_vector(f, 0);
  initialize_rational_matrix(J, 0, 0);
  initialize_rational_matrix(LU, 0, 0);

  // evaluate the function and its Jacobian
  eval_polynomial_system_rational(f, J, F, x);

  // compute an LU decomposition of the Jacobian
  retVal = LUdecomp_rational(LU, &rowswaps, J);

  if (retVal == 0)
  { // successful - compute gamma
    compute_gamma_sqr_from_LU_rational(gamma_sqr, LU, rowswaps, F, x);
  }

  // clear memory
  free(rowswaps);
  clear_rational_vector(f);
  clear_rational_matrix(J);
  clear_rational_matrix(LU);

  return retVal;
}

int compute_alpha_beta_gamma(complex_vector newX, complex_number alpha, complex_number beta, complex_number gamma, polynomial_system *F, complex_vector x, int eval_prec)
/***************************************************************\
* USAGE: compute (bound) alpha, beta, and gamma                 *
\***************************************************************/
{
  int retVal = 0, *rowswaps = NULL;
  mpf_t pivot_tol, pivot_drop_tol;
  complex_number newtonRes_sqr;
  complex_matrix LU;

  // set precision
  setPrec(eval_prec);

  // initialize
  mpf_init(pivot_tol); mpf_init(pivot_drop_tol);
  initialize_number(newtonRes_sqr);
  initialize_matrix(LU, 0, 0);

  // adjust precision on alpha, beta & gamma
  setPrec_vector(newX, eval_prec);
  setPrec_number(alpha, eval_prec);
  setPrec_number(beta, eval_prec);
  setPrec_number(gamma, eval_prec);

  // determine pivot_tol & pivot_drop_tol
  determine_pivot_tolerances(pivot_tol, pivot_drop_tol, eval_prec);

  // compute a newton iteration and an LU decomposition of the Jacobian
  retVal = newton_iteration(newX, newtonRes_sqr, LU, &rowswaps, F, x, pivot_tol, pivot_drop_tol, eval_prec);

  if (retVal == 0)
  { // successful

    // compute beta
    mpf_sqrt(beta->re, newtonRes_sqr->re);
    mpf_set_ui(beta->im, 0);

    // compute gamma
    compute_gamma_from_LU(gamma, LU, rowswaps, F, x, eval_prec);

    // compute alpha
    multiply(alpha, beta, gamma);
  }
  else if (retVal == EXACT_SOLUTION_LU_ERROR)
  { // solution with LU decomp error

    // alpha = beta = 0
    set_zero_number(alpha);
    set_zero_number(beta);

    // gamma = inf
    mpfr_set_inf(gamma->re, 1);  
    mpf_set_ui(gamma->im, 0);
  }
  else
  { // LU decomp error that is not a solution

    // alpha = beta = gamma = inf
    mpfr_set_inf(alpha->re, 1);
    mpf_set_ui(alpha->im, 0);
    mpfr_set_inf(beta->re, 1);
    mpf_set_ui(beta->im, 0);
    mpfr_set_inf(gamma->re, 1);
    mpf_set_ui(gamma->im, 0);
  }

  // clear memory
  free(rowswaps);
  mpf_clear(pivot_tol); mpf_clear(pivot_drop_tol);
  clear_number(newtonRes_sqr); 
  clear_matrix(LU);

  return retVal;
}

int compute_alpha_beta_gamma_sqr_rational(rational_complex_vector newX, rational_complex_number alpha_sqr, rational_complex_number beta_sqr, rational_complex_number gamma_sqr, polynomial_system *F, rational_complex_vector x)
/***************************************************************\
* USAGE: compute (bound) alpha^2, beta^2, and gamma^2           *
\***************************************************************/
{
  int retVal = 0, *rowswaps = NULL;
  rational_complex_matrix LU;

  // initialize
  initialize_rational_matrix(LU, 0, 0);

  // compute a newton iteration and an LU decomposition of the Jacobian
  retVal = newton_iteration_rational(newX, beta_sqr, LU, &rowswaps, F, x);

  if (retVal == 0)
  { // successful

    // compute gamma_sqr
    compute_gamma_sqr_from_LU_rational(gamma_sqr, LU, rowswaps, F, x);

    // compute alpha_sqr
    multiply_rational(alpha_sqr, beta_sqr, gamma_sqr);
  }
  else if (retVal == EXACT_SOLUTION_LU_ERROR)
  { // solution with LU decomp error

    // alpha = beta = 0
    set_zero_rational_number(alpha_sqr);
    set_zero_rational_number(beta_sqr);

    // gamma = inf
    mpq_set_ui(gamma_sqr->re, 1, 0);
    mpq_set_ui(gamma_sqr->im, 0, 1);
  }
  else
  { // LU decomp error that is not a solution

    // alpha = beta = gamma = inf
    mpq_set_ui(alpha_sqr->re, 1, 0);
    mpq_set_ui(alpha_sqr->im, 0, 1);
    mpq_set_ui(beta_sqr->re, 1, 0);
    mpq_set_ui(beta_sqr->im, 0, 1);
    mpq_set_ui(gamma_sqr->re, 1, 0);
    mpq_set_ui(gamma_sqr->im, 0, 1);
  }

  // clear memory
  free(rowswaps);
  clear_rational_matrix(LU);

  return retVal;
}






