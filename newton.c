/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   newton.c: Compute the newton residual
*/

#include "alphaCertified.h"

int newton_iteration(complex_vector Nx, complex_number newtonRes_sqr, complex_matrix LU, int **rowswaps, polynomial_system *F, complex_vector x, mpf_t pivot_tol, mpf_t pivot_drop_tol, int eval_prec)
/***************************************************************\
* USAGE: compute Nx = x - JF^-1(x)*F(x) & ||JF^-1(x)*F(x)||_2^2 *
*  along with a LU decomposiont of JF(x)                        *
\***************************************************************/
{
  int i, retVal = 0;
  complex_vector f;
  complex_matrix J;

  // initialize f & J
  initialize_vector2(f, 0, eval_prec);
  initialize_matrix2(J, 0, 0, eval_prec);

  // compute f = F(x) & J = JF(x)
  eval_polynomial_system(f, J, F, x, eval_prec);

  // compute an LU decomposition of J
  retVal = LUdecomp(LU, rowswaps, J, pivot_tol, pivot_drop_tol, eval_prec);

  if (!retVal)
  { // successful LU decomposition - compute JF^-1(x)*F(x)
    LUsolve_vector(f, LU, *rowswaps, f, eval_prec);

    // setup Nx and newtonRes_sqr
    setPrec_vector(Nx, eval_prec);
    change_size_vector(Nx, f->size);
    set_zero_number(newtonRes_sqr);

    // compute x - JF^-1(x)*F(x) & ||JF^-1(x)*F(x)||_2^2
    for (i = 0; i < f->size; i++)
    {
      subtract(Nx->coord[i], x->coord[i], f->coord[i]);

      norm_sqr_number(newtonRes_sqr->im, f->coord[i]);
      mpf_add(newtonRes_sqr->re, newtonRes_sqr->re, newtonRes_sqr->im);
    }
    mpf_set_ui(newtonRes_sqr->im, 0);    
  }
  else
  { // determine if F(x) == 0
    norm_sqr_vector(newtonRes_sqr->im, f);
    if (mpf_cmp_ui(newtonRes_sqr->im, 0) == 0)
    { // update the error code since exact solution
      retVal = EXACT_SOLUTION_LU_ERROR;

      // setup Nx and newtonRes_sqr
      setPrec_vector(Nx, eval_prec);
      copy_vector(Nx, x);
      set_zero_number(newtonRes_sqr);
    }
    mpf_set_ui(newtonRes_sqr->im, 0);    
  }

  // clear memory
  clear_vector(f);
  clear_matrix(J);

  return retVal;
}

int newton_iteration_rational(rational_complex_vector Nx, rational_complex_number newtonRes_sqr, rational_complex_matrix LU, int **rowswaps, polynomial_system *F, rational_complex_vector x)
/***************************************************************\
* USAGE: compute Nx = x - JF^-1(x)*F(x) & ||JF^-1(x)*F(x)||_2^2 *
*  along with a LU decomposiont of JF(x)                        *
\***************************************************************/
{
  int i, retVal = 0;
  rational_complex_vector f;
  rational_complex_matrix J;

  // initialize f & J
  initialize_rational_vector(f, 0);
  initialize_rational_matrix(J, 0, 0);

  // compute f = F(x) & J = JF(x)
  eval_polynomial_system_rational(f, J, F, x);

  // compute an LU decomposition of J
  retVal = LUdecomp_rational(LU, rowswaps, J);

  if (!retVal)
  { // successful LU decomposition - compute JF^-1(x)*F(x)
    LUsolve_rational_vector(f, LU, *rowswaps, f);

    // setup Nx and newtonRes_sqr
    change_size_rational_vector(Nx, f->size);
    set_zero_rational_number(newtonRes_sqr);

    // compute x - JF^-1(x)*F(x) & ||JF^-1(x)*F(x)||_2^2
    for (i = 0; i < f->size; i++)
    {
      subtract_rational(Nx->coord[i], x->coord[i], f->coord[i]);

      norm_sqr_rational_number(newtonRes_sqr->im, f->coord[i]);
      mpq_add(newtonRes_sqr->re, newtonRes_sqr->re, newtonRes_sqr->im);
    }
    mpq_set_ui(newtonRes_sqr->im, 0, 1);    
  }
  else
  { // determine if F(x) == 0
    norm_sqr_rational_vector(newtonRes_sqr->im, f);
    if (mpq_cmp_ui(newtonRes_sqr->im, 0, 1) == 0)
    { // update the error code since exact solution
      retVal = EXACT_SOLUTION_LU_ERROR;

      // setup Nx and newtonRes_sqr
      copy_rational_vector(Nx, x);
      set_zero_rational_number(newtonRes_sqr);
    }
    mpq_set_ui(newtonRes_sqr->im, 0, 1);    
  }

  // clear memory
  clear_rational_vector(f);
  clear_rational_matrix(J);

  return retVal;
}


