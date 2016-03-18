/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   refine.c: Refine the points
*/

#include "alphaCertified.h"

int need_to_refine(complex_number beta, int curr_prec, int digits)
/***************************************************************\
* USAGE: determine if we need to continue to refine             *
\***************************************************************/
{
  int rV = 0, precDigits = (int) floor(curr_prec * log10(2.0) - 2.5);

  if (precDigits <= digits)
  { // the current precision is not large enough to handle the correct number of digits!
    rV = 1;
  }
  else
  { // see if beta is small enough
    mpf_t tempMPF;
    // tempMPF = (2*10^digits)^(-1).  If beta <= tempMPF, then distance to solution is <= 2*beta <= (10^digits)^(-1)
    mpf_init2(tempMPF, curr_prec);
    mpf_set_ui(tempMPF, 10);
    mpf_pow_ui(tempMPF, tempMPF, digits);
    mpf_mul_ui(tempMPF, tempMPF, 2);
    mpf_ui_div(tempMPF, 1, tempMPF);

    if (mpf_cmp(beta->re, tempMPF) <= 0)
    { // accuracy is reached
      rV = 0;
    }
    else
    { // accuracy is not reached
      rV = 1;
    }

    mpf_clear(tempMPF);
  }

  return rV;
}

void refine_points(int numPoints, point_struct *Points, polynomial_system *F, int eval_prec, int digits)
/***************************************************************\
* USAGE: refine the points to the requested number of digits    *
\***************************************************************/
{
  int i, rV, curr_prec;

  // see if refinement is requested
  if (digits <= 0)
     return;

  // loop over the points and refine the ones that are approximate solutions
  for (i = 0; i < numPoints; i++)
    if (Points[i].isApproxSoln)
    { // refine this one
      curr_prec = eval_prec;

      // perform a newton iteration
      rV = compute_alpha_beta_gamma(Points[i].Nx, Points[i].alpha, Points[i].beta, Points[i].gamma, F, Points[i].x, curr_prec);

      // check for errors - should not occur!!
      if (rV == ERROR_LU_DECOMP)
      {
        printf("ERROR: Invalid LU decomposition!\n");
        errExit(ERROR_CONFIGURATION);
      }

      // loop until we have enough correct digits
      while (need_to_refine(Points[i].beta, curr_prec, digits))
      { // double the precision
        curr_prec *= 2;
        // set the precision and copy the new point
        setPrec_vector(Points[i].x, curr_prec);
        copy_vector(Points[i].x, Points[i].Nx);

        // update norm
        mpf_set_prec(Points[i].norm_x, curr_prec);
        norm_vector(Points[i].norm_x, Points[i].x);
  
        // perform a newton iteration
        rV = compute_alpha_beta_gamma(Points[i].Nx, Points[i].alpha, Points[i].beta, Points[i].gamma, F, Points[i].x, curr_prec);

        // check for errors - should not occur!!
        if (rV == ERROR_LU_DECOMP)
        {
          printf("ERROR: Invalid LU decomposition!\n");
          errExit(ERROR_CONFIGURATION);
        }
      }

      // set precision back
      setPrec(eval_prec);  
    }

  return;
}

int need_to_refine_rational(rational_complex_number beta_sqr, int digits)
/***************************************************************\
* USAGE: determine if we need to continue to refine             *
\***************************************************************/
{
  int rV = 0;
  mpq_t tempRat;

  mpq_init(tempRat);

  // setup the string "1/[4*10^(2*digits)]"
  int i, size = 4 + 2*digits;
  char *str = errMalloc(size * sizeof(char));
  str[0] = '1'; str[1] = '/'; str[2] = '4';
  for (i = 0; i < 2*digits; i++)
    str[3+i] = '0';
  str[3+2*digits] = '\0';
  
  // setup tempRat
  mpq_set_str(tempRat, str, 10);

  // tempRat = (2*10^digits)^(-2).  If beta^2 <= tempRat, then distance to solution is <= 2*beta <= (10^digits)^(-1)
  if (mpq_cmp(beta_sqr->re, tempRat) <= 0)
  { // accuracy is reached
    rV = 0;
  }
  else
  { // accuracy is not reached
    rV = 1;
  }

  mpq_clear(tempRat);
  free(str);

  return rV;
}

void refine_points_rational(int numPoints, rational_point_struct *Points, polynomial_system *F, int digits)
/***************************************************************\
* USAGE: refine the points to the requested number of digits    *
\***************************************************************/
{
  int i, rV;

  // see if refinement is requested
  if (digits <= 0)
     return;

  // loop over the points and refine the ones that are approximate solutions
  for (i = 0; i < numPoints; i++)
    if (Points[i].isApproxSoln)
    { // perform a newton iteration
      rV = compute_alpha_beta_gamma_sqr_rational(Points[i].Nx, Points[i].alpha_sqr, Points[i].beta_sqr, Points[i].gamma_sqr, F, Points[i].x);

      // check for errors - should not occur!!
      if (rV == ERROR_LU_DECOMP)
      {
        printf("ERROR: Invalid LU decomposition!\n");
        errExit(ERROR_CONFIGURATION);
      }

      // loop until we have enough correct digits
      while (need_to_refine_rational(Points[i].beta_sqr, digits))
      { // copy the new point
        copy_rational_vector(Points[i].x, Points[i].Nx);

        // update norm
        norm_sqr_rational_vector(Points[i].norm_sqr_x, Points[i].x);
  
        // perform a newton iteration
        rV = compute_alpha_beta_gamma_sqr_rational(Points[i].Nx, Points[i].alpha_sqr, Points[i].beta_sqr, Points[i].gamma_sqr, F, Points[i].x);

        // check for errors - should not occur!!
        if (rV == ERROR_LU_DECOMP)
        {
          printf("ERROR: Invalid LU decomposition!\n");
          errExit(ERROR_CONFIGURATION);
        }
      }
    }

  return;
}

