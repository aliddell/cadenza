/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   sqrt.c: Compute a certified rational upper bound of a square root 
*/

#include "alphaCertified_basic.h"

void refine_sqrt_upper(mpq_t sqrt_x, mpq_t x, int digits);

int sqrt_upper(mpq_t sqrt_x, mpq_t x, int digits)
/***************************************************************\
* USAGE: compute a certified rational upper bound on the sqrt(x)*
*    that is within 10^-digits of the true value 0-good,1-error *
\***************************************************************/
{
  int rV = 0, sign = mpq_sgn(x);

  // verify digits >= 0
  if (digits < 0)
  { // error
    printf("ERROR: The number of digits must be nonnegative!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // check the sign of x
  if (sign == -1)
  { // x is negative - return error
    rV = 1;
  }
  else if (sign == 0)
  { // x is zero - sqrt(x) = 0
    mpq_set_ui(sqrt_x, 0, 1);
    rV = 0;
  }
  else
  { // x is positive - compute sqrt(x)
    int prec;
    mpfr_t sqrt_x_float;
    mpf_t tempMPF;
    mpq_t sqrt_x_temp, tempMPQ;

    // determine the precision to use to get an approximation
    if (digits <= 16)
      prec = 64;
    else
    { // determine the precision needed
      prec = (int) ceil((digits + 0.5) / (32 * log10(2.0)));
      if (prec <= 2)
        prec = 64;
      else
        prec *= 32;
    }

    // initialize
    mpfr_init2(sqrt_x_float, prec);
    mpf_init2(tempMPF, prec);
    mpq_init(sqrt_x_temp);
    mpq_init(tempMPQ);
    
    // approximate sqrt(x) - round up!
    mpfr_set_q(sqrt_x_float, x, GMP_RNDU);
    mpfr_sqrt(sqrt_x_float, sqrt_x_float, GMP_RNDU);

    // convert to rational
    mpfr_get_f(tempMPF, sqrt_x_float, GMP_RNDU);
    mpq_set_f(sqrt_x_temp, tempMPF);

    // verify that sqrt_x_temp is an upper bound 
    mpq_mul(tempMPQ, sqrt_x_temp, sqrt_x_temp);
    if (mpq_cmp(tempMPQ, x) >= 0)
    { // we have an upper bound - refine & certify it
      refine_sqrt_upper(sqrt_x_temp, x, digits);
      // copy to sqrt_x
      mpq_set(sqrt_x, sqrt_x_temp);
      rV = 0;
    }
    else
    { // try again with 2*x (Newton iterations still converge quadratically!)
      mpq_add(tempMPQ, x, x);
      mpfr_set_q(sqrt_x_float, tempMPQ, GMP_RNDU);
      mpfr_sqrt(sqrt_x_float, sqrt_x_float, GMP_RNDU);

      // convert to rational
      mpfr_get_f(tempMPF, sqrt_x_float, GMP_RNDU);
      mpq_set_f(sqrt_x_temp, tempMPF);

      // verify that sqrt_x_temp is an upper bound 
      mpq_mul(tempMPQ, sqrt_x_temp, sqrt_x_temp);
      if (mpq_cmp(tempMPQ, x) >= 0)
      { // we have an upper bound - refine & certify it
        refine_sqrt_upper(sqrt_x_temp, x, digits);
        // copy to sqrt_x
        mpq_set(sqrt_x, sqrt_x_temp);
        rV = 0;
      }
      else
      { // take any upper bound
        if (mpq_cmp_ui(x, 1, 1) <= 0)
        { // 1 is an upper bound for sqrt(x)
          mpq_set_ui(sqrt_x_temp, 1, 1);
        }
        else
        { // x is an upper bound for sqrt(x)
          mpq_set(sqrt_x_temp, x);
        }
        // we have an upper bound - refine & certify it
        refine_sqrt_upper(sqrt_x_temp, x, digits);
        // copy to sqrt_x
        mpq_set(sqrt_x, sqrt_x_temp);
        rV = 0;
      }
    }

    // clear
    mpfr_clear(sqrt_x_float);
    mpf_clear(tempMPF);
    mpq_clear(sqrt_x_temp);
    mpq_clear(tempMPQ);
  }

  return rV;
}

void refine_sqrt_upper(mpq_t sqrt_x, mpq_t x, int digits)
/***************************************************************\
* USAGE: refine a rational upper bound on the sqrt(x) to the    *
*    requested tolerance: 10^-digits                            *
\***************************************************************/
{ // assume (sqrt_x)^2 >= x > 0 && sqrt_x > 0 && digits >= 0
  mpq_t new_sqrt_x, beta, tol;

  // initialize
  mpq_init(new_sqrt_x);
  mpq_init(beta);
  mpq_init(tol);

  // setup tol = (1/10)^digits
  mpq_set_ui(tol, 1, 10);
  exponentiate_mpq(tol, tol, digits);

  // loop until we have refined to the tol
  do
  { // compute new_sqrt_x & beta 
    mpq_mul(beta, sqrt_x, sqrt_x);       // (sqrt_x)^2
    mpq_sub(beta, beta, x);              // (sqrt_x)^2 - x

    mpq_add(new_sqrt_x, sqrt_x, sqrt_x); // 2*sqrt_x 
    mpq_div(beta, beta, new_sqrt_x);     // ((sqrt_x)^2 - x)/(2*sqrt_x)
    mpq_sub(new_sqrt_x, sqrt_x, beta);   // sqrt_x - ((sqrt_x)^2 - x)/(2*sqrt_x)

    // determine if 2*beta <= tol
    mpq_add(beta, beta, beta);
    if (mpq_cmp(beta, tol) <= 0)
    { // success!!
      break;
    }
    else
    { // update sqrt_x & try again
      mpq_set(sqrt_x, new_sqrt_x);
    }
  } while (1);
  
  // clear
  mpq_clear(new_sqrt_x);
  mpq_clear(beta);
  mpq_clear(tol);

  return;
}
