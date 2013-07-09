/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   eval.c: Evaluates a polynomial system and its Jacobian
*/

#include "alphaCertified.h"

void eval_exponents_derivative(complex_number y, complex_vector dy, complex_vector vars, int *exp, int numVars)
/***************************************************************\
* USAGE: compute y = vars^exp & dy = d(vars^exp)                *
\***************************************************************/
{
  int i, j, *tempExp = (int *)errMalloc(numVars * sizeof(int));
  complex_number tempNum, tempY;
  complex_vector tempVec;

  // initialize tempY, tempNum, & tempVec
  initialize_number(tempNum);
  initialize_number(tempY);
  set_one_number(tempY);
  initialize_vector(tempVec, numVars);
  set_one_vector(tempVec);

  // initialize y to 1
  set_one_number(y);

  // subtract one off of every positive exponent, compute x^tempExp & initialize dy
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    { // compute values
      tempExp[i] = exp[i] - 1;
      exponentiate(tempNum, vars->coord[i], tempExp[i]);
      multiply(tempY, tempY, tempNum);
      set_one_number(dy->coord[i]);

      for (j = 0; j < numVars; j++)
        if (j != i && exp[j] > 0)
          multiply(tempVec->coord[j], tempVec->coord[j], vars->coord[i]);
    }
    else
    { // set to 0
      tempExp[i] = 0;
      set_zero_number(dy->coord[i]);
    }

  // compute dy
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    { // multiply
      multiply(dy->coord[i], tempY, tempVec->coord[i]);
      multiply_positive_int(dy->coord[i], dy->coord[i], exp[i]);
    }

  // compute y (already set to 1 if exp == 0)
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    {
      multiply(y, tempY, tempVec->coord[i]);
      multiply(y, y, vars->coord[i]);
    }

  // clear tempY, tempNum & tempVec
  clear_number(tempNum);
  clear_number(tempY);
  clear_vector(tempVec);

  free(tempExp);
  return;
}

void eval_polynomial(complex_number func, complex_number *jac, polynomial *F, complex_vector vars)
/***************************************************************\
* USAGE: compute F(vars) & JF(vars) where F is a polynomial     *
\***************************************************************/
{
  int i, j, numVars = F->numVariables, numTerms = F->numTerms;
  complex_number coeff, tempNum;
  complex_vector tempVec;

  // error checking
  if (numVars != vars->size)
  {
    printf("ERROR: Incorrect number of variables!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // initialize coeff, tempNum & tempVec
  initialize_number(coeff);
  initialize_number(tempNum);
  initialize_vector(tempVec, numVars);

  // initialize func & jac to all 0
  set_zero_number(func);
  for (i = 0; i < numVars; i++)
    set_zero_number(jac[i]);

  // loop over the number of terms
  for (i = 0; i < numTerms; i++)
  { // compute tempNum = x^exp[i] as well as tempVec = [d/dxj(x^exp[i])]
    eval_exponents_derivative(tempNum, tempVec, vars, F->exponents[i], numVars);

    // setup coeff (convert rational number to floating point number)
    convert_rational_number(coeff, F->coeff[i]);

    // multiply by coefficient and add to func
    sum_multiply(func, tempNum, coeff);

    // multiply by coefficient and add to jac
    for (j = 0; j < numVars; j++)
      sum_multiply(jac[j], tempVec->coord[j], coeff);
  }
 
  // clear coeff, tempNum & tempVec
  clear_number(coeff);
  clear_number(tempNum);
  clear_vector(tempVec);

  return;
}

void eval_exponential(complex_number func, complex_number *jac, exponential *F, complex_vector vars)
/***************************************************************\
* USAGE: compute F(vars) & JF(vars) where F is a exponential    *
\***************************************************************/
{
  complex_number beta, tempNum;

  // initialize beta, tempNum & tempVec
  initialize_number(beta);
  initialize_number(tempNum);

  // set jac[yIndex] = 1
  set_one_number(jac[F->yIndex]);

  if (F->expFunction == 'X')
  { // F = y - exp(beta * x)
    convert_rational_number(beta, F->beta);
    multiply(func, beta, vars->coord[F->xIndex]);
    exp_number(func, func);

    // setup jax[xIndex] = - beta * exp(beta * x)
    multiply(jac[F->xIndex], func, beta);
    negate(jac[F->xIndex], jac[F->xIndex]);

    subtract(func, vars->coord[F->yIndex], func);
  }
  else if (F->expFunction == 'S')
  { // F = y - sin(beta * x) or y - sinh(beta * x)
    convert_rational_number(beta, F->beta);
    multiply(func, beta, vars->coord[F->xIndex]);

    // setup F & jax[xIndex] = - beta * cos(beta * x) or beta * cosh(beta * x)
    if (F->isHyperbolic)
    {
      cosh_number(jac[F->xIndex], func);
      sinh_number(func, func);
    }
    else
    {
      cos_number(jac[F->xIndex], func);
      sin_number(func, func);
      negate(jac[F->xIndex], jac[F->xIndex]);
    }
    // complete F
    subtract(func, vars->coord[F->yIndex], func);

    // complete JF
    multiply(jac[F->xIndex], beta, jac[F->xIndex]); 
  }
  else if (F->expFunction == 'C')
  { // F = y - cos(beta * x) or y - cosh(beta * x)
    convert_rational_number(beta, F->beta);
    multiply(func, beta, vars->coord[F->xIndex]);

    // setup F & jax[xIndex] = beta * sin(beta * x) or beta * sinh(beta * x)
    if (F->isHyperbolic)
    {
      sinh_number(jac[F->xIndex], func);
      cosh_number(func, func);
    }
    else
    {
      sin_number(jac[F->xIndex], func);
      cos_number(func, func);
    }
    // complete F
    subtract(func, vars->coord[F->yIndex], func);

    // complete JF
    multiply(jac[F->xIndex], beta, jac[F->xIndex]);
  }
  else
  { // error
    printf("ERROR: Invalid exponential function.\n\n");
    errExit(ERROR_INPUT_SYSTEM);
  }
 
  // clear beta & tempNum
  clear_number(beta);
  clear_number(tempNum);

  return;
}

void eval_polynomial_system(complex_vector func, complex_matrix jac, polynomial_system *F, complex_vector vars, int eval_prec)
/***************************************************************\
* USAGE: compute F(vars) & JF(vars) where F is a poly system    *
\***************************************************************/
{
  int i, numVars = F->numVariables, numPoly = F->numPolynomials, numExp = F->numExponentials;
  int numFuncs = numPoly + numExp;

  // error checking
  if (numVars != vars->size)
  {
    printf("ERROR: Incorrect number of variables!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // set the default precision
  setPrec(eval_prec);

  // set precision on func & jac
  setPrec_vector(func, eval_prec);
  setPrec_matrix(jac, eval_prec);  

  // setup func & jac to the correct size
  change_size_vector(func, numVars);
  change_size_matrix(jac, numVars, numFuncs);

  // loop over the polynomials - compute the function and jacobian
  for (i = 0; i < numPoly; i++)
    eval_polynomial(func->coord[i], jac->entry[i], &F->polynomials[i], vars);

  // loop over the exponentials - compute the function and jacobian
  for (i = 0; i < numExp; i++)
    eval_exponential(func->coord[i+numPoly], jac->entry[i+numPoly], &F->exponentials[i], vars);

  return;
}

void eval_exponents_derivative_rational(rational_complex_number y, rational_complex_vector dy, rational_complex_vector vars, int *exp, int numVars)
/***************************************************************\
* USAGE: compute y = vars^exp & dy = d(vars^exp)                *
\***************************************************************/
{
  int i, j, *tempExp = (int *)errMalloc(numVars * sizeof(int));
  rational_complex_number tempNum, tempY;
  rational_complex_vector tempVec;

  // initialize tempY, tempNum, & tempVec
  initialize_rational_number(tempNum);
  initialize_rational_number(tempY);
  set_one_rational_number(tempY);
  initialize_rational_vector(tempVec, numVars);
  set_one_rational_vector(tempVec);

  // initialize y to 1
  set_one_rational_number(y);

  // subtract one off of every positive exponent, compute x^tempExp & initialize dy
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    { // compute values
      tempExp[i] = exp[i] - 1;
      exponentiate_rational(tempNum, vars->coord[i], tempExp[i]);
      multiply_rational(tempY, tempY, tempNum);
      set_one_rational_number(dy->coord[i]);

      for (j = 0; j < numVars; j++)
        if (j != i && exp[j] > 0)
          multiply_rational(tempVec->coord[j], tempVec->coord[j], vars->coord[i]);
    }
    else
    { // set to 0
      tempExp[i] = 0;
      set_zero_rational_number(dy->coord[i]);
    }

  // compute dy
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    { // multiply
      multiply_rational(dy->coord[i], tempY, tempVec->coord[i]);
      multiply_rational_positive_int(dy->coord[i], dy->coord[i], exp[i]);
    }

  // compute y (already set to 1 if exp == 0)
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    {
      multiply_rational(y, tempY, tempVec->coord[i]);
      multiply_rational(y, y, vars->coord[i]);
    }

  // clear tempY, tempNum & tempVec
  clear_rational_number(tempNum);
  clear_rational_number(tempY);
  clear_rational_vector(tempVec);

  free(tempExp);
  return;
}

void eval_polynomial_rational(rational_complex_number func, rational_complex_number *jac, polynomial *F, rational_complex_vector vars)
/***************************************************************\
* USAGE: compute F(vars) & JF(vars) where F is a polynomial     *
\***************************************************************/
{
  int i, j, numVars = F->numVariables, numTerms = F->numTerms;
  rational_complex_number tempNum;
  rational_complex_vector tempVec;

  // error checking
  if (numVars != vars->size)
  {
    printf("ERROR: Incorrect number of variables!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // initialize tempNum & tempVec
  initialize_rational_number(tempNum);
  initialize_rational_vector(tempVec, numVars);

  // initialize func & jac to all 0
  set_zero_rational_number(func);
  for (i = 0; i < numVars; i++)
    set_zero_rational_number(jac[i]);
  
  // loop over the number of terms
  for (i = 0; i < numTerms; i++)
  { // compute tempNum = x^exp[i] as well as tempVec = [d/dxj(x^exp[i])]
    eval_exponents_derivative_rational(tempNum, tempVec, vars, F->exponents[i], numVars);

    // multiply by coefficient and add to func
    sum_multiply_rational(func, tempNum, F->coeff[i]);

    // multiply by coefficient and add to jac
    for (j = 0; j < numVars; j++)
      sum_multiply_rational(jac[j], tempVec->coord[j], F->coeff[i]);
  }

  // clear tempNum & tempVec
  clear_rational_number(tempNum);
  clear_rational_vector(tempVec);

  return;
}

void eval_polynomial_system_rational(rational_complex_vector func, rational_complex_matrix jac, polynomial_system *F, rational_complex_vector vars)
/***************************************************************\
* USAGE: compute F(vars) & JF(vars) where F is a poly system    *
\***************************************************************/
{
  int i, numVars = F->numVariables, numPoly = F->numPolynomials;

  // error checking
  if (numVars != vars->size)
  {
    printf("ERROR: Incorrect number of variables!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // setup func & jac to the correct size
  change_size_rational_vector(func, numVars);
  change_size_rational_matrix(jac, numVars, numPoly);

  // loop over the polynomials - compute the function and jacobian
  for (i = 0; i < numPoly; i++)
    eval_polynomial_rational(func->coord[i], jac->entry[i], &F->polynomials[i], vars);

  return;
}

void eval_exponents_rational(rational_complex_number y, rational_complex_vector vars, int *exp, int numVars)
/***************************************************************\
* USAGE: compute y = vars^exp                                   *
\***************************************************************/
{
  int i;
  rational_complex_number tempNum;

  // initialize tempNum
  initialize_rational_number(tempNum);

  // initialize y to 1
  set_one_rational_number(y);

  // compute x_i^exp[i] and multiply on to y
  for (i = 0; i < numVars; i++)
    if (exp[i] > 0)
    { // compute values
      exponentiate_rational(tempNum, vars->coord[i], exp[i]);
      multiply_rational(y, y, tempNum);
    }

  // clear tempNum
  clear_rational_number(tempNum);

  return;
}

void eval_polynomial_only_rational(rational_complex_number func, polynomial *F, rational_complex_vector vars)
/***************************************************************\
* USAGE: compute F(vars) where F is a polynomial                *
\***************************************************************/
{
  int i, numVars = F->numVariables, numTerms = F->numTerms;
  rational_complex_number tempNum;

  // error checking
  if (numVars != vars->size)
  {
    printf("ERROR: Incorrect number of variables!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // initialize tempNum
  initialize_rational_number(tempNum);

  // initialize func to 0
  set_zero_rational_number(func);

  // loop over the number of terms
  for (i = 0; i < numTerms; i++)
  { // compute tempNum = x^exp[i] 
    eval_exponents_rational(tempNum, vars, F->exponents[i], numVars);

    // multiply by coefficient and add to func
    sum_multiply_rational(func, tempNum, F->coeff[i]);
  }

  // clear tempNum
  clear_rational_number(tempNum);

  return;
}

void eval_polynomial_system_only_rational(rational_complex_vector func, polynomial_system *F, rational_complex_vector vars)
/***************************************************************\
* USAGE: compute F(vars) where F is a poly system               *
\***************************************************************/
{
  int i, numVars = F->numVariables, numPoly = F->numPolynomials;

  // error checking
  if (numVars != vars->size)
  {
    printf("ERROR: Incorrect number of variables!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // setup func to the correct size
  change_size_rational_vector(func, numVars);

  // loop over the polynomials - compute the function
  for (i = 0; i < numPoly; i++)
    eval_polynomial_only_rational(func->coord[i], &F->polynomials[i], vars);

  return;
}

