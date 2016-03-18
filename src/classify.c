/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   classify.c: Classifies the points
*/

#include "alphaCertified.h"

void classify_points(int numPoints, complex_vector *Points, polynomial_system *F, configurations *S)
/***************************************************************\
* USAGE: classify the points using the given precision          *
\***************************************************************/
{
  int i, rV, numApproxSolns = 0, numDistinctSolns = 0, numRealSolns = 0, numVars = F->numVariables;
  point_struct *Points_struct = (point_struct *)errMalloc(numPoints * sizeof(point_struct));

  // set the default precision
  setPrec(S->startingPrecision);

  // setup Points_struct and determine which ones are approximate solutions
  for (i = 0; i < numPoints; i++)
  { // initialize 
    initialize_point_struct(&Points_struct[i], numVars);

    // set to active
    Points_struct[i].isActive = 1;

    // copy point
    copy_vector(Points_struct[i].origX, Points[i]);
    copy_vector(Points_struct[i].x, Points[i]);

    // compute ||x||_2
    norm_vector(Points_struct[i].norm_x, Points_struct[i].x);

    // compute alpha, beta, & gamma (and save to original values of alpha, beta, & gamma)
    rV = compute_alpha_beta_gamma(Points_struct[i].Nx, Points_struct[i].alpha, Points_struct[i].beta, Points_struct[i].gamma, F, Points_struct[i].x, S->startingPrecision);
    set_number(Points_struct[i].origAlpha, Points_struct[i].alpha);
    set_number(Points_struct[i].origBeta, Points_struct[i].beta);
    set_number(Points_struct[i].origGamma, Points_struct[i].gamma);

    // check to see if we have successfully computed alpha, beta, & gamma
    if (rV == EXACT_SOLUTION_LU_ERROR)
    { // exact solution
      numApproxSolns += Points_struct[i].isApproxSoln = 1; 
    }
    else if (rV)
    { // error 
      Points_struct[i].isApproxSoln = 0; // unknown 
    }
    else
    { // determine if alpha is small enough to be an approximate solution
      numApproxSolns += Points_struct[i].isApproxSoln = determine_approximate_solution(Points_struct[i].alpha);
    }
  }

  if (S->algorithm >= 1)
  { // now that we have approximate solutions, isolate them
    printf("Isolating %d approximate solution%s.\n\n", numApproxSolns, numApproxSolns == 1 ? "" : "s");
    numDistinctSolns = isolate_approximate_solutions(numPoints, Points_struct, F, S->startingPrecision);

    // now that we have distinct ones, determine which ones are real, if needed
    if (S->algorithm >= 2 && F->isReal)
    { // print message and do the analysis
      printf("Classifying %d distinct approximate solution%s.\n\n", numDistinctSolns, numDistinctSolns == 1 ? "" : "s");
      if (S->realityTest)
      { // use global approach
        numRealSolns = classify_real_points_global(numPoints, Points_struct, F, S->startingPrecision);
      }
      else
      { // use local approach
        numRealSolns = classify_real_points(numPoints, Points_struct, F, S->startingPrecision);
      }
    }
  }

  // refine the solutions
  refine_points(numPoints, Points_struct, F, S->startingPrecision, S->refineDigits);

  // print the data out
  classify_output(numPoints, Points_struct, numApproxSolns, numDistinctSolns, numRealSolns, F->isReal, S, 0, F);

  // clear Points_struct
  for (i = 0; i < numPoints; i++)
    clear_point_struct(&Points_struct[i]);
  free(Points_struct);
  Points_struct = NULL;

  return;
}

void classify_points_rational(int numPoints, rational_complex_vector *Points, polynomial_system *F, configurations *S)
/***************************************************************\
* USAGE: classify the points                                    *
\***************************************************************/
{
  int i, rV, numApproxSolns = 0, numDistinctSolns = 0, numRealSolns = 0, numVars = F->numVariables;
  rational_point_struct *Points_struct = (rational_point_struct *)errMalloc(numPoints * sizeof(rational_point_struct));

  // setup Points_struct and determine which ones are approximate solutions
  for (i = 0; i < numPoints; i++)
  { // initialize 
    initialize_rational_point_struct(&Points_struct[i], numVars);

    // set to active
    Points_struct[i].isActive = 1;

    // copy rational point
    copy_rational_vector(Points_struct[i].origX, Points[i]);
    copy_rational_vector(Points_struct[i].x, Points[i]);

    // compute ||x||_2^2
    norm_sqr_rational_vector(Points_struct[i].norm_sqr_x, Points_struct[i].x);

    // compute alpha^2, beta^2, & gamma^2 (and save to original values of alpha^2, beta^2, & gamma^2)
    rV = compute_alpha_beta_gamma_sqr_rational(Points_struct[i].Nx, Points_struct[i].alpha_sqr, Points_struct[i].beta_sqr, Points_struct[i].gamma_sqr, F, Points_struct[i].x);
    set_rational_number(Points_struct[i].origAlpha_sqr, Points_struct[i].alpha_sqr);
    set_rational_number(Points_struct[i].origBeta_sqr, Points_struct[i].beta_sqr);
    set_rational_number(Points_struct[i].origGamma_sqr, Points_struct[i].gamma_sqr);

    // check to see if we have successfully computed alpha_sqr, beta_sqr, & gamma_sqr
    if (rV == EXACT_SOLUTION_LU_ERROR)
    { // exact solution
      numApproxSolns += Points_struct[i].isApproxSoln = 1;
    }
    else if (rV)
    { // error 
      Points_struct[i].isApproxSoln = 0; // unknown 
    }
    else
    { // determine if alpha is small enough to be an approximate solution
      numApproxSolns += Points_struct[i].isApproxSoln = determine_approximate_solution_rational(Points_struct[i].alpha_sqr);
    }
  }

  if (S->algorithm >= 1)
  { // now that we have approximate solutions, isolate them
    printf("Isolating %d approximate solution%s.\n\n", numApproxSolns, numApproxSolns == 1 ? "" : "s");
    numDistinctSolns = isolate_approximate_solutions_rational(numPoints, Points_struct, F);

    // now that we have distinct ones, determine which ones are real
    if (S->algorithm >= 2 && F->isReal)
    { // print message and do the analysis
      printf("Classifying %d distinct approximate solution%s.\n\n", numDistinctSolns, numDistinctSolns == 1 ? "" : "s");
      if (S->realityTest)
      { // use global approach
        numRealSolns = classify_real_points_global_rational(numPoints, Points_struct, F);
      }
      else
      { // use local approach
        numRealSolns = classify_real_points_rational(numPoints, Points_struct, F);
      }
    }
  }

  // refine the solutions
  refine_points_rational(numPoints, Points_struct, F, S->refineDigits);

  // print the data out
  classify_rational_output(numPoints, Points_struct, numApproxSolns, numDistinctSolns, numRealSolns, F->isReal, S, 0, F);

  // clear Points_struct
  for (i = 0; i < numPoints; i++)
    clear_rational_point_struct(&Points_struct[i]);
  free(Points_struct);
  Points_struct = NULL;

  return;
}

int classify_real_points(int numPoints, point_struct *Pts, polynomial_system *F, int eval_prec)
/***************************************************************\
* USAGE: classify the real points starting with given precision *
\***************************************************************/
{
  int i, curr_prec, rV = 0, numReal = 0;

  for (i = 0; i < numPoints; i++)
    if (Pts[i].isActive == 1 && Pts[i].isApproxSoln)
    { // setup curr_prec
      curr_prec = eval_prec;
    
      // loop to determine if real or not
      while ((rV = determine_real_solution(Pts[i].x, Pts[i].alpha, Pts[i].beta, Pts[i].gamma)) == 0)
      { // we need to double the precision and update
        curr_prec *= 2;
        setPrec_vector(Pts[i].x, curr_prec);
        copy_vector(Pts[i].x, Pts[i].Nx);
        rV = compute_alpha_beta_gamma(Pts[i].Nx, Pts[i].alpha, Pts[i].beta, Pts[i].gamma, F, Pts[i].x, curr_prec);

        // check for errors - should not occur!!
        if (rV == ERROR_LU_DECOMP)
        {
          printf("ERROR: Invalid LU decomposition!\n");
          errExit(ERROR_CONFIGURATION);
        }

        // update the norm
        mpf_set_prec(Pts[i].norm_x, curr_prec);
        norm_vector(Pts[i].norm_x, Pts[i].x);
      }

      if (rV > 0)
      { // real
        Pts[i].isReal = 1;
        numReal++;
      }
      else
      { // not real
        Pts[i].isReal = 0;
      }
    
      // set precision back
      setPrec(eval_prec);  
    }

  return numReal;
}

int classify_real_points_global(int numPoints, point_struct *Pts, polynomial_system *F, int eval_prec)
/***************************************************************\
* USAGE: classify the real points using the global approach     *
\***************************************************************/
{
  int i, j, numReal = 0;
  point_struct conjPt;

  // initialize
  initialize_point_struct(&conjPt, F->numVariables);
  for (i = 0; i < numPoints; i++)
    Pts[i].isReal = 1;

  // loop over the points
  for (i = 0; i < numPoints; i++)
    if (Pts[i].isActive == 1 && Pts[i].isApproxSoln && Pts[i].isReal)
    { // determine if conj(Pts[i]) does not correspond to all the other solutions below it
      setPrec_vector(conjPt.origX, eval_prec); 
      conjugate_vector(conjPt.origX, Pts[i].origX);
      setPrec_vector(conjPt.x, eval_prec); 
      conjugate_vector(conjPt.x, Pts[i].x);
      setPrec_vector(conjPt.Nx, eval_prec); 
      conjugate_vector(conjPt.Nx, Pts[i].Nx);
      mpf_set_prec(conjPt.norm_x, eval_prec);
      mpf_set(conjPt.norm_x, Pts[i].norm_x);
      conjPt.isApproxSoln = Pts[i].isApproxSoln;
      conjPt.isActive = Pts[i].isActive;
      conjPt.isReal = Pts[i].isReal;
      setPrec_number(conjPt.origAlpha, eval_prec);
      set_number(conjPt.origAlpha, Pts[i].origAlpha);
      setPrec_number(conjPt.origBeta, eval_prec);
      set_number(conjPt.origBeta, Pts[i].origBeta);
      setPrec_number(conjPt.origGamma, eval_prec);
      set_number(conjPt.origGamma, Pts[i].origGamma);
      setPrec_number(conjPt.alpha, eval_prec);
      set_number(conjPt.alpha, Pts[i].alpha);
      setPrec_number(conjPt.beta, eval_prec);
      set_number(conjPt.beta, Pts[i].beta);
      setPrec_number(conjPt.gamma, eval_prec);
      set_number(conjPt.gamma, Pts[i].gamma);

      for (j = i+1; j < numPoints && Pts[i].isReal; j++)
        if (Pts[j].isActive == 1 && Pts[j].isApproxSoln)
        { // determine if conjPt & Pts[j] correspond to the same solution
          if (is_same_solution(&conjPt, &Pts[j], F, eval_prec))
          { // cannot be real
            Pts[i].isReal = Pts[j].isReal = 0;
          }
        }

      // update the count
      numReal += Pts[i].isReal;
    }

  clear_point_struct(&conjPt);

  return numReal;
}


int classify_real_points_rational(int numPoints, rational_point_struct *Pts, polynomial_system *F)
/***************************************************************\
* USAGE: classify the real points                               *
\***************************************************************/
{
  int i, rV = 0, numReal = 0;

  for (i = 0; i < numPoints; i++)
    if (Pts[i].isActive == 1 && Pts[i].isApproxSoln)
    { // loop to determine if real or not
      while ((rV = determine_real_solution_rational(Pts[i].x, Pts[i].alpha_sqr, Pts[i].beta_sqr, Pts[i].gamma_sqr)) == 0)
      { // update
        copy_rational_vector(Pts[i].x, Pts[i].Nx);
        rV = compute_alpha_beta_gamma_sqr_rational(Pts[i].Nx, Pts[i].alpha_sqr, Pts[i].beta_sqr, Pts[i].gamma_sqr, F, Pts[i].x);

        // check for errors - should not occur!!
        if (rV == ERROR_LU_DECOMP)
        {
          printf("ERROR: Invalid LU decomposition!\n");
          errExit(ERROR_CONFIGURATION);
        }

        // update the norm
        norm_sqr_rational_vector(Pts[i].norm_sqr_x, Pts[i].x);
      }

      if (rV > 0)
      { // real
        Pts[i].isReal = 1;
        numReal++;
      }
      else
      { // not real
        Pts[i].isReal = 0;
      }
    }

  return numReal;
}

int classify_real_points_global_rational(int numPoints, rational_point_struct *Pts, polynomial_system *F)
/***************************************************************\
* USAGE: classify the real points using the global approach     *
\***************************************************************/
{
  int i, j, numReal = 0;
  rational_point_struct conjPt;

  // initialize
  initialize_rational_point_struct(&conjPt, F->numVariables);
  for (i = 0; i < numPoints; i++)
    Pts[i].isReal = 1;

  for (i = 0; i < numPoints; i++)
    if (Pts[i].isActive == 1 && Pts[i].isApproxSoln && Pts[i].isReal)
    { // determine if conj(Pts[i]) does not correspond to all the other solutions below it
      conjugate_rational_vector(conjPt.origX, Pts[i].origX);
      conjugate_rational_vector(conjPt.x, Pts[i].x);
      conjugate_rational_vector(conjPt.Nx, Pts[i].Nx);
      mpq_set(conjPt.norm_sqr_x, Pts[i].norm_sqr_x);
      conjPt.isApproxSoln = Pts[i].isApproxSoln;
      conjPt.isActive = Pts[i].isActive;
      conjPt.isReal = Pts[i].isReal;
      set_rational_number(conjPt.origAlpha_sqr, Pts[i].origAlpha_sqr);
      set_rational_number(conjPt.origBeta_sqr, Pts[i].origBeta_sqr);
      set_rational_number(conjPt.origGamma_sqr, Pts[i].origGamma_sqr);
      set_rational_number(conjPt.alpha_sqr, Pts[i].alpha_sqr);
      set_rational_number(conjPt.beta_sqr, Pts[i].beta_sqr);
      set_rational_number(conjPt.gamma_sqr, Pts[i].gamma_sqr);

      for (j = i+1; j < numPoints && Pts[i].isReal; j++)
        if (Pts[j].isActive == 1 && Pts[j].isApproxSoln)
        { // determine if conjPt & Pts[j] correspond to the same solution
          if (is_same_solution_rational(&conjPt, &Pts[j], F))
          { // cannot be real 
            Pts[i].isReal = Pts[j].isReal = 0;
          }
        }

      // update the count
      numReal += Pts[i].isReal;
    }

  clear_rational_point_struct(&conjPt);

  return numReal;
}

int is_same_solution_one_test(point_struct *Pt1, point_struct *Pt2)
/***************************************************************\
* USAGE: determine if Pt1 & Pt2 correspond to the same solution *
*  -1 - unknown, 0 - different, 1 - same                        *
\***************************************************************/
{
  int i, size = Pt1->x->size, rV = 0;
  mpf_t alpha0, twice_beta, norm_diff;
  complex_number tempNum;

  // error checking - both need to be approximate solutions
  if (!Pt1->isApproxSoln || !Pt2->isApproxSoln)
  {
    printf("ERROR: Only approximate solutions can be utilized!\n");
    errExit(ERROR_CONFIGURATION);
  }
  else if (Pt1->x->size != Pt2->x->size)
  {
    printf("ERROR: The approximate solutions need to be the same size!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // intialize memory
  mpf_init(alpha0); mpf_init(twice_beta); mpf_init(norm_diff);
  initialize_number(tempNum);

  // set alpha0 = 3/100
  mpf_set_ui(alpha0, 3);
  mpf_div_ui(alpha0, alpha0, 100);

  // compute 2*(beta1 + beta2)
  mpf_add(twice_beta, Pt1->beta->re, Pt2->beta->re);
  mpf_mul_ui(twice_beta, twice_beta, 2);

  // perform a cheap test: determine if | ||x1|| - ||x2|| | > 2*(beta1 + beta2)
  mpf_sub(norm_diff, Pt1->norm_x, Pt2->norm_x);
  mpf_abs(norm_diff, norm_diff);

  if (mpf_cmp(norm_diff, twice_beta) > 0)
  { // can not correspond to the same solution
    rV = 0;
  }
  else
  { // compute ||x1 - x2||_2
    mpf_set_ui(norm_diff, 0);
    for (i = 0; i < size; i++)
    {
      subtract(tempNum, Pt1->x->coord[i], Pt2->x->coord[i]); 
      mpf_mul(tempNum->re, tempNum->re, tempNum->re);
      mpf_mul(tempNum->im, tempNum->im, tempNum->im);
      mpf_add(tempNum->re, tempNum->re, tempNum->im);
      mpf_add(norm_diff, norm_diff, tempNum->re);
    }
    mpf_sqrt(norm_diff, norm_diff);

    // determine if ||x1 - x2|| > 2*(beta1 + beta2)
    if (mpf_cmp(norm_diff, twice_beta) > 0)
    { // can not correspond to the same solution
      rV = 0;
    }
    else
    { // see if alpha is small enough to use robust alpha theorem
      if (mpf_cmp(Pt1->alpha->re, alpha0) <= 0)
      { // determine if ||x1 - x2||_2 <= 1/(20*gamma)
        mpf_set_ui(tempNum->re, 20);
        mpf_mul(tempNum->re, tempNum->re, Pt1->gamma->re);
        mpf_ui_div(tempNum->re, 1, tempNum->re);
      
        if (mpf_cmp(norm_diff, tempNum->re) <= 0)
        { // these correspond to the same solution
         rV = 1;
        }
        else
        { // unknown
          rV = -1;
        } 
      }
      else
      { // unknown
        rV = -1;
      }

      if (rV != 1)
      { // see if alpha is small enough to use robust alpha theorem
        if (mpf_cmp(Pt2->alpha->re, alpha0) <= 0)
        { // determine if ||x1 - x2||_2 <= 1/(20*gamma)
          mpf_set_ui(tempNum->re, 20);
          mpf_mul(tempNum->re, tempNum->re, Pt2->gamma->re);
          mpf_ui_div(tempNum->re, 1, tempNum->re);
      
          if (mpf_cmp(norm_diff, tempNum->re) <= 0)
          { // these correspond to the same solution
            rV = 1;
          }
          else
          { // unknown
            rV = -1;
          }
        }
        else
        { // unknown
          rV = -1;
        }
      }
    }
  }

  // clear memory
  mpf_clear(alpha0); mpf_clear(twice_beta); mpf_clear(norm_diff);
  clear_number(tempNum);

  return rV;
}

int is_same_solution(point_struct *Pt1, point_struct *Pt2, polynomial_system *F, int eval_prec)
/***************************************************************\
* USAGE: determine if Pt1 & Pt2 correspond to the same solution *
*  0 - not the same, 1 - same                                   *
\***************************************************************/
{
  int rV = 0, curr_prec = eval_prec;

  while ((rV = is_same_solution_one_test(Pt1, Pt2)) == -1)
  { // need to perform a newton iteration and test again
    curr_prec *= 2;
    setPrec_vector(Pt1->x, curr_prec);
    setPrec_vector(Pt2->x, curr_prec);
    copy_vector(Pt1->x, Pt1->Nx);
    copy_vector(Pt2->x, Pt2->Nx);

    rV = compute_alpha_beta_gamma(Pt1->Nx, Pt1->alpha, Pt1->beta, Pt1->gamma, F, Pt1->x, curr_prec);

    // check for errors - should not occur!!
    if (rV == ERROR_LU_DECOMP)
    {
      printf("ERROR: Invalid LU decomposition!\n");
      errExit(ERROR_CONFIGURATION);
    }

    // update norm
    mpf_set_prec(Pt1->norm_x, curr_prec);
    norm_vector(Pt1->norm_x, Pt1->x);

    rV = compute_alpha_beta_gamma(Pt2->Nx, Pt2->alpha, Pt2->beta, Pt2->gamma, F, Pt2->x, curr_prec);
    // check for errors - should not occur!!
    if (rV == ERROR_LU_DECOMP)
    {
      printf("ERROR: Invalid LU decomposition!\n");
      errExit(ERROR_CONFIGURATION);
    }

    // update norm
    mpf_set_prec(Pt2->norm_x, curr_prec);
    norm_vector(Pt2->norm_x, Pt2->x);
  }

  // set precision back
  setPrec(eval_prec);  

  return rV;
}

int isolate_approximate_solutions(int numPoints, point_struct *Pts, polynomial_system *F, int eval_prec)
/***************************************************************\
* USAGE: determine which approximate solutions correspond to the*
*  same solution, starting with given precision                 *
\***************************************************************/
{
  int i, j, numDistinct = 0;

  for (i = 0; i < numPoints; i++)
    if (Pts[i].isActive == 1 && Pts[i].isApproxSoln)
    { // we have a new distinct solution
      numDistinct++;

      // we know it is distinct for j < i, so we check for j > i
      for (j = i+1; j < numPoints; j++)
        if (Pts[j].isActive == 1 && Pts[j].isApproxSoln)
        { // determine if i & j correspond to the same solution
          if (is_same_solution(&Pts[i], &Pts[j], F, eval_prec))
          { // save the correspondence between i & j
            Pts[j].isActive = -i;
          }
        }
    }

  return numDistinct;
}

int is_same_solution_one_test_rational(rational_point_struct *Pt1, rational_point_struct *Pt2)
/***************************************************************\
* USAGE: determine if Pt1 & Pt2 correspond to the same solution *
*  -1 - unknown, 0 - different, 1 - same                        *
\***************************************************************/
{
  int i, size = Pt1->x->size, rV = 0;
  mpq_t alpha0_sqr, four_beta_sqr, norm_diff_sqr;
  rational_complex_number tempNum;

  // error checking - both need to be approximate solutions
  if (!Pt1->isApproxSoln || !Pt2->isApproxSoln)
  {
    printf("ERROR: Only approximate solutions can be utilized!\n");
    errExit(ERROR_CONFIGURATION);
  }
  else if (Pt1->x->size != Pt2->x->size)
  {
    printf("ERROR: The approximate solutions need to be the same size!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // intialize memory
  mpq_init(alpha0_sqr); mpq_init(four_beta_sqr); mpq_init(norm_diff_sqr);
  initialize_rational_number(tempNum);

  // setup alpha0_sqr = (3/100)^2 = 9/10000
  mpq_set_ui(alpha0_sqr, 9, 10000);

  // compute ||x1 - x2||_2^2
  mpq_set_ui(norm_diff_sqr, 0, 1);
  for (i = 0; i < size; i++)
  {
    subtract_rational(tempNum, Pt1->x->coord[i], Pt2->x->coord[i]); 
    mpq_mul(tempNum->re, tempNum->re, tempNum->re);
    mpq_mul(tempNum->im, tempNum->im, tempNum->im);
    mpq_add(tempNum->re, tempNum->re, tempNum->im);
    mpq_add(norm_diff_sqr, norm_diff_sqr, tempNum->re);
  }
  
  // determine if ||x1 - x2||_2 > 2*(beta1 + beta2) by using two inequalities
  // (1) ||x1 - x2||_2^2 > 4*(beta1^2 + beta2^2) AND (2) (||x1 - x2||_2^2 - 4(beta1^2 + beta2^2))^2 > 64*beta1^2*beta2^2

  // compute 4*(beta1^2 + beta2^2)
  mpq_add(four_beta_sqr, Pt1->beta_sqr->re, Pt2->beta_sqr->re);
  mpq_set_ui(tempNum->re, 4, 1);
  mpq_mul(four_beta_sqr, four_beta_sqr, tempNum->re);

  // compute (||x1 - x2||_2^2 - 4*(beta1^2 + beta2^2))^2
  mpq_sub(tempNum->re, norm_diff_sqr, four_beta_sqr);
  mpq_mul(tempNum->re, tempNum->re, tempNum->re);

  // compute 64*beta1^2*beta2^2
  mpq_set_ui(tempNum->im, 64, 1);
  mpq_mul(tempNum->im, tempNum->im, Pt1->beta_sqr->re);
  mpq_mul(tempNum->im, tempNum->im, Pt2->beta_sqr->re);

  if (mpq_cmp(norm_diff_sqr, four_beta_sqr) > 0 && mpq_cmp(tempNum->re, tempNum->im) > 0)
  { // can not correspond to the same solution
    rV = 0;
  }
  else 
  { // see if alpha is small enough to use robust alpha theorem
    if (mpq_cmp(Pt1->alpha_sqr->re, alpha0_sqr) <= 0)
    { // determine if ||x1 - x2||_2^2 <= 1/(400*gamma^2)
      mpq_set_ui(tempNum->re, 400, 1);
      mpq_mul(tempNum->re, tempNum->re, Pt1->gamma_sqr->re);
      mpq_inv(tempNum->re, tempNum->re);
      
      if (mpq_cmp(norm_diff_sqr, tempNum->re) <= 0)
      { // these correspond to the same solution
        rV = 1;
      }
      else
      { // unknown
        rV = -1;
      } 
    }
    else
    { // unknown
      rV = -1;
    }
  
    if (rV != 1)
    { // see if alpha is small enough to use robust alpha theorem
      if (mpq_cmp(Pt2->alpha_sqr->re, alpha0_sqr) <= 0)
      { // determine if ||x1 - x2||_2^2 <= 1/(400*gamma^2)
        mpq_set_ui(tempNum->re, 400, 1);
        mpq_mul(tempNum->re, tempNum->re, Pt2->gamma_sqr->re);
        mpq_inv(tempNum->re, tempNum->re);
      
        if (mpq_cmp(norm_diff_sqr, tempNum->re) <= 0)
        { // these correspond to the same solution
          rV = 1;
        }
        else
        { // unknown
          rV = -1;
        }
      }
      else
      { // unknown
        rV = -1;
      }
    }
  }

  // clear memory
  mpq_clear(alpha0_sqr); mpq_clear(four_beta_sqr); mpq_clear(norm_diff_sqr);
  clear_rational_number(tempNum);

  return rV;
}

int is_same_solution_rational(rational_point_struct *Pt1, rational_point_struct *Pt2, polynomial_system *F)
/***************************************************************\
* USAGE: determine if Pt1 & Pt2 correspond to the same solution *
*  0 - not the same, 1 - same                                   *
\***************************************************************/
{
  int rV = 0;

  while ((rV = is_same_solution_one_test_rational(Pt1, Pt2)) == -1)
  { // need to perform a newton iteration and test again
    copy_rational_vector(Pt1->x, Pt1->Nx);
    copy_rational_vector(Pt2->x, Pt2->Nx);

    rV = compute_alpha_beta_gamma_sqr_rational(Pt1->Nx, Pt1->alpha_sqr, Pt1->beta_sqr, Pt1->gamma_sqr, F, Pt1->x);

    // check for errors - should not occur!!
    if (rV == ERROR_LU_DECOMP)
    {
      printf("ERROR: Invalid LU decomposition!\n");
      errExit(ERROR_CONFIGURATION);
    }

    // update norm
    norm_sqr_rational_vector(Pt1->norm_sqr_x, Pt1->x);

    rV = compute_alpha_beta_gamma_sqr_rational(Pt2->Nx, Pt2->alpha_sqr, Pt2->beta_sqr, Pt2->gamma_sqr, F, Pt2->x);

    // check for errors - should not occur!!
    if (rV == ERROR_LU_DECOMP)
    {
      printf("ERROR: Invalid LU decomposition!\n");
      errExit(ERROR_CONFIGURATION);
    }

    // update norm
    norm_sqr_rational_vector(Pt2->norm_sqr_x, Pt2->x);
  }

  return rV;
}

int isolate_approximate_solutions_rational(int numPoints, rational_point_struct *Pts, polynomial_system *F)
/***************************************************************\
* USAGE: determine which approximate solutions correspond to the*
*  same solution, starting with given precision                 *
\***************************************************************/
{
  int i, j, numDistinct = 0;

  for (i = 0; i < numPoints; i++)
    if (Pts[i].isActive == 1 && Pts[i].isApproxSoln)
    { // we have a new distinct solution
      numDistinct++;

      // we know it is distinct for j < i, so we check for j > i
      for (j = i+1; j < numPoints; j++)
        if (Pts[j].isActive == 1 && Pts[j].isApproxSoln)
        { // determine if i & j correspond to the same solution
          if (is_same_solution_rational(&Pts[i], &Pts[j], F))
          { // save the correspondence between i & j
            Pts[j].isActive = -i;
          }
        }
    }

  return numDistinct;
}

int determine_approximate_solution(complex_number alpha)
/***************************************************************\
* USAGE: determine if approximate solution: 1 - yes, 0 - no     *
\***************************************************************/
{ // we take alpha0 = (13 - 3*sqrt(17))/4
  int rV = 0;
  mpf_t alpha0;

  mpf_init(alpha0);

  // setup alpha0 = (13 - 3*sqrt(17))/4
  mpf_set_ui(alpha0, 17);
  mpf_sqrt(alpha0, alpha0);
  mpf_mul_ui(alpha0, alpha0, 3);
  mpf_ui_sub(alpha0, 13, alpha0);
  mpf_div_ui(alpha0, alpha0, 4);

  // determine if alpha <= alpha0
  rV = (mpf_cmp(alpha->re, alpha0) <= 0);

  mpf_clear(alpha0);  

  return rV;
}

int determine_approximate_solution_rational(rational_complex_number alpha_sqr)
/***************************************************************\
* USAGE: determine if approximate solution: 1 - yes, 0 - no     *
\***************************************************************/
{ // we want alpha_sqr <= alpha0^2 = ((13 - 3*sqrt(17))/4)^2 <= 1/16
  // this is equivalent to {alpha_sqr <= 1/16 && 17 <= ((161 - 8*alpha_sqr)/39)^2}
  int rV = 0;
  rational_complex_number alpha0;

  initialize_rational_number(alpha0); 

  // verify alpha_sqr <= 1/16
  mpq_set_ui(alpha0->re, 1, 16);

  // determine if alpha_sqr <= 1/16
  rV = (mpq_cmp(alpha_sqr->re, alpha0->re) <= 0);
  if (rV)
  { // now check 17 <= ((161 - 8*alpha_sqr)/39)^2
    mpq_set_ui(alpha0->re, 161, 1); 
    mpq_mul_2exp(alpha0->im, alpha_sqr->re, 3); // 2^3*alpha_sqr
    mpq_sub(alpha0->re, alpha0->re, alpha0->im);
    mpq_set_ui(alpha0->im, 1, 39);
    mpq_mul(alpha0->re, alpha0->re, alpha0->im);
    mpq_mul(alpha0->re, alpha0->re, alpha0->re);
    mpq_set_ui(alpha0->im, 17, 1);

    // determine if 17 <= ((161 - 8*alpha_sqr)/39)^2
    rV = (mpq_cmp(alpha0->im, alpha0->re) <= 0);
  }

  clear_rational_number(alpha0);    

  return rV;
}

int determine_real_solution(complex_vector x, complex_number alpha, complex_number beta, complex_number gamma)
/***************************************************************\
* USAGE: determine if real: 1 - real, 0 - unknown, -1 - not real*
\***************************************************************/
{ // we take alpha0 = 3/100
  int rV = 0;
  mpf_t tempMPF, norm_imag;

  mpf_init(tempMPF);
  mpf_init(norm_imag);

  // compute ||x - \pi(x)||_2
  norm_imag_vector(norm_imag, x);

  // compute 2*beta
  mpf_mul_ui(tempMPF, beta->re, 2);

  // determine if norm_imag > 2*beta
  if (mpf_cmp(norm_imag, tempMPF) > 0)
  { // not real
    rV = -1;
  }
  else
  { // determine if real

    // set to 3/100
    mpf_set_ui(tempMPF, 3);
    mpf_div_ui(tempMPF, tempMPF, 100);

    // determine if alpha <= 3/100
    if (mpf_cmp(alpha->re, tempMPF) <= 0)
    { // determine if gamma is finite 
      if (mpfr_inf_p(gamma->re))
      { // set tempMPF = 0
        mpf_set_ui(tempMPF, 0);
      }
      else
      { // set tempMPF = 1/(20*gamma)
        mpf_mul_ui(tempMPF, gamma->re, 20);
        mpf_ui_div(tempMPF, 1, tempMPF);
      }

      // determine if norm_imag <= tempMPF
      rV = (mpf_cmp(norm_imag, tempMPF) <= 0); // either we know it is real or unknown
    }
    else
    { // unable to determine if real
      rV = 0;
    }
  }

  // clear memory
  mpf_clear(tempMPF);  
  mpf_clear(norm_imag);

  return rV;
}

int determine_real_solution_rational(rational_complex_vector x, rational_complex_number alpha_sqr, rational_complex_number beta_sqr, rational_complex_number gamma_sqr)
/***************************************************************\
* USAGE: determine if real: 1 - real, 0 - unknown, -1 - not real*
\***************************************************************/
{ // we take alpha0 = 3/100
  int rV = 0;
  mpq_t tempRat, norm_sqr_imag;

  mpq_init(tempRat);
  mpq_init(norm_sqr_imag);

  // compute ||x - \pi(x)||_2^2
  norm_sqr_imag_rational_vector(norm_sqr_imag, x);

  // compute 4*beta_sqr
  mpq_mul_2exp(tempRat, beta_sqr->re, 2);

  // determine if norm_sqr_imag > 4*beta_sqr
  if (mpq_cmp(norm_sqr_imag, tempRat) > 0)
  { // not real
    rV = -1;
  }
  else
  { // determine if real

    // set to 9/10000
    mpq_set_ui(tempRat, 9, 10000);

    // determine if alpha_sqr <= 9/10000
    if (mpq_cmp(alpha_sqr->re, tempRat) <= 0)
    { // determine if gamma_sqr is finite
      if (mpz_cmp_ui(mpq_denref(gamma_sqr->re), 0) == 0)
      { // set tempRat = 0
        mpq_set_ui(tempRat, 0, 1);
      }
      else
      { // set tempRat = 1/(400*gamma_sqr)
        mpq_set_ui(tempRat, 400, 1);
        mpq_mul(tempRat, tempRat, gamma_sqr->re);
        mpq_inv(tempRat, tempRat);
      }

      // determine if norm_sqr_imag <= tempRat
      rV = (mpq_cmp(norm_sqr_imag, tempRat) <= 0); // either we know it is real or unknown
    }
    else
    { // unable to determine if real
      rV = 0;
    }
  }

  // clear memory
  mpq_clear(tempRat);  
  mpq_clear(norm_sqr_imag);

  return rV;
}
