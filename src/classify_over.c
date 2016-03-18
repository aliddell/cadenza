/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   classify_over.c: Classifies the points for overdetermined polynomial systems
*/

#include "alphaCertified.h"

void sort_polynomials(polynomial_system *F, int **perm)
/***************************************************************\
* USAGE: perm[i] corresponds to the poly with ith largest degree*
\***************************************************************/
{ // use a simple bubble sort the degrees from largest to smallest
  int i, j, k, sorted = 0;

  // allocate & initialize perm
  *perm = (int *)errRealloc(*perm, F->numPolynomials * sizeof(int));
  for (i = 0; i < F->numPolynomials; i++)
    (*perm)[i] = i;

  // sort
  for (i = 0; i < F->numPolynomials; i++)
  { // initialize sorted
    sorted = 1;
    for (j = i+1; j < F->numPolynomials; j++)
      if (F->polynomials[(*perm)[j-1]].degree < F->polynomials[(*perm)[j]].degree)
      { // swap
        k = (*perm)[j-1];
        (*perm)[j-1] = (*perm)[j];
        (*perm)[j] = k;
        sorted = 0;
      }
    
    // end loop if sorted
    if (sorted)
      i = F->numPolynomials; 
  }

  return;
}

int compare_exponents(int **expStruct, int numRows, int numCols, int *newExp)
/***************************************************************\
* USAGE: determine if newExp appears in expStruct and its loc   *
\***************************************************************/
{
  int i, j, rV = -1, same = 0;

  for (i = 0; i < numRows && rV == -1; i++)
  { // determine if the same as newExp
    same = 1;    
    for (j = 0; j < numCols && same; j++)
      if (expStruct[i][j] != newExp[j])
        same = 0;

    // see if they are the same
    if (same)
      rV = i;
  }

  return rV;
}

void randomize_polynomials(polynomial *randPoly, polynomial_system *F, int *perm, rational_complex_matrix A)
/***************************************************************\
* USAGE: perform the actual randomization                       *
\***************************************************************/
{
  int i, j, k, l, rV, curr = 0, numTerms_extra = 0, numVars = F->numVariables, numPoly = F->numPolynomials;

  // count the number of terms in the polynomials at the bottom that are randomized up
  for (i = F->numVariables; i < numPoly; i++)
    numTerms_extra += F->polynomials[perm[i]].numTerms;

  for (i = 0; i < F->numVariables; i++)
  { // set the number of variables & terms in the ith polynomial of F_rand
    randPoly[i].numVariables = numVars;
    randPoly[i].numTerms = F->polynomials[perm[i]].numTerms + numTerms_extra;
     
    // allocate memory for ith polynomial in F_rand
    mpq_init(randPoly[i].norm_sqr);
    randPoly[i].coeff = (rational_complex_number *)errMalloc(randPoly[i].numTerms * sizeof(rational_complex_number));
    randPoly[i].exponents = (int **)errMalloc(randPoly[i].numTerms * sizeof(int *));

    // setup coefficients & exponents
    curr = 0;
    for (j = 0; j < F->polynomials[perm[i]].numTerms; j++)
      if ((rV = compare_exponents(randPoly[i].exponents, curr, numVars, F->polynomials[perm[i]].exponents[j])) >= 0)
      { // this one corresponds one already found - add on to coefficient
        add_rational(randPoly[i].coeff[rV], randPoly[i].coeff[rV], F->polynomials[perm[i]].coeff[j]);
      }
      else // new monomial
      { // allocate & initialize memory
        randPoly[i].exponents[curr] = (int *)errMalloc(numVars * sizeof(int));
        initialize_rational_number(randPoly[i].coeff[curr]);
 
        // copy exponents
        for (k = 0; k < numVars; k++)
          randPoly[i].exponents[curr][k] = F->polynomials[perm[i]].exponents[j][k];

        // copy coeff
        set_rational_number(randPoly[i].coeff[curr], F->polynomials[perm[i]].coeff[j]);

        curr++;
      }

    // randomize in the other ones
    for (l = numVars; l < numPoly; l++)
      for (j = 0; j < F->polynomials[perm[l]].numTerms; j++)
        if ((rV = compare_exponents(randPoly[i].exponents, curr, numVars, F->polynomials[perm[l]].exponents[j])) >= 0)
        { // this one corresponds one already found - multiply and add on to coefficient
          sum_multiply_rational(randPoly[i].coeff[rV], A->entry[i][l-numVars], F->polynomials[perm[l]].coeff[j]);
        }
        else // new monomial
        { // allocate & initialize memory
          randPoly[i].exponents[curr] = (int *)errMalloc(numVars * sizeof(int));
          initialize_rational_number(randPoly[i].coeff[curr]);

          // copy exponents
          for (k = 0; k < numVars; k++)
            randPoly[i].exponents[curr][k] = F->polynomials[perm[l]].exponents[j][k];

          // setup coeff
          multiply_rational(randPoly[i].coeff[curr], A->entry[i][l-numVars], F->polynomials[perm[l]].coeff[j]);

          curr++;
        }    

    // adjust the number of terms and size of coeff & exponents
    randPoly[i].numTerms = curr;
    randPoly[i].coeff = (rational_complex_number *)errRealloc(randPoly[i].coeff, randPoly[i].numTerms * sizeof(rational_complex_number));
    randPoly[i].exponents = (int **)errRealloc(randPoly[i].exponents, randPoly[i].numTerms * sizeof(int *));

    // setup degree - same as corresponding one in F
    randPoly[i].degree = F->polynomials[perm[i]].degree;

    // compute norm_sqr
    norm_sqr_polynomial(randPoly[i].norm_sqr, &randPoly[i]);
  }

  return;
}

void setup_randomized_polynomials(polynomial_system *F_rand, polynomial_system *F, int *perm)
/***************************************************************\
* USAGE: setup the randomized polynomials in F_rand             *
\***************************************************************/
{
  int A_rows, A_cols;
  rational_complex_matrix A;

  // find the size of A (we are randomizing with [I A])
  A_rows = F->numVariables;
  A_cols = F->numPolynomials - F->numVariables;

  initialize_rational_matrix(A, A_rows, A_cols);

  // compute a random matrix 
  if (F->isReal)
    random_real_rational_matrix(A, A_rows, A_cols);  
  else
    random_rational_matrix(A, A_rows, A_cols);

  // perform the randomization now that everything is setup
  randomize_polynomials(F_rand->polynomials, F, perm, A);

  clear_rational_matrix(A);

  return;
}

void setup_randomized_systems(int numSystems, polynomial_system *F_rand, polynomial_system *F)
/***************************************************************\
* USAGE: create 'numSystems' number of randomized poly systems  *
\***************************************************************/
{
  int i, j, *perm = NULL;

  // setup perm to sort the polynomials
  sort_polynomials(F, &perm);

  for (i = 0; i < numSystems; i++)
  { // initialize F_rand[i]
    initialize_polynomial_system(&F_rand[i]);

    // setup number of variables & functions == number of variables in original system
    F_rand[i].numPolynomials = F_rand[i].numVariables = F->numVariables;

    // copy maximum degree & isReal
    F_rand[i].maximumDegree = F->maximumDegree;
    F_rand[i].isReal = F->isReal;

    // set norm_sqr to 0
    mpq_set_ui(F_rand[i].norm_sqr, 0, 1);

    // allocate the number of polynomials
    F_rand[i].polynomials = (polynomial *)errMalloc(F_rand[i].numPolynomials * sizeof(polynomial));

    // setup the polynomials
    setup_randomized_polynomials(&F_rand[i], F, perm);

    // compute the norm
    mpq_set_ui(F_rand[i].norm_sqr, 0, 1);
    for (j = 0; j < F_rand[i].numPolynomials; j++)
      mpq_add(F_rand[i].norm_sqr, F_rand[i].norm_sqr, F_rand[i].polynomials[j].norm_sqr);
  }

  free(perm);

  return;
}

void classify_points_over(int numPoints, complex_vector *Points, polynomial_system *F, configurations *S)
/***************************************************************\
* USAGE: classify the points using the given precision          *
\***************************************************************/
{
  int i, numApproxSolns = 0, numDistinctSolns = 0, numRealSolns = 0, numVars = F->numVariables;
  point_struct *Points_struct = (point_struct *)errMalloc(numPoints * sizeof(point_struct));
  polynomial_system *F_rand = (polynomial_system *)errMalloc(S->numRandomSystems * sizeof(polynomial_system));

  // setup randomized systems
  setup_randomized_systems(S->numRandomSystems, F_rand, F);

  // set the default precision
  setPrec(S->startingPrecision);

  // setup Points_struct and determine which ones are approximate solutions to every randomized system
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

    // initialize isApproxSoln
    Points_struct[i].isApproxSoln = 1;

    // determine if it is an approx soln for every randomized system
    determine_over_solution(&Points_struct[i].isApproxSoln, &Points_struct[i], F_rand, S->numRandomSystems, S->randomDigits, S->startingPrecision);

    if (Points_struct[i].isApproxSoln)
    { // update number of approx solutions
      numApproxSolns++;
    }
  }

  if (S->algorithm >= 1)
  { // now that we have approximate solutions, isolate them
    printf("Isolating %d approximate solution%s.\n\n", numApproxSolns, numApproxSolns == 1 ? "" : "s");
    numDistinctSolns = isolate_approximate_solutions(numPoints, Points_struct, &F_rand[0], S->startingPrecision);

    // now that we have distinct ones, determine which ones are real
    if (S->algorithm >= 2 && F->isReal)
    { // print message and do the analysis
      printf("Classifying %d distinct approximate solution%s.\n\n", numDistinctSolns, numDistinctSolns == 1 ? "" : "s");
      if (S->realityTest)
      { // use global approach
        numRealSolns = classify_real_points_global(numPoints, Points_struct, &F_rand[0], S->startingPrecision);
      }
      else
      { // use local approach
        numRealSolns = classify_real_points(numPoints, Points_struct, &F_rand[0], S->startingPrecision);
      }
    }
  }

  // refine the solutions
  refine_points(numPoints, Points_struct, &F_rand[0], S->startingPrecision, S->refineDigits);

  // print the data out
  classify_output(numPoints, Points_struct, numApproxSolns, numDistinctSolns, numRealSolns, F->isReal, S, 1, F);

  // clear Points_struct
  for (i = 0; i < numPoints; i++)
    clear_point_struct(&Points_struct[i]);
  free(Points_struct);
  Points_struct = NULL;

  return;
}

void classify_points_over_rational(int numPoints, rational_complex_vector *Points, polynomial_system *F, configurations *S)
/***************************************************************\
* USAGE: classify the points                                    *
\***************************************************************/
{
  int i, numApproxSolns = 0, numDistinctSolns = 0, numRealSolns = 0, numVars = F->numVariables;
  rational_point_struct *Points_struct = (rational_point_struct *)errMalloc(numPoints * sizeof(rational_point_struct));
  polynomial_system *F_rand = (polynomial_system *)errMalloc(S->numRandomSystems * sizeof(polynomial_system));

  // setup randomized systems
  setup_randomized_systems(S->numRandomSystems, F_rand, F);

  // setup Points_struct and determine which ones are approximate solutions to every randomized system
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

    // initialize isApproxSoln
    Points_struct[i].isApproxSoln = 1;

    // determine if it is an approx soln for every randomized system
    determine_over_solution_rational(&Points_struct[i].isApproxSoln, &Points_struct[i], F_rand, S->numRandomSystems, S->randomDigits);

    if (Points_struct[i].isApproxSoln)
    { // update number of approx solutions
      numApproxSolns++;
    }
  }

  if (S->algorithm >= 1)
  { // now that we have approximate solutions, isolate them
    printf("Isolating %d approximate solution%s.\n\n", numApproxSolns, numApproxSolns == 1 ? "" : "s");
    numDistinctSolns = isolate_approximate_solutions_rational(numPoints, Points_struct, &F_rand[0]);

    // now that we have distinct ones, determine which ones are real
    if (S->algorithm >= 2 && F->isReal)
    { // print message and do the analysis
      printf("Classifying %d distinct approximate solution%s.\n\n", numDistinctSolns, numDistinctSolns == 1 ? "" : "s");
      if (S->realityTest)
      { // use global approach
        numRealSolns = classify_real_points_global_rational(numPoints, Points_struct, &F_rand[0]);
      }
      else
      { // use local approach
        numRealSolns = classify_real_points_rational(numPoints, Points_struct, &F_rand[0]);
      }
    }
  }

  // refine the solutions
  refine_points_rational(numPoints, Points_struct, &F_rand[0], S->refineDigits);

  // print the data out
  classify_rational_output(numPoints, Points_struct, numApproxSolns, numDistinctSolns, numRealSolns, F->isReal, S, 1, F);

  // clear Points_struct
  for (i = 0; i < numPoints; i++)
    clear_rational_point_struct(&Points_struct[i]);
  free(Points_struct);
  Points_struct = NULL;

  return;
}

int determine_over_associated_solution(point_struct *Pts, int numPts, int digits, int eval_prec)
/***************************************************************\
* USAGE: determine if associated solutions are close enough     *
*   1 - yes, 0 - unknown, -1 - no                               *
\***************************************************************/
{
  int i, j, k, numVars, rV = 0;
  mpf_t twice_beta, norm_diff, max_dist, tol;
  complex_number tempNum;

  // verify numPts >= 2 && digits >= 0
  if (numPts < 2)
  {
    printf("ERROR: Need at least 2 points!\n");
    errExit(ERROR_CONFIGURATION);
  }
  else if (digits < 0)
  {
    printf("ERROR: The number of digits must be nonnegative!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // set the precision
  setPrec(eval_prec);

  // intialize memory
  mpf_init(twice_beta); 
  mpf_init(norm_diff);
  mpf_init(max_dist);
  mpf_init(tol);
  initialize_number(tempNum);

  // initialize maximum
  mpf_set_ui(max_dist, 0);

  // setup tol = 10^-digits
  mpf_set_si(tol, -digits);
  mpfr_exp10(tol, tol, __gmp_default_rounding_mode);

  // find the number of variables
  numVars = Pts[0].x->size;

  // determine if they can correspond to the same point
  for (i = 0; i < numPts && !rV; i++)
    for (j = i+1; j < numPts && !rV; j++)
    { // compute ||xi - xj||_2
      mpf_set_ui(norm_diff, 0);
      for (k = 0; k < numVars; k++)
      { 
        subtract(tempNum, Pts[i].x->coord[k], Pts[j].x->coord[k]); 
        mpf_mul(tempNum->re, tempNum->re, tempNum->re);
        mpf_mul(tempNum->im, tempNum->im, tempNum->im);
        mpf_add(tempNum->re, tempNum->re, tempNum->im);
        mpf_add(norm_diff, norm_diff, tempNum->re);
      }
      mpf_sqrt(norm_diff, norm_diff);

      // determine if ||xi - xj|| > 2*(betai + betaj) 

      // compute 2*(betai + betaj)
      mpf_add(twice_beta, Pts[i].beta->re, Pts[j].beta->re);
      mpf_mul_ui(twice_beta, twice_beta, 2);

      if (mpf_cmp(norm_diff, twice_beta) > 0)
      { // can not correspond to the same solution
        rV = -1;
      }
      else
      { // update maximum of ||pi - pj||_2 <= 2*betai + 2*betaj + ||xi - xj||_2
        mpf_add(twice_beta, twice_beta, norm_diff);

        if (mpf_cmp(twice_beta, max_dist) > 0)
        { // update maximum distance
          mpf_set(max_dist, twice_beta);
        }
      }
    }

  if (!rV)
  { // check to see if the maximum is small enough
    if (mpf_cmp(max_dist, tol) <= 0)
    { // all of them are close enough
      rV = 1;
    }
    else
    { // unknown
      rV = 0;
    }
  }

  // clear memory
  mpf_clear(twice_beta); 
  mpf_clear(norm_diff);
  mpf_clear(max_dist);
  mpf_clear(tol);
  clear_number(tempNum);

  return rV;
}

int determine_over_associated_solution_rational(rational_point_struct *Pts, int numPts, int digits)
/***************************************************************\
* USAGE: determine if associated solutions are close enough     *
*   1 - yes, 0 - unknown, -1 - no                               *
\***************************************************************/
{
  int i, j, k, numVars, rV = 0, allSmall = 1;
  mpq_t beta, four_beta_sqr, norm_diff_sqr, tol, tol_sqr, q, q_sqr;
  rational_complex_number tempNum;

  // verify numPts >= 2 && digits >= 0
  if (numPts < 2)
  {
    printf("ERROR: Need at least 2 points!\n");
    errExit(ERROR_CONFIGURATION);
  }
  else if (digits < 0)
  {
    printf("ERROR: The number of digits must be nonnegative!\n");
    errExit(ERROR_CONFIGURATION);
  }

  // intialize memory
  mpq_init(beta);
  mpq_init(four_beta_sqr); 
  mpq_init(norm_diff_sqr);
  mpq_init(tol);
  mpq_init(tol_sqr);
  mpq_init(q);
  mpq_init(q_sqr);
  initialize_rational_number(tempNum);

  // setup tol = (1/10)^digits
  mpq_set_ui(tol, 1, 10);
  exponentiate_mpq(tol, tol, digits);
  mpq_mul(tol_sqr, tol, tol);

  // find the number of variables
  numVars = Pts[0].x->size;

  // determine if they can not correspond to the same solution
  for (i = 0; i < numPts && !rV; i++)
    for (j = i+1; j < numPts && !rV; j++)
    { // compute ||xi - xj||_2^2
      mpq_set_ui(norm_diff_sqr, 0, 1);
      for (k = 0; k < numVars; k++)
      { 
        subtract_rational(tempNum, Pts[i].x->coord[k], Pts[j].x->coord[k]); 
        mpq_mul(tempNum->re, tempNum->re, tempNum->re);
        mpq_mul(tempNum->im, tempNum->im, tempNum->im);
        mpq_add(tempNum->re, tempNum->re, tempNum->im);
        mpq_add(norm_diff_sqr, norm_diff_sqr, tempNum->re);
      }

      // determine if ||xi - xj||_2 > 2*(betai + betaj) by using two inequalities
      // (1) ||xi - xj||_2^2 > 4*(betai^2 + betaj^2)
      // (2) (||xi - xj||_2^2 - 4(betai^2 + betaj^2))^2 > 64*betai^2*betaj^2

      // compute 4*(betai^2 + betaj^2)
      mpq_add(four_beta_sqr, Pts[i].beta_sqr->re, Pts[j].beta_sqr->re);
      mpq_set_ui(tempNum->re, 4, 1);
      mpq_mul(four_beta_sqr, four_beta_sqr, tempNum->re);

      // verify (1)
      if (mpq_cmp(norm_diff_sqr, four_beta_sqr) > 0)
      { // compute (||xi - xj||_2^2 - 4*(betai^2 + betaj^2))^2
        mpq_sub(tempNum->re, norm_diff_sqr, four_beta_sqr);
        mpq_mul(tempNum->re, tempNum->re, tempNum->re);

        // compute 64*betai^2*betaj^2
        mpq_set_ui(tempNum->im, 64, 1);
        mpq_mul(tempNum->im, tempNum->im, Pts[i].beta_sqr->re);
        mpq_mul(tempNum->im, tempNum->im, Pts[j].beta_sqr->re);

        // verify (2)
        if (mpq_cmp(tempNum->re, tempNum->im) > 0)
        { // can not correspond to the same solution
          rV = -1;
        }
      }

      if (!rV && allSmall)
      { // determine if ||xi - xj||_2 + 2*(betai + betaj) < tol by using the four inequalities
        // (1) tol^2 > ||xi - xj||_2^2
        // Define q = (tol^2 + ||xi - xj||_2^2 - 4(betai^2 + betaj^2))/2
        // (2) q > 0
        // (3) q^2 > 16*betai^2*betaj^2 + ||xi - xj||_2^2*tol^2
        // (4) ((q^2 - 16*beta^2*betaj^2 - ||xi - xj||_2^2*tol^2)/2)^2 > 16*betai^2*betaj^2*||xi - xj||_2^2*tol^2

        // verify (1)
        if (mpq_cmp(tol_sqr, norm_diff_sqr) > 0)
        { // compute q
          mpq_add(q, tol_sqr, norm_diff_sqr);
          mpq_sub(q, q, four_beta_sqr);
          mpq_set_ui(tempNum->re, 1, 2);
          mpq_mul(q, q, tempNum->re);

          // verify (2)
          if (mpq_sgn(q) > 0)
          { // compute q^2
            mpq_mul(q_sqr, q, q);

            // compute 16*betai^2*betaj^2
            mpq_set_ui(tempNum->re, 16, 1);
            mpq_mul(tempNum->re, tempNum->re, Pts[i].beta_sqr->re);
            mpq_mul(tempNum->re, tempNum->re, Pts[j].beta_sqr->re);

            // compute 16*betai^2*betaj^2 + ||xi - xj||_2^2*tol^2
            mpq_mul(tempNum->im, norm_diff_sqr, tol_sqr);
            mpq_add(tempNum->im, tempNum->im, tempNum->re);
 
            // verify (3)
            if (mpq_cmp(q_sqr, tempNum->im) > 0)
            { // compute ((q^2 - 16*beta^2*betaj^2 - ||xi - xj||_2^2*tol^2)/2)^2
              mpq_sub(q_sqr, q_sqr, tempNum->im);
              mpq_set_ui(tempNum->im, 1, 2);
              mpq_mul(q_sqr, q_sqr, tempNum->im);
              mpq_mul(q_sqr, q_sqr, q_sqr);

              // compute 16*betai^2*betaj^2*||xi - xj||_2^2*tol^2
              mpq_mul(tempNum->re, tempNum->re, norm_diff_sqr);
              mpq_mul(tempNum->re, tempNum->re, tol_sqr);

              // verify (4) 
              if (mpq_cmp(q_sqr, tempNum->re) <= 0)
              { // not small enough
                allSmall = 0;
              }
            }
            else
            { // not small enough
              allSmall = 0;
            }
          }
          else
          { // not small enough
            allSmall = 0;
          }
        }
        else
        { // not small enough
          allSmall = 0;
        }
      }
    }

  if (!rV)
  { // check to see if all are small enough
    if (allSmall)
    { // all of them are small enough
      rV = 1;
    }
    else
    { // unknown
      rV = 0;
    }
  }

  // clear memory
  mpq_clear(beta);
  mpq_clear(four_beta_sqr); 
  mpq_clear(norm_diff_sqr);
  mpq_clear(tol);
  mpq_clear(tol_sqr);
  mpq_clear(q);
  mpq_clear(q_sqr);
  clear_rational_number(tempNum);

  return rV;
}

void determine_over_solution(int *isApproxSoln, point_struct *Pt, polynomial_system *F_rand, int numRandomSystems, int randomDigits, int eval_prec)
/***************************************************************\
* USAGE: determine if overdetermined solution                   *
\***************************************************************/
{ // we want the associated solutions to be within '10^-randomDigits'
  int i, rV = 0;
  point_struct *randPts = (point_struct *)errMalloc(numRandomSystems * sizeof(point_struct));

  // initialize
  *isApproxSoln = 1;
  for (i = 0; i < numRandomSystems; i++)
    initialize_point_struct(&randPts[i], F_rand[0].numVariables);

  // determine if approx soln
  for (i = 0; i < numRandomSystems && *isApproxSoln; i++)
  { // copy point
    copy_vector(randPts[i].origX, Pt->origX);
    copy_vector(randPts[i].x, Pt->x);

    // copy ||x||_2^2
    mpf_set(randPts[i].norm_x, Pt->norm_x);

    // compute alpha, beta, & gamma for the ith randomized system
    rV = compute_alpha_beta_gamma(randPts[i].Nx, randPts[i].alpha, randPts[i].beta, randPts[i].gamma, &F_rand[i], randPts[i].x, eval_prec);

    // save to original values
    set_number(randPts[i].origAlpha, randPts[i].alpha);
    set_number(randPts[i].origBeta, randPts[i].beta);
    set_number(randPts[i].origGamma, randPts[i].gamma);

    // check to see if we have successfully computed alpha_sqr, beta_sqr, & gamma_sqr
    if (rV == EXACT_SOLUTION_LU_ERROR)
    { // exact solution
      *isApproxSoln = 1;
    }
    else if (rV)
    { // error
      *isApproxSoln = 0; // unknown 
    }
    else
    { // determine if alpha is small enough to be an approximate solution
      *isApproxSoln = determine_approximate_solution(randPts[i].alpha);
    }
  }

  if (*isApproxSoln)
  { // determine if the associated solutions correspond to 10^-randomDigits or are certified not the same
    while ((rV = determine_over_associated_solution(randPts, numRandomSystems, randomDigits, eval_prec)) == 0)
    { // update
      for (i = 0; i < numRandomSystems; i++)
      { // copy and compute newton iteration
        copy_vector(randPts[i].x, randPts[i].Nx);
        rV = compute_alpha_beta_gamma(randPts[i].Nx, randPts[i].alpha, randPts[i].beta, randPts[i].gamma, &F_rand[i], randPts[i].x, eval_prec);

        // check for errors - should not occur!!
        if (rV == ERROR_LU_DECOMP)
        {
          printf("ERROR: Invalid LU decomposition!\n");
          errExit(ERROR_CONFIGURATION);
        }

        // update the norm
        norm_vector(randPts[i].norm_x, randPts[i].x);
      }
    }

    if (rV > 0)
    { // overdetermined solution
      *isApproxSoln = 1;
    }
    else
    { // not overdetermined solution
      *isApproxSoln = 0;
    }
  }

  // use randPts[0] to setup data in Pt
  copy_vector(Pt->x, randPts[0].x);
  copy_vector(Pt->Nx, randPts[0].Nx);
  mpf_set(Pt->norm_x, randPts[0].norm_x);
  set_number(Pt->origAlpha, randPts[0].origAlpha);
  set_number(Pt->origBeta, randPts[0].origBeta);
  set_number(Pt->origGamma, randPts[0].origGamma);
  set_number(Pt->alpha, randPts[0].alpha);
  set_number(Pt->beta, randPts[0].beta);
  set_number(Pt->gamma, randPts[0].gamma);

  // clear randPts
  for (i = 0; i < numRandomSystems; i++)
    clear_point_struct(&randPts[i]);
  free(randPts);

  return;
}

void determine_over_solution_rational(int *isApproxSoln, rational_point_struct *Pt, polynomial_system *F_rand, int numRandomSystems, int randomDigits)
/***************************************************************\
* USAGE: determine if overdetermined solution                   *
\***************************************************************/
{ // we want the associated solutions to be within '10^-randomDigits'
  int i, rV = 0;
  rational_point_struct *randPts = (rational_point_struct *)errMalloc(numRandomSystems * sizeof(rational_point_struct));

  // initialize
  *isApproxSoln = 1;
  for (i = 0; i < numRandomSystems; i++)
    initialize_rational_point_struct(&randPts[i], F_rand[0].numVariables);

  // determine if approx soln
  for (i = 0; i < numRandomSystems && *isApproxSoln; i++)
  { // copy rational point
    copy_rational_vector(randPts[i].origX, Pt->origX);
    copy_rational_vector(randPts[i].x, Pt->x);

    // copy ||x||_2^2
    mpq_set(randPts[i].norm_sqr_x, Pt->norm_sqr_x);

    // compute alpha_sqr, beta_sqr, & gamma_sqr for the ith randomized system
    rV = compute_alpha_beta_gamma_sqr_rational(randPts[i].Nx, randPts[i].alpha_sqr, randPts[i].beta_sqr, randPts[i].gamma_sqr, &F_rand[i], randPts[i].x);

    // save to original values
    set_rational_number(randPts[i].origAlpha_sqr, randPts[i].alpha_sqr);
    set_rational_number(randPts[i].origBeta_sqr, randPts[i].beta_sqr);
    set_rational_number(randPts[i].origGamma_sqr, randPts[i].gamma_sqr);

    // check to see if we have successfully computed alpha_sqr, beta_sqr, & gamma_sqr
    if (rV == EXACT_SOLUTION_LU_ERROR)
    { // exact solution
      *isApproxSoln = 1;
    }
    else if (rV)
    { // error
      *isApproxSoln = 0; // unknown 
    }
    else
    { // determine if alpha is small enough to be an approximate solution
      *isApproxSoln = determine_approximate_solution_rational(randPts[i].alpha_sqr);
    }
  }

  if (*isApproxSoln)
  { // determine if the associated solutions correspond to 10^-randomDigits or are certified not the same
    while ((rV = determine_over_associated_solution_rational(randPts, numRandomSystems, randomDigits)) == 0)
    { // update
      for (i = 0; i < numRandomSystems; i++)
      { // copy and compute newton iteration
        copy_rational_vector(randPts[i].x, randPts[i].Nx);
        rV = compute_alpha_beta_gamma_sqr_rational(randPts[i].Nx, randPts[i].alpha_sqr, randPts[i].beta_sqr, randPts[i].gamma_sqr, &F_rand[i], randPts[i].x);

        // check for errors - should not occur!!
        if (rV == ERROR_LU_DECOMP)
        {
          printf("ERROR: Invalid LU decomposition!\n");
          errExit(ERROR_CONFIGURATION);
        }

        // update the norm
        norm_sqr_rational_vector(randPts[i].norm_sqr_x, randPts[i].x);
      }
    }

    if (rV > 0)
    { // overdetermined solutions
      *isApproxSoln = 1;
    }
    else
    { // not overdetermined solution
      *isApproxSoln = 0;
    }
  }

  // use randPts[0] to setup data in Pt
  copy_rational_vector(Pt->x, randPts[0].x);
  copy_rational_vector(Pt->Nx, randPts[0].Nx);
  mpq_set(Pt->norm_sqr_x, randPts[0].norm_sqr_x);
  set_rational_number(Pt->origAlpha_sqr, randPts[0].origAlpha_sqr);
  set_rational_number(Pt->origBeta_sqr, randPts[0].origBeta_sqr);
  set_rational_number(Pt->origGamma_sqr, randPts[0].origGamma_sqr);
  set_rational_number(Pt->alpha_sqr, randPts[0].alpha_sqr);
  set_rational_number(Pt->beta_sqr, randPts[0].beta_sqr);
  set_rational_number(Pt->gamma_sqr, randPts[0].gamma_sqr);

  // clear randPts
  for (i = 0; i < numRandomSystems; i++)
    clear_rational_point_struct(&randPts[i]);
  free(randPts);

  return;
}
