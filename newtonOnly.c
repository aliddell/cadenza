/*
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   newtonOnly.c: Only perform newton iterations on the points
*/

#include "alphaCertified.h"

void create_newton_only_summary(int numPoints, point_struct *Pts, configurations *S, polynomial_system *F)
/***************************************************************\
* USAGE: create summary file                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  FILE *OUT = fopen("summary", "w");

  fprintf(OUT, "Floating point (%d bits) summary:\n\n", S->startingPrecision);

  // loop over the points
  for (i = 0; i < numPoints; i++)
  { // print the data for the ith point
    fprintf(OUT, "-------------------------\nPoint %d\n", i);
    print_vector_coordinate(OUT, 0, Pts[i].origX);

    // print original values of alpha, beta, and gamma
    fprintf(OUT, "Original values:\n");
    fprintf(OUT, " alpha < "); mpf_out_str(OUT, base, digits, Pts[i].origAlpha->re); fprintf(OUT, "\n");
    fprintf(OUT, " beta ~= "); mpf_out_str(OUT, base, digits, Pts[i].origBeta->re); fprintf(OUT, "\n");
    fprintf(OUT, " gamma < "); mpf_out_str(OUT, base, digits, Pts[i].origGamma->re); fprintf(OUT, "\n");

    // print final values of alpha, beta, and gamma
    fprintf(OUT, "Final values:\n");
    fprintf(OUT, " alpha < "); mpf_out_str(OUT, base, digits, Pts[i].alpha->re); fprintf(OUT, "\n");
    fprintf(OUT, " beta ~= "); mpf_out_str(OUT, base, digits, Pts[i].beta->re); fprintf(OUT, "\n");
    fprintf(OUT, " gamma < "); mpf_out_str(OUT, base, digits, Pts[i].gamma->re); fprintf(OUT, "\n");
  }
  fprintf(OUT, "\n");

  // print configurations and message about alphaCertified
  configuration_summary(OUT, S, F);

  fclose(OUT);

  return;
}

void newton_only_output(int numPoints, point_struct *Pts, configurations *S, polynomial_system *F)
/***************************************************************\
* USAGE: create output files                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  FILE *F1 = NULL;

  // create the best approximations file
  F1 = fopen("refinedPoints", "w");
  fprintf(F1, "%d\n\n", numPoints);
  for (i = 0; i < numPoints; i++)
  { // print to F1
    print_vector_coordinate(F1, 0, Pts[i].x);
    fprintf(F1, "\n");
  }
  fclose(F1);

  // create certified and unknown files
  F1 = fopen("constantValues", "w");
  fprintf(F1, "%d\n\n", numPoints);
  for (i = 0; i < numPoints; i++)
  { // print alpha, beta & gamma for origSoln
    mpf_out_str(F1, base, digits, Pts[i].origAlpha->re);
    fprintf(F1, "\n");
    mpf_out_str(F1, base, digits, Pts[i].origBeta->re);
    fprintf(F1, "\n");
    mpf_out_str(F1, base, digits, Pts[i].origGamma->re);
    fprintf(F1, "\n\n");
  }
  fclose(F1);

  // create summary file
  create_newton_only_summary(numPoints, Pts, S, F);

  // print to screen the files that were created
  printf("\n------------------------------------------------------------------------------------------------------------\n");
  printf("The following files have been created:\n\n");
  printf("constantValues:       A list of the values of alpha, beta, and gamma for the points.\n");
  printf("refinedPoints:        A list of points after the requested Newton iterations.\n");
  printf("summary:              A human-readable summary for each point - main output file.\n");
  printf("------------------------------------------------------------------------------------------------------------\n");

  return;
}

void newton_only(int numPoints, complex_vector *Points, polynomial_system *F, configurations *S)
/***************************************************************\
* USAGE: perform newton iterations only on the points           *
\***************************************************************/
{
  int i, j, rV, curr_prec, numVars = F->numVariables;
  point_struct *Points_struct = (point_struct *)errMalloc(numPoints * sizeof(point_struct));

  // setup Points_struct and perform the newton iterations
  for (i = 0; i < numPoints; i++)
  { // initialize
    curr_prec = S->startingPrecision;
    setPrec(curr_prec);
    initialize_point_struct(&Points_struct[i], numVars);

    // copy point
    copy_vector(Points_struct[i].origX, Points[i]);
    copy_vector(Points_struct[i].x, Points[i]);

    // compute alpha, beta, & gamma (and save to original values of alpha, beta, & gamma)
    rV = compute_alpha_beta_gamma(Points_struct[i].Nx, Points_struct[i].alpha, Points_struct[i].beta, Points_struct[i].gamma, F, Points_struct[i].x, S->startingPrecision);
    set_number(Points_struct[i].origAlpha, Points_struct[i].alpha);
    set_number(Points_struct[i].origBeta, Points_struct[i].beta);
    set_number(Points_struct[i].origGamma, Points_struct[i].gamma);

    // perform the requested number of Newton iterations
    for (j = 0; j < S->newtonIts && !rV; j++)
    { // set x to Nx
      curr_prec *= 2;
      setPrec_vector(Points_struct[i].x, curr_prec);
      copy_vector(Points_struct[i].x, Points_struct[i].Nx);

      // compute alpha, beta & gamma for x
      rV = compute_alpha_beta_gamma(Points_struct[i].Nx, Points_struct[i].alpha, Points_struct[i].beta, Points_struct[i].gamma, F, Points_struct[i].x, curr_prec);
    }
  }

  // print the data out
  newton_only_output(numPoints, Points_struct, S, F);

  // clear Points_struct
  for (i = 0; i < numPoints; i++)
    clear_point_struct(&Points_struct[i]);
  free(Points_struct);
  Points_struct = NULL;

  return;
}

void create_newton_only_summary_rational(int numPoints, rational_point_struct *Pts, configurations *S, polynomial_system *F)
/***************************************************************\
* USAGE: create summary file                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  mpf_t tempMPF;
  FILE *OUT = fopen("summary", "w");

  mpf_init2(tempMPF, 1024);

  fprintf(OUT, "Rational summary:\n\n");

  // loop over the points
  for (i = 0; i < numPoints; i++)
  { // print the data for the ith point
    fprintf(OUT, "-------------------------\nPoint %d\n", i);
    print_rational_vector_coordinate(OUT, Pts[i].origX);

    // print original approximate values of alpha, beta, and gamma
    fprintf(OUT, "Original values:\n");
    mpf_set_q(tempMPF, Pts[i].origAlpha_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    fprintf(OUT, " alpha < "); mpf_out_str(OUT, base, digits, tempMPF); fprintf(OUT, "\n");
    mpf_set_q(tempMPF, Pts[i].origBeta_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    fprintf(OUT, " beta ~= "); mpf_out_str(OUT, base, digits, tempMPF); fprintf(OUT, "\n");
    mpf_set_q(tempMPF, Pts[i].origGamma_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    fprintf(OUT, " gamma < "); mpf_out_str(OUT, base, digits, tempMPF); fprintf(OUT, "\n");


    // print approximate values of alpha, beta, and gamma
    fprintf(OUT, "Final values:\n");
    mpf_set_q(tempMPF, Pts[i].alpha_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    fprintf(OUT, " alpha < "); mpf_out_str(OUT, base, digits, tempMPF); fprintf(OUT, "\n");
    mpf_set_q(tempMPF, Pts[i].beta_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    fprintf(OUT, " beta ~= "); mpf_out_str(OUT, base, digits, tempMPF); fprintf(OUT, "\n");
    mpf_set_q(tempMPF, Pts[i].gamma_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    fprintf(OUT, " gamma < "); mpf_out_str(OUT, base, digits, tempMPF); fprintf(OUT, "\n");
  }

  fprintf(OUT, "\n");

  // print configurations and message about alphaCertified
  configuration_summary(OUT, S, F);

  fclose(OUT);
  mpf_clear(tempMPF);

  return;
}

void newton_only_output_rational(int numPoints, rational_point_struct *Pts, configurations *S, polynomial_system *F)
/***************************************************************\
* USAGE: create output files                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  mpf_t tempMPF;
  FILE *F1 = NULL;

  mpf_init2(tempMPF, 1024);

  // create the best approximations file
  F1 = fopen("refinedPoints", "w");
  fprintf(F1, "%d\n\n", numPoints);
  for (i = 0; i < numPoints; i++)
  { // print to F1
    print_rational_vector_coordinate(F1, Pts[i].x);
    fprintf(F1, "\n");
  }
  fclose(F1);

  // create certified and unknown files
  F1 = fopen("constantValues", "w");
  fprintf(F1, "%d\n\n", numPoints);
  for (i = 0; i < numPoints; i++)
  { // print alpha, beta & gamma for origSoln
    mpf_set_q(tempMPF, Pts[i].origAlpha_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    mpf_out_str(F1, base, digits, tempMPF);
    fprintf(F1, "\n");
    mpf_set_q(tempMPF, Pts[i].origBeta_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    mpf_out_str(F1, base, digits, tempMPF);
    fprintf(F1, "\n");
    mpf_set_q(tempMPF, Pts[i].origGamma_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    mpf_out_str(F1, base, digits, tempMPF);
    fprintf(F1, "\n\n");
  }
  fclose(F1);

  // create summary file
  create_newton_only_summary_rational(numPoints, Pts, S, F);

  // print to screen the files that were created
  printf("\n------------------------------------------------------------------------------------------------------------\n");
  printf("The following files have been created:\n\n");
  printf("constantValues:       A list of the values of alpha, beta, and gamma for the points.\n");
  printf("refinedPoints:        A list of points after the requested Newton iterations.\n");
  printf("summary:              A human-readable summary for each point - main output file.\n");
  printf("------------------------------------------------------------------------------------------------------------\n");

  mpf_clear(tempMPF);

  return;
}

void newton_only_rational(int numPoints, rational_complex_vector *Points, polynomial_system *F, configurations *S)
/***************************************************************\
* USAGE: perform newton iterations only on the points           *
\***************************************************************/
{
  int i, j, rV, numVars = F->numVariables;
  rational_point_struct *Points_struct = (rational_point_struct *)errMalloc(numPoints * sizeof(rational_point_struct));

  // setup Points_struct and perform the newton iterations
  for (i = 0; i < numPoints; i++)
  { // initialize
    initialize_rational_point_struct(&Points_struct[i], numVars);

    // copy point
    copy_rational_vector(Points_struct[i].origX, Points[i]);
    copy_rational_vector(Points_struct[i].x, Points[i]);

    // compute alpha, beta, & gamma (and save to original values of alpha, beta, & gamma)
    rV = compute_alpha_beta_gamma_sqr_rational(Points_struct[i].Nx, Points_struct[i].alpha_sqr, Points_struct[i].beta_sqr, Points_struct[i].gamma_sqr, F, Points_struct[i].x);
    set_rational_number(Points_struct[i].origAlpha_sqr, Points_struct[i].alpha_sqr);
    set_rational_number(Points_struct[i].origBeta_sqr, Points_struct[i].beta_sqr);
    set_rational_number(Points_struct[i].origGamma_sqr, Points_struct[i].gamma_sqr);

    // perform the requested number of Newton iterations
    for (j = 0; j < S->newtonIts && !rV; j++)
    { // set x to Nx
      copy_rational_vector(Points_struct[i].x, Points_struct[i].Nx);

      // compute alpha, beta & gamma for x
      rV = compute_alpha_beta_gamma_sqr_rational(Points_struct[i].Nx, Points_struct[i].alpha_sqr, Points_struct[i].beta_sqr, Points_struct[i].gamma_sqr, F, Points_struct[i].x);
    }
  }

  // print the data out
  newton_only_output_rational(numPoints, Points_struct, S, F);

  // clear Points_struct
  for (i = 0; i < numPoints; i++)
    clear_rational_point_struct(&Points_struct[i]);
  free(Points_struct);
  Points_struct = NULL;

  return;
}


