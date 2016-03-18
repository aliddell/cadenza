/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   output.c: Create output - tables to screen and data to files
*/

#include "alphaCertified.h"

void print_configurations(FILE *OUT, configurations *S)
/***************************************************************\
* USAGE: prints summary of configurations                       *
\***************************************************************/
{
  configurations S_default;
  load_default_settings(&S_default);

  if (S_default.algorithm != S->algorithm)
    fprintf(OUT, "ALGORITHM: %d;\n", S->algorithm);
  if (S_default.arithmeticType != S->arithmeticType)
    fprintf(OUT, "ARITHMETICTYPE: %d;\n", S->arithmeticType);
  if (S->arithmeticType == 1 && S_default.startingPrecision != S->startingPrecision)
    fprintf(OUT, "PRECISION: %d;\n", S->startingPrecision);
  if (S_default.refineDigits != S->refineDigits)
    fprintf(OUT, "REFINEDIGITS: %d;\n", S->refineDigits);
  if (S_default.numRandomSystems != S->numRandomSystems)
    fprintf(OUT, "NUMRANDOMSYSTEMS: %d;\n", S->numRandomSystems);
  if (S_default.randomDigits != S->randomDigits)
    fprintf(OUT, "RANDOMDIGITS: %d;\n", S->randomDigits);
  if (S_default.newtonOnly != S->newtonOnly)
    fprintf(OUT, "NEWTONONLY: %d;\n", S->newtonOnly);
  if (S_default.newtonIts != S->newtonIts)
    fprintf(OUT, "NUMITERATIONS: %d;\n", S->newtonIts);
  if (S_default.realityCheck != S->realityCheck)
    fprintf(OUT, "REALITYCHECK: %d;\n", S->realityCheck);
  if (S_default.realityTest != S->realityTest)
    fprintf(OUT, "REALITYTEST: %d;\n", S->realityTest);
  fprintf(OUT, "RANDOMSEED: %u;\n\n", S->randomSeed);

  return;
}

void configuration_summary(FILE *OUT, configurations *S, polynomial_system *F)
/***************************************************************\
* USAGE: prints summary of configurations and version           *
\***************************************************************/
{
  fprintf(OUT, "*************** Configurations ****************\n\n");
  print_configurations(OUT, S);

  fprintf(OUT, "************* Version information *************\n");
  print_welcome_message(OUT);

  fprintf(OUT, "******************* System ********************\n");
  print_polynomial_system(OUT, F);

  return;
}

void display_output_files(int isReal, int algorithm)
/***************************************************************\
* USAGE: displays the output files created                      *
\***************************************************************/
{
  // print to screen the files that were created
  printf("\n------------------------------------------------------------------------------------------------------------\n");
  printf("The following files have been created:\n\n");
  printf("approxSolns:          A list of points that are certified approximate solutions.\n");
  printf("constantValues:       A list of the values of alpha, beta, and gamma for the points.\n");
  if (algorithm >= 1)
    printf("distinctSolns:        A list of points that correspond to distinct solutions.\n");
  printf("isApproxSoln:         A list which describes if the ith point is an approximate solution.\n");
  if (algorithm >= 1)
    printf("isDistinctSoln:       A list which describes if the ith point is listed in 'distinctSolns'.\n");
  if (algorithm >= 2 &&isReal)
  {
    printf("isRealSoln:           A list which describes if the ith point corresponds to a real solution.\n");
    printf("nonrealDistinctSolns: A list of points that correspond to distinct nonreal solutions.\n");
    printf("realDistinctSolns:    A list of points that correspond to distinct real solutions.\n");
  }
  if (algorithm >= 1)
    printf("redundantSolns:       A list of points that correspond to the same solution as one in 'distinctSolns'.\n");
  printf("refinedPoints:        A list of points that are the best internally computed approximation of each solution.\n");
  printf("summary:              A human-readable summary for each point - main output file.\n");
  printf("unknownPoints:        A list of points that which cannot be certified as approximate solutions.\n");
  printf("------------------------------------------------------------------------------------------------------------\n");

  return;
}

void display_summary_table(FILE *OUT, int numPoints, int numApproxSolns, int numDistinctSolns, int numRealSolns, int isReal, int algorithm, int realityCheck, int realityTest)
/***************************************************************\
* USAGE: create summary table                                   *
\***************************************************************/
{
  fprintf(OUT, "Number of points tested:           %d\n", numPoints);
  fprintf(OUT, "Certified approximate solutions:   %d\n", numApproxSolns);
  if (algorithm >= 1)
    fprintf(OUT, "Certified distinct solutions:      %d\n", numDistinctSolns);
  if (algorithm >= 2 && isReal)
  {
    fprintf(OUT, "Certified real distinct solutions: %d", numRealSolns);
    if (realityTest)
    { // describe what was actually computed
      fprintf(OUT, "**\n\n");
      fprintf(OUT, "** alphaCertified has found that the conjugate\n   of each of these points do not correspond to\n   the same solution as any other approximate\n   solution (REALITYTEST: %d).\n", realityTest);
    }
    else if (realityCheck == -1)
    { // describe what was actually computed
      fprintf(OUT, "** \n\n");
      fprintf(OUT, "** alphaCertified has found a real point which\n   is an approximate solution corresponding\n   to each of these solutions (REALITYCHECK: %d).\n", realityCheck);
    }
    else
      fprintf(OUT, "\n");
  }

  return;
}

void create_summary(int numPoints, point_struct *Pts, int numApproxSolns, int numDistinctSolns, int numRealSolns, int isReal, configurations *S, int overDet, polynomial_system *F)
/***************************************************************\
* USAGE: create summary file                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  FILE *OUT = fopen("summary", "w");

  if (overDet)
  {
    fprintf(OUT, "Heuristic floating point (%d bits) soft certification results for an overdetermined system.\n", S->startingPrecision);
    fprintf(OUT, "(alphaCertified used %d random systems where each was required to have\na solution in a ball of radius 1e-%d)\n\n", S->numRandomSystems, S->randomDigits);
  }
  else
  { // square system
    fprintf(OUT, "Floating point (%d bits) soft certification results:\n\n", S->startingPrecision);
  }
  display_summary_table(OUT, numPoints, numApproxSolns, numDistinctSolns, numRealSolns, isReal, S->algorithm, S->realityCheck, S->realityTest);

  // loop over the points
  for (i = 0; i < numPoints; i++)
  { // print the data for the ith point
    fprintf(OUT, "-------------------------\nPoint %d\n", i);
    print_vector_coordinate(OUT, 0, Pts[i].origX);

    // results
    fprintf(OUT, "Approx solution: %s\n", Pts[i].isApproxSoln ? "Yes" : "No");
    if (S->algorithm >= 1 && Pts[i].isApproxSoln)
    { // determine if distinct from above solutions
      fprintf(OUT, "Distinct from the solutions above: %s", Pts[i].isActive == 1 ? "Yes" : "No");
      if (Pts[i].isActive <= 0)
        fprintf(OUT, " (Point %d)", -Pts[i].isActive);
      fprintf(OUT, "\n");
     
      if (S->algorithm >= 2 && isReal)
      {
        if (Pts[i].isActive == 1)
          fprintf(OUT, "Real solution: %s\n", Pts[i].isReal ? "Yes" : "No");
        else
          fprintf(OUT, "Real solution: %s\n", Pts[-Pts[i].isActive].isReal ? "Yes" : "No");
      }
    }

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

  // print configurations, message about alphaCertified, and F
  configuration_summary(OUT, S, F);

  fclose(OUT);

  return;
}

void create_rational_summary(int numPoints, rational_point_struct *Pts, int numApproxSolns, int numDistinctSolns, int numRealSolns, int isReal, configurations *S, int overDet, polynomial_system *F)
/***************************************************************\
* USAGE: create summary file                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  mpf_t tempMPF;
  FILE *OUT = fopen("summary", "w");

  mpf_init2(tempMPF, 1024);

  if (overDet)
  {
    fprintf(OUT, "Heuristic rational certification results for an overdetermined system.\n");
    fprintf(OUT, "(alphaCertified used %d random systems where each was required to have\na solution in a ball of radius 1e-%d)\n\n", S->numRandomSystems, S->randomDigits);
  }
  else
  { // square system
    fprintf(OUT, "Rational certification results:\n\n");
  }
  display_summary_table(OUT, numPoints, numApproxSolns, numDistinctSolns, numRealSolns, isReal, S->algorithm, S->realityCheck, S->realityTest);

  // loop over the points
  for (i = 0; i < numPoints; i++)
  { // print the data for the ith point
    fprintf(OUT, "-------------------------\nPoint %d\n", i);
    print_rational_vector_coordinate(OUT, Pts[i].origX);

    // results
    fprintf(OUT, "Approx solution: %s\n", Pts[i].isApproxSoln ? "Yes" : "No");
    if (S->algorithm >= 1 && Pts[i].isApproxSoln)
    { // determine if distinct from above solutions
      fprintf(OUT, "Distinct from the solutions above: %s", Pts[i].isActive == 1 ? "Yes" : "No");
      if (Pts[i].isActive <= 0)
        fprintf(OUT, " (Point %d)", -Pts[i].isActive);
      fprintf(OUT, "\n");
     
      if (S->algorithm >= 2 && isReal)
      {
        if (Pts[i].isActive == 1)
          fprintf(OUT, "Real solution: %s\n", Pts[i].isReal ? "Yes" : "No");
        else
          fprintf(OUT, "Real solution: %s\n", Pts[-Pts[i].isActive].isReal ? "Yes" : "No");
      }
    }

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

void classify_output(int numPoints, point_struct *Pts, int numApproxSolns, int numDistinctSolns, int numRealSolns, int isReal, configurations *S, int overDet, polynomial_system *F)
/***************************************************************\
* USAGE: create output files                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  FILE *F1 = NULL, *F2 = NULL, *F3 = NULL, *F4 = NULL;

  // print summary
  if (overDet)
  {
    printf("Heuristic floating point (%d bits) soft certification results for an overdetermined system.\n", S->startingPrecision);
    printf("(alphaCertified used %d random systems where each was required to have\na solution in a ball of radius 1e-%d)\n\n", S->numRandomSystems, S->randomDigits);
  }
  else
  { // square system
    printf("Floating point (%d bits) soft certification results:\n\n", S->startingPrecision);
  }

  display_summary_table(stdout, numPoints, numApproxSolns, numDistinctSolns, numRealSolns, isReal, S->algorithm, S->realityCheck, S->realityTest);

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
  F1 = fopen("approxSolns", "w");
  fprintf(F1, "%d\n\n", numApproxSolns);
  F2 = fopen("unknownPoints", "w");
  fprintf(F2, "%d\n\n", numPoints - numApproxSolns);
  F3 = fopen("isApproxSoln", "w");
  fprintf(F3, "%d\n\n", numPoints);
  F4 = fopen("constantValues", "w");
  fprintf(F4, "%d\n\n", numPoints);
  for (i = 0; i < numPoints; i++)
  { // print isApproxSoln
    fprintf(F3, "%d\n\n", Pts[i].isApproxSoln);

    // print alpha, beta & gamma for origSoln
    mpf_out_str(F4, base, digits, Pts[i].origAlpha->re); 
    fprintf(F4, "\n");
    mpf_out_str(F4, base, digits, Pts[i].origBeta->re); 
    fprintf(F4, "\n");
    mpf_out_str(F4, base, digits, Pts[i].origGamma->re); 
    fprintf(F4, "\n\n");

    if (Pts[i].isApproxSoln)
    { // print to F1
      print_vector_coordinate(F1, 0, Pts[i].origX);
      fprintf(F1, "\n");
    }
    else
    { // print to F2
      print_vector_coordinate(F2, 0, Pts[i].origX);
      fprintf(F2, "\n");
    }
  }
  fclose(F1);
  fclose(F2);
  fclose(F3);
  fclose(F4);

  if (S->algorithm >= 1)
  { // create distinct and redundant files
    F1 = fopen("distinctSolns", "w");
    fprintf(F1, "%d\n\n", numDistinctSolns);
    F2 = fopen("redundantSolns", "w");
    fprintf(F2, "%d\n\n", numApproxSolns - numDistinctSolns);
    F3 = fopen("isDistinctSoln", "w");
    fprintf(F3, "%d\n\n", numPoints);
    for (i = 0; i < numPoints; i++)
      if (!Pts[i].isApproxSoln)
      { // print -2 - unknown
        fprintf(F3, "-2\n\n");
      }
      else 
      { // certified solution
        if (Pts[i].isActive == 1)
        { // print to F1
          print_vector_coordinate(F1, 0, Pts[i].origX);
          fprintf(F1, "\n");
          // print to F3 - certified distinct among list of other distinct ones
          fprintf(F3, "-1\n\n");
        }
        else
        { // print to F2
          print_vector_coordinate(F2, 0, Pts[i].origX);
          fprintf(F2, "\n");
          // print to F3 - print number of point that it corresponds with
          fprintf(F3, "%d\n\n", -Pts[i].isActive);
        }
      }
    fclose(F1);
    fclose(F2);
    fclose(F3);
  }

  // create real distinct and nonreal distinct files, if needed
  if (S->algorithm >= 2 && isReal)
  {
    F1 = fopen("realDistinctSolns", "w");
    fprintf(F1, "%d\n\n", numRealSolns);
    F2 = fopen("nonrealDistinctSolns", "w");
    fprintf(F2, "%d\n\n", numDistinctSolns - numRealSolns);
    F3 = fopen("isRealSoln", "w");
    fprintf(F3, "%d\n\n", numPoints);
    for (i = 0; i < numPoints; i++)
      if (!Pts[i].isApproxSoln)
      { // print -2 - unknown
        fprintf(F3, "-2\n\n");
      }
      else if (Pts[i].isActive <= 0)
      { // print if it corresponds to a real solution
        fprintf(F3, "%d\n\n", Pts[-Pts[i].isActive].isReal);
      }
      else
      { // certified distinct solution
        fprintf(F3, "%d\n\n", Pts[i].isReal);
        if (Pts[i].isReal)
        { // print to F1
          print_vector_coordinate(F1, 0, Pts[i].origX);
          fprintf(F1, "\n");
        }
        else if (Pts[i].isApproxSoln && Pts[i].isActive == 1)
        { // print to F2
          print_vector_coordinate(F2, 0, Pts[i].origX);
          fprintf(F2, "\n");
        }
      }
    fclose(F1);
    fclose(F2);
    fclose(F3);
  }

  // create summary file
  create_summary(numPoints, Pts, numApproxSolns, numDistinctSolns, numRealSolns, isReal, S, overDet, F);

  // print to screen the files that were created
  display_output_files(isReal, S->algorithm);

  return;
}

void classify_rational_output(int numPoints, rational_point_struct *Pts, int numApproxSolns, int numDistinctSolns, int numRealSolns, int isReal, configurations *S, int overDet, polynomial_system *F)
/***************************************************************\
* USAGE: create output files                                    *
\***************************************************************/
{
  int i, base = 10, digits = 16;
  mpf_t tempMPF;
  FILE *F1 = NULL, *F2 = NULL, *F3 = NULL, *F4 = NULL;

  mpf_init2(tempMPF, 1024);

  // print summary
  if (overDet)
  {
    printf("Heuristic rational certification results for an overdetermined system.\n");
    printf("(alphaCertified used %d random systems where each was required to have\na solution in a ball of radius 1e-%d)\n\n", S->numRandomSystems, S->randomDigits);
  }
  else
  { // square system
    printf("Rational certification results:\n\n");
  }
  display_summary_table(stdout, numPoints, numApproxSolns, numDistinctSolns, numRealSolns, isReal, S->algorithm, S->realityCheck, S->realityTest);

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
  F1 = fopen("approxSolns", "w");
  fprintf(F1, "%d\n\n", numApproxSolns);
  F2 = fopen("unknownPoints", "w");
  fprintf(F2, "%d\n\n", numPoints - numApproxSolns);
  F3 = fopen("isApproxSoln", "w");
  fprintf(F3, "%d\n\n", numPoints);
  F4 = fopen("constantValues", "w");
  fprintf(F4, "%d\n\n", numPoints);
  for (i = 0; i < numPoints; i++)
  { // print isApproxSoln
    fprintf(F3, "%d\n\n", Pts[i].isApproxSoln);

    mpf_set_q(tempMPF, Pts[i].origAlpha_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    mpf_out_str(F4, base, digits, tempMPF); 
    fprintf(F4, "\n");
    mpf_set_q(tempMPF, Pts[i].origBeta_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    mpf_out_str(F4, base, digits, tempMPF); 
    fprintf(F4, "\n");
    mpf_set_q(tempMPF, Pts[i].origGamma_sqr->re);
    mpf_sqrt(tempMPF, tempMPF);
    mpf_out_str(F4, base, digits, tempMPF); 
    fprintf(F4, "\n\n");

    if (Pts[i].isApproxSoln)
    { // print to F1
      print_rational_vector_coordinate(F1, Pts[i].origX);
      fprintf(F1, "\n");
    }
    else
    { // print to F2
      print_rational_vector_coordinate(F2, Pts[i].origX);
      fprintf(F2, "\n");
    }
  }
  fclose(F1);
  fclose(F2);
  fclose(F3);
  fclose(F4);

  if (S->algorithm >= 1)
  { // create distinct and redundant files
    F1 = fopen("distinctSolns", "w");
    fprintf(F1, "%d\n\n", numDistinctSolns);
    F2 = fopen("redundantSolns", "w");
    fprintf(F2, "%d\n\n", numApproxSolns - numDistinctSolns);
    F3 = fopen("isDistinctSoln", "w");
    fprintf(F3, "%d\n\n", numPoints);
    for (i = 0; i < numPoints; i++)
      if (!Pts[i].isApproxSoln)
      { // print -2 - unknown
        fprintf(F3, "-2\n\n");
      }
      else 
      { // certified solution
        if (Pts[i].isActive == 1)
        { // print to F1
          print_rational_vector_coordinate(F1, Pts[i].origX);
          fprintf(F1, "\n");
          // print to F3 - certified distinct among list of other distinct ones
          fprintf(F3, "-1\n\n");
        }
        else
        { // print to F2
          print_rational_vector_coordinate(F2, Pts[i].origX);
          fprintf(F2, "\n");
          // print to F3 - print number of point that it corresponds with
          fprintf(F3, "%d\n\n", -Pts[i].isActive);
        }
      }
    fclose(F1);
    fclose(F2);
    fclose(F3);
  }

  // create real distinct and nonreal distinct files, if needed
  if (S->algorithm >= 2 && isReal)
  {
    F1 = fopen("realDistinctSolns", "w");
    fprintf(F1, "%d\n\n", numRealSolns);
    F2 = fopen("nonrealDistinctSolns", "w");
    fprintf(F2, "%d\n\n", numDistinctSolns - numRealSolns);
    F3 = fopen("isRealSoln", "w");
    fprintf(F3, "%d\n\n", numPoints);
    for (i = 0; i < numPoints; i++)
      if (!Pts[i].isApproxSoln)
      { // print -2 - unknown
        fprintf(F3, "-2\n\n");
      }
      else if (Pts[i].isActive <= 0)
      { // print if it corresponds to a real solution
        fprintf(F3, "%d\n\n", Pts[-Pts[i].isActive].isReal);
      }
      else
      { // certified distinct solution
        fprintf(F3, "%d\n\n", Pts[i].isReal);
        if (Pts[i].isReal)
        { // print to F1
          print_rational_vector_coordinate(F1, Pts[i].origX);
          fprintf(F1, "\n");
        }
        else if (Pts[i].isApproxSoln && Pts[i].isActive == 1)
        { // print to F2
          print_rational_vector_coordinate(F2, Pts[i].origX);
          fprintf(F2, "\n");
        }
      }
    fclose(F1);
    fclose(F2);
    fclose(F3);
  }

  // create summary file
  create_rational_summary(numPoints, Pts, numApproxSolns, numDistinctSolns, numRealSolns, isReal, S, overDet, F);

  // print to screen the files that were created
  display_output_files(isReal, S->algorithm);

  mpf_clear(tempMPF);

  return;
}



