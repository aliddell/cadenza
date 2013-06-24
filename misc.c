/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   misc.c: Contains miscellaneous functions used throughout the program
*/

#include "alphaCertified.h"
#include <ctype.h>

void setPrec(int prec)
/***************************************************************\
* USAGE: set default precision                                  *
\***************************************************************/
{
  mpf_set_default_prec(prec);

  return;
}

void errExit(int errorCode)
/***************************************************************\
* USAGE: exit alphaCertified                                    *
\***************************************************************/
{
  if (errorCode == 0)
    errorCode = ERROR_OTHER;

  printf("%s\n", ERROR_MESSAGE);

  exit(errorCode);

  return;
}

void *errMalloc(size_t size)
/***************************************************************\
* USAGE: allocate memory with error checking                    *
\***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate
    return NULL;
  }
  else
  { // try to allocate memory
    void *x = malloc(size);
    if (x == NULL)
    {
      printf("ERROR: errMalloc was unable to allocate memory (%d)!\n", (int) size);
      errExit(ERROR_MEMORY_ALLOCATION);
    }
    return x;
  }
}

void *errRealloc(void *ptr, size_t size)
/***************************************************************\
* USAGE: reallocate memory with error checking                  *
\***************************************************************/
{
  if (size <= 0)
  { // nothing to reallocate - free memory and return NULL
    free(ptr);
    ptr = NULL;
  }
  else
  { // try to reallocate memory
    ptr = realloc(ptr, size);
    if (ptr == NULL)
    {
      printf("ERROR: errRealloc was unable to allocate memory (%d)!\n", (int) size);
      errExit(ERROR_MEMORY_ALLOCATION);
    }
  }
  return ptr;
}

void setup_polynomial(polynomial *F, int numVars, FILE *IN, int polyNumber)
/***************************************************************\
* USAGE: setup the next polynomial described in IN              *
\***************************************************************/
{
  int i, j, rV, max, base = 10;

  // initialize the number of variables, degree & isReal
  F->numVariables = numVars;
  F->degree = 0;
  F->isReal = 1;

  // read in the number of terms
  fscanf(IN, "%d\n", &F->numTerms);

  // error checking - want number of terms >= 0
  if (F->numTerms <= 0)
  { // error
    printf("ERROR: The number of terms (%d) must be positive.\n", F->numTerms);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  }

  // allocate memory
  mpq_init(F->norm_sqr);
  F->coeff = (rational_complex_number *)errMalloc(F->numTerms * sizeof(rational_complex_number));
  F->exponents = (int **)errMalloc(F->numTerms * sizeof(int *));

  // setup the terms
  for (i = 0; i < F->numTerms; i++)
  { // allocate & initialize memory
    F->exponents[i] = (int *)errMalloc(numVars * sizeof(int));
    initialize_rational_number(F->coeff[i]);

    // read in exponents, compute degree, and perform error checking - want exponents >= 0
    max = 0;
    for (j = 0; j < numVars; j++)
    {
      fscanf(IN, "%d", &F->exponents[i][j]);
      if (F->exponents[i][j] < 0)
      { // error
        printf("ERROR: The exponent for variable %d in monomial %d of polynomial %d must be nonnegative (%d).\n", j+1, i+1, polyNumber+1, F->exponents[i][j]);
        // close file and exit
        fclose(IN);
        errExit(ERROR_INPUT_SYSTEM);
      }
      max += F->exponents[i][j];
    }

    // update degree, if needed
    if (max > F->degree)
      F->degree = max;

    // read in real & imaginary part of coefficient
    rV = mpq_inp_str(F->coeff[i]->re, IN, base);
    if (rV == 0)
    { // error in reading the coefficient
      printf("ERROR: There appears to be an error when reading in the real part of the\n       coefficient for monomial %d of polynomial %d.\n", i+1, polyNumber+1);
      // close file and exit
      fclose(IN);
      errExit(ERROR_INPUT_SYSTEM);
    } 
    rV = mpq_inp_str(F->coeff[i]->im, IN, base);
    if (rV == 0)
    { // error in reading the coefficient
      printf("ERROR: There appears to be an error when reading in the imaginary part of the\n       coefficient for monomial %d of polynomial %d.\n", i+1, polyNumber+1);
      // close file and exit
      fclose(IN);
      errExit(ERROR_INPUT_SYSTEM);
    } 

    // setup in canonical form
    mpq_canonicalize(F->coeff[i]->re);
    mpq_canonicalize(F->coeff[i]->im);

    // update isReal, if needed
    if (F->isReal && mpq_cmp_ui(F->coeff[i]->im, 0, 1) != 0)
      F->isReal = 0;
  }

  // compute norm_sqr
  norm_sqr_polynomial(F->norm_sqr, F);
  
  return;
}

int readInExpFunction(FILE *IN, char *expFunction, int *isHyperbolic)
/***************************************************************\
* USAGE: read in exponential function - skip white space        *
\***************************************************************/
{
  int rV = 0;  // 0 - error, 1 - read in a char
  char c;

  do
  { // read in char
    c = fgetc(IN);
  } while (c != EOF && isspace((int) c));

  if (c == EOF)
  { // return error
    rV = 0;
  }
  else
  { // save char
    rV = 1;
    *expFunction = c;

    // initialize isHyperbolic
    *isHyperbolic = 0;

    // read in next char and determine what to do
    c = fgetc(IN);
    if (c == 'H' || c == 'h')
    { // using hyperbolic
      *isHyperbolic = 1;
    }
    else
    { // put char back
      ungetc(c, IN);
    }
  }

  return rV;
}

void setup_exponential(exponential *F, int numVars, FILE *IN, int expNumber, int yIndex)
/***************************************************************\
* USAGE: setup the next exponential described in IN             *
\***************************************************************/
{
  int rV, base = 10;

  // read in the x variable index
  rV = fscanf(IN, "%d", &F->xIndex);
  if (rV == 0)
  { // error in reading the x variable index
    printf("ERROR: There appears to be an error when reading in the variable index\n       for exponential %d.\n", expNumber+1);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  } 

  // error checking - want variable index to be in 1,2,..,numVars
  if (F->xIndex < 1 || F->xIndex > numVars)
  { // error
    printf("ERROR: The variable index for exponential %d must be between 1 and %d.\n", expNumber+1, numVars);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  }

  // setup x & y variable index
  F->xIndex--;
  F->yIndex = yIndex;

  // read in the type of exponential function: exp, sin, or cos
  rV = readInExpFunction(IN, &F->expFunction, &F->isHyperbolic);
  if (rV == 0)
  { // error in reading the exponential function
    printf("ERROR: There appears to be an error when reading in the exponential function\n       for exponential %d.\n", expNumber+1);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  } 

  // error checking
  if (F->expFunction != 'X' && F->expFunction != 'S' && F->expFunction != 'C')
  { // error
    printf("ERROR: The exponential function for exponential %d appears to be incorrect.\n", expNumber+1);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  }

  // read in real & imaginary part of beta
  initialize_rational_number(F->beta);
  rV = mpq_inp_str(F->beta->re, IN, base);
  if (rV == 0)
  { // error in reading the coefficient
    printf("ERROR: There appears to be an error when reading in the real part of the\n       exponential constant for exponential %d.\n", expNumber+1);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  } 
  rV = mpq_inp_str(F->beta->im, IN, base);
  if (rV == 0)
  { // error in reading the coefficient
    printf("ERROR: There appears to be an error when reading in the imaginary part of the\n       exponential constant exponential %d.\n", expNumber+1);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  } 

  // setup in canonical form
  mpq_canonicalize(F->beta->re);
  mpq_canonicalize(F->beta->im);

  // update isReal
  F->isReal = !mpq_cmp_ui(F->beta->im, 0, 1);

  return;
}

void polynomial_system_check(FILE *IN)
/***************************************************************\
* USAGE: check the polynomial system for decimal points, e, & E *
\***************************************************************/
{
  char c;

  while ((c = fgetc(IN)) != EOF)
  { // check c
    if (c == '.')
    {
      printf("ERROR: To prevent floating point coefficients, no decimal points are allowed in the polynomial system file.\n");
      // close file and exit
      fclose(IN);
      errExit(ERROR_INPUT_SYSTEM);
    }
    else if (c == 'e')
    {
      printf("ERROR: To prevent floating point coefficients, 'e' is allowed in the polynomial system file.\n");
      // close file and exit
      fclose(IN);
      errExit(ERROR_INPUT_SYSTEM);
    }    
    else if (c == 'E')
    {
      printf("ERROR: To prevent floating point coefficients, 'E' is allowed in the polynomial system file.\n");
      // close file and exit
      fclose(IN);
      errExit(ERROR_INPUT_SYSTEM);
    }    
  }
 
  // rewind the file
  rewind(IN);

  return;
}

void setup_polynomial_system(polynomial_system *F, char *fileName, int realityCheck, int arithmeticType)
/***************************************************************\
* USAGE: setup the polynomial system described in fileName      *
\***************************************************************/
{
  int i, hasRealCoeff = 1;
  FILE *IN = fopen(fileName, "r");

  // verify file exists
  if (IN == NULL)
  {
    printf("\nERROR: '%s' does not exist!\n", fileName);
    errExit(ERROR_FILE_NOT_EXIST); 
  }

  // print message about rational coefficients
  printf("Please note that all coefficients must be complex rational numbers.\n\n");

  // check for decimal points, e, and E in the polynomial system - rewinds file if okay, closes file if error
  polynomial_system_check(IN);

  // initialize polynomial system
  initialize_polynomial_system(F);

  // read in the number of variables & number of polynomials
  fscanf(IN, "%d%d\n", &F->numVariables, &F->numPolynomials);

  // want 0 < variables && 0 < polynomials
  if (F->numVariables <= 0)
  { // error
    printf("ERROR: The number of variables (%d) must be positive.\n", F->numVariables);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  }
  else if (F->numPolynomials <= 0)
  { // error
    printf("ERROR: The number of polynomials (%d) must be positive.\n", F->numPolynomials);
    // close file and exit
    fclose(IN);
    errExit(ERROR_INPUT_SYSTEM);
  }

  // determine if using exponentials - if so, verify arithmetic type and allocate memory
  if (F->numVariables <= F->numPolynomials)
  { // polynomial only system
    printf("alphaCertified is using the polynomial certification algorithms.\n\n");
    F->numExponentials = 0;
    F->exponentials = NULL;
  }
  else // numVariables > numPolynomials
  { // polynomial exponential system
    if (arithmeticType)
    { // floating point - allocate memory
      printf("alphaCertified is using the polynomial-exponential certification algorithms.\n\n");
      F->numExponentials = F->numVariables - F->numPolynomials;
      F->exponentials = (exponential *)errMalloc(F->numExponentials * sizeof(exponential));
    }
    else
    { // error
      printf("ERROR: Exponentials may only be used with floating point arithmetic.\n\n");
      // close file and exit
      fclose(IN);
      errExit(ERROR_INPUT_SYSTEM);
    }
  }

  // allocate memory in F for polynomials
  F->polynomials = (polynomial *)errMalloc(F->numPolynomials * sizeof(polynomial));

  // read in the polynomials, compute maximum degree, hasRealCoeff, and norm_sqr
  mpq_set_ui(F->norm_sqr, 0, 1);
  for (i = 0; i < F->numPolynomials; i++)
  {
    setup_polynomial(&F->polynomials[i], F->numVariables, IN, i);
    if (F->polynomials[i].degree > F->maximumDegree)
      F->maximumDegree = F->polynomials[i].degree;
    hasRealCoeff = hasRealCoeff && F->polynomials[i].isReal;
    mpq_add(F->norm_sqr, F->norm_sqr, F->polynomials[i].norm_sqr);
  }

  // read in the exponentials and update hasRealCoeff
  for (i = 0; i < F->numExponentials; i++)
  {
    setup_exponential(&F->exponentials[i], F->numPolynomials, IN, i, F->numPolynomials + i);
    hasRealCoeff = hasRealCoeff && F->exponentials[i].isReal;    
  }  

  // setup F->isReal
  F->isReal = 0; // initialize to 0
  if (F->numPolynomials != F->numVariables)
  { // not square system - reality only depends upon coefficients
    F->isReal = hasRealCoeff;
  }
  else if (realityCheck < 0)
  { // square system that user says is real
    F->isReal = 1;
  }
  else // realityCheck >= 0
  { // square system - perform checks based on the settings
    if (hasRealCoeff)
    { // system has real coefficients
      F->isReal = 1;
    }
    else if (realityCheck > 0)
    { // test for conjugate map
      if (square_conj_map_test(F))
      { // system is invariant under conjugation
        F->isReal = 1;
      }
    }
  }

  // close IN
  fclose(IN);

  return;
}

void initialize_polynomial_system(polynomial_system *F)
/***************************************************************\
* USAGE: initialize the memory in F                             *
\***************************************************************/
{
  F->numVariables = F->numPolynomials = F->maximumDegree = F->isReal = F->numExponentials = 0;
  mpq_init(F->norm_sqr);
  F->polynomials = NULL;
  F->exponentials = NULL;

  return;
}

void clear_polynomial(polynomial *F)
/***************************************************************\
* USAGE: release the memory in F                                *
\***************************************************************/
{
  int i;

  mpq_clear(F->norm_sqr);
  for (i = 0; i < F->numTerms; i++)
  {
    clear_rational_number(F->coeff[i]);
    free(F->exponents[i]);
  }
  free(F->coeff);
  free(F->exponents);
  F->coeff = NULL;
  F->exponents = NULL;
  F->numVariables = F->numTerms = F->degree = F->isReal = 0;

  return;
}

void clear_exponential(exponential *F)
/***************************************************************\
* USAGE: release the memory in F                                *
\***************************************************************/
{
  clear_rational_number(F->beta);
  F->expFunction = '\n';
  F->xIndex = F->yIndex = F->isHyperbolic = F->isReal = 0;

  return;
}

void clear_polynomial_system(polynomial_system *F)
/***************************************************************\
* USAGE: release the memory in F                                *
\***************************************************************/
{
  int i;

  for (i = 0; i < F->numPolynomials; i++)
    clear_polynomial(&F->polynomials[i]);
  free(F->polynomials);
  F->polynomials = NULL;
  mpq_clear(F->norm_sqr);

  for (i = 0; i < F->numExponentials; i++)
    clear_exponential(&F->exponentials[i]);
  free(F->exponentials);
  F->exponentials = NULL;

  F->numVariables = F->numPolynomials = F->maximumDegree = F->isReal = F->numExponentials = 0;

  return;
}

void print_number(FILE *OUT, int digits, complex_number z)
/***************************************************************\
* USAGE: prints z to OUT                                        *
\***************************************************************/
{
  int base = 10;

  // verify that z is a valid number
  if (mpfr_number_p(z->re) && mpfr_number_p(z->im))
  { // print to OUT
    mpf_out_str(OUT, base, digits, z->re);
    if (mpfr_sgn(z->im) >= 0)
      fprintf(OUT, "+");
    mpf_out_str(OUT, base, digits, z->im);
    fprintf(OUT, "*I");
  }
  else
    fprintf(OUT, "NaN+NaN*I");

  return;
}

void print_rational_number(FILE *OUT, rational_complex_number z)
/***************************************************************\
* USAGE: prints z to OUT                                        *
\***************************************************************/
{
  int base = 10;

  mpq_out_str(OUT, base, z->re);
  if (mpq_sgn(z->im) >= 0)
    fprintf(OUT, "+");
  mpq_out_str(OUT, base, z->im);
  fprintf(OUT, "*I");

  return;
}

void print_vector_coordinate(FILE *OUT, int digits, complex_vector z)
/***************************************************************\
* USAGE: prints z                                               *
\***************************************************************/
{
  int i, size = z->size, base = 10;

  if (size > 0)
  { // print each coordinate
    for (i = 0; i < size; i++)
    {
      mpf_out_str(OUT, base, digits, z->coord[i]->re);
      fprintf(OUT, " ");
      mpf_out_str(OUT, base, digits, z->coord[i]->im);
      fprintf(OUT, "\n");
    }
  } 

  return;
}

void print_vector(FILE *OUT, int digits, complex_vector z)
/***************************************************************\
* USAGE: prints z                                               *
\***************************************************************/
{
  int i, size = z->size;

  if (size <= 0)
  { // nothing to print
    fprintf(OUT, "[];\n");
  }
  else
  { // size > 0
    for (i = 0; i < size; i++)
    {
      if (i == 0)
        fprintf(OUT, "[");

      fprintf(OUT, "[");
      print_number(OUT, digits, z->coord[i]);
      fprintf(OUT, "]");

      if (i == size - 1)
        fprintf(OUT, "];\n\n");
      else
        fprintf(OUT, ";\n");
    }
  } 

  return;
}

void print_rational_vector_coordinate(FILE *OUT, rational_complex_vector z)
/***************************************************************\
* USAGE: prints z                                               *
\***************************************************************/
{
  int i, size = z->size, base = 10;

  if (size > 0)
  { // print each coordinate
    for (i = 0; i < size; i++)
    {
      mpq_out_str(OUT, base, z->coord[i]->re);
      fprintf(OUT, " ");
      mpq_out_str(OUT, base, z->coord[i]->im);
      fprintf(OUT, "\n");
    }
  } 

  return;
}

void print_rational_vector(FILE *OUT, rational_complex_vector z)
/***************************************************************\
* USAGE: prints z                                               *
\***************************************************************/
{
  int i, size = z->size;

  if (size <= 0)
  { // nothing to print
    fprintf(OUT, "[];\n");
  }
  else
  { // size > 0
    for (i = 0; i < size; i++)
    {
      if (i == 0)
        fprintf(OUT, "[");

      fprintf(OUT, "[");
      print_rational_number(OUT, z->coord[i]);
      fprintf(OUT, "]");

      if (i == size - 1)
        fprintf(OUT, "];\n\n");
      else
        fprintf(OUT, ";\n");
    }
  } 

  return;
}

void print_matrix(FILE *OUT, int digits, complex_matrix z)
/***************************************************************\
* USAGE: prints z                                               *
\***************************************************************/
{
  int i, j, rows = z->rows, cols = z->cols;

  if (rows <= 0 || cols <= 0)
  { // nothing to print
    fprintf(OUT, "[];\n");
  }
  else
  { // rows > 0 && cols > 0
    for (i = 0; i < rows; i++)
    {
      if (i == 0)
        fprintf(OUT, "[");

      fprintf(OUT, "[");
      for (j = 0; j < cols; j++)
      {
        print_number(OUT, digits, z->entry[i][j]);
        if (j + 1 < cols)
          fprintf(OUT, ", ");
        else
          fprintf(OUT, "]");
      }

      if (i == rows - 1)
        fprintf(OUT, "];\n\n");
      else
        fprintf(OUT, ";\n");
    }
  } 

  return;
}

void print_rational_matrix(FILE *OUT, rational_complex_matrix z)
/***************************************************************\
* USAGE: prints z                                               *
\***************************************************************/
{
  int i, j, rows = z->rows, cols = z->cols;

  if (rows <= 0 || cols <= 0)
  { // nothing to print
    fprintf(OUT, "[];\n");
  }
  else
  { // rows > 0 && cols > 0
    for (i = 0; i < rows; i++)
    {
      if (i == 0)
        fprintf(OUT, "[");

      fprintf(OUT, "[");
      for (j = 0; j < cols; j++)
      {
        print_rational_number(OUT, z->entry[i][j]);
        if (j + 1 < cols)
          fprintf(OUT, ", ");
        else
          fprintf(OUT, "]");
      }

      if (i == rows - 1)
        fprintf(OUT, "];\n\n");
      else
        fprintf(OUT, ";\n");
    }
  } 

  return;
}

void determine_pivot_tolerances(mpf_t pivot_tol, mpf_t pivot_drop_tol, int prec)
/***************************************************************\
* USAGE: compute tolerances based on prec                       *
\***************************************************************/
{
  int num_digits = (int) floor(prec * log10(2.0) - 2.5);

  // setup pivot_tol
  mpf_set_prec(pivot_tol, prec);
  mpf_set_ui(pivot_tol, 10);
  mpf_pow_ui(pivot_tol, pivot_tol, num_digits);
  mpf_ui_div(pivot_tol, 1, pivot_tol);
  
  // setup pivot_drop_tol
  mpf_set_prec(pivot_drop_tol, prec);
  mpf_set_ui(pivot_drop_tol, 10);
  mpf_pow_ui(pivot_drop_tol, pivot_drop_tol, num_digits - 3);
  
  return;
}

void load_floating_points(int *numPoints, complex_vector **points, int numVars, char *PtsFile)
/***************************************************************\
* USAGE: load points from PtsFile                               *
\***************************************************************/
{
  int i, j, rV, base = 10;
  FILE *IN = fopen(PtsFile, "r");

  // error checking - file must exist
  if (IN == NULL)
  {
    printf("\nERROR: '%s' does not exist!\n", PtsFile);
    errExit(ERROR_FILE_NOT_EXIST); 
  }

  // read in the number of points
  rV = fscanf(IN, "%d", numPoints);

  // error checking
  if (rV != 1)
  { 
    printf("\nERROR: Unable to read the number of points stored in '%s'.\n", PtsFile);
    errExit(ERROR_FILE_NOT_EXIST); 
  }
  else if (*numPoints <= 0)
  {
    printf("ERROR: The number of points in '%s' must be positive!\n", PtsFile);
    errExit(ERROR_CONFIGURATION);
  }

  // allocate points
  *points = (complex_vector *)errMalloc((*numPoints) * sizeof(complex_vector));
  
  // read in the points
  for (i = 0; i < *numPoints; i++)
  { // setup points[i]
    initialize_vector((*points)[i], numVars);
  
    for (j = 0; j < numVars; j++)
    { // setup real part
      rV = mpf_inp_str((*points)[i]->coord[j]->re, IN, base);

      // error checking
      if (rV == 0)
      { 
        printf("\nERROR: Unable to read in %d floating point vectors from '%s'.\n", *numPoints, PtsFile);
        errExit(ERROR_FILE_NOT_EXIST); 
      }

      // setup imag part
      rV = mpf_inp_str((*points)[i]->coord[j]->im, IN, base);

      // error checking
      if (rV == 0)
      { 
        printf("\nERROR: Unable to read in %d floating point vectors from '%s'.\n", *numPoints, PtsFile);
        errExit(ERROR_FILE_NOT_EXIST); 
      }
    }  
  }

  // close file
  fclose(IN);

  return;
}

void load_rational_points(int *numPoints, rational_complex_vector **points, int numVars, char *PtsFile)
/***************************************************************\
* USAGE: load rational points from PtsFile                      *
\***************************************************************/
{
  int i, j, rV, base = 10;
  FILE *IN = fopen(PtsFile, "r");

  // error checking - file must exist
  if (IN == NULL)
  {
    printf("\nERROR: '%s' does not exist!\n", PtsFile);
    errExit(ERROR_FILE_NOT_EXIST); 
  }

  // read in the number of points
  rV = fscanf(IN, "%d", numPoints);

  // error checking
  if (rV != 1)
  { 
    printf("\nERROR: Unable to read the number of points stored in '%s'.\n", PtsFile);
    errExit(ERROR_FILE_NOT_EXIST); 
  }
  else if (*numPoints <= 0)
  {
    printf("ERROR: The number of points in '%s' must be positive!\n", PtsFile);
    errExit(ERROR_CONFIGURATION);
  }

  // allocate points
  *points = (rational_complex_vector *)errMalloc((*numPoints) * sizeof(rational_complex_vector));
  
  // read in the points
  for (i = 0; i < *numPoints; i++)
  { // setup points[i]
    initialize_rational_vector((*points)[i], numVars);
  
    for (j = 0; j < numVars; j++)
    { // setup real part
      rV = mpq_inp_str((*points)[i]->coord[j]->re, IN, base);

      // error checking
      if (rV == 0)
      { 
        printf("\nERROR: Unable to read in %d rational vectors from '%s'.\n", *numPoints, PtsFile);
        errExit(ERROR_FILE_NOT_EXIST); 
      }

      // setup imag part
      rV = mpq_inp_str((*points)[i]->coord[j]->im, IN, base);

      // error checking
      if (rV == 0)
      { 
        printf("\nERROR: Unable to read in %d rational vectors from '%s'.\n", *numPoints, PtsFile);
        errExit(ERROR_FILE_NOT_EXIST); 
      }

      // canonicalize rational numbers
      mpq_canonicalize((*points)[i]->coord[j]->re);
      mpq_canonicalize((*points)[i]->coord[j]->im);
    }  
  }

  // close file
  fclose(IN);

  return;
}

void initialize_point_struct(point_struct *Pt, int numVariables)
/***************************************************************\
* USAGE: initialize the memory                                  *
\***************************************************************/
{ // initialize the memory
  initialize_vector(Pt->origX, numVariables);
  initialize_vector(Pt->x, numVariables);
  initialize_vector(Pt->Nx, numVariables);
  mpf_init(Pt->norm_x);
  initialize_number(Pt->origAlpha);
  initialize_number(Pt->origBeta);
  initialize_number(Pt->origGamma);
  initialize_number(Pt->alpha);
  initialize_number(Pt->beta);
  initialize_number(Pt->gamma);
  Pt->isApproxSoln = Pt->isActive = Pt->isReal = 0;

  return;
}

void initialize_rational_point_struct(rational_point_struct *Pt, int numVariables)
/***************************************************************\
* USAGE: initialize the memory                                  *
\***************************************************************/
{ // initialize the memory
  initialize_rational_vector(Pt->origX, numVariables);
  initialize_rational_vector(Pt->x, numVariables);
  initialize_rational_vector(Pt->Nx, numVariables);
  mpq_init(Pt->norm_sqr_x);
  initialize_rational_number(Pt->origAlpha_sqr);
  initialize_rational_number(Pt->origBeta_sqr);
  initialize_rational_number(Pt->origGamma_sqr);
  initialize_rational_number(Pt->alpha_sqr);
  initialize_rational_number(Pt->beta_sqr);
  initialize_rational_number(Pt->gamma_sqr);
  Pt->isApproxSoln = Pt->isActive = Pt->isReal = 0;

  return;
}

void clear_point_struct(point_struct *Pt)
/***************************************************************\
* USAGE: clear the memory                                       *
\***************************************************************/
{ // clear the memory
  clear_vector(Pt->origX);
  clear_vector(Pt->x);
  clear_vector(Pt->Nx);
  mpf_clear(Pt->norm_x);
  clear_number(Pt->origAlpha);
  clear_number(Pt->origBeta);
  clear_number(Pt->origGamma);
  clear_number(Pt->alpha);
  clear_number(Pt->beta);
  clear_number(Pt->gamma);
  Pt->isApproxSoln = Pt->isActive = Pt->isReal = 0;

  return;
}

void clear_rational_point_struct(rational_point_struct *Pt)
/***************************************************************\
* USAGE: clear the memory                                       *
\***************************************************************/
{ // initialize the memory
  clear_rational_vector(Pt->origX);
  clear_rational_vector(Pt->x);
  clear_rational_vector(Pt->Nx);
  mpq_clear(Pt->norm_sqr_x);
  clear_rational_number(Pt->origAlpha_sqr);
  clear_rational_number(Pt->origBeta_sqr);
  clear_rational_number(Pt->origGamma_sqr);
  clear_rational_number(Pt->alpha_sqr);
  clear_rational_number(Pt->beta_sqr);
  clear_rational_number(Pt->gamma_sqr);
  Pt->isApproxSoln = Pt->isActive = Pt->isReal = 0;

  return;
}

void create_random_number_str(char **str)
/***************************************************************\
* USAGE: creates a random number in [-1,1]                      *
\***************************************************************/
{
  int i, counter = 0, numDigits = RATIONALDIGITLENGTH;

  // allocate str to size
  *str = (char *)errRealloc(*str, (2*numDigits + 4) * sizeof(char));

  // random sign
  if (rand() % 2)
  {
    (*str)[counter] = '-';
    counter++;
  }

  for (i = 0; i < numDigits; i++)
  {
    (*str)[counter] = 48 + (rand() % 10); // ASCII for the digits
    counter++;
  }
  (*str)[counter] = '/';
  counter++;
  (*str)[counter] = '1';
  counter++;
  for (i = 0; i < numDigits; i++)
  {
    (*str)[counter] = 48 + (rand() % 10); // ASCII for the digits
    counter++;
  }
  (*str)[counter] = '\0';

  return;
}

void generate_random_rational(rational_complex_number x)
/***************************************************************\
* USAGE: generate a random rational number of unit modulus      *
\***************************************************************/
{
  int base = 10;
  char *str = NULL;
  mpq_t t, t_sqr;

  mpq_init(t);
  mpq_init(t_sqr);  

  // generate a random number
  create_random_number_str(&str);
  mpq_set_str(t, str, base);
  mpq_canonicalize(t);

  // compute t_sqr = t^2
  mpq_mul(t_sqr, t, t);

  // compute denominator
  mpq_set_ui(x->im, 1, 1);
  mpq_add(x->im, x->im, t_sqr); // 1 + t^2

  // compute numerator
  mpq_set_ui(x->re, 1, 1);
  mpq_sub(x->re, x->re, t_sqr); // 1 - t^2

  // divide to compute real part
  mpq_div(x->re, x->re, x->im);

  // compute imaginary part
  mpq_set_ui(x->im, 1, 1);
  mpq_add(x->im, x->re, x->im); // 1 + x->re
  mpq_mul(x->im, x->im, t);     // t*(1+x->re)

  // clear memory
  mpq_clear(t);
  mpq_clear(t_sqr);
  free(str);

  return;
}

void generate_random_real_rational_vector(rational_complex_vector x, int size)
/***************************************************************\
* USAGE: generate a random real rational vector of unit length  *
\***************************************************************/
{ 
  // set x to the correct size & set to zero
  change_size_rational_vector(x, size);
  set_zero_rational_vector(x);

  if (size == 1)
  { // setup x to either 1 or -1
    mpq_set_si(x->coord[0]->re, rand() % 2 ? 1 : -1, 1);
  }
  else if (size > 1)
  { // setup x to real unit vector
    int i, base = 10;
    char *str = NULL;
    mpq_t sum_sqr, *t = (mpq_t *)errMalloc((size - 1) * sizeof(mpq_t)), *t_sqr = (mpq_t *)errMalloc((size - 1) * sizeof(mpq_t));

    // initialize and set to random numbers
    mpq_init(sum_sqr);
    mpq_set_ui(sum_sqr, 0, 1);
    for (i = 0; i < size - 1; i++)
    {
      mpq_init(t[i]);
      mpq_init(t_sqr[i]);
    
      create_random_number_str(&str);
      mpq_set_str(t[i], str, base);
      mpq_canonicalize(t[i]);

      mpq_mul(t_sqr[i], t[i], t[i]);
      mpq_add(sum_sqr, sum_sqr, t_sqr[i]);
    }

    // setup the first entry
    mpq_set_ui(x->coord[0]->re, 1, 1);
    mpq_add(x->coord[0]->im, x->coord[0]->re, sum_sqr); // 1 + sum(t^2)
    mpq_sub(sum_sqr, x->coord[0]->re, sum_sqr); // 1 - sum(t^2)
    mpq_div(x->coord[0]->re, sum_sqr, x->coord[0]->im); // (1 - sum(t^2)) / (1 + sum(t^2))
    mpq_set_ui(x->coord[0]->im, 0, 1);

    // compute x[0] + 1
    mpq_set_ui(sum_sqr, 1, 1);
    mpq_add(sum_sqr, sum_sqr, x->coord[0]->re);

    // setup other entries: t[i] * (x[0] + 1)
    for (i = 1; i < size; i++)
      mpq_mul(x->coord[i]->re, t[i-1], sum_sqr);

    // clear memory
    mpq_clear(sum_sqr);
    for (i = 0; i < size - 1; i++)
    {
      mpq_clear(t[i]);
      mpq_clear(t_sqr[i]);
    }
    free(t);
    free(t_sqr);
    free(str);
  }

  return;
}

void generate_random_rational_vector(rational_complex_vector x, int size)
/***************************************************************\
* USAGE: generate a random rational vector of unit length       *
\***************************************************************/
{
  int i;
  rational_complex_vector xR;
  initialize_rational_vector(xR, 2 * size);

  // compute a real unit vector of twice the size
  generate_random_real_rational_vector(xR, 2 * size);

  // setup x from xR
  change_size_rational_vector(x, size);
  for (i = 0; i < size; i++)
  {
    mpq_set(x->coord[i]->re, xR->coord[2*i]->re);
    mpq_set(x->coord[i]->im, xR->coord[2*i + 1]->re);
  }

  // clear
  clear_rational_vector(xR);

  return;
}

void random_rational_long_matrix(rational_complex_matrix A, int rows, int cols)
/***************************************************************\
* USAGE: create a random long rational matrix                   *
\***************************************************************/
{
  int i, j, k;
  rational_complex_number gamma;
  rational_complex_vector x, u;

  // verify that rows >= cols
  if (rows < cols)
  {
    printf("ERROR: The number of rows can not be less than the number of columns!!\n");
    errExit(ERROR_CONFIGURATION);
  }

  initialize_rational_number(gamma);
  initialize_rational_vector(x, rows);
  initialize_rational_vector(u, rows);

  // find a random gamma
  generate_random_rational(gamma);

  // set A to a square matrix
  change_size_rational_matrix(A, rows, rows);

  // initialize A to 'gamma * Id' matrix for random gamma of unit modulus
  for (i = 0; i < rows; i++)
    for (j = 0; j < rows; j++)
      if (i == j)
      {
        set_rational_number(A->entry[i][j], gamma);
      }
      else
      {
        set_zero_rational_number(A->entry[i][j]);
      }

  // generate random vectors of unit length and apply Householder transformation to A
  for (j = 0; j < cols; j++)
  { // generate random vector of correct size
    generate_random_rational_vector(u, rows - j);

    // compute x = A*u taking into account the size of u    
    for (i = 0; i < rows; i++)
    {
      set_zero_rational_number(x->coord[i]);
      for (k = j; k < rows; k++)
      {
        sum_multiply_rational(x->coord[i], A->entry[i][k], u->coord[k-j]);
      }
    }

    // update A := A - 2*(A*u)*u^* taking into account the size of u
    for (i = 0; i < rows; i++)
      for (k = j; k < rows; k++)
      {
        householder_rational_multiplication(A->entry[i][k], A->entry[i][k], x->coord[i], u->coord[k-j]);
      }
  }

  // set A to the correct size
  A->cols = cols;

  clear_rational_number(gamma);
  clear_rational_vector(x);
  clear_rational_vector(u);

  return;
}

void random_rational_matrix(rational_complex_matrix A, int rows, int cols)
/***************************************************************\
* USAGE: create a random rational matrix                        *
\***************************************************************/
{
  if (rows >= cols)
  { // create a long matrix
    random_rational_long_matrix(A, rows, cols);
  }
  else
  { // create a long matrix and then transpose
    rational_complex_matrix B;
    initialize_rational_matrix(B, cols, rows);
    random_rational_long_matrix(B, cols, rows);
    transpose_rational_matrix(A, B);
    clear_rational_matrix(B);
  }

  return;
}

void random_real_rational_long_matrix(rational_complex_matrix A, int rows, int cols)
/***************************************************************\
* USAGE: create a real random long rational matrix              *
\***************************************************************/
{
  int i, j, k;
  rational_complex_number gamma;
  rational_complex_vector x, u;

  // verify that rows >= cols
  if (rows < cols)
  {
    printf("ERROR: The number of rows can not be less than the number of columns!!\n");
    errExit(ERROR_CONFIGURATION);
  }

  initialize_rational_number(gamma);
  initialize_rational_vector(x, rows);
  initialize_rational_vector(u, rows);

  // find a real random gamma
  generate_random_rational(gamma);
  mpq_set_ui(gamma->im, 0, 1);

  // set A to a square matrix
  change_size_rational_matrix(A, rows, rows);

  // initialize A to 'Id' matrix
  for (i = 0; i < rows; i++)
    for (j = 0; j < rows; j++)
      if (i == j)
      {
        set_rational_number(A->entry[i][j], gamma);
      }
      else
      {
        set_zero_rational_number(A->entry[i][j]);
      }

  // generate random vectors of unit length and apply Householder transformation to A
  for (j = 0; j < cols; j++)
  { // generate random vector of correct size
    generate_random_real_rational_vector(u, rows - j);

    // compute x = A*u taking into account the size of u    
    for (i = 0; i < rows; i++)
    {
      set_zero_rational_number(x->coord[i]);
      for (k = j; k < rows; k++)
      {
        sum_multiply_rational(x->coord[i], A->entry[i][k], u->coord[k-j]);
      }
    }

    // update A := A - 2*(A*u)*u^* taking into account the size of u
    for (i = 0; i < rows; i++)
      for (k = j; k < rows; k++)
      {
        householder_rational_multiplication(A->entry[i][k], A->entry[i][k], x->coord[i], u->coord[k-j]);
      }
  }

  // set A to the correct size
  A->cols = cols;

  clear_rational_number(gamma);
  clear_rational_vector(x);
  clear_rational_vector(u);

  return;
}

void random_real_rational_matrix(rational_complex_matrix A, int rows, int cols)
/***************************************************************\
* USAGE: create a random real rational matrix                   *
\***************************************************************/
{
  if (rows >= cols)
  { // create a long matrix
    random_real_rational_long_matrix(A, rows, cols);
  }
  else
  { // create a long matrix and then transpose
    rational_complex_matrix B;
    initialize_rational_matrix(B, cols, rows);
    random_real_rational_long_matrix(B, cols, rows);
    transpose_rational_matrix(A, B);
    clear_rational_matrix(B);
  }

  return;
}

void transpose_rational_matrix(rational_complex_matrix A, rational_complex_matrix B)  
/***************************************************************\
* USAGE: A = B^* (conjugate transpose)                          *
\***************************************************************/
{
  int i, j, rows = B->rows, cols = B->cols;

  if (A != B)
  { // setup A
    change_size_rational_matrix(A, cols, rows);

    for (i = 0; i < cols; i++)
      for (j = 0; j < rows; j++)
      {
        conjugate_rational(A->entry[i][j], B->entry[j][i]);
      }
  }
  else // A == B
  { // need to use a temporary matrix
    rational_complex_matrix tempMat;

    // copy B to tempMat
    initialize_rational_matrix(tempMat, rows, cols);
    copy_rational_matrix(tempMat, B);

    // setup A
    change_size_rational_matrix(A, cols, rows);
 
    for (i = 0; i < cols; i++)
      for (j = 0; j < rows; j++)
      {
        conjugate_rational(A->entry[i][j], tempMat->entry[j][i]);
      }

    // clear tempMat
    clear_rational_matrix(tempMat);
  }
 
  return;
}

void print_welcome_message(FILE *OUT)
/***************************************************************\
* USAGE: print welcome message                                  *
\***************************************************************/
{
  int i, max, count[3] = {0, 0, 0};

  // count the number of characeters in each message
  count[0] = snprintf(NULL, 0, "alphaCertified v%s (%s)\n", VERSION_STRING, DATE_STRING);
  count[1] = snprintf(NULL, 0, "Jonathan D. Hauenstein and Frank Sottile\n");
  count[2] = snprintf(NULL, 0, "GMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());

  // find the maximum
  max = count[0];
  if (max < count[1])
    max = count[1];
  if (max < count[2])
    max = count[2];

  fprintf(OUT, "\n");

  // pad and then print each line
  count[0] = (max - count[0]) / 2 + 2;
  for (i = 0; i < count[0]; i++)
    fprintf(OUT, " ");
  fprintf(OUT, "alphaCertified v%s (%s)\n", VERSION_STRING, DATE_STRING);
  count[1] = (max - count[1]) / 2 + 2;
  for (i = 0; i < count[1]; i++)
    fprintf(OUT, " ");
  fprintf(OUT,   "Jonathan D. Hauenstein and Frank Sottile\n");
  count[2] = (max - count[2]) / 2 + 2;
  for (i = 0; i < count[2]; i++)
    fprintf(OUT, " ");
  fprintf(OUT,   "GMP v%d.%d.%d & MPFR v%s\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());

  return;
}

int square_conj_map_test(polynomial_system *F)
/***************************************************************\
* USAGE: determine if F is invariant under conjugation          *
\***************************************************************/
{
  int i, j, found = 0, isReal = 1, numVars = F->numVariables, base = 10;
  char *str = NULL;
  mpq_t scale;
  rational_complex_vector x, func1, func2;

  if (F->numVariables != F->numPolynomials)
  {
    printf("\nERROR: The polynomial must be square!\n");
    errExit(ERROR_INVALID_SIZE);
  }

  mpq_init(scale);
  initialize_rational_vector(x, numVars);
  initialize_rational_vector(func1, numVars);
  initialize_rational_vector(func2, numVars);

  // generate a random rational real vector x
  generate_random_real_rational_vector(x, numVars);

  // generate a random rational number in [-1,1]
  create_random_number_str(&str);
  mpq_set_str(scale, str, base);
  mpq_canonicalize(scale);

  // scale x
  for (i = 0; i < numVars; i++)
   multiply_rational_number(x->coord[i], x->coord[i], scale);

  // evaluate F at x
  eval_polynomial_system_only_rational(func1, F, x);

  // conjugate func1
  conjugate_rational_vector(func2, func1);
 
  // determine if, as sets, func1 == func2
  isReal = 1;
  for (i = 0; i < numVars && isReal; i++)
  { // determine if func1[i] is in func2
    found = 0;
    for (j = 0; j < numVars && !found; j++)
      if (mpq_equal(func1->coord[i]->re, func2->coord[j]->re) && mpq_equal(func1->coord[i]->im, func2->coord[j]->im))
        found = 1;

    if (!found)
      isReal = 0;
  }

  free(str);
  str = NULL;
  mpq_clear(scale);
  clear_rational_vector(x);
  clear_rational_vector(func1);
  clear_rational_vector(func2);

  return isReal;
}

int square_newton_real_map_test(polynomial_system *F)
/***************************************************************\
* USAGE: determine if N_F is a real map where F is square       *
\***************************************************************/
{
  int i, retVal, isReal = 1, numVars = F->numVariables, base = 10, *rowswaps = NULL;
  char *str = NULL;
  mpq_t scale;
  rational_complex_number r;
  rational_complex_vector x, y;
  rational_complex_matrix LU;

  if (F->numVariables != F->numPolynomials)
  {
    printf("\nERROR: The polynomial must be square!\n");
    errExit(ERROR_INVALID_SIZE);
  }

  mpq_init(scale);
  initialize_rational_number(r);
  initialize_rational_vector(x, numVars);
  initialize_rational_vector(y, numVars);
  initialize_rational_matrix(LU, numVars, numVars);

  // generate a random rational real vector x
  generate_random_real_rational_vector(x, numVars);

  // generate a random rational number in [-1,1]
  create_random_number_str(&str);
  mpq_set_str(scale, str, base);
  mpq_canonicalize(scale);

  // scale x
  for (i = 0; i < numVars; i++)
   multiply_rational_number(x->coord[i], x->coord[i], scale);
 
  // compute a newton iteration for x
  retVal = newton_iteration_rational(y, r, LU, &rowswaps, F, x);

  if (retVal)
  { // Jacobian is rank deficient for random point
    printf("\nERROR: The Jacobian is rank deficient at a random point!\n");
    errExit(ERROR_INPUT_SYSTEM);
  }

  // determine the imaginary parts are zero
  mpq_set_ui(scale, 0, 1);
  for (i = 0; i < numVars && isReal; i++)
  {
    if (!mpq_equal(y->coord[i]->im, scale))
      isReal = 0;
  }

  // clear memory
  free(str);
  str = NULL;
  free(rowswaps);
  rowswaps = NULL;
  mpq_clear(scale);
  clear_rational_number(r);
  clear_rational_vector(x);
  clear_rational_vector(y);
  clear_rational_matrix(LU);

  return isReal;
}

int print_coeff(FILE *OUT, rational_complex_number z, int somethingPrinted)
/***************************************************************\
* USAGE: prints z to OUT in user-friendly way                 *
\***************************************************************/
{
  int rV = -2; // 0 if z is zero (nothing printed), 1 if z is one (nothing printed), -1 if z is -1 (nothing printed), -2 otherwise (something printed)
  int base = 10;

  if (mpq_cmp_ui(z->re, 0, 1) == 0)
  { // real part is zero
    if (mpq_cmp_ui(z->im, 0, 1) == 0)
    { // imag part is zero
      rV = 0;
    }
    else
    { // imag part is nonzero
      if (somethingPrinted && mpq_sgn(z->im) >= 0)
        fprintf(OUT, "+");

      if (mpq_cmp_ui(z->im, 1, 1) == 0)
        fprintf(OUT, "I");
      else if (mpq_cmp_si(z->im, -1, 1) == 0)
        fprintf(OUT, "-I");
      else
      {
        mpq_out_str(OUT, base, z->im);
        fprintf(OUT, "*I");
      }
    }
  }
  else
  { // real part is nonzero
    if (mpq_cmp_ui(z->im, 0, 1) == 0)
    { // imag part is zero
      if (mpq_cmp_ui(z->re, 1, 1) == 0)
      { // value is 1
        rV = 1;
      }
      else if (mpq_cmp_si(z->re, -1, 1) == 0)
      { // value is -1
        rV = -1;
      }
      else
      {
        if (somethingPrinted && mpq_sgn(z->re) >= 0)
          fprintf(OUT, "+");
        mpq_out_str(OUT, base, z->re);
      }
    }
    else
    { // imag part is nonzero
      if (somethingPrinted)
        fprintf(OUT, "+");
      fprintf(OUT, "(");
      mpq_out_str(OUT, base, z->re);
      if (mpq_sgn(z->im) >= 0)
        fprintf(OUT, "+");
      mpq_out_str(OUT, base, z->im);
      fprintf(OUT, "*I)");
    }
  }

  return rV;
}

void print_monomial(FILE *OUT, int *exponents, int numVars, int numFuncs, int coeff_rV, int somethingPrinted)
/***************************************************************\
* USAGE: prints monomial to OUT in user-friendly way            *
\***************************************************************/
{
  int i, numNonZero = 0;

  // compute totalDegree of monomial
  for (i = 0; i < numVars; i++)
    if (exponents[i] != 0)
      numNonZero++;

  // determine if totalDegree is 0 or not
  if (numNonZero == 0)
  { // monomial is 1
    if (coeff_rV == 1)
    { // print 1
      if (somethingPrinted)
        fprintf(OUT, "+");
      fprintf(OUT, "1");
    }
    else if (coeff_rV == -1)
    { // print -1
      fprintf(OUT, "-1");
    }
  }
  else if (numNonZero == 1)
  { // print monomial
    if (coeff_rV == 1)
    { // coefficient is 1
      if (somethingPrinted)
        fprintf(OUT, "+");
    }
    else if (coeff_rV == -1)
    { // coefficient is -1
      fprintf(OUT, "-");
    }
    else
    { // coefficient is something else
      fprintf(OUT, "*");
    }
    
    for (i = 0; i < numVars; i++)
      if (exponents[i] == 1)
        fprintf(OUT, "%c%d", i < numFuncs ? 'x' : 'y', i < numFuncs ? i+1 : i+1-numFuncs);
      else if (exponents[i] > 1)
        fprintf(OUT, "%c%d^%d", i < numFuncs ? 'x' : 'y', i < numFuncs ? i+1 : i+1-numFuncs, exponents[i]);
  }
  else
  { // print monomial
    numNonZero = 1;
    if (coeff_rV == 1)
    { // coefficient is 1
      numNonZero = 0;
      if (somethingPrinted)
        fprintf(OUT, "+");
    }
    else if (coeff_rV == -1)
    { // coefficient is -1
      fprintf(OUT, "-");
      numNonZero = 0;
    }

    for (i = 0; i < numVars; i++)
      if (exponents[i] == 1)
      { // print variable
        if (numNonZero == 0)
          fprintf(OUT, "%c%d", i < numFuncs ? 'x' : 'y', i < numFuncs ? i+1 : i+1-numFuncs);
        else
          fprintf(OUT, "*%c%d", i < numFuncs ? 'x' : 'y', i < numFuncs ? i+1 : i+1-numFuncs);
        numNonZero++;
      }
      else if (exponents[i] > 1)
      { // print variable
        if (numNonZero == 0)
          fprintf(OUT, "%c%d", i < numFuncs ? 'x' : 'y', i < numFuncs ? i+1 : i+1-numFuncs);
        else
          fprintf(OUT, "*%c%d", i < numFuncs ? 'x' : 'y', i < numFuncs ? i+1 : i+1-numFuncs);
        numNonZero++;
      }
  }

  return;
}

void print_polynomial(FILE *OUT, polynomial *F, int numFuncs)
/***************************************************************\
* USAGE: user-friendly printing of polynomial F to a file       *
\***************************************************************/
{
  int i, rV, somethingPrinted = 0, numTerms = F->numTerms, numVars = F->numVariables;

  for (i = 0; i < numTerms; i++)
  { 
    // print coefficient
    rV = print_coeff(OUT, F->coeff[i], somethingPrinted);
    if (rV)
    { // print monomial
      print_monomial(OUT, F->exponents[i], numVars, numFuncs, rV, somethingPrinted);

      // something has been printed
      somethingPrinted = 1;
    }
  }

  // make sure something is printed
  if (!somethingPrinted)
    fprintf(OUT, "0");

  return;
}

void print_exponential(FILE *OUT, exponential *F, int numFuncs)
/***************************************************************\
* USAGE: user-friendly printing of exponential F to a file      *
\***************************************************************/
{
  int rV;

  // print y variable
  fprintf(OUT, "y%d-", F->yIndex+1-numFuncs);

  // print function
  if (F->expFunction == 'X')
    fprintf(OUT, "exp(");
  else
  { // print either sin or cos
    if (F->expFunction == 'S')
      fprintf(OUT, "sin");
    else if (F->expFunction == 'C')
      fprintf(OUT, "cos");

    // determine if using hyperbolic or circular
    if (F->isHyperbolic)
      fprintf(OUT, "h");

    fprintf(OUT, "(");
  }

  // print beta
  rV = print_coeff(OUT, F->beta, 0);

  if (rV == 0)
  { // beta is zero
    fprintf(OUT, "0)");
  }
  else
  { // print sign and then x variable
    if (rV == -1)
      fprintf(OUT, "-");
    else if (rV != 1)
      fprintf(OUT, "*");

    // print x variable
    fprintf(OUT, "x%d)", F->xIndex+1);
  }

  return;
}

void print_polynomial_system(FILE *OUT, polynomial_system *F)
/***************************************************************\
* USAGE: user-friendly printing of polynomial system F to a file*
\***************************************************************/
{
  int i, numVars = F->numVariables, numPolys = F->numPolynomials, numExps = F->numExponentials;
  int numFuncs = numPolys + numExps;

  // print variables
  fprintf(OUT, "\nvariable ");
  for (i = 0; i < numVars; i++)
  {
    fprintf(OUT, "%c%d", i < numPolys ? 'x' : 'y', i < numPolys ? i + 1 : i + 1 - numPolys);
    if (i+1 < numVars)
      fprintf(OUT, ",");
    else
      fprintf(OUT, ";\n");
  }

  // print functions
  fprintf(OUT, "function ");
  for (i = 0; i < numFuncs; i++)
  {
    fprintf(OUT, "F%d", i + 1);
    if (i + 1 < numVars)
      fprintf(OUT, ",");
    else
      fprintf(OUT, ";\n\n");
  }

  // print the polynomials
  for (i = 0; i < numPolys; i++)
  {
    fprintf(OUT, "F%d = ", i + 1);
    print_polynomial(OUT, &F->polynomials[i], F->numPolynomials);
    fprintf(OUT, ";\n");
  }

  // print the exponentials
  for (i = 0; i < numExps; i++)
  {
    fprintf(OUT, "F%d = ", numPolys + i + 1);
    print_exponential(OUT, &F->exponentials[i], F->numPolynomials);
    fprintf(OUT, ";\n");
  }

  fprintf(OUT, "\n");

  return;
}


