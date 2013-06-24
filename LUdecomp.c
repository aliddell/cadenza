/* 
   alphaCertified
   Jonathan Hauenstein & Frank Sottile
   May 7, 2010
   Copyright 2010

   LUdecomp.c: Compute an LU decomposition of a matrix using scaled partial pivoting
*/

#include "alphaCertified.h"

int LUdecomp(complex_matrix LU, int **rowswaps, complex_matrix A, mpf_t pivot_tol, mpf_t pivot_drop_tol, int LU_prec)
/***************************************************************\
* USAGE: compute an LU decomposition of A, use tolerances to    *
* decide if the LU decomp should be computed or return error    *
\***************************************************************/
{
  int i, j, k, l, pivot, n = A->rows, *rownum = NULL;
  mpf_t max, tempMPF, prevNorm, *scale = NULL;
  complex_number tempComp;

  // error checking - need a square matrix
  if (A->rows != A->cols)
  {
    printf("ERROR: In order to compute an LU decomposition, the matrix must be square!!\n");
    errExit(ERROR_INVALID_SIZE);
  }

  // set the default precision
  setPrec(LU_prec);

  // initialize memory
  mpf_init(max); mpf_init(tempMPF); mpf_init(prevNorm);
  initialize_number(tempComp);

  // setup rowswaps & scale
  rownum = *rowswaps = (int *)errRealloc(*rowswaps, n * sizeof(int));
  scale = (mpf_t *)errMalloc(n * sizeof(mpf_t));
  for (i = 0; i < n; i++)
  { // initialize rownum[i] & scale[i]
    rownum[i] = i;
    mpf_init(scale[i]);
    mpf_set_ui(scale[i], 0);
  }

  // setup LU
  setPrec_matrix(LU, LU_prec);
  change_size_matrix(LU, n, n);
  
  // find scale factors for each row and copy A to LU
  mpf_set_ui(prevNorm, 0);
  for (i = 0; i < n; i++)
  { 
    set_number(LU->entry[i][0], A->entry[i][0]);
    norm_sqr_number(max, LU->entry[i][0]);
    for (j = 1; j < n; j++)
    {
      set_number(LU->entry[i][j], A->entry[i][j]);
      norm_sqr_number(tempMPF, LU->entry[i][j]);
      if (mpf_cmp(max, tempMPF) < 0)
        mpf_set(max, tempMPF);
    }
    // compute sqrt of max
    mpf_sqrt(max, max);

    // compare with pivot tolerance
    if (mpf_cmp(max, pivot_tol) < 0)
    { // give up!
      mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
      clear_number(tempComp);
      for (i = 0; i < n; i++)
        mpf_clear(scale[i]);
      free(scale);
      scale = NULL;
      rownum = NULL;

      return ERROR_LU_DECOMP;
    }

    // reciprocate scale[i] - turn / into *
    mpf_ui_div(scale[i], 1, max); 
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find scaled maximum element in column k
    norm_number(max, LU->entry[rownum[k]][k]);
    mpf_mul(max, max, scale[rownum[k]]);
    for (i = k + 1; i < n; i++)
    {
      norm_number(tempMPF, LU->entry[rownum[i]][k]);
      mpf_mul(tempMPF, tempMPF, scale[rownum[i]]);
      if (mpf_cmp(max, tempMPF) < 0)
      {
        pivot = i;
        mpf_set(max, tempMPF);
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change
    mpf_mul(tempMPF, max, pivot_drop_tol);
    if (mpf_cmp(max, pivot_tol) < 0 || (k > 0 && mpf_cmp(prevNorm, tempMPF) > 0))
    { // clear
      mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
      clear_number(tempComp);
      for (i = 0; i < n; i++)
        mpf_clear(scale[i]);
      free(scale);
      scale = NULL;
      rownum = NULL;

      return ERROR_LU_DECOMP;
    }
    else
    { // update prevNorm
      mpf_set(prevNorm, max);
    }

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
    }

    // update elements in "bottom" block of LU
    pivot = rownum[k];
    recip_number(tempComp, LU->entry[pivot][k]);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];

      // compute muliplier
      multiply(LU->entry[l][k], LU->entry[l][k], tempComp);

      // update elements of row i
      for (j = k + 1; j < n; j++)
        subtract_multiply(LU->entry[l][j], LU->entry[l][k], LU->entry[pivot][j]);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change
  l = n - 1;
  pivot = rownum[l];
  norm_number(max, LU->entry[pivot][l]);
  mpf_mul(max, max, scale[pivot]);
  mpf_mul(tempMPF, max, pivot_drop_tol);
  if (mpf_cmp(max, pivot_tol) < 0 || mpf_cmp(prevNorm, tempMPF) > 0)
  { // clear
    mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
    clear_number(tempComp);
    for (i = 0; i < n; i++)
      mpf_clear(scale[i]);
    free(scale);
    scale = NULL;
    rownum = NULL;

    return ERROR_LU_DECOMP;
  }

  mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
  clear_number(tempComp);
  for (i = 0; i < n; i++)
    mpf_clear(scale[i]);
  free(scale);
  scale = NULL;
  rownum = NULL;

  return 0;
}

int LUdecomp_rational(rational_complex_matrix LU, int **rowswaps, rational_complex_matrix A)
/***************************************************************\
* USAGE: compute an LU decomposition of A                       *
\***************************************************************/
{
  int i, j, k, l, pivot, n = A->rows, *rownum = NULL;
  mpf_t max, tempMPF, *scale = NULL;
  complex_number tempComp;
  rational_complex_number tempRat;

  // error checking - need a square matrix
  if (A->rows != A->cols)
  {
    printf("ERROR: In order to compute an LU decomposition, the matrix must be square!!\n");
    errExit(ERROR_INVALID_SIZE);
  }

  // initialize memory
  mpf_init(max); mpf_init(tempMPF); 
  initialize_number(tempComp);
  initialize_rational_number(tempRat);

  // setup rowswaps & scale
  rownum = *rowswaps = (int *)errRealloc(*rowswaps, n * sizeof(int));
  scale = (mpf_t *)errMalloc(n * sizeof(mpf_t));
  for (i = 0; i < n; i++)
  { // initialize rownum[i] & scale[i]
    rownum[i] = i;
    mpf_init(scale[i]);
    mpf_set_ui(scale[i], 0);
  }

  // setup LU
  change_size_rational_matrix(LU, n, n);
  
  // find scale factors for each row and copy A to LU
  for (i = 0; i < n; i++)
  { 
    set_rational_number(LU->entry[i][0], A->entry[i][0]);
    convert_rational_number(tempComp, LU->entry[i][0]);
    norm_sqr_number(max, tempComp);
    for (j = 1; j < n; j++)
    {
      set_rational_number(LU->entry[i][j], A->entry[i][j]);
      convert_rational_number(tempComp, LU->entry[i][j]);
      norm_sqr_number(tempMPF, tempComp);
      if (mpf_cmp(max, tempMPF) < 0)
        mpf_set(max, tempMPF);
    }
    // compute sqrt of max
    mpf_sqrt(max, max);

    // reciprocate scale[i] - turn / into *
    mpf_ui_div(scale[i], 1, max); 
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find scaled maximum element in column k
    convert_rational_number(tempComp, LU->entry[rownum[k]][k]);
    norm_number(max, tempComp);
    mpf_mul(max, max, scale[rownum[k]]);
    for (i = k + 1; i < n; i++)
    {
      convert_rational_number(tempComp, LU->entry[rownum[i]][k]);
      norm_number(tempMPF, tempComp);
      mpf_mul(tempMPF, tempMPF, scale[rownum[i]]);
      if (mpf_cmp(max, tempMPF) < 0)
      {
        pivot = i;
        mpf_set(max, tempMPF);
      }
    }

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
    }

    // verify the element is nonzero
    if (mpq_cmp_ui(LU->entry[rownum[k]][k]->re, 0, 1) == 0 && mpq_cmp_ui(LU->entry[rownum[k]][k]->im, 0, 1) == 0)
    { // pivot element is zero!
      mpf_clear(max); mpf_clear(tempMPF);
      clear_number(tempComp);
      clear_rational_number(tempRat);
      for (i = 0; i < n; i++)
        mpf_clear(scale[i]);
      free(scale);
      scale = NULL;
      rownum = NULL;

      return ERROR_LU_DECOMP;
    }

    // update elements in "bottom" block of LU
    pivot = rownum[k];
    recip_rational_number(tempRat, LU->entry[pivot][k]);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];

      // compute muliplier
      multiply_rational(LU->entry[l][k], LU->entry[l][k], tempRat);

      // update elements of row i
      for (j = k + 1; j < n; j++)
        subtract_multiply_rational(LU->entry[l][j], LU->entry[l][k], LU->entry[pivot][j]);
    }
  }

  // check to make sure the last scaled pivot element is not zero
  l = n - 1;
  if (mpq_cmp_ui(LU->entry[rownum[l]][l]->re, 0, 1) == 0 && mpq_cmp_ui(LU->entry[rownum[l]][l]->im, 0, 1) == 0)
  { // pivot element is zero!
    mpf_clear(max); mpf_clear(tempMPF);
    clear_number(tempComp);
    clear_rational_number(tempRat);
    for (i = 0; i < n; i++)
      mpf_clear(scale[i]);
    free(scale);
    scale = NULL;
    rownum = NULL;

    return ERROR_LU_DECOMP;
  }

  mpf_clear(max); mpf_clear(tempMPF); 
  clear_number(tempComp);
  clear_rational_number(tempRat);
  for (i = 0; i < n; i++)
    mpf_clear(scale[i]);
  free(scale);
  scale = NULL;
  rownum = NULL;

  return 0;
}

void LUsolve(complex_matrix X, complex_matrix LU, int *rownum, complex_matrix B, int LU_prec)
/***************************************************************\
* USAGE: solve A*X = B, where A is stored in LU & rownum        *
\***************************************************************/
{
  int i, j, k, l, m, pivot, n = LU->rows, p = B->cols;
  complex_matrix Y;
  complex_number tempComp;

  if (LU->rows != LU->cols)
  {
    printf("ERROR: In order to solve using an LU decomposition, the matrix must be square!!\n");
    errExit(ERROR_INVALID_SIZE);
  }
  else if (LU->rows != B->rows)
  {
    printf("ERROR: The number of rows must match!!\n");
    errExit(ERROR_INVALID_SIZE);
  }

  // set the default precision
  setPrec(LU_prec);

  // initialize tempComp & Y
  initialize_number(tempComp);
  initialize_matrix(Y, n, p);

  // setup X - same size as B
  setPrec_matrix(X, LU_prec);
  change_size_matrix(X, n, p);

  // first compute L*Y = P*B, where P is defined in rownum
  for (k = 0; k < p; k++)
  { // calculate the kth column of Y
    set_number(Y->entry[0][k], B->entry[rownum[0]][k]);

    for (i = 1; i < n; i++)
    { // calculate Y[i][k]
      l = rownum[i];
      negate(Y->entry[i][k], B->entry[l][k]);
      for (j = 0; j < i; j++)
        sum_multiply(Y->entry[i][k], LU->entry[l][j], Y->entry[j][k]);
      negate(Y->entry[i][k], Y->entry[i][k]);
    }
  }

  // then compute U*X = Y
  l = n - 1;
  pivot = rownum[l];
  for (k = 0; k < p; k++)
  { // calculate the kth column of X
    divide(X->entry[l][k], Y->entry[l][k], LU->entry[pivot][l]);
    for (i = n - 2; i >= 0; i--)
    {
      set_number(tempComp, Y->entry[i][k]);
      m = rownum[i];
      for (j = i + 1; j < n; j++)
        subtract_multiply(tempComp, LU->entry[m][j], X->entry[j][k]);
      divide(X->entry[i][k], tempComp, LU->entry[m][i]);
    }
  }

  clear_matrix(Y);
  clear_number(tempComp);

  return;
}

void LUsolve_rational(rational_complex_matrix X, rational_complex_matrix LU, int *rownum, rational_complex_matrix B)
/***************************************************************\
* USAGE: solve A*X = B, where A is stored in LU & rownum        *
\***************************************************************/
{
  int i, j, k, l, m, pivot, n = LU->rows, p = B->cols;
  rational_complex_matrix Y;
  rational_complex_number tempComp;

  if (LU->rows != LU->cols)
  {
    printf("ERROR: In order to solve using an LU decomposition, the matrix must be square!!\n");
    errExit(ERROR_INVALID_SIZE);
  }
  else if (LU->rows != B->rows)
  {
    printf("ERROR: The number of rows must match!!\n");
    errExit(ERROR_INVALID_SIZE);
  }

  // initialize tempComp & Y
  initialize_rational_number(tempComp);
  initialize_rational_matrix(Y, n, p);

  // setup X - same size as B
  change_size_rational_matrix(X, n, p);

  // first compute L*Y = P*B, where P is defined in rownum
  for (k = 0; k < p; k++)
  { // calculate the kth column of Y
    set_rational_number(Y->entry[0][k], B->entry[rownum[0]][k]);

    for (i = 1; i < n; i++)
    { // calculate Y[i][k]
      l = rownum[i];
      negate_rational(Y->entry[i][k], B->entry[l][k]);
      for (j = 0; j < i; j++)
        sum_multiply_rational(Y->entry[i][k], LU->entry[l][j], Y->entry[j][k]);
      negate_rational(Y->entry[i][k], Y->entry[i][k]);
    }
  }


  // then compute U*X = Y
  l = n-1;
  pivot = rownum[l];
  for (k = 0; k < p; k++)
  { // calculate the kth column of X
    divide_rational(X->entry[l][k], Y->entry[l][k], LU->entry[pivot][l]);
    for (i = n - 2; i >= 0; i--)
    {
      set_rational_number(tempComp, Y->entry[i][k]);
      m = rownum[i];
      for (j = i + 1; j < n; j++)
        subtract_multiply_rational(tempComp, LU->entry[m][j], X->entry[j][k]);
      divide_rational(X->entry[i][k], tempComp, LU->entry[m][i]);
    }
  }

  clear_rational_matrix(Y);
  clear_rational_number(tempComp);

  return;
}

void LUsolve_vector(complex_vector X, complex_matrix LU, int *rownum, complex_vector B, int LU_prec)
/***************************************************************\
* USAGE: solve A*X = B, where A is stored in LU & rownum        *
\***************************************************************/
{
  int i, size = B->size;
  complex_matrix Bmat, Xmat;

  // initialize
  initialize_matrix2(Bmat, size, 1, LU_prec);
  initialize_matrix2(Xmat, size, 1, LU_prec);

  // copy B to Bmat
  for (i = 0; i < size; i++)
    set_number(Bmat->entry[i][0], B->coord[i]);

  // solve for Xmat
  LUsolve(Xmat, LU, rownum, Bmat, LU_prec);

  // copy Xmat to X
  setPrec_vector(X, LU_prec);
  change_size_vector(X, size);
  for (i = 0; i < size; i++)
    set_number(X->coord[i], Xmat->entry[i][0]);

  clear_matrix(Bmat);
  clear_matrix(Xmat);

  return;
}

void LUsolve_rational_vector(rational_complex_vector X, rational_complex_matrix LU, int *rownum, rational_complex_vector B)
/***************************************************************\
* USAGE: solve A*X = B, where A is stored in LU & rownum        *
\***************************************************************/
{
  int i, size = B->size;
  rational_complex_matrix Bmat, Xmat;

  // initialize
  initialize_rational_matrix(Bmat, size, 1);
  initialize_rational_matrix(Xmat, size, 1);

  // copy B to Bmat
  for (i = 0; i < size; i++)
    set_rational_number(Bmat->entry[i][0], B->coord[i]);

  // solve for Xmat
  LUsolve_rational(Xmat, LU, rownum, Bmat);

  // copy Xmat to X
  change_size_rational_vector(X, size);
  for (i = 0; i < size; i++)
    set_rational_number(X->coord[i], Xmat->entry[i][0]);

  clear_rational_matrix(Bmat);
  clear_rational_matrix(Xmat);

  return;
}


