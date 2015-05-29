/* file:        svd.c
** author:      Andrea Vedaldi (based on Dianne Cook's code).
** description: Definition of function dsvd.
**/

/** @file svd.h
 ** @brief SVD computation.
 **/

/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  this code and associated documentation,  provided       *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *     Copyright (c) held by Dianne Cook                    *
 *  All Rights Reserved.                                    *
 *                                                          *
 *  Questions and comments are welcome, and I request       *
 *  that you share any modifications with me.               *
 *                                                          *
 *                Dianne Cook                               *
 *             dicook@iastate.edu                           *
 *                                                          *
 ************************************************************/

#include"svd.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <vector>
#include "Matrix.h"

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double PYTHAG(double a, double b)
{
  double at = fabs(a), bt = fabs(b), ct, result;
  
  if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
  else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
  else result = 0.0;
  return(result);
}

/*
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 */

int 
_svd(float **a, int m, int n, float *w, float **v)
{
  int flag, i, its, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  double *rv1;
  
  if (m < n) 
    {
      /* fprintf(stderr, "#rows must be >=  #cols \n"); */
      return(0);
    }
  
  rv1 = (double *)malloc((unsigned int) n*sizeof(double));
  
  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++) 
    {
      /* left-hand reduction */
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = 0.0;
      if (i < m) 
        {
          for (k = i; k < m; k++) 
            scale += fabs((double)a[k][i]);
          if (scale) 
            {
              for (k = i; k < m; k++) 
                {
                  a[k][i] = (float)((double)a[k][i]/scale);
                  s += ((double)a[k][i] * (double)a[k][i]);
                }
              f = (double)a[i][i];
              g = -SIGN(sqrt(s), f);
              h = f * g - s;
              a[i][i] = (float)(f - g);
              if (i != n - 1) 
                {
                  for (j = l; j < n; j++) 
                    {
                      for (s = 0.0, k = i; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                      f = s / h;
                      for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
              for (k = i; k < m; k++) 
                a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
      w[i] = (float)(scale * g);
      
      /* right-hand reduction */
      g = s = scale = 0.0;
      if (i < m && i != n - 1) 
        {
          for (k = l; k < n; k++) 
            scale += fabs((double)a[i][k]);
          if (scale) 
            {
              for (k = l; k < n; k++) 
                {
                  a[i][k] = (float)((double)a[i][k]/scale);
                  s += ((double)a[i][k] * (double)a[i][k]);
                }
              f = (double)a[i][l];
              g = -SIGN(sqrt(s), f);
              h = f * g - s;
              a[i][l] = (float)(f - g);
              for (k = l; k < n; k++) 
                rv1[k] = (double)a[i][k] / h;
              if (i != m - 1) 
                {
                  for (j = l; j < m; j++) 
                    {
                      for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[j][k] * (double)a[i][k]);
                      for (k = l; k < n; k++) 
                        a[j][k] += (float)(s * rv1[k]);
                    }
                }
              for (k = l; k < n; k++) 
                a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
      anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
  /* accumulate the right-hand transformation */
  for (i = n - 1; i >= 0; i--) 
    {
      if (i < n - 1) 
        {
          if (g) 
            {
              for (j = l; j < n; j++)
                v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
              /* double division to avoid underflow */
              for (j = l; j < n; j++) 
                {
                  for (s = 0.0, k = l; k < n; k++) 
                    s += ((double)a[i][k] * (double)v[k][j]);
                  for (k = l; k < n; k++) 
                    v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
          for (j = l; j < n; j++) 
            v[i][j] = v[j][i] = 0.0;
        }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }
  
  /* accumulate the left-hand transformation */
  for (i = n - 1; i >= 0; i--) 
    {
      l = i + 1;
      g = (double)w[i];
        if (i < n - 1) 
          for (j = l; j < n; j++) 
            a[i][j] = 0.0;
        if (g) 
          {
            g = 1.0 / g;
            if (i != n - 1) 
            {
              for (j = l; j < n; j++) 
                {
                  for (s = 0.0, k = l; k < m; k++) 
                    s += ((double)a[k][i] * (double)a[k][j]);
                  f = (s / (double)a[i][i]) * g;
                  for (k = i; k < m; k++) 
                    a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
              a[j][i] = (float)((double)a[j][i]*g);
          }
        else 
          {
            for (j = i; j < m; j++) 
              a[j][i] = 0.0;
          }
        ++a[i][i];
    }
  
  /* diagonalize the bidiagonal form */
  for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
      for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
          flag = 1;
          for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                  {
                    flag = 0;
                    break;
                  }
                if (fabs((double)w[nm]) + anorm == anorm) 
                  break;
            }
          if (flag) 
            {
              c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                  {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                      {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                          {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                          }
                      }
                  }
            }
          z = (double)w[k];
          if (l == k) 
            {                  /* convergence */
              if (z < 0.0) 
                {              /* make singular value nonnegative */
                  w[k] = (float)(-z);
                  for (j = 0; j < n; j++) 
                    v[j][k] = (-v[j][k]);
                }
              break;
            }
          if (its >= 30) {
                free((void*) rv1);
                /*      fprintf(stderr, "No convergence after 30,000! iterations \n");*/
                return(0);
          }
          
          /* shift from bottom 2 x 2 minor */
          x = (double)w[l];
          nm = k - 1;
          y = (double)w[nm];
          g = rv1[nm];
          h = rv1[k];
          f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
          g = PYTHAG(f, 1.0);
          f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
          /* next QR transformation */
          c = s = 1.0;
          for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                  {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                  }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                  {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                  }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                  {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                  }
            }
          rv1[l] = 0.0;
          rv1[k] = f;
          w[k] = (float)x;
        }
    }
  free((void*) rv1);
  return(1);
}

float**
_alloc_matrix(int M, int N)
{
  int i;
  float** result = (float**) malloc(sizeof(float*) * M) ;
  for(i = 0 ; i < M ; ++i) 
    result[i] = (float*) malloc(sizeof(float) * N) ;
  return result ;
}

void
_free_matrix(float** A, int M)
{
  int i;
  for(i = 0 ; i < M ; ++i) 
    free(A[i]) ;
  free(A) ;
}

void
_copy_fortran_matrix_to_matrix(float** A, const float* B, int M, int N)
{
  int i, j ;
  for(i = 0 ; i < M ; ++i) {
    for(j = 0 ; j < N ; ++j) {
      A[i][j] = *(B + i + j*M) ;
     }
  }
}

void
_copy_matrix_to_fortran_matrix(float*B, float** A, int M, int N)
{
  int i, j ;  
  for(i = 0 ; i < M ; ++i) {
    for(j = 0 ; j < N ; ++j) {
      *(B + i + j*M) = A[i][j] ;
    }
  }
}

#define _max(a,b) (((a)>(b))?(a):(b))
#define _min(a,b) (((a)<(b))?(a):(b))

/** @brief SVD decomposition routine.
 **
 ** The function  takes the MxN matrix  A, (M >= N)  and decomposes it
 ** into  U  D  V,  where  U,V  are  the  left  and  right  orthogonal
 ** transformation matrices,  and D is  a diagonal matrix  of singular
 ** values.
 **
 ** Matrix U has the same dimensions as A and its colums span only the
 ** range of A.  Matrix U can  be made MxM square by adding arbitrarly
 ** other orghogonal columns.
 **
 ** The algorithm  can only handle  the case M  >= N.  If you  need to
 ** deal with the  other case, you can compute the  SVD of A transpose
 ** and then adjust the result.
 **
 ** @param A MxN matrix to be decomposed
 ** @param M row dimension of A
 ** @param N column dimension of A
 ** @param U MxN matrix of the left orthogonal transformation U
 ** @param S N vector of the singular values of A
 ** @param V NxN matrix of the right orthogonal transformation V
 **
 ** @return 1 if succes, 0 otherwise.
 **/
int 
svd(const Matrix& A, int M, int N, Matrix& U, Matrix& S, Matrix& V)
{
  float* A_array = (float*) malloc(sizeof(float) * M * N);
  float* U_array = (float*) malloc(sizeof(float) * M * N);
  float* S_array = (float*) malloc(sizeof(float) * N);
  float* V_array = (float*) malloc(sizeof(float) * N * N);

   // Copy data from our A matrix into array
   const std::vector<double>& A_data = A.getData();
   for ( int i = 0; i < M * N; ++i )
   {
     A_array[i] = (float) A_data[i];
   }
        
   float** AA = _alloc_matrix(M,N) ;
   float** VV = _alloc_matrix(N,N) ;
   int success ;
   _copy_fortran_matrix_to_matrix(AA, A_array, M, N) ;

   success = _svd(AA,M,N,S_array,VV) ;

   _copy_matrix_to_fortran_matrix(U_array,AA,M,N) ;
   _copy_matrix_to_fortran_matrix(V_array,VV,N,N) ;

   _free_matrix(AA,M) ;
   _free_matrix(VV,N) ;
   
   // Copy data to our matrices
   for ( int i = 0; i < M * N; ++i )
   {
     U.getData(i) = U_array[i];
   }
   for ( int i = 0; i < N; ++i )
   {
     S.getData(i) = S_array[i];
   }
   for ( int i = 0; i < N * N; ++i )
   {
     V.getData(i) = V_array[i];
   }
      
   free(A_array);
   free(U_array);
   free(S_array);
   free(V_array);
   
   return success ;
}
