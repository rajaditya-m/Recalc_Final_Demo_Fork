/* file:        inverse.cpp
** author:      Andrea Vedaldi
** description: Example of usage of the math routines.
**/

#include "inverse.h"
#include "svd.h"
#include <stdlib.h>
#include "Matrix.h"

/** @file inverse.h
 ** @brief Matrix inverse computation.
 **/

const double epsilon = 1E-16 ;

/** @brief Computes the inverse of a matrix.
 **
 ** The results is written into @c B.
 **
 ** @param A A square matrix.
 ** @param N Number of rows/columns.
 ** @param B A square matrix of the same dimension of A.
 ** @return 0 if the matrix is found singular; 1 otherwise.
 **/
int
inverse(const Matrix& A, int N,  Matrix& B)
{
  Matrix U = B;
  Matrix S(1, N);
  Matrix V(N, N);
  Matrix T(N, N);

  /*float* U = B ;
  float* S = malloc(N*sizeof(float)) ;
  float* V = malloc(N*N*sizeof(float)) ;
  float* T = malloc(N*N*sizeof(float)) ;*/
  
  int singular = 0 ;
  int i,j,k ;

  /* Takes the SVD. */
  svd(A, N, N, U, S, V) ;

  /* Check wether A is singular. */
  for(i = 0 ; i < N ; ++i) 
    singular |= 
      (-epsilon < S(0,i)) || 
      ( S(0,i) < epsilon) ;

  /* T = inv(S)*U'. */
  for(j = 0 ; j < N ; ++j)
    for(i = 0 ; i < N ; ++i) 
      T(i,j) = U(j,i)/S(0,i) ;

  /* B = U = V*T .*/
  for(j = 0 ; j < N ; ++j) {
    for(i = 0 ; i < N ; ++i) {
      U(i,j) = 0.0 ;
      for(k = 0 ; k < N ; ++k) {
        U(i,j) +=  V(i,k) * T(k,j) ;     
      }
    }
  }
  B = U;

  //free(V) ;
  //free(S) ;
  //free(T) ;
  return ! singular ;
}
