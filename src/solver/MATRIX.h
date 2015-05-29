///////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2002 - 2014, Huamin Wang
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////
//  MATRIX
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef __WHMIN_MATRIX_H__
#define __WHMIN_MATRIX_H__
#include <math.h>
#include "SOLVER_MATH.h"
#include <Eigen/Dense>
#include "vector_io.h"
#include "print_macro.h"
namespace solver {
template <class T>
class MATRIX {
public:
  int ni;
  int nj;
  int size;
  T A[4096];  //The matrix cannot be larger than 4096

  //**************************************************************************************
  // Initialize the matrix to certain values
  //**************************************************************************************
  MATRIX() {}
  MATRIX(T* v, int _ni, int _nj) {
    ni = _ni;
    nj = _nj;
    size = ni * nj;
    memcpy(A, v, sizeof(double) * size);
  }

  void Set(T value, int _ni, int _nj) {
    ni = _ni;
    nj = _nj;
    size = ni * nj;
    for (int id = 0; id < ni * nj; id++) A[id] = value;
  }

  void SetZero() {
    memset(A, 0, sizeof(T) * size);
  }

  void Set(int _ni, int _nj) {
    ni = _ni;
    nj = _nj;
    size = ni * nj;
    for (int id = 0; id < ni * nj; id++) A[id] = RandomFloat() * 2 - 1;
  }

  void Set(T value0, T value1, int _nj) {
    ni = _nj;
    nj = _nj;
    size = ni * nj;
    for (int i = 0; i < ni; i++) {
      A[i * nj + i] = value0;
      for (int j = i + 1; j < nj; j++)
        A[i * nj + j] = A[j * nj + i] = value1;
    }
  }

  //**************************************************************************************
  // Get the matrix from another one
  //**************************************************************************************

  void Get(MATRIX<T> &m) {
    ni = m.ni;
    nj = m.nj;
    size = ni * nj;
    memcpy(A, m.A, sizeof(T)*ni * nj);
  }

  //**************************************************************************************
  // Print the matrix for debugging
  //**************************************************************************************

  void Print() {
    for (int i = 0; i < ni; i++) {
      for (int j = 0; j < nj; j++) printf("%f; ", A[i * nj + j]);
      printf("\n");
    }
  }

  //**************************************************************************************
  // Cholesky_Decomposition

  bool Cholesky_Decomposition() {
    //    double tmp[4096];
    //    memcpy(tmp, A, sizeof(double) * ni * nj);
    if (ni != nj)  printf("ERROR: cannot decompose a non-square matrix.\n");
    for (int i = 0; i < ni; i++) {
      for (int j = i; j < nj; j++) {
        T sum = A[i * nj + j];
        for (int k = 0; k < i; k++) sum -= A[k * nj + i] * A[k * nj + j];

        if (i == j) {
          if (sum <= 0) {
            //            P(sum);
            //            typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
            //            Eigen::Map<Mat> mat(tmp, ni, nj);
            //            WriteEigenMatrixToMatlab(mat, "/tmp/mat");
            //            exit(0);
            return false;
          }
          A[i * nj + i] = sqrt(sum);
        } else A[i * nj + j] = sum / A[i * nj + i];
      }
    }
    //Set the lower part to zeros
    for (int i = 0; i < ni; i++)
      for (int j = i + 1; j < nj; j++)
        A[j * nj + i] = 0;
    return true;
  }

  //**************************************************************************************
  // Square_Tranpose
  //**************************************************************************************
  void Square_Tranpose() {
    if (ni != nj)  printf("ERROR: cannot handle non-square matrix in Square_Tranpose.\n");
    for (int i = 0; i < ni; i++)
      for (int j = i + 1; j < nj; j++) {
        T temp = A[i * ni + j];
        A[i * ni + j] = A[j * ni + i];
        A[j * ni + i] = temp;
      }
  }

  //**************************************************************************************
  // Upper triangular Matrix Inverse
  // This function computes the inverse of the matrix, if it is an upper triangular matrix
  //**************************************************************************************

  void Solve(T x[]) {
    int n = ni;
    for (int i = n - 1; i >= 0; i--) {
      for (int j = i + 1; j < n; j++) x[i] -= A[i * n + j] * x[j];
      x[i] /= A[i * n + i];
    }
  }

  void Upper_Inverse(MATRIX<T> &mr) {
    if (ni != nj)  printf("ERROR: cannot handle non-square matrix in Upper_Inverse.\n");
    mr.ni = ni;
    mr.nj = nj;
    mr.Set(T(1.0), T(0.0), ni); //set mr to identity
    for (int i = 0; i < ni; i++)
      Solve(&mr.A[i * ni]);
    mr.Square_Tranpose();
  }

  //**************************************************************************************
  // Basic arithematic operations
  //**************************************************************************************

  static void Add(MATRIX<T> &m0, MATRIX<T> &m1, MATRIX<T> &mr) {      //mr=m0+m1
    if (m0.ni != m1.ni || m0.nj != m1.nj)  printf("ERROR: Dimensions are not match in Add.\n");
    mr.ni = m0.ni;
    mr.nj = m0.nj;
    for (int id = 0; id < m0.ni * m0.nj; id++)
      mr.A[id] = m0.A[id] + m1.A[id];
  }

  static void Subtract(MATRIX<T> &m0, MATRIX<T> &m1, MATRIX<T> &mr) {   //mr=m0-m1
    if (m0.ni != m1.ni || m0.nj != m1.nj)  printf("ERROR: Dimensions are not match in Subtract.\n");
    mr.ni = m0.ni;
    mr.nj = m0.nj;
    for (int id = 0; id < m0.ni * m0.nj; id++)
      mr.A[id] = m0.A[id] - m1.A[id];
  }

  static void Multiply(MATRIX<T> &m0, MATRIX<T> &m1, MATRIX<T> &mr) {   //mr=m0*m1
    if (m0.nj != m1.ni)  printf("ERROR: Dimensions are not match in Multiply.\n");
    mr.ni = m0.ni;
    mr.nj = m1.nj;
    memset(mr.A, 0, sizeof(T)*mr.ni * mr.nj);
    for (int i = 0; i < mr.ni; i++)
      for (int j = 0; j < mr.nj; j++)
        for (int k = 0; k < m0.nj; k++)
          mr.A[i * mr.nj + j] += m0.A[i * m0.nj + k] * m1.A[k * m1.nj + j];
  }

  static void Transpose_Multiply(MATRIX<T> &m0, MATRIX<T> &m1, MATRIX<T> &mr) { //mr=m0'*m1
    if (m0.ni != m1.ni)  printf("ERROR: Dimensions are not match in Multiply.\n");
    mr.ni = m0.nj;
    mr.nj = m1.nj;
    memset(mr.A, 0, sizeof(T)*mr.ni * mr.nj);
    for (int k = 0; k < m0.ni; k++)
      for (int i = 0; i < mr.ni; i++)
        for (int j = 0; j < mr.nj; j++)
          mr.A[i * mr.nj + j] += m0.A[k * m0.nj + i] * m1.A[k * m1.nj + j];
  }

  static void Multiply(MATRIX<T> &m0, T x[], T r[]) {     //r=m0*x
    memset(r, 0, sizeof(T)*m0.ni);
    for (int i = 0; i < m0.ni; i++)
      for (int j = 0; j < m0.nj; j++)
        r[i] += m0.A[i * m0.nj + j] * x[j];
  }

  static void Transpose_Multiply(MATRIX<T> &m0, T x[], T r[]) { //r=m0'*x
    memset(r, 0, sizeof(T)*m0.nj);
    for (int j = 0; j < m0.nj; j++)
      for (int i = 0; i < m0.ni; i++)
        r[j] += m0.A[i * m0.nj + j] * x[i];
  }

};
}

#endif
