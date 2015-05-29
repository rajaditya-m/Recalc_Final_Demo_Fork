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
//  BLOCK_MATRIX_GRAPH
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef __WHMIN_BLOCK_MATRIX_GRAPH_H__
#define __WHMIN_BLOCK_MATRIX_GRAPH_H__
#include <vector>
#include <iostream>
#include <math.h>
#include "MATRIX.h"
#include <assert.h>
#include <climits>
//#include "print_macro.h"

namespace solver {
template <class T>
class BLOCK_MATRIX_GRAPH {
public:
  enum {
    kSuccess = -1
  };
  int max_number;
  int max_e_number;

  MATRIX<T>  *Av;   //the matrices defined at the vertices    (diagonal)
  MATRIX<T>  *Ae;   //the upper matrices defined at the edges (off-diagonal)
  int  *E;      //edge list. (low vertex goes first)
  int   (*VE)[64];  //the adjacent edges of the vertices
  int  *ve_number;

  int size_;
  std::vector<int> vv2edge_;
  int e_number0;

  int   number;   //vertex number
  int   e_number; //edge number

  BLOCK_MATRIX_GRAPH(std::vector<std::vector<int> >& topology, std::vector<int>& block_size, int additional_entry = 0);
  BLOCK_MATRIX_GRAPH(int block_num, std::vector<int> edges, std::vector<int> block_size);

  BLOCK_MATRIX_GRAPH() {
    max_number = 32;
    max_e_number = 128;
    Av  = new MATRIX<T>[max_number];
    Ae  = new MATRIX<T>[max_e_number];
    E = new int[max_e_number * 2];
    VE  = new int[max_number][64];
    ve_number = new int[max_number];

    //A simple test case
    // 1 ---- 0 ---- 4
    // |      |
    // 2 ---- 3
    number  = 5;
    e_number = 5;
    E[0] = 0;
    E[1] = 1;
    E[2] = 1;
    E[3] = 2;
    E[4] = 2;
    E[5] = 3;
    E[6] = 0;
    E[7] = 3;
    E[8] = 0;
    E[9] = 4;

    Av[0].Set(8.0 + RandomFloat(), 1.0, 2);
    Av[1].Set(8.0 + RandomFloat(), 1.0, 3);
    Av[2].Set(8.0 + RandomFloat(), 1.0, 2);
    Av[3].Set(8.0 + RandomFloat(), 1.0, 3);
    Av[4].Set(8.0 + RandomFloat(), 1.0, 3);

    //Ae[0].Set(1.0, 2, 3);
    //Ae[1].Set(2.0, 3, 2);
    //Ae[2].Set(3.0, 2, 3);
    //Ae[3].Set(2.0, 2, 3);
    //Ae[4].Set(3.0, 2, 3);

    Ae[0].Set(2, 3);
    Ae[1].Set(3, 2);
    Ae[2].Set(2, 3);
    Ae[3].Set(2, 3);
    Ae[4].Set(2, 3);

    BuildVE();
  }

  void BuildVE() {
    //Build VE automatically
    for (int i = 0; i < number; i++)
      ve_number[i] = 0;
    for (int e = 0; e < e_number; e++) {
      int vi = E[e * 2 + 0];
      int vj = E[e * 2 + 1];
      VE[vi][ve_number[vi]++] = e;
      VE[vj][ve_number[vj]++] = e;
    }
  }

  ~BLOCK_MATRIX_GRAPH() {
    delete[] Av;
    delete[] Ae;
    delete[] E;
    delete[] VE;
    delete[] ve_number;
  }

  void ResetTopology() {
    e_number = e_number0;
    BuildVE();
  }

  void SetMatrixZero() {
    for (int i = 0; i < number; ++i) {
      Av[i].SetZero();
    }
    for (int i = 0; i < e_number0; ++i) {
      Ae[i].SetZero();
    }
  }

  void Get(BLOCK_MATRIX_GRAPH<T> &g) {
    number  = g.number;
    e_number = g.e_number;

    //copy the matrices
    for (int i = 0; i < number; i++)
      Av[i].Get(g.Av[i]);
    for (int e = 0; e < e_number; e++)
      Ae[e].Get(g.Ae[e]);

    //copy the edge list
    memcpy(E, g.E, sizeof(int)*e_number * 2);

    BuildVE();
  }


  //**************************************************************************************
  //  Cholesky_Decomposition
  //  The vertex matrices will store the inverse of its decomposed upper triangular matrix
  //**************************************************************************************
  int Cholesky_Decomposition();
  int SolveWithDecomposedMatrix(T b[], T x[]);

  //**************************************************************************************
  //  Solve the matrix system x=A\b
  //**************************************************************************************
  int Solve(T b[], T x[]);


  //**************************************************************************************
  //  Multiply: perform matrix-vector product
  //**************************************************************************************
  void Multiply(T x[], T r[]);

  //**************************************************************************************
  //  Solve the matrix system A\x by Conjugate Gradient (without conditioning)
  //**************************************************************************************
  std::pair<int, T> CG_Solve(T b[], T x[], int max_iteration = 300, T residual_threshold = 1e-6) {
    int length = 0;
    for (int i = 0; i < number; i++) length += Av[i].ni;

    T  *p = new T[length];
    T  *r = new T[length];
    T *ap = new T[length];

    Multiply(x, r);
    for (int i = 0; i < length; i++) {
      r[i] = b[i] - r[i];
      p[i] = r[i];
    }

    T rr_dot, pap_dot;
    rr_dot = 0;
    for (int i = 0; i < length; i++) rr_dot += r[i] * r[i];

    T avg_error = 0;
    int k;
    for (k = 0; k < max_iteration; k++) {
      //Get ap
      Multiply(p, ap);
      //Get alpha
      pap_dot = 0;
      for (int i = 0; i < length; i++) pap_dot += p[i] * ap[i];
      T alpha = rr_dot / pap_dot;
      //Update x and r
      for (int i = 0; i < length; i++) {
        x[i] += alpha * p[i];
        r[i] -= alpha * ap[i];
      }
      //Update rr_dot
      T new_rr_dot = 0;
      for (int i = 0; i < length; i++) new_rr_dot += r[i] * r[i];
      avg_error = (new_rr_dot / length);
      if (avg_error < residual_threshold) {
        std::cout << "CG: iteration " << k << " with residual " << new_rr_dot << std::endl;
        break;
      }
      //Get beta and update p
      T beta = new_rr_dot / rr_dot;
      rr_dot = new_rr_dot;
      for (int i = 0; i < length; i++) p[i] = r[i] + beta * p[i];
    }
    if (k == max_iteration)  std::cout << "ERROR: CG_Solver fails to converge.\n" << std::endl;

    delete[] ap;
    delete[] p;
    delete[] r;
    return std::make_pair(k, avg_error);
  }

  void BuildEdgeIndex();

  inline int GetEdgeIdx(int v0, int v1) {
    return vv2edge_[v0 * number + v1];
  }

};


}

#endif
