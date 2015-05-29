#ifndef CONJUGATE_GRADIENT_SOLVER_H_
#define CONJUGATE_GRADIENT_SOLVER_H_
#include <iostream>
#include <functional>
#include <iostream>
#include <algorithm>
#include "blas_wrapper.h"
#include "buffer_manager.h"
#include "print_macro.h" // CURRENT_LINE

#include "blas_wrapper.h"
//#include "conjugate_gradient_solver.h"
//#include "macro.h"
// typedef float Float;

static const int kDefaultMaxIteration = 200;
static const float kDefaultResidualThreadshold = float(1e-16);

template <class Float>
class ConjugateGradientSolver
{
public:
  typedef std::function<Float (Float, Float)> Transformer;
  typedef std::function<bool (int)> Filter;
  typedef std::function<void (Float*, Float*)> Matrix;
  typedef std::pair<int, Float> Info; // <iteration to solve, residual>
private:
  Float *residual_, *search_direction_, *temp_ ;
  Float normalizer_;
  Filter filter_;
  int size_ ;
  BufferManager buffer_manager;
  // Disallow assignment of one solver to another
  ConjugateGradientSolver();
  ConjugateGradientSolver(ConjugateGradientSolver&);
  const ConjugateGradientSolver& operator=(ConjugateGradientSolver&);
  inline void ForEach(Float* v1, Float* v2, Float* out, Transformer transform)
  {
    for (int i = 0; i < size_; i++) {
      if (filter_(i)) {
        out[i] = transform(v1[i], v2[i]);
      }
    }
  }
  inline Float Dot(Float* v1, Float* v2)
  {
    Float dot = 0;
    for (int i = 0; i < size_; i++) {
      if (filter_(i)) {
        dot += v1[i] * v2[i];
      }
    }
    return dot;
  }
public:

  ConjugateGradientSolver(int size) : size_(size)
  {
#ifdef MKL
    mkl_set_num_threads(4);
#endif
    normalizer_ = Float(1.0) / size;
    residual_ = buffer_manager.Malloc<Float>(size_);
    search_direction_ = buffer_manager.Malloc<Float>(size_);
    temp_ = buffer_manager.Malloc<Float>(size_);
  }

  ~ConjugateGradientSolver() {}

  void Resize(int size)
  {
    size_ = size;
  }

  Float set_normalizer(int n)
  {
    normalizer_ = Float(1.0) / n;
    return normalizer_;
  }

  // Compute residual = rhs - A * x
  inline void ComputeResidual(Float* x, Float* rhs, Matrix A, Float* residual)
  {
    A(x, residual);
    for (int i = 0; i < size_; i++) {
      residual[i] = rhs[i] - residual[i];
    }
  }

  Float ComputeResidual2Norm(Float* x, Float* rhs, Matrix A)
  {
    ComputeResidual(x, rhs, A, temp_);
    Float residual_2_norm;
    blas::dot(size_, temp_, temp_, residual_2_norm);
    return residual_2_norm;
  }

  //	A * result = rhs
  Info Solve(Float* rhs,
             Float* result,
             Matrix A,
             Filter filter,
             int max_iteration = kDefaultMaxIteration,
             Float residual_threshold = kDefaultResidualThreadshold)
  {
    Float current_residual = 0;
    max_iteration = (max_iteration == 0) ? kDefaultMaxIteration : max_iteration;
    filter_ = filter;
    memset(temp_, 0, sizeof(Float) * size_);
    memset(result, 0, sizeof(Float) * size_); // result = 0
    memcpy(residual_, rhs, sizeof(Float) * size_); // r = rhs - A*0 = rhs
    memcpy(search_direction_, residual_, sizeof(Float) * size_);  // p = rhs
    Float denumerator;
    denumerator = Dot(residual_, residual_);
    Float residual_2_norm = denumerator;
    int iteration = 0 ;
    if (denumerator < -0.0000001 || denumerator > 0.00000001) {
      for (iteration = 1; iteration <= max_iteration; iteration++) {
        A(search_direction_, temp_) ; // temp = A * p
        Float alpha = Dot(search_direction_, temp_);
        alpha = denumerator / alpha ;  // alpha = r*r/p*A*p
        ForEach(search_direction_, result, result, [&](Float a, Float b) { return alpha * a + b; } ); // result += alpha * search_direction_;
        ForEach(temp_, residual_, residual_, [&](Float a, Float b) { return -alpha * a + b; } ); // residual_ -= alpha * A * search_direction_;
        residual_2_norm = Dot(residual_, residual_);
        if (residual_2_norm != residual_2_norm) {
          std::cerr << CURRENT_LINE << " => " << iteration << " iterations with residual " << residual_2_norm << std::endl;
          getchar();
        }
        current_residual = residual_2_norm * normalizer_;
        if (current_residual < residual_threshold) {
//          std::cout << CURRENT_LINE << " => " << iteration << " iterations with residual " << residual_2_norm <<std::endl;
          break ;
        }
        Float beta = residual_2_norm / denumerator;
        ForEach(residual_, search_direction_, search_direction_, [&](Float a, Float b) { return a + b * beta;} );  // p = r+beta*p
        denumerator = residual_2_norm;
      }
      if (iteration > max_iteration) {
        std::cerr << CURRENT_LINE << " => warning: maximum " << max_iteration << " iteration reached with residual " << residual_2_norm << std::endl;
      }
    }
    return Info(iteration, current_residual) ;
  }
  //	A * result = rhs
  Info Solve(Float* rhs,
             Float* result,
             Matrix A,
             int max_iteration = kDefaultMaxIteration,
             Float residual_threshold = kDefaultResidualThreadshold);
};
extern template class ConjugateGradientSolver<float>;
extern template class ConjugateGradientSolver<double>;

#endif // ConjugateGradientSolver
