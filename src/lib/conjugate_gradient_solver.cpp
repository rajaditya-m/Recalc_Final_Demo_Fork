#include "conjugate_gradient_solver.h"

template <class Float>
typename ConjugateGradientSolver<Float>::Info ConjugateGradientSolver<Float>::Solve(Float *rhs,
                                                                                    Float *result,
                                                                                    typename ConjugateGradientSolver::Matrix A,
                                                                                    int max_iteration, Float residual_threshold)
{
  Float current_residual = 0;
  max_iteration = (max_iteration == 0) ? kDefaultMaxIteration : max_iteration;
  //  memset(result, 0, sizeof(Float) * size_); // result = 0
  A(result, temp_);
  memcpy(residual_, rhs, sizeof(Float) * size_);  // p = rhs
  blas::axpy(size_, -1.0f, temp_, residual_);
  memcpy(search_direction_, residual_, sizeof(Float) * size_);  // p = rhs
  Float denumerator;
  blas::dot(size_, residual_, residual_, denumerator);
  Float residual_2_norm = denumerator;
  int iteration = 0 ;
  ASSERT(denumerator == denumerator);
//  P(denumerator * normalizer_);
  if (denumerator * normalizer_ > residual_threshold) {
    for (iteration = 1; iteration <= max_iteration; iteration++) {
      A(search_direction_, temp_) ; // temp = A * p
      Float alpha;
      blas::dot(size_, search_direction_, temp_, alpha);
      alpha = denumerator / alpha ;  // alpha = r*r/p*A*p
      blas::axpy(size_, alpha, search_direction_, result); // x = x+alpha*p
      blas::axpy(size_, -alpha, temp_, residual_);  // r = r-alpha*A*p
      blas::dot(size_, residual_, residual_, residual_2_norm) ;  // residual = r*r
      if (residual_2_norm != residual_2_norm) {
        std::cerr << CURRENT_LINE << " => " << iteration << " iterations with residual " << residual_2_norm << std::endl;
        getchar();
      }
      current_residual = residual_2_norm * normalizer_;
      if (current_residual < residual_threshold) {
//        std::cout << CURRENT_LINE << " => " << iteration << " iterations with residual " << current_residual << std::endl;
        break ;
      }
      Float beta = residual_2_norm / denumerator;
      blas::scal(size_, beta, search_direction_) ; // p = r+beta*p
      blas::axpy(size_, 1.0 , residual_, search_direction_);
      denumerator = residual_2_norm;
    }
    if (iteration > max_iteration) {
      std::cerr << CURRENT_LINE << " => warning: maximum " << max_iteration << " iteration reached with residual " << residual_2_norm << std::endl;
    }
  } else {
//    std::cout << CURRENT_LINE << " => " << " 0 iterations with residual " << current_residual << std::endl;
  }
  return Info(iteration, current_residual) ;
}

template class ConjugateGradientSolver<float>;
template class ConjugateGradientSolver<double>;
