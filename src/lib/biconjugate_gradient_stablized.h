#pragma once
#include <Eigen/Dense>
#include <functional>

template <class Float>
class BiconjugateGradientStablized
{
public:
  typedef typename Eigen::Matrix<Float, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Map<Vec> MapVec;
  typedef std::function<void (Float*, Float*)> Matrix;
  typedef std::pair<int, Float> Info; // <iteration to solve, residual>

  BiconjugateGradientStablized(int size)
    : r0(size)
    , r1(size)
    , p0(size)
    , p1(size)
    , v0(size)
    , v1(size)
    , s(size)
    , t(size)
    , r0_hat(size)
    , residual(size)
    , size_(size)
  {
  }

  ~BiconjugateGradientStablized() {}

  int size() {
    return size_;
  }

  void Resize(int size) {
    if (size != size_) {
      size_ = size;
      r0 = Vec(size_);
      r0_hat = Vec(size_);
      residual = Vec(size_);
      s = Vec(size_);
      t = Vec(size_);
      r1 = Vec(size_);
      p0 = Vec(size_);
      p1 = Vec(size_);
      v0 = Vec(size_);
      v1 = Vec(size_);
    }
  }

/*
 From: http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
  r0 = b − Ax0
  Choose an arbitrary vector r̂0 such that (r̂0, r0) ≠ 0, e.g., r̂0 = r0
  ρ0 = α = ω0 = 1
  v0 = p0 = 0
  For i = 1, 2, 3, …
  ρi = (r̂0, ri−1)
  β = (ρi/ρi−1)(α/ωi−1)
  pi = ri−1 + β(pi−1 − ωi−1vi−1)
  vi = Api
  α = ρi/(r̂0, vi)
  s = ri−1 − αvi
  t = As
  ωi = (t, s)/(t, t)
  xi = xi−1 + αpi + ωis
  If xi is accurate enough then quit
  ri = s − ωit
*/
  Info Solve(Float* b, Float* x, Matrix A, int max_iteration, Float residual_threshold) {
    MapVec map_b(b, size_);
    MapVec map_x(x, size_);
    A(x, &r0[0]);
    r0 = map_b - r0;
    Float error = r0.dot(r0) / size_;
    if (error < residual_threshold) {
      return Info(0, 0);
    }
    r0_hat = r0;
    Float rho0 = 1, alpha = 1, omega0 = 1;
    v0.setZero();
    p0.setZero();
    int iter = 1;
    for (; iter <= max_iteration; ++iter) {
      Float rho1 = r0_hat.dot(r0);
      Float beta = (rho1 / rho0) * (alpha / omega0);
      p1 = r0 + beta * (p0 - omega0 * v0);
      A(&p1[0], &v1[0]);
      alpha = rho1 / (r0_hat.dot(v1));
      s = r0 - alpha * v1;
      A(&s[0], &t[0]);
      Float omega1 = t.dot(s) / t.dot(t);
      map_x = map_x + alpha * p1 + omega1 * s;
      A(&map_x[0], &residual[0]);
      residual -= map_b;
      error = residual.dot(residual) / size_;
      if (error < residual_threshold) {
        return Info(iter, error);
      }
      r1 = s - omega1 * t;
      rho0 = rho1;
      omega0 = omega1;
      v0 = v1;
      p0 = p1;
      r0 = r1;
    }
    return Info(iter, error);
  }

private:
  Vec r0, r1, p0, p1, v0, v1, s, t;
  Vec r0_hat, residual;
  int size_;
};
