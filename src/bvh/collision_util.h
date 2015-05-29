#ifndef COLLISION_UTIL_H
#define COLLISION_UTIL_H
#include <iostream>
#include <math.h>
namespace collision_util
{

template <class T>
T Abs(T val) {
  return (val >= 0) ? val : -val;
}

template <class T>
inline void SubVec3(const T* a, const T* b, T* result) {
  result[0] = a[0] - b[0];
  result[1] = a[1] - b[1];
  result[2] = a[2] - b[2];
}

/// Compute v1 X v2
template <class T>
inline void Cross3(const T* v1, const T* v2, T* cross) {
  cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

template <class T>
inline T Dot3(const T* v1, const T* v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/**
 * @brief Checks if a line segment intersects with a triangle
 *        based on paper "Fast, Minimum Storage Ray/Triangle Intersection"
 * @param segment_start starting position of the segment
 * @param seent_end ending position of the segment
 * @param tri_v0 triangle vert 0
 * @param tri_v1 triangle vert 1
 * @param tri_v2 triangle vert 2
 * @param threshold error tolerance for intersection
 * @param u tri_v0 * (1 - u - v) + tri_v1 * u + tri_v2 * v == intersection point
 * @param v
 * @param t segment_start * (1 - t) + segment_start * t == intersection point
 * @return true if there is intersection, false otherwise
 */
template <class T>
bool SegmentTriangleIntersection(const T* o, const T* e, const T* v0, const T* v1, const T* v2, T threshold, T& u, T& v, T& t )
{
  const T kEpsilon = 1e-10;
  T d[3];
  SubVec3(e, o, d);
  T e1[3], e2[3], pvec[3];
  SubVec3(v1, v0, e1);
  SubVec3(v2, v0, e2);
  Cross3(d, e2, pvec);
  const T det = Dot3(e1, pvec);
  if (fabs(det) <= kEpsilon)
    return false;

  T invdet = T(1) / det;
  T tvec[3];
  SubVec3(o, v0, tvec);
  u = Dot3(tvec, pvec) * invdet;
  if (u < -threshold || 1 + threshold < u)
    return false;

  T qvec[3];
  Cross3(tvec, e1, qvec);
  v = Dot3(d, qvec) * invdet;
  if (-threshold > v || 1 + threshold < u + v )
    return false;
  t = Dot3(e2, qvec) * invdet;
  if (t < -threshold || t >= 1 + threshold)
    return false;

  return true;
}

template <class T>
inline T Determinant4(const T* p0, const T* p1, const T* p2, const T* p3)
{
  T det = p1[2] * p2[1] * p3[0] - p0[2] * p2[1] * p3[0] -
          p1[1] * p2[2] * p3[0] + p0[1] * p2[2] * p3[0] +

          p0[2] * p1[1] * p3[0] - p0[1] * p1[2] * p3[0] -
          p1[2] * p2[0] * p3[1] + p0[2] * p2[0] * p3[1] +

          p1[0] * p2[2] * p3[1] - p0[0] * p2[2] * p3[1] -
          p0[2] * p1[0] * p3[1] + p0[0] * p1[2] * p3[1] +

          p1[1] * p2[0] * p3[2] - p0[1] * p2[0] * p3[2] -
          p1[0] * p2[1] * p3[2] + p0[0] * p2[1] * p3[2] +

          p0[1] * p1[0] * p3[2] - p0[0] * p1[1] * p3[2] -
          p0[2] * p1[1] * p2[0] + p0[1] * p1[2] * p2[0] +

          p0[2] * p1[0] * p2[1] - p0[0] * p1[2] * p2[1] -
          p0[1] * p1[0] * p2[2] + p0[0] * p1[1] * p2[2];
  return det;
}

/**
 * @param p_i: the tetrahedron vertex position
 * @param p: the point barycentric coordinate of which this function computes
 * @param barycentricCoord: result of barycentric coordinate computation
 * @return true if p is inside the tet, false o/w
 */
template <class T>
inline bool IsInsideTet(const T* p0, const T* p1, const T* p2, const T* p3,
                        const T* p, T threshold, T* barycentricCoord)
{
  const T kEpsilon = T(1e-12);
  const T det = Determinant4(p0, p1, p2, p3);
  //  if (det < EPSILON && det > -EPSILON) {
  if (det < kEpsilon && det > -kEpsilon) {
    std::cerr << "GetBarycentricCoordinate() => Tet volume (" << det << ") is almost flat!!" << std::endl;
    return false;
  }
  const T det_inv = T(1.0) / det;
  barycentricCoord[0] = Determinant4(p, p1, p2, p3) * det_inv;
  barycentricCoord[1] = Determinant4(p0, p, p2, p3) * det_inv;
  barycentricCoord[2] = Determinant4(p0, p1, p, p3) * det_inv;
  barycentricCoord[3] = Determinant4(p0, p1, p2, p) * det_inv;

  if (barycentricCoord[0] < threshold
      || barycentricCoord[1] < threshold
      || barycentricCoord[2] < threshold
      || barycentricCoord[3] < threshold) {
    return false;
  } else {
    return true;
  }
}

}
#endif // COLLISION_UTIL_H
