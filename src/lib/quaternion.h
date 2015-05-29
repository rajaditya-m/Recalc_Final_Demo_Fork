/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "quaternion" library , Copyright (C) 2007 CMU                         *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef _QUATERNION_H_
#define _QUATERNION_H_

/*
  Quaternion C++ class.
  This class implements quaternions and the commonly used 
  algebraic operations on quaternions. The class is templated:
  you can use either float or double precision.
  Supports using quaternions to represent/manipulate rotations.

  q = s + x * i + y * j + z * k

  If using quaternions to represent rotations, q must be a unit quaternion. 
    Then the following relates q to the corresponding rotation:
      s = cos(angle/2)
      (x,y,z) = sin(angle/2) * axis_of_rotation
      (axis_of_rotation is unit length)

  See also example.cpp .
*/

#include <stdio.h>
#include <math.h>


// forward declarations for the friend template
template <typename Real> class Quaternion;
template <typename Real> Quaternion<Real> operator* (Real alpha, Quaternion<Real> q2);

template <typename Real>
class Quaternion
{
public:

  inline Quaternion(); // q = 0
  inline Quaternion(Real s, Real x, Real y, Real z); // q = s + x * i + y * j + z * k
  inline Quaternion(Real s); // q = s + 0 * i + 0 * j + 0 * k;

  // Makes the unit quaternion corresponding to a rotation around axis 'unitAxis' of angle 'angle'.
  // Angle is in radians; unitAxis must be a unit 3D vector.
  inline Quaternion(Real angle, Real unitAxis[3]);

  inline void Set(Real s, Real x, Real y, Real z); // sets quaternion to the new value

  inline Real Gets() const;
  inline Real Getx() const;
  inline Real Gety() const;
  inline Real Getz() const;

  inline Quaternion operator+ (const Quaternion q2) const; // q3 = q1+q2
  inline Quaternion operator- (const Quaternion q2) const; // q3 = q1-q2
  inline Quaternion operator* (const Quaternion q2) const; // q3 = q1 * q2
  inline Quaternion operator/ (const Quaternion q2) const; // q3 = q1 / q2

  // Multiply quaternion with a scalar; e.g. q1 = alpha * q2;
  friend Quaternion<Real> operator* (Real alpha, const Quaternion<Real> q2)
  {
    return Quaternion<Real>(alpha * q2.s, alpha * q2.x, alpha * q2.y, alpha * q2.z);
  }

  inline Quaternion conj(); // q2 = q1.conj()

  inline Quaternion & operator= (const Quaternion rhs); // q2 = q1;
  inline Quaternion & operator= (Real s); // sets quaternion equal to the scalar quaternion s
  inline int operator== (const Quaternion rhs) const; // q2 == q1
  inline int operator!= (const Quaternion rhs) const; // q2 != q1

  void Normalize(); // q.Normalize() scales q such that it is unit size

  inline void MoveToRightHalfSphere(); //  if scalar part (that is, 's') is negative, this will multiply the quaternion by -1
  
  inline Real Norm2() const; // returns the squared norm of the quaternion, i.e. s*s + x*x + y*y + z*z
  inline Real Norm() const { return sqrt(Norm2()); }

  // Transforms the quaternion to the corresponding rotation matrix.
  // Quaternion is assumed to be a unit quaternion.
  // R is a 3x3 orthogonal matrix and will be returned in row-major order.
  inline void Quaternion2Matrix(Real * R) const;

  // Transforms the given matrix (assumed orthogonal) into one of the two corresponding quaternions.
  // Matrix is assumed to be in row-major order.
  // There are two quaternions corresponding to a rotation (and they have opposite signs). You can't directly control which one you get, but you can force the Real part to be non-negative by a subsequent call to MoveToRightHalfSphere() .
  // This implementation follows David Baraff's SIGGRAPH course notes:
  // http://www.cs.cmu.edu/~baraff/pbm/pbm.html
  static Quaternion Matrix2Quaternion(Real * R);
  
  // Returns the angle of rotation (in radians), and the unit rotation axis corresponding to the quaternion.
  // Assumes a unit quaternion (use Normalize() to remove any noise due to floating point errors).
  // If s >= 0, the angle will be on the interval [0,pi] .
  // If s < 0, the angle will be on the interval (pi,2pi]. To get a representation where the angle is on [0, pi], you can manually flip the sign of the returned unitAxis, and use the angle of 2pi-angle. Alternatively, you can use MoveToRightHalfSphere before calling GetRotation to ensure that you are always in the s >= 0 case.
  inline void GetRotation(Real * angle, Real unitAxis[3]) const;

  // Returns (x,y,z) = sin(theta/2) * axis, where
  //   theta is the angle of rotation, theta is on [-pi,pi), and axis is the unit axis of rotation.
  // Assumes a unit quaternion.
  // Note: this routine is a bit exotic; I expect it to be not so widely used.
  inline void GetSinExponential(Real * x, Real * y, Real * z) const;

  // Prints the quaternion to stdout.
  inline void Print() const;

protected:
  Real s,x,y,z;
};

template <typename Real>
inline Quaternion<Real>::Quaternion()
{
  s = x = y = z = 0;
}

template <typename Real>
inline Quaternion<Real>::Quaternion(Real s_)
{
  s = s_;
  x = y = z = 0;
}

template <typename Real>
inline Quaternion<Real>::Quaternion(Real angle, Real unitAxis[3])
{
  s = cos(angle/2.0);
  Real sin2 = sin(angle/2.0);
  x = sin2 * unitAxis[0];
  y = sin2 * unitAxis[1];
  z = sin2 * unitAxis[2];
}

template <typename Real>
inline Quaternion<Real>::Quaternion(Real s_, Real x_, Real y_, Real z_)
{
  s = s_;
  x = x_;
  y = y_;
  z = z_;
}

template <typename Real>
inline void Quaternion<Real>::Set(Real s_g, Real x_g, Real y_g, Real z_g) // sets quaternion to the new value
{
  s = s_g;
  x = x_g;
  y = y_g;
  z = z_g;
}

template <typename Real>
inline Real Quaternion<Real>::Gets() const { return s; }

template <typename Real>
inline Real Quaternion<Real>::Getx() const { return x; }

template <typename Real>
inline Real Quaternion<Real>::Gety() const { return y; }

template <typename Real>
inline Real Quaternion<Real>::Getz() const { return z; }

template <typename Real>
inline Quaternion<Real> & Quaternion<Real>::operator= (const Quaternion<Real> rhs)
{
  s = rhs.s;
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;

  return *this;
}

template <typename Real>
inline Quaternion<Real> & Quaternion<Real>::operator= (Real s_g)
{
  s = s_g;
  x = 0;
  y = 0;
  z = 0;

  return *this;
}

template <typename Real>
inline int Quaternion<Real>::operator== (const Quaternion<Real> rhs) const
{
  return ((s == rhs.s) && (x == rhs.x) &&
          (y == rhs.y) && (z == rhs.z));
}

template <typename Real>
inline int Quaternion<Real>::operator!= (const Quaternion<Real> rhs) const
{
  return ((s != rhs.s) || (x != rhs.x) ||
          (y != rhs.y) || (z != rhs.z));
}

template <typename Real>
inline void Quaternion<Real>::Normalize()
{
  Real invNorm;
  invNorm = (Real)1.0 / (Real)sqrt(Norm2());

  s *= invNorm;
  x *= invNorm;
  y *= invNorm;
  z *= invNorm;

}

template <typename Real>
inline Real Quaternion<Real>::Norm2() const
{ 
  return (s*s + x*x + y*y + z*z); 
}

template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator+ (const Quaternion<Real> q2) const
{
  Quaternion<Real> w(s + q2.s, x + q2.x, y + q2.y, z + q2.z);

  return w;  
}

template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator- (const Quaternion<Real> q2) const
{
  Quaternion<Real> w(s - q2.s, x - q2.x, y - q2.y, z - q2.z);
  return w;  
}

template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator* (const Quaternion<Real> q2) const
{
  Quaternion<Real> w(
        s * q2.s - x * q2.x - y    * q2.y - z * q2.z,
        s * q2.x + q2.s * x + y    * q2.z - q2.y * z,
        s * q2.y + q2.s * y + q2.x * z    - x    * q2.z,
        s * q2.z + q2.s * z + x    * q2.y - q2.x * y);

  return w;  
}

template <typename Real>
inline Quaternion<Real> Quaternion<Real>::operator/ (const Quaternion<Real> q2) const
{
  // compute invQ2 = q2^{-1}
  Quaternion<Real> invQ2;
  Real invNorm2 = 1.0 / q2.Norm2();
  invQ2.s = q2.s * invNorm2;
  invQ2.x = -q2.x * invNorm2;
  invQ2.y = -q2.y * invNorm2;
  invQ2.z = -q2.z * invNorm2;

  // result = *this * invQ2
  return (*this * invQ2); 

}

template <typename Real>
inline Quaternion<Real> Quaternion<Real>::conj()
{
  Quaternion<Real> w(s,-x,-y,-z);
  return w;
}

// Transforms the quaternion to the corresponding rotation matrix.
// Quaternion is assumed to be a unit quaternion.
// R is a 3x3 orthogonal matrix and will be returned in row-major order.
template <typename Real>
inline void Quaternion<Real>::Quaternion2Matrix(Real * R) const
{
  R[0] = 1 - 2*y*y - 2*z*z; R[1] = 2*x*y - 2*s*z;     R[2] = 2*x*z + 2*s*y;
  R[3] = 2*x*y + 2*s*z;     R[4] = 1 - 2*x*x - 2*z*z; R[5] = 2*y*z - 2*s*x;
  R[6] = 2*x*z - 2*s*y;     R[7] = 2*y*z + 2*s*x;     R[8] = 1 - 2*x*x - 2*y*y;
}

// Returns (x,y,z) = sin(theta/2) * axis, where
//   theta is the angle of rotation, theta\in\{-pi,pi\}, and
//   axis is the unit axis of rotation.
template <typename Real>
inline void Quaternion<Real>::GetSinExponential(Real * sex, Real * sey, Real * sez)const
{
  if (s<0)
  {
    *sex = -x;
    *sey = -y;
    *sez = -z;
  }
  else
  {
    *sex = x;
    *sey = y;
    *sez = z;
  }
}

template <typename Real>
inline void Quaternion<Real>::GetRotation(Real * angle, Real unitAxis[3]) const
{
  if ((s >= ((Real)1)) || (s <= (Real)(-1)))
  {
    // identity; this check is necessary to avoid problems with acos if s is 1 + eps
    *angle = 0;
    unitAxis[0] = 1;
    unitAxis[0] = 0;
    unitAxis[0] = 0;
    return;
  }

  *angle = 2.0 * acos(s);
  Real sin2 = x*x + y*y + z*z; //sin^2(*angle / 2.0)

  if (sin2 == 0)
  {
    // identity rotation; angle is zero, any axis is equally good
    unitAxis[0] = 1;
    unitAxis[0] = 0;
    unitAxis[0] = 0;
  }
  else
  {
    Real inv = 1.0 / sqrt(sin2); // note: *angle / 2.0 is on [0,pi], so sin(*angle / 2.0) >= 0, and therefore the sign of sqrt can be safely taken positive
    unitAxis[0] = x * inv;
    unitAxis[1] = y * inv;
    unitAxis[2] = z * inv;
  }
}

template <typename Real>
inline void Quaternion<Real>::MoveToRightHalfSphere()
{
  if (s<0)
  {
    s *= -1;
    x *= -1;
    y *= -1;
    z *= -1;
  }
}

template <typename Real>
inline void Quaternion<Real>::Print() const
{
  printf("%f + %fi + %fj + %fk\n",s,x,y,z);
}
#endif

