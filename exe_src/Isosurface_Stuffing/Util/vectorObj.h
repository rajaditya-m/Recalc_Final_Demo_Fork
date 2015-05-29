/**************************************************
Copyright 2005 by Ari Shapiro and Petros Faloutsos

DANCE
Dynamic ANimation and Control Environment

 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

This software is distributed for noncommercial use in the hope that it will 
be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
to anyone for the consequences	of using it or for whether it serves any 
particular purpose or works at all. No warranty is made about the software 
or its performance. Commercial use is prohibited. 

Any plugin code written for DANCE belongs to the developer of that plugin,
who is free to license that code in any manner desired.

Content and code development by third parties (such as FLTK, Python, 
ImageMagick, ODE) may be governed by different licenses.
You may modify and distribute this software as long as you give credit 
to the original authors by including the following text in every file 
that is distributed:

/*********************************************************
	Copyright 2005 by Ari Shapiro and Petros Faloutsos

	DANCE
	Dynamic ANimation and Control Environment
	-----------------------------------------
	AUTHOR:
		Ari Shapiro (ashapiro@cs.ucla.edu)
	ORIGINAL AUTHORS: 
		Victor Ng (victorng@dgp.toronto.edu)
		Petros Faloutsos (pfal@cs.ucla.edu)
	CONTRIBUTORS:
		Yong Cao (abingcao@cs.ucla.edu)
		Paco Abad (fjabad@dsic.upv.es)
**********************************************************/

#ifndef _VECTOROBJ_
#define _VECTOROBJ_

#include "mathdefs.h"
#include "Quaternion.h"


class VectorObj {
 protected:
  double vec[3];
  
 public:
  
  VectorObj(double a, double b, double c) {
    vec[0] = a; vec[1] = b; vec[2] = c;}
  VectorObj(double *a) { 
    vec[0] = a[0]; vec[1] = a[1]; vec[2] = a[2];}

  VectorObj() {
    vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;}
  ~VectorObj() {}
  
  inline void clear() { vec[0] = vec[1] = vec[2] = 0;}

  inline double *data() { return vec;}
  inline const double* data() const { return vec;}

  inline double &x() { return vec[0];}
  inline double &y() { return vec[1];}
  inline double &z() { return vec[2];}
  
  inline double x() const { return vec[0];}
  inline double y() const { return vec[1];}
  inline double z() const { return vec[2];}
  
  inline void assign(double a, double b, double c) {
    vec[0] = a; vec[1] = b; vec[2] = c;
  }
  inline double &operator[](int i) { return vec[i];}
  inline const double &operator[](int i) const { return vec[i];}

  inline void toState(double *state);
  inline void fromState(double *state);

  inline VectorObj operator*(double a) const;
  inline VectorObj operator/(double a) const;

  inline void operator+=(double a[3]);
  inline void operator+=(const VectorObj &a);
  inline void operator-=(double a[3]);
  inline void operator-=(const VectorObj &a);
  inline void operator+=(const double a);

  inline void operator/=(double a);
  inline void operator*=(double a);

  inline bool operator==(const VectorObj &a) const;
  inline bool operator!=(const VectorObj &a) const;

  inline VectorObj& operator=(const VectorObj& a);
  inline VectorObj& operator=(double a);
  inline VectorObj& operator=(double *a);
  inline VectorObj& operator=(const double *a);

  inline void writeTo(double *a);
  inline void addTo(double *a);

  // this returns a vector which is the product of the x's, y's, and z's
  inline VectorObj mult(const VectorObj &a);

  // returns a vector which is the division of x's, y's, and z's
  inline VectorObj div(const VectorObj &a);

  inline double dot(const VectorObj &a) const;
  inline double length() const;
  
  inline VectorObj operator+(const VectorObj& a) const;
  inline VectorObj operator-(const VectorObj& a) const;
  
  inline VectorObj &normalize();

  inline friend VectorObj operator-(const VectorObj &a) {
    VectorObj result;
    result[0] = -a[0];
    result[1] = -a[1];
    result[2] = -a[2];
    
    return result;
  }

  inline friend VectorObj operator*(double a, const VectorObj &v) {
    return v*a;
  }

  inline friend double norm(const VectorObj &v) {
	  return VecNorm(v.vec);
  }

  inline friend double normSquared(const VectorObj &v) {
	  return VecNormSquared(v.vec);
  }

  void lerp(const VectorObj &from, const VectorObj &to, double t);


 
  Quaternion operator*(Quaternion &q);

  VectorObj cross(const VectorObj &vec) const;

  void printInterp();

  bool isNan();
  bool isHuge(double log10Val);

};

inline void VectorObj::toState(double *state) {
  state[0] = vec[0];
  state[1] = vec[1];
  state[2] = vec[2];
}

inline void VectorObj::fromState(double *state) {
  vec[0] = state[0];
  vec[1] = state[1];
  vec[2] = state[2];
}

inline VectorObj VectorObj::operator*(double a) const {
  VectorObj result(*this);
  
  result[0] *= a;
  result[1] *= a;
  result[2] *= a;
  
  return result;
}

inline VectorObj VectorObj::operator/(double a) const {
  VectorObj result(*this);
  
  result[0] /= a;
  result[1] /= a;
  result[2] /= a;
  
  return result;
}

inline void VectorObj::operator+=(const double a) {
  vec[0] += a;
  vec[1] += a;
  vec[2] += a;
}

inline void VectorObj::operator+=(double a[3]) {
  vec[0] += a[0];
  vec[1] += a[1];
  vec[2] += a[2];
}

inline void VectorObj::operator+=(const VectorObj &a) {
  vec[0] += a[0];
  vec[1] += a[1];
  vec[2] += a[2];
}

inline void VectorObj::operator-=(double a[3]) {
  vec[0] -= a[0];
  vec[1] -= a[1];
  vec[2] -= a[2];
}

inline void VectorObj::operator-=(const VectorObj &a) {
  vec[0] -= a[0];
  vec[1] -= a[1];
  vec[2] -= a[2];
}

inline void VectorObj::operator/=(double a) {
  vec[0] /= a;
  vec[1] /= a;
  vec[2] /= a;
}

inline void VectorObj::operator*=(double a) {
  vec[0] *= a;
  vec[1] *= a;
  vec[2] *= a;
}

inline bool VectorObj::operator==(const VectorObj &a) const {
  
  return vec[0]==a[0] && vec[1]==a[1] && vec[2]==a[2];
}

inline bool VectorObj::operator!=(const VectorObj &a) const {
  return vec[0]!=a[0] || vec[1]!=a[1] || vec[2]!=a[2];
}

inline VectorObj& VectorObj::operator=(const VectorObj& a) {
  vec[0] = a[0];
  vec[1] = a[1];
  vec[2] = a[2];
  
  return *this;
}

inline VectorObj& VectorObj::operator=(const double *a) {
  
  vec[0] = a[0];
  vec[1] = a[1];
  vec[2] = a[2];

  return *this;
}

inline VectorObj& VectorObj::operator=(double a) {
  vec[0] = vec[1] = vec[2] = a;
  return *this;
}

inline VectorObj& VectorObj::operator=(double *a) {
  vec[0] = a[0];
  vec[1] = a[1];
  vec[2] = a[2];
  
  return *this;
}

inline VectorObj VectorObj::mult(const VectorObj &a) {
  VectorObj result;
  
  result[0] = vec[0]*a[0];
  result[1] = vec[1]*a[1];
  result[2] = vec[2]*a[2];
  
  return result;
}

inline VectorObj VectorObj::div(const VectorObj &a) {
  VectorObj result;
  
  result[0] = vec[0]/a[0];
  result[1] = vec[1]/a[1];
  result[2] = vec[2]/a[2];
  
  return result;
}

inline double VectorObj::dot(const VectorObj &a) const {
  return vec[0]*a[0] + vec[1]*a[1] + vec[2]*a[2];
}

inline double VectorObj::length() const {
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

inline VectorObj VectorObj::operator+(const VectorObj& a) const {
  VectorObj result;
  
  result[0] = vec[0]+a[0];
  result[1] = vec[1]+a[1];
  result[2] = vec[2]+a[2];
  
  return result;
}

inline VectorObj VectorObj::operator-(const VectorObj& a) const {
  VectorObj result;
  
  result[0] = vec[0]-a[0];
  result[1] = vec[1]-a[1];
  result[2] = vec[2]-a[2];
  
  return result;
}

inline VectorObj &VectorObj::normalize() {
  *this /= length();
  
  return *this;
}

inline void VectorObj::writeTo(double *a) {
  a[0] = vec[0];
  a[1] = vec[1];
  a[2] = vec[2];
}

inline void VectorObj::addTo(double *a) {
  a[0] += vec[0];
  a[1] += vec[1];
  a[2] += vec[2];
}

#endif
