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

#ifndef	QUATERNION_H
#define	QUATERNION_H

#include <cmath>
#include <cstdio>
#include "vector.h"

class Quaternion {

public:

  // initializes itself to zero
  Quaternion();
  // copies another quaternion
  Quaternion(Quaternion* q);
  // this constructor expects a normalized axis
  Quaternion(double axis[3], double angle);
  // this constructor expects an axis whose length is the angle of
  // rotation in radians
  Quaternion(double axis[3]);
  Quaternion(double a, double b, double c, double d);

  inline double* data() { return e;}
  inline double getS() const { return e[3];}
  inline void getVector(double *a) const { 
    a[0] = e[0]; a[1] = e[1]; a[2] = e[2];}

  void Slerp(Quaternion *targetQuat, double t, Quaternion *interpolated);

  // axes be normalized
  void setAxisAngle(double axis[3], double angle);

  // set quaternion from axis/magnitude notation
  void setAxisAngle(double axis[3]);

  void identity();

  void rotateAxis(double from[3], double to[3]);
  void getAxisMagnitude(double axis[3]) const;
  void getAxisAngle(double axis[3], double *angle) const;

  int toAzimElevTwist(double *azim, double *elev, double *twist);
  double m_Twist;
  double m_Azimuth;
  double m_Elevation;
  void set(double azim, double elev, double twist);

  inline double	get(int	i) { return (e[i]); };
  inline double& operator[](int i) { return e[i];}
  inline double operator[](int i) const { return e[i];}



  double *toMatrix(double m[4][4]) ;
  double *toMatrix(double m[3][3]) ;
  float  *toMatrix(float m[4][4]);
  void fromMatrix(double m[3][3]) ;
  void fromMatrix(double m[4][4]) ;
  void fromMatrix(float m[4][4]) ;
  void fromVector(double v[4]);
  inline void toVector(double v[4]) {v[0] = e[0]; v[1] = e[1]; v[2] = e[2]; v[3] = e[3];};
  inline void set(double x, double y, double z,	double w)
	{ e[0] = x; e[1] = y; e[2] = z;	e[3] = w; calcAzimElevTwist();};
  void set(double rad, double axis[3]) ;
  void set(int index, double value);
  void set(int index, double value, bool batch);
  void print(FILE *fp) { fprintf(fp,"%f	%f %f %f\n",e[0],e[1],e[2],e[3]) ; } ;
  void multiply(Quaternion *q1,	Quaternion *q2)	;
  void multiply(Quaternion *q) ;
  void copy(Quaternion *q) { e[0] = q->e[0]; e[1] = q->e[1] ; e[2] = q->e[2] ;
  e[3] = q->e[3] ; calcAzimElevTwist(); } ;

  void operator*=(const double &a);
  double length() const;
  double lengthSq() const;

  Quaternion& normalize();
  void invert();

  Quaternion operator+(const Quaternion &q) const;
  Quaternion operator*(const Quaternion &q) const;
  Quaternion operator-(const Quaternion &q) const;

  void operator+=(const Quaternion &q);

  // rotates the given point without resorting to 3x3 matrix
  void rotatePoint(double point[3]) const;

  // this takes the derivative of the quaternion with respect to the
  // 3x1 axis angle defining it. Taking the derivative of a 1x4 quat w.r.t.
  // a 3x1 vector produces a 3x4 tensor array of 3 quaternions
  // returns a 1 if the derivative is singular, zero otherwise
  int getDeriv(Quaternion derivTensor[3]) const;

  // this takes the derivative of the quaternion w.r.t. time. 
  void getVelocity(double angularVelocity[3]);

  // this takes the 2nd derivative of the quaternion w.r.t. time. 
  void getAcceleration(double angularVelocity[3], 
		       double angularAcceleration[3]);

  // this takes the derivative of the quaternion w.r.t. time, and then
  // takes the derivative of the quaternion w.r.t. omega. Taking the 
  // time, then vector derivative of a quaternion takes less computation
  // using this function. Interestingly enough, omega is not an input
  // to this function

  void getTimeVecDeriv(Quaternion derivTensor[3]) const;
  
  // returns the square roots for this quaternion. There are two roots
  void getSquareRoot(Quaternion roots[2]) const;

  void printInterp();

private:
  double e[4] ;
  void calcAzimElevTwist();

} ;


#endif
