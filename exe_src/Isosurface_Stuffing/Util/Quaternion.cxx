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

#include <cstdio>
#include "Quaternion.h"
#include "vectorObj.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Quaternion::Quaternion() {
  
  static double ax[] = {
    1.0, 0.0, 0.0
  };

  setAxisAngle(ax, 0);
}

Quaternion::Quaternion(Quaternion* q)
{
	double* data = q->data();
	for (int x = 0; x < 4; x++)
		e[x] = data[x];
	this->calcAzimElevTwist();
}

Quaternion::Quaternion(double a, double b, double c, double d) {
  
  e[0] = a;
  e[1] = b;
  e[2] = c;
  e[3] = d;
  this->calcAzimElevTwist();

}

Quaternion::Quaternion(double axis[3], double angle) {
  
  setAxisAngle(axis, angle);
  this->calcAzimElevTwist();
}

Quaternion::Quaternion(double axis[3]) {
  
  VectorObj vec = axis;
  double angle = vec.length();
  
  if (angle==0) {
    vec.assign(1, 0, 0);
  } else {
    vec.normalize();
  }
  
  setAxisAngle(vec.data(), angle);
  this->calcAzimElevTwist();

}

double *Quaternion::toMatrix(double m[3][3])
{

  double mat[4][4]; toMatrix(mat);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      m[i][j] = mat[i][j];
  return &m[0][0];
}


// creates a rotation matrix 4x4 from a	qauternion
// To matrix follows the column	vector convention
// v' =	m*v
double * Quaternion::toMatrix(double m[4][4])
{

	double w,x,y,z,s,xs,ys,zs,wx,wy,wz,xx,yy,zz,xy,yz,xz;

	x = e[0];
	y = e[1];
	z = e[2];
	w = e[3];

	s = 2.0/(x*x + y*y + z*z + w*w);
	xs = x*s; ys = y*s; zs = z*s;

	wx = w*xs; wy =	w*ys; wz = w*zs;
	xx = x*xs; xy =	x*ys; xz = x*zs;
	yy = y*ys; yz =	y*zs; zz = z*zs;

	m[0][0]	= 1.0 -	(yy+zz);
	m[0][1]	= xy+wz;
	m[0][2]	= xz-wy;
	m[0][3]	= 0.0;

	m[1][0]	= xy - wz;
	m[1][1]	= 1.0 -	(xx+zz);
	m[1][2]	= yz + wx;
	m[1][3]	= 0.0;

	m[2][0]	= xz+wy;
	m[2][1]	= yz-wx;
	m[2][2]	= 1.0 -	(xx+yy);
	m[2][3]	= 0.0;

	m[3][0]	= 0.0;
	m[3][1]	= 0.0;
	m[3][2]	= 0.0;
	m[3][3]	= 1.0;

    return &m[0][0];
}

// creates a rotation matrix 4x4 from a	qauternion
// To matrix follows the column	vector convention
// v' =	m*v
float * Quaternion::toMatrix(float m[4][4])
{

	double mat[4][4];
	toMatrix(mat);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			m[i][j] = (float)mat[i][j];

	return &m[0][0];
}

void Quaternion::fromMatrix(double m[3][3])
{
	double mat[4][4];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			mat[i][j] = m[i][j];
		}
	mat[0][3] = 0.0;
	mat[1][3] = 0.0;
	mat[2][3] = 0.0;

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;
	fromMatrix(mat);
}

void Quaternion::fromMatrix(float m[4][4])
{
	double mat[4][4];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			mat[i][j] = m[i][j];
		}
	mat[0][3] = 0.0;
	mat[1][3] = 0.0;
	mat[2][3] = 0.0;

	mat[3][0] = 0.0;
	mat[3][1] = 0.0;
	mat[3][2] = 0.0;
	mat[3][3] = 1.0;
	fromMatrix(mat);
}


// Matrix is expected to be used v'=m*v	where v	is
// of a	column vector
void Quaternion::fromMatrix(double m[4][4])
{
    double tr,s;
    int	i,j,k;

    tr = m[0][0]+m[1][1]+m[2][2];
    if (tr + 1.0 > 1e-10) {
	s = sqrt(tr+1.0);
	e[3] = s*0.5;
	s = 0.5/s;

	e[0] = (m[1][2]-m[2][1])*s;
	e[1] = (m[2][0]-m[0][2])*s;
	e[2] = (m[0][1]-m[1][0])*s;
    }
    else {
	i = 0;
	if (m[1][1] > m[0][0]) i = 1;
	if (m[2][2] > m[i][i]) i = 2;
	j = (i+1) % 3; k = (j+1) % 3;

	s = sqrt((m[i][i]-(m[j][j]+m[k][k])) + 1.0);
	e[i] = s*0.5;
	s = 0.5/s;
	e[j] = (m[i][j]	+ m[j][i])*s;
	e[k] = (m[i][k]	+ m[k][i])*s;
	e[3] = (m[j][k]	- m[k][j])*s;
    }

	calcAzimElevTwist();
    return ;
}

void Quaternion::fromVector(double v[4]) 
{
	e[0] = v[0]; e[1] = v[1]; e[2] = v[2]; e[3] = v[3];
	calcAzimElevTwist();
}


// Does not set Azim, Elev, Twist (called from Set(Azim,Elev,Twist)
void Quaternion::set(double rad, double axis[3])
{

    double angle2 = rad*0.5 ;
    double sin2	= sin(angle2) ;

    // Normalize axis.
    double mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    axis[0] /= mag;
    axis[1] /= mag;
    axis[2] /= mag;

    e[0] = sin2*axis[0]	;
    e[1] = sin2*axis[1]	;
    e[2] = sin2*axis[2]	;
    e[3] = cos(angle2) ;

}

// multiplies q1*q2 and	stores the result in the calling quaternion.
void Quaternion::multiply(Quaternion *q1, Quaternion *q2)
{
    Vector t1, t2 ;

    VecCrossProd(t2, q1->e, q2->e) ;
    VecNumMul(t1, q2->e, q1->e[3]) ;
    VecAdd(t2, t2, t1) ;
    VecNumMul(t1, q1->e, q2->e[3]) ;
    VecAdd(e, t1, t2) ;
    e[3] = q1->e[3] * q2->e[3] - VecDotProd(q1->e, q2->e) ;
    calcAzimElevTwist();
}

void Quaternion::multiply(Quaternion *q)
{
    Quaternion q1 ;
    q1.copy(this) ;
    multiply(&q1, q) ;
}

void Quaternion::set(double azim, double elev, double twist)
{

  // Convert angles to radians
  m_Azimuth = M_PI/180.0*azim;
  // We only allow positive elevation values. If a negative one is
  // entered, we convert it to a positive value.
  if (elev < 0.0) elev *= -1.0;
  m_Elevation = M_PI/180.0*elev;
  m_Twist = M_PI/180.0*twist;
  
  // Adjust to zero if very small.
  if (fabs(m_Azimuth) < 0.000001) m_Azimuth = 0.0;
  if (fabs(m_Elevation) < 0.000001) m_Elevation = 0.0;
  if (fabs(m_Twist) < 0.000001) m_Twist = 0.0;
  
  // Convert to spherical coordinates with radius 1
  Vector axis;
  axis[0] = sin(m_Elevation)*cos(m_Azimuth);
  axis[1] = cos(m_Elevation);
  axis[2] = sin(m_Elevation)*sin(m_Azimuth);
  
  set(m_Twist,axis);
}

int Quaternion::toAzimElevTwist(double * azim, double * elev, double * twist)
{
  double conv = 180.0 / M_PI ;
  *azim = m_Azimuth*conv ;
  *elev = m_Elevation*conv ;
  *twist = m_Twist*conv ;
  return 1;
}

// given an axis, and an angle to rotate through, set the quaterion

void Quaternion::rotateAxis(double f[3], double t[3]) {

  VectorObj from, to;
  VectorObj rotAxis;
  double amt;

  from = f;
  to = t;
  from.normalize() ;
  to.normalize() ;
  rotAxis = from.cross(to);
  amt = acos(from.dot(to));
  
  if (fabs(amt)<0.0001 ||
      fabs(amt-M_PI*2)<0.0001) {
    rotAxis.assign(1, 0, 0);
  }

  setAxisAngle(rotAxis.data(), amt);

}

void Quaternion::identity() {
  
  e[0] = 0;
  e[1] = 0;
  e[2] = 0;
  e[3] = 1;
  
}

void Quaternion::setAxisAngle(double a[3]) {
  
  VectorObj axis;
  double theta;

  axis = a;
  
  theta = axis.length();
  if (theta>0) {
    axis /= theta;
  }

  setAxisAngle(axis.data(), theta);

}

void Quaternion::setAxisAngle(double axis[3], double angle) {
  
	Vector v ;
 VecCopy(v, axis) ;
 VecNormalize(v) ;
  double sin2 = sin(angle/2.0);

  e[0] = v[0]*sin2;
  e[1] = v[1]*sin2;
  e[2] = v[2]*sin2;
  e[3] = cos(angle/2.0);

  
}

void Quaternion::getAxisMagnitude(double axis[3]) const {
  
  double angle;
  
  getAxisAngle(axis, &angle);
  
  axis[0] *= angle;
  axis[1] *= angle;
  axis[2] *= angle;

}

// Returns the axis of rotation and angle in the quaternion.
void Quaternion::getAxisAngle(double axis [3], double *angle) const {
  double half_angle = acos(e[3]);
  *angle = half_angle*2.0;
  
  double sin2	= sin(half_angle);
  if (fabs(sin2) < 0.0001) { // No twist rotation
    axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0;
    return;
  }
  
  // Axis of rotation.
  axis[0] = e[0]/sin2;
  axis[1] = e[1]/sin2;
  axis[2] = e[2]/sin2;
}

// calcAzimElevTwist:
// Converts stored quaternion into azimuth, elevation and twist angles.
// This should be called in cases where the m_Twist, m_Azimuth and
// m_Elevation have not be initialized or provided. DO NOT call this if
// we are setting the quaternion using set(azim,elev,twist).
//
void Quaternion::calcAzimElevTwist() {

  // We use m_Twist to distinguish between positive and negative twist angles.
  // An acos operation would always return a positive angle which would lose the
  // negative angle information and incorrectly produce a positive sin instead of
  // a negative one.
  double half_angle = acos(e[3]);
  m_Twist = half_angle*2.0;
  
  double sin2	= sin(half_angle);
  if (fabs(sin2) < 0.0001) { // No twist rotation
    m_Elevation = 0.0; m_Azimuth = 0.0;
    // Note, this may cause some discontinuities in the rotation
    // axis, as a quaternion with a zero twist angle can have an infinite 
    // number of rotation axes.
    return;
  }
  
  // Axis of rotation.
  double rotaxis[3];
  rotaxis[0] = e[0]/sin2;
  rotaxis[1] = e[1]/sin2;
  rotaxis[2] = e[2]/sin2;
  // Renormalize to account for truncation errors.
  VecNormalize(rotaxis);
  
  
  // Inverse operation to find spherical coordinates.
  double projXZ = sqrt(rotaxis[0]*rotaxis[0] + rotaxis[2]*rotaxis[2]);
  
  // m_Elevation has 0 at the positive y axis and Pi at the negative y axis
  // arctan gives Pi/2 at the positive y axis and -Pi/2 at the negative y axis
  // Therefore, we perform a linear mapping to go to our correct elevation space
  //
  m_Elevation = -1.0*atan2(rotaxis[1],projXZ)+M_PI*0.5;
  m_Azimuth = atan2(rotaxis[2],rotaxis[0]);
}


double acosSafe( double x ) {
	if( x < -1 ) {
		x = -1;
	} else if( x > 1 ) {
		x = 1;
	}

	return( acos(x) );
}

void Quaternion::Slerp(Quaternion *targetQuat, double t, 
		       Quaternion *interpolated) {

double dotProd = 
		e[0] * targetQuat->e[0] + e[1] * targetQuat->e[1] +
		e[2] * targetQuat->e[2] + e[3] * targetQuat->e[3];

	double theta;

	if( dotProd < 0 ) {
		theta = acosSafe(-dotProd);
	} else {
		theta = acosSafe(dotProd);
	}

	if( theta < 1e-10) {
		memcpy( interpolated->e, this->e, 4*sizeof(double) );
		return;
	}

	double coeff1 = sin((1.0-t)*theta) / sin( theta );
	double coeff2 = sin(t*theta) / sin( theta );

	if( dotProd < 0 ) {
		for( int i = 0; i < 4; i++ ) {
			interpolated->e[i] = -coeff1*this->e[i] + coeff2*targetQuat->e[i];
		}
	} else {
		for( int i = 0; i < 4; i++ ) {
			interpolated->e[i] = coeff1*this->e[i] + coeff2*targetQuat->e[i];
		}
	}
}

void Quaternion::operator*=(const double &a) {
  
  e[0] *= a;
  e[1] *= a;
  e[2] *= a;
  e[3] *= a;
}

double Quaternion::length() const {
  
  return sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]+e[3]*e[3]);

}

double Quaternion::lengthSq() const {
  
  return e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3];
}

Quaternion& Quaternion::normalize() {
  
  double l;
  
  l = length();

  if (l!=0) {
    *this *= 1.0/l;
  } else {
    setAxisAngle(VectorObj(1, 0, 0).data(), 0);
  }
  
  return *this;
}

void Quaternion::invert() {
  
  double l2;
  
  l2 = e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3];
  
  e[0] /= -l2;
  e[1] /= -l2;
  e[2] /= -l2;
  e[3] /= l2;
  
}

Quaternion Quaternion::operator-(const Quaternion &q) const {
  Quaternion result;
  
  result[0] = e[0]-q[0];
  result[1] = e[1]-q[1];
  result[2] = e[2]-q[2];
  result[3] = e[3]-q[3];

  return result;
}

void Quaternion::operator+=(const Quaternion &q) {
  
  e[0] += q[0];
  e[1] += q[1];
  e[2] += q[2];
  e[3] += q[3];
}

Quaternion Quaternion::operator+(const Quaternion &q) const {
  
  Quaternion result;
  
  result[0] = e[0]+q[0];
  result[1] = e[1]+q[1];
  result[2] = e[2]+q[2];
  result[3] = e[3]+q[3];

  return result;
}

Quaternion Quaternion::operator*(const Quaternion &Q) const {
  
  double s1, s2;
  VectorObj v1, v2, vr;
  Quaternion result;
  
  s1 = getS();
  s2 = Q.getS();
  
  getVector(v1.data());
  Q.getVector(v2.data());
  
  result[3] = s1*s2-v1.dot(v2);
  
  vr = v2*s1 + v1*s2 + v1.cross(v2);
  
  result[0] = vr[0];
  result[1] = vr[1];
  result[2] = vr[2];

  return result;
}


void Quaternion::getVelocity(double w[3]) {
  
  Quaternion om, tmp;
  
  om[3] = 0;
  om[0] = w[0];
  om[1] = w[1];
  om[2] = w[2];

  tmp = *this;

  *this = om*tmp;
  *this *= 0.5;

}

void Quaternion::getAcceleration(double om[3], 
				 double omDot[3]) {
  
  Quaternion t0, t1, t3;
  
  t0[3] = 0;
  t0[0] = om[0];
  t0[1] = om[1];
  t0[2] = om[2];
  
  t1 = *this;
  t1.getVelocity(om);

  t0 = t0 * t1;
  
  t1[3] = 0;
  t1[0] = omDot[0];
  t1[1] = omDot[1];
  t1[2] = omDot[2];

  t3 = t1 * *this;
  
  t3 = t3 + t0;
  t3 *= 0.5;
  
  *this = t3;
}


void Quaternion::getTimeVecDeriv(Quaternion tensor[3]) const {
  
  int i;
  VectorObj v, I, tmp;
  double s;
  
  v[0] = e[0];
  v[1] = e[1];
  v[2] = e[2];
  s = e[3];

  for(i=0 ; i < 3 ; i++) {
    tensor[i][3] = -0.5*v[i];
    
    I.clear();
    I[i] = 0.5;
    
    tmp = I*s;
    tmp += I.cross(v);
    
    tensor[i][0] = v[0];
    tensor[i][1] = v[1];
    tensor[i][2] = v[2];

  }
}

void Quaternion::rotatePoint(double point[3]) const {
  
  Quaternion Q, P, Qinv;
  
  Q = Qinv = *this;
  Qinv.invert();
  
  P[0] = point[0];
  P[1] = point[1];
  P[2] = point[2];
  P[3] = 0;
  
  Q = Q * P * Qinv;

  point[0] = Q[0];
  point[1] = Q[1];
  point[2] = Q[2];

}

void Quaternion::getSquareRoot(Quaternion roots[2]) const {
  
  double w, x, y, z;
  double len;
  double a[2], r;
  double b, c, d;
  int i;

  x = e[0];
  y = e[1];
  z = e[2];
  w = e[3];

  len = w*w + x*x + y*y + z*z;

  r = (w+sqrt(len))/2.0;

  a[0] = sqrt(r);
  a[1] = -sqrt(r);

  for (i=0 ; i < 2 ; i++) {
    b = x/(2*a[i]);
    c = y/(2*a[i]);
    d = z/(2*a[i]);
    
    roots[i][3] = a[i];
    roots[i][0] = b;
    roots[i][1] = c;
    roots[i][2] = d;
  }
}

void Quaternion::set(int index, double value)
{
	e[index] = value; 
	calcAzimElevTwist();
}

void Quaternion::set(int index, double value, bool batch)
{
	if (!batch)
		set(index, value);
	else 
		e[index] =	value;
}

