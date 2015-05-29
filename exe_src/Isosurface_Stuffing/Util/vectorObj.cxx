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

#include "vectorObj.h"
#include <cfloat>

#define _VLERP(a, b, c) ((1.0-c)*a+c*b)

void VectorObj::lerp(const VectorObj &from, const VectorObj &to, double t) {
  
  
  vec[0] = _VLERP(from[0], to[0], t);
  vec[1] = _VLERP(from[1], to[1], t);
  vec[2] = _VLERP(from[2], to[2], t);

}





Quaternion VectorObj::operator*(Quaternion &q) {
  
  Quaternion result, tmp;
  
  tmp[0] = vec[0];
  tmp[1] = vec[1];
  tmp[2] = vec[2];
  tmp[3] = 0;
  
  result.multiply(&tmp, &q);
  
  return result;
}

VectorObj VectorObj::cross(const VectorObj &a) const {
  
  VectorObj result;

  result[0] = vec[1]*a[2] - vec[2]*a[1];
  result[1] = -vec[0]*a[2] + vec[2]*a[0];
  result[2] = vec[0]*a[1] - vec[1]*a[0];

  return result;
}

bool VectorObj::isNan() {
  
#ifdef WIN32
  return (!_finite(vec[0]) || !_finite(vec[1]) || !_finite(vec[2]));
#endif
#ifdef __linux
  return (isnan(vec[0]) || isnan(vec[1]) || isnan(vec[2]));
#endif
}

bool VectorObj::isHuge(double log10Val) {
  
  double val;

  if (isNan()) {
    return true;
  }

  val = pow((double) 10, (double)log10Val);
  
  if (fabs(vec[0])>val || fabs(vec[1])>val || fabs(vec[2])>val) {
    return true;
  }
  
  return false;
}
