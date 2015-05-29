/**************************************************************************
	D.A.N.C.E.
	Dynamic Animation aNd Control Environment
	----------------------------------------------
	ORIGINAL AUTHORS: 
		Victor Ng (victorng@dgp.toronto.edu)
		Petros Faloutsos (pfal@cs.ucla.edu)
	CONTRIBUTORS:
		Ari Shapiro (ashapiro@cs.ucla.edu)
		Yong Cao (abingcao@cs.ucla.edu)
-----------------------------------------------
	
 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

 This software is distributed for noncommercial use in the hope that it will 
 be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
 to anyone for the consequences	of using it or for whether it serves any 
 particular purpose or works at all. No warranty is made about the software 
 or its performance.
***************************************************************************/

#ifndef	myVECTOR_H
#define	myVECTOR_H
typedef	double Vector[3] ;

 void printVector(const Vector v) ;

 void VecAdd(Vector c, const Vector a, const Vector b) ;

 void setVector(Vector v, double a, double b, double c)	;

 void VecSubtract(Vector c, const Vector a, const Vector b)	;

 int VecEq(const Vector a, const Vector b)	;

 void zeroVector(Vector	a) ;

 void VecCopy(Vector c, const Vector a) ;

 void VecSwap(Vector a, Vector b);

 double	VecDotProd(const Vector a, const Vector b) ;

 void VecCrossProd(Vector c, const Vector a, const Vector b) ;

 void VecInter(Vector c, const Vector a, const Vector b, double t) ;

 void VecNumMul(Vector c, const Vector a, double n) ;

 void VecScale(Vector c, double	n) ;

 double	VecLength(const Vector v) ;

 double VecNorm(const Vector v) ;

 double VecNormSquared(const Vector v) ;

 void VecNormalize(Vector v) ;

 double	*aVecNormalize(double *v, int n) ;

 double	aVecLength(const double *v, int n) ;

#endif























