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

#ifndef	myMATH_H
#define	myMATH_H

#include "vector.h"
#include "CoordSystem.h"

#ifndef PI
#define PI 3.1415926535897
#endif

 double DNCRandom(double max, double min) ;

/*** Point ***/

 double	inter(double a,	double b, double t) ;


/**** Array ****/
 void D2ArrayCopy( int n,int m,double *c,double *a) ;

 void I2ArrayCopy( int n,int m,int *c,int *a) ;

 void transpArray(double to[4][4], double from[4][4]);

 double * transpArray(double *at, double	*a, int	n, int m) ;
 float *  transpArrayFloat(float *at, float *a, int n, int m) ;

 double * addArray(double *c, double *a,	double *b, int n, int m) ;

 double * subtractArray(double *c, double *a, double *b,	int n, int m) ;

 double * multNumArray( double num, double * a, int n, int m) ;

 double * MultNumArray( double *c, double num, double * a, int n, int m)	;

 double * multArray(double *c, double *a, double	*b, int	n1, int	m, int n2) ;

 double qT_M_q(double *m, double	*q, int	n) ;

 void printArray( double *a, int	n, int m) ;

 int isDiagonal(double *a,int n)	;

/*************** end array ****************/

 double Norm(double q[],	int n) ;

 /*** other ***/
 double Bernst3( double x, int i) ;

 double dBernst3( double	x, int i) ;

 double ddBernst3( double x, int	i) ;

 double * EulerIntegr(double dq[], double q[], int n, double dt)	;


/**********************************************
 Rotation using	matrices
 **********************************************/

 void rotateCSorigin_mat(CoordSystem *cs, double	rot[3][3]) ;
 void rotatePoint_mat(Vector point, double rot[3][3]) ;
 void transformPoint_mat(Vector point, double rot[4][4])	;
 void transformPoints_mat(Vector *points, int n, Vector	*newpoints, double transmat[4][4]) ;
 void rotPoint_mat4(Vector point, double transmat[4][4]) ;
 void relativeToFrame(double mat[4][4], double wmat[4][4], double rmat[4][4]);

 void XRotatePoints(Vector *points, int npoints, double degrees) ;
 void YRotatePoints(Vector *points, int npoints, double degrees) ;
 void ZRotatePoints(Vector *points, int npoints, double degrees) ;


void embedCsInCs(CoordSystem * cs_guest, CoordSystem * cs_host)	;

double *constrTransfMatrixFromCs(double	m[4][4], CoordSystem *cs) ;

double * quatToMat(double q[4],double m[4][4]) ;

void setIdentMat(double *m, int n) ;

void compRotMat4(double c[4][4], double m1[4][4], double m2[4][4]) ;


double *lineIntersection(Vector v1, Vector v2, Vector w1, Vector w2, 
			 Vector intersection, int infinite) ;
#endif


