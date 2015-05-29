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

#include <stdlib.h>
#include <math.h>
#include "myMath.h"

#define	EPSI	 (1.0e-6)
#define	MAX_ITER  100

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


int converg(double *x,double *x1, int n) ;


double VecLength(double *v, int size)
{
    double sum = 0 ;
    for( int i = 0 ; i < size ; i++)
    {
	sum += v[i]*v[i] ;
    }
    return sqrt(sum) ;
}


void
printPoints(double * geom, int size, int dim)
{
	int	i,j ;

	for( i = 0 ; i < size ;	i++)
	{
		for( j = 0 ; j < dim ; j++ )
			printf("%f ", *(geom+i*dim+j)) ;
		printf("\n") ;
	}
}

double
inter(double a,	double b, double t)
{
	return a*(1.0-t)+b*t ;
}

/**** Array Functions ***/

void D2ArrayCopy( int n,int m,double *c,double *a)
{
	int i, j ;

	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < m ; j++)
			*(c+i*m+j) = *(a+i*m+j)	;
	return ;
}

void I2ArrayCopy( int n, int m,	int *c,	int *a)
{
    int	i, j ;

    for( i = 0 ; i < n ; i++)
	for( j = 0 ; j < m ; j++)
	    *(c+i*m+j) = *(a+i*m+j) ;
	return ;
}

void transpArray(double to[4][4], double from[4][4]) {
  
  int i, j;
  
  for (i=0 ; i < 4 ; i++) {
    for (j=0 ; j < 4 ; j++) {
      to[i][j] = from[j][i];
    }
  }
}

double * transpArray(double *at, double	*a, int	n, int m)
{
    int	i, j ;

	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < m ; j++)
			*(at+j*n+i) = *(a+i*m+j) ;
	return at ;
}

float  * transpArrayFloat(float	*at, float *a, int n, int m)
{
    int	i, j ;

	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < m ; j++)
			*(at+j*n+i) = *(a+i*m+j) ;
	return at ;
}


double * addArray(double *c, double *a,	double *b, int n, int m)
{
    int	i, j ;

	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < m ; j++)
			*(c+i*m+j) = *(a+i*m+j)	+ *(b+i*m+j) ;
	return c ;
}

double * subtractArray(double *c, double *a, double *b,	int n, int m)
{
    int	i, j ;

	for( i = 0 ; i < n ; i++)
		for( j = 0 ; j < m ; j++)
			*(c+i*m+j) = *(a+i*m+j)	- *(b+i*m+j) ;
	return c ;
}


double * multNumArray( double num, double * a, int n, int m)
{
    int	i, j ;

	for( i = 0 ; i < n ; i++)
	for( j = 0 ; j < m ; j++)
			*(a+i*m+j) *= num ;
	return (double *) a ;
}


double * MultNumArray( double *c, double num, double * a, int n, int m)
{
    int	i, j ;

	for( i = 0 ; i < n ; i++)
	for( j = 0 ; j < m ; j++)
			*(c+i*m+j) = *(a+i*m+j)	* num ;
	return (double *) c ;
}

// multArray:
//
// returns C = AB, matrix product of n1 x m and m x n2 matrix
double * multArray(double *c, double *a, double	*b, int	n1, int	m, int n2)
{
    int	i, j, k	;

	for( i = 0 ; i < n1 ; i++)
		for ( j	= 0 ; j	< n2 ; j++)
		{
			*(c+i*n2+j) = 0.0 ;
			for( k = 0 ; k < m ; k++)
				*(c+i*n2+j) += *(a+i*m+k) * (*(b+k*n2+j))  ;
				/** c[i][j] = Sa[i][k]*b[k][j] **/
		}
	return (double *) c ;
}



double qT_M_q(double *m, double	*q, int	n)
{
	double rel, rel1 ;
	int i,j	;

	rel = 0.0 ;
	for( i = 0 ; i < n ; i++)
	{
		rel1 = 0.0 ;
		for( j = 0 ; j < n ; j++)
			rel1 = rel1 + *(m+i*n+j) * (*(q+j)) ;
		rel = rel+rel1*(*(q+i))	;
	}
	return rel ;
}


void printArray( double	*a, int	n, int m)
{
	int i,j	;

	for( i = 0 ; i < n ; i++)
	{
	    for( j = 0 ; j < m ; j++)
		printf("%lf ", *(a+i*m+j)) ;
	    printf("\n") ;
	}
	return ;
}


double Norm(double q[],	int n)
{
	double x = 0.0 ;
	int	i ;

	for ( i	= 0 ; i	< n ; i++)
		x += q[i]*q[i] ;
	return (double)	sqrt(x)	;
}


/*******************************************
  PROC:	ddBerstn3()
  DOES:	returns	the second derivative with respect to x
********************************************/
double ddBernst3( double x, int	i)
{
    double  t ;

    switch(i)
    {
    case 0:
	t = 6*(1-x) ;
	break ;
    case 1:
	t = 18*x-12 ;
	break ;
    case 2:
	t = 6-18*x ;
	break ;
    case 3:
	t = 6*x	;
	break ;
    }
    return t ;
}

/*******************************************
  PROC:	dBerstn3()
  DOES:	returns	the first derivative with respect to x
********************************************/
double dBernst3( double	x, int i)
{
    double  t ;

    switch(i)
    {
    case 0:
	t = -3*(1-x)*(1-x) ;
	break ;
    case 1:
	t = (1-x)*(3-9*x) ;
	break ;
    case 2:
	t = x*(6-9*x) ;
	break ;
    case 3:
	t = 3*x*x ;
	break ;
    }
    return t ;
}

/** BernsteinPolynomials degree	3 **/
double Bernst3(	double x, int i)
{
    double  t ;

    switch(i)
    {
    case 0:
	t = (1-x) ;
	t = t*t*t ;
	break ;
    case 1:
	t = 1-x	;
	t = 3*x*t*t ;
	break ;
    case 2:
	t = 3*x*x*(1-x)	;
	break ;
    case 3:
	t = x*x*x ;
	break ;
    default:
	fprintf(stderr,"Bernst3, illegal exponent!\n") ;
	break ;
    }
    return t ;
}


double * EulerIntegr(double dq[], double q[], int n, double dt)
{
	int i ;

	for( i = 0; i <	n ; i++)
		q[i] = q[i] + dt*dq[i] ;
	return &q[0] ;
}




/***************************************************
 *	      Rotation functions		   *
 *	    ----------------------		   *
 * All of them accept 2	quaternions  an	angle in   *
 * in rad and an axis. If  quaternions are given   *
 * then	they are				   *
 * used. Otherwise quaternions are constructed.	   *
 * rotatePoint2D is an exception (2D rotation.	   *
 * ----------------------------------------------- *
 * ATTENTION !!	: In order for the functions to	   *
 * work	the quaternions	should be of unit mag.	   *
 * Currently for effeciency I have left	the norma- *
 * lization fot	the calling function.		   *
 * However the functions construct unit	quaternions*
 ***************************************************/

void rotatePoint_mat(Vector point, double rot[3][3])
{
Vector temp ;

multArray(&temp[0], &rot[0][0], &point[0], 3,3,1) ;
point[0] = temp[0] ; point[1] = temp[1] ; point[2] =temp[2] ;
return ;
}

void rotatePoint2D(Vector point, double	theta)
{
	double costh = cos(theta) ;
	double sinth = sin(theta) ;
	double x,y ;

	x = point[0]*costh - point[1]*sinth ;
	y = point[0]*sinth + point[1]*costh ;
	point[0] = x ;
	point[1] = y ;
	point[2] = 0.0 ;
	return ;
}


void rotateCSorigin_mat(CoordSystem *cs, double	rot[3][3])
{
    rotatePoint_mat(cs->x, rot)	;
    rotatePoint_mat(cs->y, rot)	;
    rotatePoint_mat(cs->z, rot)	;
    return ;
}

void rotateVector_mat(Vector point, double rot[3][3])
{
    Vector temp	;

    multArray(&temp[0],	&rot[0][0], &point[0], 3,3,1) ;
    point[0] = temp[0] ; point[1] = temp[1] ; point[2] =temp[2]	;
    return ;
}

//  transformPoint_mat:
//	transform given	point with transformation matrix by premultiplying
//	the matrix which is in column-major form
//
//	point	: Vector of size 3 doubles
//	transmat: Transformation matrix	in column-major	form.
//
void transformPoint_mat(Vector point, double transmat[4][4])
{
    double transpoint[4];
    transpoint[0] = point[0];
    transpoint[1] = point[1];
    transpoint[2] = point[2];
    transpoint[3] = 1.0;

    double newpoint[4]={0.0,0.0,0.0,0.0};

    // transmat	is in column-major order and is	pre-multiplied with
    // point (a	row vector)
    for	(int i=0; i < 4; i++)
	for (int j=0; j	< 4; j++)
		newpoint[i] += transmat[j][i]*transpoint[j];

    double n1 =	1.0 / newpoint[3];
    point[0] = newpoint[0]*n1 ;
    point[1] = newpoint[1]*n1 ;
    point[2] = newpoint[2]*n1 ;
    return ;
}


//  transformPoints_mat:
//	transform given	points with transformation matrix by premultiplying
//	the matrix which is in column-major form
//
//	point	: Array	of vectors of size 3 doubles
//	transmat: Transformation matrix	in column-major	form.
//	n	: number of points
//
void transformPoints_mat(Vector	*points, int n,	Vector *newpoints,
			 double	transmat[4][4])
{
    double transpoint[4];
    double newpoint[4] ;

    for( int np	= 0 ; np < n ; np++ )
    {
	transpoint[0] =	points[np][0];
	transpoint[1] =	points[np][1];
	transpoint[2] =	points[np][2];
	transpoint[3] =	1.0;

	newpoint[0] = newpoint[1] = newpoint[2]	= newpoint[3] =	0 ;


	// transmat is in column-major order and is pre-multiplied with
	// point (a row	vector)
	for (int i=0; i	< 4; i++)
	    for	(int j=0; j < 4; j++)
		newpoint[i] += transmat[j][i]*transpoint[j];

	double n1 = 1.0	/ newpoint[3];
	newpoints[np][0] = newpoint[0]*n1 ;
	newpoints[np][1] = newpoint[1]*n1 ;
	newpoints[np][2] = newpoint[2]*n1 ;
    }
    return ;
}

// as the above	but does only the rotation
void rotVector_mat4(Vector point, double	transmat[4][4])
{
    double newpoint[3] = {0.0,0.0,0.0};

    // transmat	is in column-major order and is	pre-multiplied with
    // point (a	row vector)
    for	(int i=0; i < 3; i++)
	for (int j=0; j	< 3; j++)
		newpoint[i] += transmat[j][i]*point[j];

    point[0] = newpoint[0];
    point[1] = newpoint[1];
    point[2] = newpoint[2] ;
    return ;
}





/************************************************************
	PROC: converg()
	DOES:	checks how close to successive solutions are.
		If they	are close then returns true.
*************************************************************/

int converg(double *x,double *x1, int n)
{
	int i ;
	double s ;

	s = 0 ;
	for( i = 0 ; i < n ; i++)
		s = s +	fabs(x[i]-x1[i]) ;

	if ( s/n < EPSI)
		return 1 ;
	else return 0 ;
}


/*****************************************************
  PROC:	embedCsInCs()
  DOES:	changes	the cs_guest such that it is given with
	respect	to the cs_host.	Both systems must
	initially be with respect to the same one
*******************************************************/
void embedCsInCs(CoordSystem *cs_guest,	CoordSystem *cs_host)
{
    double x,y,z ;

    VecSubtract(cs_guest->origin, cs_guest->origin, cs_host->origin) ;

    /* Now I have the guest origin in the reference frame of the */
    /* host frame. I have to transform it with respect to the host */
    /* frame itself which may be rotated with respect to its reference	frame */
    x =	VecDotProd(cs_guest->origin,cs_host->x)	;
    y =	VecDotProd(cs_guest->origin,cs_host->y)	;
    z =	VecDotProd(cs_guest->origin,cs_host->z)	;
    cs_guest->origin[0]	= x ;
    cs_guest->origin[1]	= y ;
    cs_guest->origin[2]	= z ;

    x =	VecDotProd(cs_guest->x,cs_host->x) ;
    y =	VecDotProd(cs_guest->x,cs_host->y) ;
    z =	VecDotProd(cs_guest->x,cs_host->z) ;
    cs_guest->x[0] = x ;
    cs_guest->x[1] = y ;
    cs_guest->x[2] = z ;

    x =	VecDotProd(cs_guest->y,cs_host->x) ;
    y =	VecDotProd(cs_guest->y,cs_host->y) ;
    z =	VecDotProd(cs_guest->y,cs_host->z) ;

    cs_guest->y[0] = x ;
    cs_guest->y[1] = y ;
    cs_guest->y[2] = z ;

    x =	VecDotProd(cs_guest->z,cs_host->x) ;
    y =	VecDotProd(cs_guest->z,cs_host->y) ;
    z =	VecDotProd(cs_guest->z,cs_host->z) ;

    cs_guest->z[0] = x ;
    cs_guest->z[1] = y ;
    cs_guest->z[2] = z ;

    /* Normalize unit vectors */
    VecNormalize(cs_guest->x) ;
    VecNormalize(cs_guest->y) ;
    VecNormalize(cs_guest->z) ;

    return ;
}

/**************************************************************
  PROC:	constrTransfMatrixFromCs()
  DOES:
   constructs the  transpose (following	the OpenGL convention
   of the transformation matrix	that transforms	points from
   cs coordinate to the	reference frame	of cs
****************************************************************/
double *constrTransfMatrixFromCs(double	m[4][4], CoordSystem *cs)
{
/*    This also	works
    m[0] = cs->x[0] ;
    m[1] = cs->x[1] ;
    m[2] = cs->x[2] ;
    m[3] = 0.0 ;
    m[4] = cs->y[0] ;
    m[5] = cs->y[1] ;
    m[6] = cs->y[2] ;
    m[7] = 0.0 ;
    m[8] = cs->z[0] ;
    m[9] = cs->z[1] ;
    m[10] = cs->z[2] ;
    m[11]= 0.0 ;
    m[12] = cs->origin[0] ;
    m[13] = cs->origin[1] ;
    m[14] = cs->origin[2] ;
    m[15] = 1.0	;
*/
    /* construct the transpose following the OpenGl convention */
    /* for storing matrices */
    VecCopy(m[0],cs->x)	;
    VecCopy(m[1],cs->y)	;
    VecCopy(m[2],cs->z)	;
    m[3][3] = 1.0 ;
    m[3][0] = cs->origin[0] ;
    m[3][1] = cs->origin[1] ;
    m[3][2] = cs->origin[2] ;
    m[0][3] = 0.0 ;
    m[1][3] = 0.0 ;
    m[2][3] = 0.0 ;
    return &m[0][0] ;
}

/**************************************************
  PROC:	invSmart4()
  DOES:	inverts	a 4x4 transformation matrix that consists
	of rotation and	translation only. The 3x3 part is
	orthonormal so there is	a smart	way to do it.
	The matrix is in column	order (following the
	OpenGL transpose convention. Thus I have to multiply
	from the left side.
	The inversion is done as follows:
	a = R*T	-> a^(-1) = T^(-1)*R^(-1) etc
	To avoid coputations with zeroes I do the
	matrix multiplication myself.
***************************************************/

double *invSmart4old(double inv[4][4],double a[4][4])
{
    int	i,j ;
    double t[3]	;

    for( i = 0 ; i < 3 ; i++)
	for( j = 0 ; j < 3 ; j++)
	    inv[i][j] =	a[j][i]	;
    inv[0][3] =	0.0 ;
    inv[1][3] =	0.0 ;
    inv[2][3] =	0.0 ;
    inv[3][3] =	1.0 ;

    t[0] = a[3][0] ;
    t[1] = a[3][1] ;
    t[2] = a[3][2] ;
    inv[3][0] =	-(t[0]*inv[0][0]+t[1]*inv[1][0]+t[2]*inv[2][0])	;
    inv[3][1] =	-(t[0]*inv[0][1]+t[1]*inv[1][1]+t[2]*inv[2][1])	;
    inv[3][2] =	-(t[0]*inv[0][2]+t[1]*inv[1][2]+t[2]*inv[2][2])	;

    
    return &inv[0][0] ;
}


// PROC: invSmart4()
// DOES: inverts a matrix of the form:
//	 A = [ a00, a01, a02, 0,
//	       a10, a11, a12, 0,
//	       a20, a21, a22, 0,
//	       tx, ty, tz, 1]

double *invSmart4(double inv[4][4],double a[4][4])
{
    double t4,t6,t8,t10,t12,t14,t17,t27,
      t28, t31,	t32, t35, t36,t49,t51,t59 ;
    double tx =	a[3][0]	;
    double ty =	a[3][1]	;
    double tz =	a[3][2]	;

    double a00 = a[0][0] ;
    double a01 = a[0][1] ;
    double a02 = a[0][2] ;
    double a10 = a[1][0] ;
    double a11 = a[1][1] ;
    double a12 = a[1][2] ;
    double a20 = a[2][0] ;
    double a21 = a[2][1] ;
    double a22 = a[2][2] ;

    t4 = a00*a11;
    t6 = a00*a21;
    t8 = a10*a01;
    t10	= a10*a21;
    t12	= a20*a01;
    t14	= a20*a11;
    t17	= 1.0/(t4*a22-t6*a12-t8*a22+t10*a02+t12*a12-t14*a02);
    t27	= a10*a22;
    t28	= a20*a12;
    t31	= a00*a22;
    t32	= a20*a02;
    t35	= a00*a12;
    t36	= a10*a02;
    t49	= tx*a11;
    t51	= tx*a21;
    t59	= tx*a01;
    inv[0][0] =	(a11*a22-a21*a12)*t17;
    inv[0][1] =	-(a01*a22-a21*a02)*t17;
    inv[0][2] =	(a01*a12-a11*a02)*t17;
    inv[0][3] =	0.0;
    inv[1][0] =	-(t27-t28)*t17;
    inv[1][1] =	(t31-t32)*t17;
    inv[1][2] =	-(t35-t36)*t17;
    inv[1][3] =	0.0;
    inv[2][0] =	(t10-t14)*t17;
    inv[2][1] =	-(t6-t12)*t17;
    inv[2][2] =	(t4-t8)*t17;
    inv[2][3] =	0.0;
    inv[3][0] =	-(t10*tz-t27*ty-t14*tz+t28*ty+t49*a22-t51*a12)*t17;
    inv[3][1] =	(t6*tz-t31*ty-t12*tz+t32*ty+t59*a22-t51*a02)*t17;
    inv[3][2] =	-(t4*tz-t35*ty-t8*tz+t36*ty+t59*a12-t49*a02)*t17;
    inv[3][3] =	1.0;

    return &inv[0][0] ;
}

// Matrix in row-major order!! [ Rot T]
//                             [  0  1]  
double * quatToMat(double q[4],double m[4][4])
{
    double w,x,y,z;
    double xs,ys,zs,wx,wy,wz,xx,yy,zz,xy,yz,xz;

    w =	q[0];
    x =	q[1];
    y =	q[2];
    z =	q[3];

    xs = x + x;
    ys = y + y;
    zs = z + z;
    wx = w*xs;
    wy = w*ys;
    wz = w*zs;
    xx = x*xs;
    yy = y*ys;
    zz = z*zs;
    xy = x*ys;
    yz = y*zs;
    xz = x*zs;
    m[0][0] = 1.0 - (yy+zz);
    m[0][1] = xy-wz;
    m[0][2] = xz+wy;
    m[0][3] = 0.0;
    m[1][0] = xy + wz;
    m[1][1] = 1.0 - (xx+zz);
    m[1][2] = yz - wx;
    m[1][3] = 0.0;
    m[2][0] = xz-wy;
    m[2][1] = yz + wx;
    m[2][2] = 1.0 - (xx+yy);
    m[2][3] = 0.0;
    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = 0.0;
    m[3][3] = 1.0;

    return &m[0][0] ;
}

//------------------------------------------
void setIdentMat(double	*mat,int n)
{
    int	i, j ;

    for( i = 0 ; i < n ; i++)
	for( j = 0 ; j < n ; j++ )
	    if(i == j)
		*(mat+i*n+j) = 1.0 ;
	    else
		*(mat+i*n+j) = 0.0 ;
}

int isDiagonal(double *a,int n)
{
    int	i,j ;
    int	res = TRUE ;

    for( i = 0 ; i < n ; i++ )
	for( j = 0 ; j < n ; j++ )
	    if(	i != j)
		if(*(a+i*n+j) != (double) 0.0)
		    res	= FALSE	;
    return res ;
}

// multipies two 4x4 matrices using only the 3x3 part
// (rotational one). The translation one is left untouched!!
void compRotMat4(double	c[4][4], double	m1[4][4], double m2[4][4])
{
    int	i,j,k ;

    for( i = 0 ; i < 3 ; i++)
    {

	for( j = 0 ; j < 3 ; j++ )
	{
	    c[i][j] = 0.0 ;
	    for( k = 0 ; k < 3 ; k++)
		c[i][j]	+= m1[i][k]*m2[k][j] ;
	}
    }
}

// relativeToFrame:
//
// Returns the new coordinate frame axes of wmat relative to the coordinate
// frame rmat. 
//
// NOTE: does not take into account translation. Assumes 0,0,0 translation.
void relativeToFrame(double mat[4][4], double wmat[4][4], double rmat[4][4])
{
	// Take	transpose of reference frame.
	double tmat[4][4];
	transpArray(&tmat[0][0], &rmat[0][0], 4,4);
	multArray(&mat[0][0], &tmat[0][0], &wmat[0][0],4,4,4);
}

// Calculates the intersection of the line defined by points v1,v2 with the
// line defined by w1,w2. If infinite is 1 then the lines are considered infinite.
// Otherwise the lines are treated as line segments.
// If there is no intersection NULL is returned.
double *lineIntersection(Vector v1, Vector v2, Vector w1, Vector w2, 
			 Vector intersection, int infinite )
{
    int i ;
    double t2, x2x1, ay,py,az,pz ;

    x2x1 = v2[0] - v1[0] ;

    py = (w2[0] - w1[0])*(v2[1] - v1[1]) +
         (v2[0] - v1[0])*(w1[1] - w2[1]) ;
    ay = w1[1]*(v2[0] - v1[0]) +
         v2[1]*(v1[0] - w1[0]) +
         v1[1]*(w1[0] - v2[0]) ;
    if ( (fabs(py) < EPSI) && (fabs(ay) < EPSI) )
    {  /** means same projection on xy-plane
         try to see the xz plane
         **/
        pz = (w2[0] - w1[0])*(v2[2] - v1[2]) +
             (v2[0] - v1[0])*(w1[2] - w2[2]) ;
        az = w1[2]*(v2[0] - v1[0]) +
             v2[2]*(v1[0] - w1[0]) +
             v1[2]*(w1[0] - v2[0]) ;
        if ( (fabs(pz) < EPSI) && (fabs(az) < EPSI) )
        {
            fprintf(stderr, "interscpLines: lines coincide\n") ;
            return NULL ;
        }
        else if ( (fabs(pz) < EPSI) && (fabs(az) > EPSI) )
        {
            return NULL ;
        }
        else
        {
	    //printf("pz = %lf\n", pz) ;
            t2 = az / pz ;
            for( i = 0 ; i < 3 ; i++)
                intersection[i] = (1-t2)*w1[i] + t2*w2[i] ;
        }
    }
    else if ( (fabs(py) < EPSI) && (fabs(ay) > EPSI) )
    {
        return NULL ;
    }
    else
    {
	//printf("py = %lf\n", py) ;
        t2 = ay / py ;
        for( i = 0 ; i < 3 ; i++)
            intersection[i] = (1-t2)*w1[i] + t2*w2[i] ;
    }

    // if infinites is not 1
    // check to see if the intersection point is in the line segments
    if( infinite == 0 )
    {
	if( (v1[0] - intersection[0])*(v2[0] - intersection[0]) > 0.0 )
	    return NULL ;
	else if( (v1[1] - intersection[1])*(v2[1] - intersection[1]) > 0.0 )
	    return NULL ;
	else if( (v1[2] - intersection[2])*(v2[2] - intersection[2]) > 0.0 )
	    return NULL ;
	else if( (w1[0] - intersection[0])*(w2[0] - intersection[0]) > 0.0 )
	    return NULL ;
	else if( (w1[1] - intersection[1])*(w2[1] - intersection[1]) > 0.0 )
	    return NULL ;
	else if( (w1[2] - intersection[2])*(w2[2] - intersection[2]) > 0.0 )
	    return NULL ;
    }
   
    return (double *) &intersection[0] ;

}

void XRotatePoints(Vector *points, int npoints, double degrees)
{
    // transform degrees into rad
    double angle = degrees *0.01745329252 ;  // M_PI / 180.0 ;
   
    double c = cos(angle) ;
    double s = sin(angle) ;
   
    for( int i = 0 ; i < npoints ; i++ )
    {
	Vector point ;
	VecCopy(point, points[i]) ;
	points[i][1] = point[1]*c - point[2]*s ;
	points[i][2] = point[1]*s + point[2]*c ;
    }
	
}

void YRotatePoints(Vector *points, int npoints, double degrees)
{
    // transform degrees into rad
    double angle = degrees *0.01745329252 ; // M_PI / 180.0 ;
   
    double c = cos(angle) ;
    double s = sin(angle) ;
   
    for( int i = 0 ; i < npoints ; i++ )
    {
	Vector point ;
	VecCopy(point, points[i]) ;
	points[i][0] = point[0]*c + point[2]*s ;
	points[i][2] = -point[0]*s + point[2]*c ;
    }
	
}

void ZRotatePoints(Vector *points, int npoints, double degrees)
{
    // transform degrees into rad
    double angle = degrees *0.01745329252 ;  // M_PI / 180.0 ;
   
    double c = cos(angle) ;
    double s = sin(angle) ;
   
    for( int i = 0 ; i < npoints ; i++ )
    {
	Vector point ;
	VecCopy(point, points[i]) ;
	points[i][0] = point[0]*c - point[1]*s ;
	points[i][1] = point[0]*s + point[1]*c ;
    }
	
}



double DNCRandom(double max, double min)
{
    double x;
//    int sgn = 1 ;
	
	
#ifdef WIN32 
    static double factor = 1.0 / (double) RAND_MAX ;
    x = (double) rand()*factor ;
#else
    x = drand48() ;
    //  printf(" x = %lf ", x) ;
#endif
    // x is in [0.0, 1.0) so scale it to the proper values
	
	return x*(max - min) + min ;

}
