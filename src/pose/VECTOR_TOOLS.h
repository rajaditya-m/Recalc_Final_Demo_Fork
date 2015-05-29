#ifndef __WHMIN_VECTOR_TOOLS_H__
#define __WHMIN_VECTOR_TOOLS_H__
#include <math.h>
#include "MY_MATH.h"
#include "IO_FUNC.h"
//#include "TIMER.h"
#include "timer.h"

#ifndef FORCEINLINE
#define FORCEINLINE inline
#endif



template <class T> FORCEINLINE
T Dot(T *input0, T *input1, int number = 3) {
  double ret = 0;
  for (int i = 0; i < number; i++)	ret += input0[i] * input1[i];
  return ret;
}


template <class T> FORCEINLINE
void Cross(T *input0, T *input1, T *output) {
  output[0] = input0[1] * input1[2] - input1[1] * input0[2];
  output[1] = -input0[0] * input1[2] + input1[0] * input0[2];
  output[2] = input0[0] * input1[1] - input1[0] * input0[1];
}

FORCEINLINE
double Norm(double *input, int number = 3) {
  double ret = 0;
  for (int i = 0; i < number; i++)	if (ret < fabs(input[i]))	ret = fabs(input[i]);
  return ret;
}

template <class T> FORCEINLINE
T Magnitude_Squared(T *input) {
  return input[0] * input[0] + input[1] * input[1] + input[2] * input[2];
}

template <class T> FORCEINLINE
T Magnitude(T *input) {
  return sqrtf(input[0] * input[0] + input[1] * input[1] + input[2] * input[2]);
}



template <class T> FORCEINLINE
T Normalize(T *input) {
  double m = Magnitude(input);
  if (m < 1e-10f)	{
    printf("ERROR: vector cannot be normalized. Length of vector is %f\n", m);
    return 0;
  }
  double inv_m = 1 / m;
  input[0] *= inv_m;
  input[1] *= inv_m;
  input[2] *= inv_m;
  return m;
}

FORCEINLINE
void Matrix_Product_4(double *a, double *b, double *r) {
  //r=a*b
  double temp[16];
  memset(temp, 0, sizeof(double) * 16);
  for (int i = 0; i < 4; i++)	for (int j = 0; j < 4; j++)	for (int k = 0; k < 4; k++)
        temp[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];

  memcpy(r, temp, sizeof(double) * 16);
}

template <class T> FORCEINLINE
T Normal(T *p0, T *p1, T *p2, T *normal) {
  T e0[3], e1[3];
  for (int i = 0; i < 3; i++) {
    e0[i] = p1[i] - p0[i];
    e1[i] = p2[i] - p0[i];
  }
  Cross(e0, e1, normal);
  T normal_length = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
  if (normal_length == 0) return normal_length;
  normal_length = T(sqrt(normal_length));
  for (int i = 0; i < 3; i++)
    normal[i] /= normal_length;
  return normal_length;
}

template <class T> FORCEINLINE
T Det_Matrix(T *x) {
  return x[0] * (x[4] * x[8] - x[7] * x[5]) + x[3] * (x[7] * x[2] - x[1] * x[8]) + x[6] * (x[1] * x[5] - x[4] * x[2]);
}

template <class T> FORCEINLINE
T Matrix_Inverse(T *x, T *inv_m) {
  inv_m[0] = x[4] * x[8] - x[7] * x[5];
  inv_m[1] = x[7] * x[2] - x[1] * x[8];
  inv_m[2] = x[1] * x[5] - x[4] * x[2];

  inv_m[3] = x[5] * x[6] - x[3] * x[8];
  inv_m[4] = x[0] * x[8] - x[2] * x[6];
  inv_m[5] = x[2] * x[3] - x[0] * x[5];

  inv_m[6] = x[3] * x[7] - x[4] * x[6];
  inv_m[7] = x[1] * x[6] - x[0] * x[7];
  inv_m[8] = x[0] * x[4] - x[1] * x[3];

  T det = (x[0] * inv_m[0] + x[3] * inv_m[1] + x[6] * inv_m[2]);
  T inv_det = 1 / det;
  for (int i = 0; i < 9; i++)
    inv_m[i] *= inv_det;
  return det;
}

FORCEINLINE
void Matrix_Inverse_2(double *m, double *inv_m) {
  double inv_det = 1 / (m[0] * m[3] - m[1] * m[2]);
  inv_m[0] = m[3] * inv_det;
  inv_m[1] = -m[1] * inv_det;
  inv_m[2] = -m[2] * inv_det;
  inv_m[3] = m[0] * inv_det;
}

template <class T> FORCEINLINE
void Matrix_Transpose(T *m, T *x) {
  memcpy(x, m, sizeof(T) * 9);
  Swap(x[1], x[3]);
  Swap(x[2], x[6]);
  Swap(x[5], x[7]);
}

template <class T> FORCEINLINE
void Matrix_Transpose_3(T *m, T *x) {
  memcpy(x, m, sizeof(T) * 9);
  Swap(x[1], x[3]);
  Swap(x[2], x[6]);
  Swap(x[5], x[7]);
}

FORCEINLINE
void Matrix_Add(double *A, double a, double *B, double b, double *c) {
  for (int i = 0; i < 9; i++)
    c[i] = A[i] * a + B[i] * b;
}



FORCEINLINE
void Chelosky_Factorization(double *x, double *m) {
  m[0] = sqrt(x[0]);
  m[3] = x[3] / m[0];
  m[6] = x[6] / m[0];


  m[1] = 0;
  m[4] = sqrt(x[4] - m[3] * m[3]);
  m[7] = (x[7] - m[6] * m[3]) / m[4];

  m[2] = 0;
  m[5] = 0;
  m[8] = sqrt(x[8] - m[6] * m[6] - m[7] * m[7]);
}

template <class T>
FORCEINLINE
void Matrix_Product_2(T *M0, T *M1, T *R) {
  R[0] = M0[0] * M1[0] + M0[1] * M1[2];
  R[1] = M0[0] * M1[1] + M0[1] * M1[3];
  R[2] = M0[2] * M1[0] + M0[3] * M1[2];
  R[3] = M0[2] * M1[1] + M0[3] * M1[3];
}
template <class T> FORCEINLINE
void Matrix_Product_3(T *M0, T *M1, T *R) {
  R[0] = M0[0] * M1[0] + M0[1] * M1[3] + M0[2] * M1[6];
  R[1] = M0[0] * M1[1] + M0[1] * M1[4] + M0[2] * M1[7];
  R[2] = M0[0] * M1[2] + M0[1] * M1[5] + M0[2] * M1[8];
  R[3] = M0[3] * M1[0] + M0[4] * M1[3] + M0[5] * M1[6];
  R[4] = M0[3] * M1[1] + M0[4] * M1[4] + M0[5] * M1[7];
  R[5] = M0[3] * M1[2] + M0[4] * M1[5] + M0[5] * M1[8];
  R[6] = M0[6] * M1[0] + M0[7] * M1[3] + M0[8] * M1[6];
  R[7] = M0[6] * M1[1] + M0[7] * M1[4] + M0[8] * M1[7];
  R[8] = M0[6] * M1[2] + M0[7] * M1[5] + M0[8] * M1[8];
}

/*
XWU_FORCEINLINE
void Matrix_Product(double *m0, double *m1, double *x, int nx, int ny, int nz)
{
  for(int i=0; i<nx; i++)
  for(int j=0; j<nz; j++)
  {
    x[i*nz+j]=0;
    for(int k=0; k<ny; k++)
      x[i*nz+j]+=m0[i*ny+k]*m1[k*nz+j];
  }
}*/


FORCEINLINE
void Matrix_Assign(double *a, double *b) {
  memcpy(b, a, sizeof(double) * 9);
}

FORCEINLINE
void Matrix_Substract(double *a, double *b, double *r) {
  for (int i = 0; i < 9; i++)
    r[i] = a[i] - b[i];
}

template <class T> FORCEINLINE
T Distance(T *a, T *b) {
  return T(sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
}


FORCEINLINE double cbrt_5f(double f) {
  unsigned int* p = (unsigned int *) &f;
  *p = *p / 3 + 709921077;
  return f;
}

FORCEINLINE double cbrta_newtonf(const double a, const double x) {
  return a - (1.0f / 3.0f) * (a - x / (a * a));
}


FORCEINLINE double my_cbrt(double d) {
  double a = cbrt_5f(d);
  a = cbrta_newtonf(a, d);
  a = cbrta_newtonf(a, d);
  //a = cbrta_newtonf(a, d);
  return cbrta_newtonf(a, d);
}

template <class T>
FORCEINLINE
T Newton_Cubic_Root(T a, T b, T c, T d, T x) {
  //Find a root from start
  T  f, g;

  //start from zero
  f = ((a * x + b) * x + c) * x + d;
  g = (3 * a * x + 2 * b) * x + c;

  //printf("enter\n");

  for (int i = 0; i < 16; i++) {
    if (fabs(f) < 1e-14f)	return x;
    T step = f / g;
    if (step > 0.1f)	step = 0.1f;
    if (step < -0.1f)	step = -0.1f;
    x -= step;
    //printf("use this %.14f; %f\n", f_value, x);
    f = ((a * x + b) * x + c) * x + d;
    g = 3 * a * x * x + 2 * b * x + c;
  }
  //printf("error: %f\n", f);
  return x;
}

template <class T>
FORCEINLINE
T Newton_Cubic_Root(T a, T b, T c, T d) {
  //Find a root between 0 and 1
  T x, f, g;

  //start from zero
  x = 0;
  f = ((a * x + b) * x + c) * x + d;
  g = (3 * a * x + 2 * b) * x + c;

  if ((f > 0 && g < 0) || (f < 0 && g > 0)) {
    for (int i = 0; i < 16; i++) {
      if (fabs(f) < 1e-14f)	return x;
      T step = f / g;
      if (step > 0.1f)	step = 0.1f;
      if (step < -0.1f)	step = -0.1f;
      x -= step;
      //printf("use this %.14f; %f\n", f_value, x);
      f = ((a * x + b) * x + c) * x + d;
      g = 3 * a * x * x + 2 * b * x + c;
    }
  }
  // save 0's result
  T x0_result = x;
  T f0_result = fabs(f);


  //start from 1
  x = 1;
  f = ((a * x + b) * x + c) * x + d;
  g = (3 * a * x + 2 * b) * x + c;
  if ((f > 0 && g > 0) || (f < 0 && g < 0)) {
    //printf("use that %.14f, %.14f\n", f_value, g_value);
    for (int i = 0; i < 16; i++) {
      if (fabs(f) < 1e-14f)	return x;
      //printf("use that %.14f; %f\n", f_value, x);
      T step = f / g;
      if (step > 0.1f)	step = 0.1f;
      if (step < -0.1f)	step = -0.1f;
      x -= step;
      f = ((a * x + b) * x + c) * x + d;
      g = 3 * a * x * x + 2 * b * x + c;
    }
  }
  // save 1's result
  T x1_result = x;
  T f1_result = fabs(f);


  //use bipart
  T start_x	= 0;
  T end_x		= 1;
  T start_f	= ((a * start_x + b) * start_x + c) * start_x + d;
  T end_f		= ((a * end_x + b) * end_x + c) * end_x + d;
  T mid_x, mid_f;
  for (int i = 0; i < 16; i++) {
    mid_x = (start_x + end_x) * 0.5f;
    mid_f = ((a * mid_x + b) * mid_x + c) * mid_x + d;

    if ((start_f > 0 && mid_f < 0) || (start_f < 0 && mid_f > 0)) {
      end_x = mid_x;
      end_f = mid_f;
    } else {
      start_x = mid_x;
      start_f = mid_f;
    }
  }
  T	fm_result = fabs(mid_f);

  if ((f0_result < f1_result) && (f0_result < fm_result))		return x0_result;
  else if (f1_result < fm_result)						return x1_result;
  return mid_x;
}

FORCEINLINE
int Cubic_Root(double a, double b, double c, double *root) {
  double delta = 4 * a * a - 12 * b;
  if (delta < 0)	delta = 0;
  //if(delta<0)	{printf("ERROR: only 1/2 root. %f (%f, %f, %f)\n", delta, a, b, c); getchar();}

  delta = sqrt(delta);
  double x0 = (-2 * a - delta) * 0.16666666666f;
  double x1 = (-2 * a + delta) * 0.16666666666f;

  root[0] = Newton_Cubic_Root<double>(1, a, b, c, x0 - 0.1f);
  root[1] = Newton_Cubic_Root<double>(1, a, b, c, (x0 + x1) * 0.5f);
  root[2] = Newton_Cubic_Root<double>(1, a, b, c, x1 + 0.1f);

  //printf("here: %f, %f, %f\n", root[0], root[1], root[2]);

  return 3;


  /*	double q=(3*b-a*a)*0.1111111111111f;
    double r=(9*a*b-27*c-2*a*a*a)*0.0185185185185f;
    double q3=q*q*q;
    double delta=q3+r*r;
    if(delta>1e-9f)
    {
      delta=sqrtf(delta);

      double s, t;
      if(delta+r>0)	s= my_cbrt(delta+r);
      else			s=-my_cbrt(-delta-r);
      if(r-delta>0)	t= my_cbrt(r-delta);
      else			t=-my_cbrt(delta-r);
      root[0]=s+t-a*0.333333333f;
      return 1;
    }
    else
    {
      double rho=sqrtf(fabsf(-q3));


      double length=my_cbrt(rho);
      double cos_theta;
      double r_over_rho=r/rho;

      if(rho<1e-9f)		r_over_rho= 1;
      if(r_over_rho> 1)	r_over_rho= 1;
      if(r_over_rho<-1)	r_over_rho=-1;

      cos_theta=cosf(acosf(r_over_rho)*0.333333333f);

    printf("ABC: %f, %f, %f; %f, %f, %f\n", a, b, c, q, r, delta);
    printf("cos: %f; %f; %f, %f\n", cos_theta, r/rho, r, rho);

      double l_cos=length*cos_theta;
      double l_sin=length*sqrtf(1-cos_theta*cos_theta)*1.73205081f;
      double a_over_3=a*0.333333333f;
      root[0]=2*l_cos-a_over_3;
      root[1]= -l_cos-a_over_3-l_sin;
      root[2]= -l_cos-a_over_3+l_sin;
      return 3;
    }*/
  return 0;
}

inline
void Eigenvalue_Decomposition_32(double *A, double *ev, double *V, int id = 0) {
  (void) id;
  double a = -A[0] - A[4] - A[8];
  double b = A[0] * A[4] + A[0] * A[8] + A[4] * A[8] - A[2] * A[6] - A[5] * A[7] - A[1] * A[3];
  double c = -A[0] * A[4] * A[8] - A[2] * A[3] * A[7] - A[1] * A[5] * A[6] + A[2] * A[4] * A[6] + A[0] * A[5] * A[7] + A[1] * A[3] * A[8];

  //if(id==203)	printf("abc:%f, %f, %f\n", a, b, c);

  int root_number = Cubic_Root(a, b, c, ev);
  if (root_number == 1) {
    printf("WARNING: only found 1 root.\n");
    return;
  }
  if (V == 0)	return;

  //printf("ABC: %f, %f, %f (%f, %f, %f)\n", a, b, c, ev[0], ev[1], ev[2]);


  double sub0 = fabs(A[0] - ev[0]) + fabs(A[4] - ev[0]) + fabs(A[8] - ev[0]);
  double sub1 = fabs(A[0] - ev[1]) + fabs(A[4] - ev[1]) + fabs(A[8] - ev[1]);
  double sub2 = fabs(A[0] - ev[2]) + fabs(A[4] - ev[2]) + fabs(A[8] - ev[2]);

  if (sub1 > sub0)	{
    Swap(ev[0], ev[1]);
    Swap(sub0, sub1);
  }
  if (sub2 > sub0)	{
    Swap(ev[0], ev[2]);
    Swap(sub0, sub2);
  }
  if (sub1 > sub2)	{
    Swap(ev[2], ev[1]);
    Swap(sub2, sub1);
  }


  double TA[9];
  TA[0] = A[0] - ev[0];
  TA[1] = A[1];
  TA[2] = A[2];
  TA[3] = A[3];
  TA[4] = A[4] - ev[0];
  TA[5] = A[5];
  TA[6] = A[6];
  TA[7] = A[7];
  TA[8] = A[8] - ev[0];

  double det_0 = TA[4] * TA[8] - TA[5] * TA[7];
  double det_1 = TA[3] * TA[8] - TA[5] * TA[6];
  double det_2 = TA[3] * TA[7] - TA[4] * TA[6];
  double det_4 = TA[0] * TA[8] - TA[2] * TA[6];
  double det_5 = TA[0] * TA[7] - TA[1] * TA[6];
  double det_8 = TA[0] * TA[4] - TA[1] * TA[3];
  double	max_det = 0;
  int		det_id;
  if (fabs(det_0) > max_det)	{
    max_det = fabs(det_0);
    det_id = 0;
  }
  if (fabs(det_1) > max_det)	{
    max_det = fabs(det_1);
    det_id = 1;
  }
  if (fabs(det_2) > max_det)	{
    max_det = fabs(det_2);
    det_id = 2;
  }
  if (fabs(det_4) > max_det)	{
    max_det = fabs(det_4);
    det_id = 4;
  }
  if (fabs(det_5) > max_det)	{
    max_det = fabs(det_5);
    det_id = 5;
  }
  if (fabs(det_8) > max_det)	{
    max_det = fabs(det_8);
    det_id = 8;
  }

  if (max_det < 1e-6f) {
    double	max_ele = 0;
    int		max_id;

    if (fabs(TA[0]) > max_ele)	{
      max_ele = fabs(TA[0]);
      max_id = 0;
    }
    if (fabs(TA[1]) > max_ele)	{
      max_ele = fabs(TA[1]);
      max_id = 1;
    }
    if (fabs(TA[2]) > max_ele)	{
      max_ele = fabs(TA[2]);
      max_id = 2;
    }
    if (fabs(TA[3]) > max_ele)	{
      max_ele = fabs(TA[3]);
      max_id = 3;
    }
    if (fabs(TA[4]) > max_ele)	{
      max_ele = fabs(TA[4]);
      max_id = 4;
    }
    if (fabs(TA[5]) > max_ele)	{
      max_ele = fabs(TA[5]);
      max_id = 5;
    }
    if (fabs(TA[6]) > max_ele)	{
      max_ele = fabs(TA[6]);
      max_id = 6;
    }
    if (fabs(TA[7]) > max_ele)	{
      max_ele = fabs(TA[7]);
      max_id = 7;
    }
    if (fabs(TA[8]) > max_ele)	{
      max_ele = fabs(TA[8]);
      max_id = 8;
    }

    if (max_ele < 1e-10f)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = 0;
    } else if (max_id == 0)	{
      V[1] = 1;
      V[2] = 0;
      V[0] = -TA[1] / TA[0];
    } else if (max_id == 1)	{
      V[0] = 1;
      V[2] = 0;
      V[1] = -TA[0] / TA[1];
    } else if (max_id == 2)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = -TA[0] / TA[2];
    } else if (max_id == 3)	{
      V[1] = 1;
      V[2] = 0;
      V[0] = -TA[4] / TA[3];
    } else if (max_id == 4)	{
      V[0] = 1;
      V[2] = 0;
      V[1] = -TA[3] / TA[4];
    } else if (max_id == 5)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = -TA[3] / TA[5];
    } else if (max_id == 6)	{
      V[1] = 1;
      V[2] = 0;
      V[0] = -TA[7] / TA[6];
    } else if (max_id == 7)	{
      V[0] = 1;
      V[2] = 0;
      V[1] = -TA[6] / TA[7];
    } else if (max_id == 8)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = -TA[6] / TA[8];
    }
  } else if (det_id == 0)		{
    double inv_det = 1.0f / det_0;
    V[0] = 1;
    V[1] = -det_1 * inv_det;
    V[2] = det_2 * inv_det;
  } else if (det_id == 1)		{
    double inv_det = 1.0f / det_1;
    V[1] = 1;
    V[0] = -det_0 * inv_det;
    V[2] = -det_2 * inv_det;
  } else if (det_id == 2)		{
    double inv_det = 1.0f / det_2;
    V[2] = 1;
    V[1] = -det_1 * inv_det;
    V[0] = det_0 * inv_det;
  } else if (det_id == 4)		{
    double inv_det = 1.0f / det_4;
    V[1] = 1;
    V[0] = -det_1 * inv_det;
    V[2] = -det_5 * inv_det;
  } else if (det_id == 5)		{
    double inv_det = 1.0f / det_5;
    V[2] = 1;
    V[1] = -det_4 * inv_det;
    V[0] = det_1 * inv_det;
  } else if (det_id == 8)		{
    double inv_det = 1.0f / det_8;
    V[2] = 1;
    V[1] = -det_5 * inv_det;
    V[0] = det_2 * inv_det;
  }

  //Normalize eigenvector 0
  double inv_length = 1.0f / sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
  V[0] *= inv_length;
  V[1] *= inv_length;
  V[2] *= inv_length;


  //Get new TA
  TA[0] = A[0] - ev[1];
  TA[4] = A[4] - ev[1];
  TA[8] = A[8] - ev[1];

  double det[3][3];
  det[0][0] = V[1] * TA[2] - V[2] * TA[1];
  det[0][1] = V[0] * TA[2] - V[2] * TA[0];
  det[0][2] = V[0] * TA[1] - V[1] * TA[0];
  det[1][0] = V[1] * TA[5] - V[2] * TA[4];
  det[1][1] = V[0] * TA[5] - V[2] * TA[3];
  det[1][2] = V[0] * TA[4] - V[1] * TA[3];
  det[2][0] = V[1] * TA[8] - V[2] * TA[7];
  det[2][1] = V[0] * TA[8] - V[2] * TA[6];
  det[2][2] = V[0] * TA[7] - V[1] * TA[6];

  max_det = 0;
  int r_id, c_id;
  if (fabs(det[0][0]) > max_det)	{
    max_det = fabs(det[0][0]);
    r_id = 0;
    c_id = 0;
  }
  if (fabs(det[0][1]) > max_det)	{
    max_det = fabs(det[0][1]);
    r_id = 0;
    c_id = 1;
  }
  if (fabs(det[0][2]) > max_det)	{
    max_det = fabs(det[0][2]);
    r_id = 0;
    c_id = 2;
  }
  if (fabs(det[1][0]) > max_det)	{
    max_det = fabs(det[1][0]);
    r_id = 1;
    c_id = 0;
  }
  if (fabs(det[1][1]) > max_det)	{
    max_det = fabs(det[1][1]);
    r_id = 1;
    c_id = 1;
  }
  if (fabs(det[1][2]) > max_det)	{
    max_det = fabs(det[1][2]);
    r_id = 1;
    c_id = 2;
  }
  if (fabs(det[2][0]) > max_det)	{
    max_det = fabs(det[2][0]);
    r_id = 2;
    c_id = 0;
  }
  if (fabs(det[2][1]) > max_det)	{
    max_det = fabs(det[2][1]);
    r_id = 2;
    c_id = 1;
  }
  if (fabs(det[2][2]) > max_det)	{
    max_det = fabs(det[2][2]);
    r_id = 2;
    c_id = 2;
  }


  if (max_det < 1e-6f) {
    if (fabs(V[2]) < fabs(V[0]) && fabs(V[2]) < fabs(V[1])) {
      V[3] = V[1];
      V[4] = -V[0];
      V[5] = 0;
    } else if (fabs(V[0]) < fabs(V[1])) {
      V[3] = 0;
      V[4] = V[2];
      V[5] = -V[1];
    } else {
      V[3] = V[2];
      V[4] = 0;
      V[5] = -V[0];
    }
  } else if (c_id == 0) {
    double inv_det = 1.0f / det[r_id][c_id];
    V[3] = 1;
    V[4] = -det[r_id][1] * inv_det;
    V[5] = det[r_id][2] * inv_det;
  } else if (c_id == 1) {
    double inv_det = 1.0f / det[r_id][c_id];
    V[4] = 1;
    V[3] = -det[r_id][0] * inv_det;
    V[5] = -det[r_id][2] * inv_det;

  } else if (c_id == 2) {
    double inv_det = 1.0f / det[r_id][c_id];
    V[5] = 1;
    V[4] = -det[r_id][1] * inv_det;
    V[3] = det[r_id][0] * inv_det;
  }

  inv_length = 1.0f / sqrt(V[3] * V[3] + V[4] * V[4] + V[5] * V[5]);
  V[3] *= inv_length;
  V[4] *= inv_length;
  V[5] *= inv_length;

  Cross(&V[0], &V[3], &V[6]);


  //printf("TA: %f, %f, %f; %f, %f, %f; %f, %f, %f\n", TA[0], TA[1], TA[2], TA[3], TA[4], TA[5], TA[6], TA[7], TA[8]);
  //printf("here: %d (%f)\n", det_id, max_det);
  //printf("v: %f, %f, %f \n", V[3], V[4], V[5]);


  //if( fabsf((A[0]-ev[0])*V[0]+A[1]*V[1]+A[2]*V[2])>1e-6f &&
  //	fabsf(A[3]*V[0]+(A[4]-ev[0])*V[1]+A[5]*V[2])>1e-6f &&
  //	fabsf(A[6]*V[0]+A[7]*V[1]+(A[8]-ev[0])*V[2])>1e-6f)
  //		printf("T0: %f, %f, %f\n", (A[0]-ev[0])*V[0]+A[1]*V[1]+A[2]*V[2], A[3]*V[0]+(A[4]-ev[0])*V[1]+A[5]*V[2], A[6]*V[0]+A[7]*V[1]+(A[8]-ev[0])*V[2]);

  //if( fabsf((A[0]-ev[1])*V[3]+A[1]*V[4]+A[2]*V[5])>1e-6f &&
  //	fabsf(A[3]*V[3]+(A[4]-ev[1])*V[4]+A[5]*V[5])>1e-6f &&
  //	fabsf(A[6]*V[3]+A[7]*V[4]+(A[8]-ev[1])*V[5])>1e-6f)
  //		printf("T1: %f, %f, %f\n", (A[0]-ev[1])*V[3]+A[1]*V[4]+A[2]*V[5], A[3]*V[3]+(A[4]-ev[1])*V[4]+A[5]*V[5], A[6]*V[3]+A[7]*V[4]+(A[8]-ev[1])*V[5]);
  //printf("T2: %f, %f, %f\n", (A[0]-ev[2])*V[6]+A[1]*V[7]+A[2]*V[8], A[3]*V[6]+(A[4]-ev[2])*V[7]+A[5]*V[8], A[6]*V[6]+A[7]*V[7]+(A[8]-ev[2])*V[8]);


  //getchar();
}

template <class T> FORCEINLINE
void Eigenvalue_Decomposition_3(T *A, T *ev, T *V, int id = 0) {
  (void) id;
  T a = -A[0] - A[4] - A[8];
  T b = A[0] * A[4] + A[0] * A[8] + A[4] * A[8] - A[2] * A[6] - A[5] * A[7] - A[1] * A[3];
  T c = -A[0] * A[4] * A[8] - A[2] * A[3] * A[7] - A[1] * A[5] * A[6] + A[2] * A[4] * A[6] + A[0] * A[5] * A[7] + A[1] * A[3] * A[8];

  //if(id==203)	printf("abc:%f, %f, %f\n", a, b, c);

  int root_number = Cubic_Root(a, b, c, ev);
  if (root_number == 1) {
    printf("WARNING: only found 1 root.\n");
    return;
  }
  if (V == 0)	return;

  //printf("ABC: %f, %f, %f (%f, %f, %f)\n", a, b, c, ev[0], ev[1], ev[2]);


  T sub0 = fabsf(A[0] - ev[0]) + fabsf(A[4] - ev[0]) + fabsf(A[8] - ev[0]);
  T sub1 = fabsf(A[0] - ev[1]) + fabsf(A[4] - ev[1]) + fabsf(A[8] - ev[1]);
  T sub2 = fabsf(A[0] - ev[2]) + fabsf(A[4] - ev[2]) + fabsf(A[8] - ev[2]);

  if (sub1 > sub0)	{
    Swap(ev[0], ev[1]);
    Swap(sub0, sub1);
  }
  if (sub2 > sub0)	{
    Swap(ev[0], ev[2]);
    Swap(sub0, sub2);
  }
  if (sub1 > sub2)	{
    Swap(ev[2], ev[1]);
    Swap(sub2, sub1);
  }


  T TA[9];
  TA[0] = A[0] - ev[0];
  TA[1] = A[1];
  TA[2] = A[2];
  TA[3] = A[3];
  TA[4] = A[4] - ev[0];
  TA[5] = A[5];
  TA[6] = A[6];
  TA[7] = A[7];
  TA[8] = A[8] - ev[0];

  T det_0 = TA[4] * TA[8] - TA[5] * TA[7];
  T det_1 = TA[3] * TA[8] - TA[5] * TA[6];
  T det_2 = TA[3] * TA[7] - TA[4] * TA[6];
  T det_4 = TA[0] * TA[8] - TA[2] * TA[6];
  T det_5 = TA[0] * TA[7] - TA[1] * TA[6];
  T det_8 = TA[0] * TA[4] - TA[1] * TA[3];
  T max_det = 0;
  int		det_id;
  if (fabsf(det_0) > max_det)	{
    max_det = fabsf(det_0);
    det_id = 0;
  }
  if (fabsf(det_1) > max_det)	{
    max_det = fabsf(det_1);
    det_id = 1;
  }
  if (fabsf(det_2) > max_det)	{
    max_det = fabsf(det_2);
    det_id = 2;
  }
  if (fabsf(det_4) > max_det)	{
    max_det = fabsf(det_4);
    det_id = 4;
  }
  if (fabsf(det_5) > max_det)	{
    max_det = fabsf(det_5);
    det_id = 5;
  }
  if (fabsf(det_8) > max_det)	{
    max_det = fabsf(det_8);
    det_id = 8;
  }

  if (max_det < 1e-6f) {
    T max_ele = 0;
    int		max_id;

    if (fabsf(TA[0]) > max_ele)	{
      max_ele = fabsf(TA[0]);
      max_id = 0;
    }
    if (fabsf(TA[1]) > max_ele)	{
      max_ele = fabsf(TA[1]);
      max_id = 1;
    }
    if (fabsf(TA[2]) > max_ele)	{
      max_ele = fabsf(TA[2]);
      max_id = 2;
    }
    if (fabsf(TA[3]) > max_ele)	{
      max_ele = fabsf(TA[3]);
      max_id = 3;
    }
    if (fabsf(TA[4]) > max_ele)	{
      max_ele = fabsf(TA[4]);
      max_id = 4;
    }
    if (fabsf(TA[5]) > max_ele)	{
      max_ele = fabsf(TA[5]);
      max_id = 5;
    }
    if (fabsf(TA[6]) > max_ele)	{
      max_ele = fabsf(TA[6]);
      max_id = 6;
    }
    if (fabsf(TA[7]) > max_ele)	{
      max_ele = fabsf(TA[7]);
      max_id = 7;
    }
    if (fabsf(TA[8]) > max_ele)	{
      max_ele = fabsf(TA[8]);
      max_id = 8;
    }

    if (max_ele < 1e-10f)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = 0;
    } else if (max_id == 0)	{
      V[1] = 1;
      V[2] = 0;
      V[0] = -TA[1] / TA[0];
    } else if (max_id == 1)	{
      V[0] = 1;
      V[2] = 0;
      V[1] = -TA[0] / TA[1];
    } else if (max_id == 2)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = -TA[0] / TA[2];
    } else if (max_id == 3)	{
      V[1] = 1;
      V[2] = 0;
      V[0] = -TA[4] / TA[3];
    } else if (max_id == 4)	{
      V[0] = 1;
      V[2] = 0;
      V[1] = -TA[3] / TA[4];
    } else if (max_id == 5)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = -TA[3] / TA[5];
    } else if (max_id == 6)	{
      V[1] = 1;
      V[2] = 0;
      V[0] = -TA[7] / TA[6];
    } else if (max_id == 7)	{
      V[0] = 1;
      V[2] = 0;
      V[1] = -TA[6] / TA[7];
    } else if (max_id == 8)	{
      V[0] = 1;
      V[1] = 0;
      V[2] = -TA[6] / TA[8];
    }
  } else if (det_id == 0)		{
    T inv_det = 1.0f / det_0;
    V[0] = 1;
    V[1] = -det_1 * inv_det;
    V[2] = det_2 * inv_det;
  } else if (det_id == 1)		{
    T inv_det = 1.0f / det_1;
    V[1] = 1;
    V[0] = -det_0 * inv_det;
    V[2] = -det_2 * inv_det;
  } else if (det_id == 2)		{
    T inv_det = 1.0f / det_2;
    V[2] = 1;
    V[1] = -det_1 * inv_det;
    V[0] = det_0 * inv_det;
  } else if (det_id == 4)		{
    T inv_det = 1.0f / det_4;
    V[1] = 1;
    V[0] = -det_1 * inv_det;
    V[2] = -det_5 * inv_det;
  } else if (det_id == 5)		{
    T inv_det = 1.0f / det_5;
    V[2] = 1;
    V[1] = -det_4 * inv_det;
    V[0] = det_1 * inv_det;
  } else if (det_id == 8)		{
    T inv_det = 1.0f / det_8;
    V[2] = 1;
    V[1] = -det_5 * inv_det;
    V[0] = det_2 * inv_det;
  }

  //Normalize eigenvector 0
  T inv_length = 1.0f / sqrtf(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
  V[0] *= inv_length;
  V[1] *= inv_length;
  V[2] *= inv_length;


  //Get new TA
  TA[0] = A[0] - ev[1];
  TA[4] = A[4] - ev[1];
  TA[8] = A[8] - ev[1];

  T det[3][3];
  det[0][0] = V[1] * TA[2] - V[2] * TA[1];
  det[0][1] = V[0] * TA[2] - V[2] * TA[0];
  det[0][2] = V[0] * TA[1] - V[1] * TA[0];
  det[1][0] = V[1] * TA[5] - V[2] * TA[4];
  det[1][1] = V[0] * TA[5] - V[2] * TA[3];
  det[1][2] = V[0] * TA[4] - V[1] * TA[3];
  det[2][0] = V[1] * TA[8] - V[2] * TA[7];
  det[2][1] = V[0] * TA[8] - V[2] * TA[6];
  det[2][2] = V[0] * TA[7] - V[1] * TA[6];

  max_det = 0;
  int r_id, c_id;
  if (fabsf(det[0][0]) > max_det)	{
    max_det = fabsf(det[0][0]);
    r_id = 0;
    c_id = 0;
  }
  if (fabsf(det[0][1]) > max_det)	{
    max_det = fabsf(det[0][1]);
    r_id = 0;
    c_id = 1;
  }
  if (fabsf(det[0][2]) > max_det)	{
    max_det = fabsf(det[0][2]);
    r_id = 0;
    c_id = 2;
  }
  if (fabsf(det[1][0]) > max_det)	{
    max_det = fabsf(det[1][0]);
    r_id = 1;
    c_id = 0;
  }
  if (fabsf(det[1][1]) > max_det)	{
    max_det = fabsf(det[1][1]);
    r_id = 1;
    c_id = 1;
  }
  if (fabsf(det[1][2]) > max_det)	{
    max_det = fabsf(det[1][2]);
    r_id = 1;
    c_id = 2;
  }
  if (fabsf(det[2][0]) > max_det)	{
    max_det = fabsf(det[2][0]);
    r_id = 2;
    c_id = 0;
  }
  if (fabsf(det[2][1]) > max_det)	{
    max_det = fabsf(det[2][1]);
    r_id = 2;
    c_id = 1;
  }
  if (fabsf(det[2][2]) > max_det)	{
    max_det = fabsf(det[2][2]);
    r_id = 2;
    c_id = 2;
  }


  if (max_det < 1e-6f) {
    if (fabsf(V[2]) < fabsf(V[0]) && fabsf(V[2]) < fabsf(V[1])) {
      V[3] = V[1];
      V[4] = -V[0];
      V[5] = 0;
    } else if (fabsf(V[0]) < fabsf(V[1])) {
      V[3] = 0;
      V[4] = V[2];
      V[5] = -V[1];
    } else {
      V[3] = V[2];
      V[4] = 0;
      V[5] = -V[0];
    }
  } else if (c_id == 0) {
    T inv_det = 1.0f / det[r_id][c_id];
    V[3] = 1;
    V[4] = -det[r_id][1] * inv_det;
    V[5] = det[r_id][2] * inv_det;
  } else if (c_id == 1) {
    T inv_det = 1.0f / det[r_id][c_id];
    V[4] = 1;
    V[3] = -det[r_id][0] * inv_det;
    V[5] = -det[r_id][2] * inv_det;

  } else if (c_id == 2) {
    T inv_det = 1.0f / det[r_id][c_id];
    V[5] = 1;
    V[4] = -det[r_id][1] * inv_det;
    V[3] = det[r_id][0] * inv_det;
  }

  inv_length = 1.0f / sqrtf(V[3] * V[3] + V[4] * V[4] + V[5] * V[5]);
  V[3] *= inv_length;
  V[4] *= inv_length;
  V[5] *= inv_length;

  Cross(&V[0], &V[3], &V[6]);


  //printf("TA: %f, %f, %f; %f, %f, %f; %f, %f, %f\n", TA[0], TA[1], TA[2], TA[3], TA[4], TA[5], TA[6], TA[7], TA[8]);
  //printf("here: %d (%f)\n", det_id, max_det);
  //printf("v: %f, %f, %f \n", V[3], V[4], V[5]);


  //if( fabsf((A[0]-ev[0])*V[0]+A[1]*V[1]+A[2]*V[2])>1e-6f &&
  //	fabsf(A[3]*V[0]+(A[4]-ev[0])*V[1]+A[5]*V[2])>1e-6f &&
  //	fabsf(A[6]*V[0]+A[7]*V[1]+(A[8]-ev[0])*V[2])>1e-6f)
  //		printf("T0: %f, %f, %f\n", (A[0]-ev[0])*V[0]+A[1]*V[1]+A[2]*V[2], A[3]*V[0]+(A[4]-ev[0])*V[1]+A[5]*V[2], A[6]*V[0]+A[7]*V[1]+(A[8]-ev[0])*V[2]);

  //if( fabsf((A[0]-ev[1])*V[3]+A[1]*V[4]+A[2]*V[5])>1e-6f &&
  //	fabsf(A[3]*V[3]+(A[4]-ev[1])*V[4]+A[5]*V[5])>1e-6f &&
  //	fabsf(A[6]*V[3]+A[7]*V[4]+(A[8]-ev[1])*V[5])>1e-6f)
  //		printf("T1: %f, %f, %f\n", (A[0]-ev[1])*V[3]+A[1]*V[4]+A[2]*V[5], A[3]*V[3]+(A[4]-ev[1])*V[4]+A[5]*V[5], A[6]*V[3]+A[7]*V[4]+(A[8]-ev[1])*V[5]);
  //printf("T2: %f, %f, %f\n", (A[0]-ev[2])*V[6]+A[1]*V[7]+A[2]*V[8], A[3]*V[6]+(A[4]-ev[2])*V[7]+A[5]*V[8], A[6]*V[6]+A[7]*V[7]+(A[8]-ev[2])*V[8]);


  //getchar();
}

FORCEINLINE
void tred2(double a[4][4], int n, double d[], double e[]) {
  int l, k, j, i;
  double scale, hh, h, g, f;
  for (i = n; i >= 2; i--) {
    l = i - 1;
    h = scale = 0.0;
    if (l > 1) {
      for (k = 1; k <= l; k++)	scale += fabs(a[i][k]);
      if (scale == 0.0) 	e[i] = a[i][l];
      else {
        for (k = 1; k <= l; k++) {
          a[i][k] /= scale;
          h += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale * g;
        h -= f * g;		// Now h is equation (11.2.4).
        a[i][l] = f - g;	// Store u in the ith row of a.
        f = 0.0;
        for (j = 1; j <= l; j++) {
          /* Next statement can be omitted if eigenvectors not wanted */
          a[j][i] = a[i][j] / h;	// Store u=H in ith column of a.
          g = 0.0;				// Form an element of A  u in g.
          for (k = 1; k <= j; k++)
            g += a[j][k] * a[i][k];
          for (k = j + 1; k <= l; k++)
            g += a[k][j] * a[i][k];
          e[j] = g / h;

          f += e[j] * a[i][j];
        }
        hh = f / (h + h); //Form K, equation (11.2.11).
        for (j = 1; j <= l; j++) {
          //Form q and store in e overwriting p.
          f = a[i][j];
          e[j] = g = e[j] - hh * f;
          for (k = 1; k <= j; k++) //Reduce a, equation (11.2.13).
            a[j][k] -= (f * e[k] + g * a[i][k]);
        }
      }
    } else
      e[i] = a[i][l];
    d[i] = h;

    //printf("i %d: %f\n", i, h);
  }

  /* Next statement can be omitted if eigenvectors not wanted */
  d[1] = 0.0;
  e[1] = 0.0;
  /* Contents of this loop can be omitted if eigenvectors not
  wanted except for statement d[i]=a[i][i]; */
  for (i = 1; i <= n; i++) {
    //Begin accumulation of transformationmal=
    l = i - 1; //trices.
    if (d[i]) {
      //This block skipped when i=1.
      for (j = 1; j <= l; j++) {
        g = 0.0;
        for (k = 1; k <= l; k++) //Use u and u=H stored in a to form PQ.
          g += a[i][k] * a[k][j];
        for (k = 1; k <= l; k++)
          a[k][j] -= g * a[k][i];
      }
    }
    d[i] = a[i][i];	//This statement remains.
    a[i][i] = 1.0;	//Reset row and column of a to identity
    for (j = 1; j <= l; j++)
      a[j][i] = a[i][j] = 0.0;	//matrix for next iteration.
  }
}

FORCEINLINE double eigenSign(double a, double b) {
  return ((b) >= 0.0 ? fabs(a) : -fabs(a));
}


FORCEINLINE double eigenPythag(double a, double b) {
  double absa, absb;
  absa = fabs(a);
  absb = fabs(b);

  if (absa > absb) {
    double value = absb / absa;
    return absa * sqrt(1.0 + value * value);
  } else {
    double value = absa / absb;
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + value * value));
  }
}

FORCEINLINE
void tqli(double d[], double e[], int n, double z[4][4]) {
  double pythag(double a, double b);
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  for (i = 2; i <= n; i++) e[i - 1] = e[i]; //Convenient to renumber the ele[n]=0.0; ements of e.
  for (l = 1; l <= n; l++) {
    iter = 0;
    do {
      for (m = l; m <= n - 1; m++) { //Look for a single small subdiagonal
        //element to split
        //the matrix.
        dd = fabs(d[m]) + fabs(d[m + 1]);
        if ((double)(fabs(e[m]) + dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) printf("Too many iterations in tqli");
        g = (d[l + 1] - d[l]) / (2.0 * e[l]); //Form shift.
        r = eigenPythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + eigenSign(r, g)); //This is dm − ks.
        s = c = 1.0;
        p = 0.0;

        for (i = m - 1; i >= l; i--) { //A plane rotation as in the original
          //QL, followed by Givens
          //rotations to restore tridiagonal
          //form.
          f = s * e[i];
          b = c * e[i];
          e[i + 1] = (r = eigenPythag(f, g));
          if (r == 0.0) { //Recover from underflow.
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          d[i + 1] = g + (p = s * r);
          g = c * r - b;
          /* Next loop can be omitted if eigenvectors not wanted*/
          for (k = 1; k <= n; k++) { //Form eigenvectors.
            f = z[k][i + 1];
            z[k][i + 1] = s * z[k][i] + c * f;
            z[k][i] = c * z[k][i] - s * f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
}

inline void  New_Eigenvalue_Decomposition_3(double *A, double *ev, double *V, int id = 0) {
  (void) id;
  double a[4][4];
  a[1][1] = A[0];
  a[1][2] = A[1];
  a[1][3] = A[2];
  a[2][1] = A[3];
  a[2][2] = A[4];
  a[2][3] = A[5];
  a[3][1] = A[6];
  a[3][2] = A[7];
  a[3][3] = A[8];

  double z[4][4];
  (void) z;
  double d[4], e[4];

  tred2(a, 3, d, e);
  //	printf("a: %f, %f, %f\n", a[1][1], a[1][2], a[1][3]);
  //	printf("a: %f, %f, %f\n", a[2][1], a[2][2], a[2][3]);
  //	printf("a: %f, %f, %f\n", a[3][1], a[3][2], a[3][3]);

  tqli(d, e, 3, a);

  ev[0] = d[1];
  ev[1] = d[2];
  ev[2] = d[3];


  V[0] = a[1][1];
  V[1] = a[2][1];
  V[2] = a[3][1];
  V[3] = a[1][2];
  V[4] = a[2][2];
  V[5] = a[3][2];
  V[6] = a[1][3];
  V[7] = a[2][3];
  V[8] = a[3][3];

  //	printf("a: %f, %f, %f\n", a[1][1], a[1][2], a[1][3]);
  //	printf("a: %f, %f, %f\n", a[2][1], a[2][2], a[2][3]);
  //	printf("a: %f, %f, %f\n", a[3][1], a[3][2], a[3][3]);

  /*
    printf("a: %f, %f, %f\n", z[1][1], z[1][2], z[1][3]);
    printf("a: %f, %f, %f\n", z[2][1], z[2][2], z[2][3]);
    printf("a: %f, %f, %f\n", z[3][1], z[3][2], z[3][3]);

    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      V[i*3+j]=0;
      for(int k=0; k<3; k++)
        V[i*3+j]+=a[i+1][k+1]*z[k+1][j+1];
    }*/

}

#define MDot(x, y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
template <class T> FORCEINLINE
void Gaussian_Elimination(T *a, int n, T *b)
{
  int* indxc=new int[n];
  int* indxr=new int[n];
  int* ipiv =new int[n];
  int i,icol,irow,j,k,l,ll;
  T big,dum,pivinv,temp;

  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++)
  {
    big=0.0;
    for (j=0;j<n;j++)
      if (ipiv[j] != 1)
        for (k=0;k<n;k++)
        {
          if (ipiv[k] ==0)
          {
            if (fabs(a[j*n+k]) >= big)
            {
              big=fabs(a[j*n+k]);
              irow=j;
              icol=k;
            }
          }
        }
    ++(ipiv[icol]);

    if (irow != icol)
    {
      for (l=0;l<n;l++) {temp=a[irow*n+l]; a[irow*n+l]=a[icol*n+l]; a[icol*n+l]=temp;}
      temp=b[irow]; b[irow]=b[icol]; b[icol]=temp;
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol*n+icol] == 0.0) printf("Error: Singular Matrix in Gaussian_Elimination.");
    pivinv=1.0/a[icol*n+icol];
    a[icol*n+icol]=1.0;
    for (l=0;l<n;l++) a[icol*n+l] *= pivinv;
    b[icol] *= pivinv;

    for (ll=0;ll<n;ll++)
      if (ll != icol)
      {
        dum=a[ll*n+icol];
        a[ll*n+icol]=0.0;
        for (l=0;l<n;l++) a[ll*n+l] -= a[icol*n+l]*dum;
        b[ll] -= b[icol]*dum;
      }
  }

  for (l=n-1;l>1;l--)
  {
    if (indxr[l] != indxc[l])
    for (k=0;k<n;k++)
    {
      temp=a[k*n+indxr[l]];
      a[k*n+indxr[l]]=a[k*n+indxc[l]];
      a[k*n+indxc[l]]=temp;
    }
  }
  delete []ipiv;
  delete []indxr;
  delete []indxc;
}




#endif
