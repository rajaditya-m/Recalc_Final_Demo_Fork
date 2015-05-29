//**************************************************************************************
//  Copyright (C) 2002 - 2004, Huamin Wang
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//**************************************************************************************
// MY_MATH library.
// Contents:
// Advanced math functions, Min, Max functions, floating value to integer functions,
// Random functions, and miscs.
//**************************************************************************************

#ifndef __MY_MATH_H__
#define __MY_MATH_H__

#define MY_PI		3.14159265358979323846f
#define PI_180		0.0174532925f
#define INV_PI		0.31830988618379067154f
#define INV_TWOPI	0.15915494309189533577f
#define MY_INFINITY 999999999
#define MY_INFINITE 999999999
#define INV_255		0.00392156862745098039f
#ifndef EPSILON
#define EPSILON		1e-7f
#endif


#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef FORCEINLINE
#define FORCEINLINE inline
#endif

//**************************************************************************************
// Math Functions.
//**************************************************************************************
template<class T> FORCEINLINE
T Sqr(const T a)
{
    return a*a;
}

template<class T> FORCEINLINE
T Cube(const T a)
{
    return a*a*a;
}

FORCEINLINE
int Mod(int a, int b)
{
    int n = int(a/b);
    a -= n*b;
    if (a<0) a += b;
    return a;
}

FORCEINLINE
int Period(const int x, const int l, const int u)
{
    return Mod(x-l, u-l+1)+l;
}

FORCEINLINE
float Log_2(float x)
{
    static float invLog2 = 1.f / logf(2.f);
    return logf(x) * invLog2;
}

FORCEINLINE
bool Is_Power2(int v)
{
    return (v & (v - 1)) == 0;
}

FORCEINLINE
unsigned int Round_Up_Power2(unsigned int v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v+1;
}

//**************************************************************************************
// Min, Max Functions.
//**************************************************************************************

#define MMin(a, b)	((a)<(b)?(a):(b))

#define MMax(a, b)	((a)>(b)?(a):(b))

#define Clamp(a, low, high)  (((a)>(high))?(high):(((a)<(low))?(low):(a)))

template<class T> FORCEINLINE
T Min(const T a,const T b)
{
    return a<b?a:b;
}

template<class T> FORCEINLINE
T Min(const T a,const T b,const T c)
{
    if(a<b && a<c) return a;
    if(b<c) return b;
    return c;
}

template<class T> FORCEINLINE
T Min(const T a,const T b,const T c, const T d)
{
    T r=a;
    if(b<r) r=b;
    if(c<r) r=c;
    if(d<r) r=d;
    return r;
}

template<class T> FORCEINLINE
T Min(const T a,const T b,const T c, const T d, const T e)
{
    T r=a;
    if(b<r) r=b;
    if(c<r) r=c;
    if(d<r) r=d;
    if(e<r)	r=e;
    return r;
}

template<class T> FORCEINLINE
T Max(const T a,const T b)
{
    if(a>b) return a;
    return b;
}

template<class T> FORCEINLINE
T Max(const T a,const T b,const T c)
{
    if(a>b && a>c) return a;
    if(b>c) return b;
    return c;
}

template<class T> FORCEINLINE
T Max(const T a,const T b,const T c,const T d)
{
    T r=a;
    if(b>r) r=b;
    if(c>r) r=c;
    if(d>r) r=d;
    return r;
}

template<class T> FORCEINLINE
T Max(const T a,const T b,const T c,const T d,const T e)
{
    T r=a;
    if(b>r) r=b;
    if(c>r) r=c;
    if(d>r) r=d;
    if(e>r) r=e;
    return r;
}

template<class T> FORCEINLINE
T Max_Abs(const T a, const T b, const T c, const T d, const T e)
{
    return Max(fabs(a),fabs(b),fabs(c),fabs(d),fabs(e));
}

template<class T> FORCEINLINE
T Max_By_Abs(const T a, const T b)
{
    return fabsf(a)>fabsf(b)?a:b;
}

template<class T> FORCEINLINE
T Min_By_Abs(const T a, const T b)
{
    return fabsf(a)<fabsf(b)?a:b;
}

template<class T> FORCEINLINE
T Min_By_Abs(const T a, const T b, const T c)
{
    if(fabs(a)<fabs(b) && fabs(a)<fabs(c)) return a;
    if(fabs(b)<fabs(c)) return b;
    return c;
}

//**************************************************************************************
// Integer Functions.
//**************************************************************************************
#if (defined(__linux__) && defined(__i386__)) || defined(WIN32)
#define FAST_INT 1
#endif
#define _doublemagicroundeps		(.5-1.4e-11)	//almost .5f = .5f - 1e^(number of exp bit)

FORCEINLINE
int Round(double val)
{
#ifdef FAST_INT
#define _doublemagic				double(6755399441055744.0)	//2^52 * 1.5,  uses limited precision to floor
    val		= val + _doublemagic;
    return ((long*)&val)[0];
#undef _doublemagic
#else
    return int (val+_doublemagicroundeps);
#endif
}

FORCEINLINE
int Float_to_Int(double val)
{
#ifdef FAST_INT
    return (val<0) ?  Round(val+_doublemagicroundeps) :
           Round(val-_doublemagicroundeps);
#else
    return (int)val;
#endif
}

FORCEINLINE
int Floor(double val)
{
#ifdef FAST_INT
    return Round(val - _doublemagicroundeps);
#else
    return (int)floorf(val);
#endif
}

FORCEINLINE
int Ceiling(double val)
{
#ifdef FAST_INT
    return Round(val + _doublemagicroundeps);
#else
    return (int)ceilf(val);
#endif
}

//**************************************************************************************
// Math Functions.
//**************************************************************************************

/*template<class T> FORCEINLINE
T Sign(const T a)
{
	if(a<0) return -1;
	return 1;
}*/

template<class T> FORCEINLINE
void Swap(T &a, T &b)
{
    T c=a;
    a=b;
    b=c;
}

template<class T> FORCEINLINE
T Radian(T degree)
{
    return ((T)MY_PI/180.f) * degree;
}

template<class T> FORCEINLINE
T Degree(T radian)
{
    return (180.f/(T)MY_PI) * radian;
}

//**************************************************************************************
// Interpolation Functions.
//**************************************************************************************
template<class T> FORCEINLINE
T Lerp(const T c0, const T c1, const T a)
{
    return (1-a)*c0+a*c1;
}

template<class T> FORCEINLINE
T Lerp(const T c00, const T c10, const T c01, const T c11, const T a, const T b)
{
    return (1-b)*((1-a)*c00+a*c10) + b*((1-a)*c01+a*c11);
}

template<class T> FORCEINLINE
T Lerp(	const T c000, const T c100, const T c010, const T c110,
        const T c001, const T c101, const T c011, const T c111,
        const T a, const T b, const T c)
{
    return (1-c)*( (1-b)*((1-a)*c000+a*c100) + b*((1-a)*c010+a*c110) )+
           c*( (1-b)*((1-a)*c001+a*c101) + b*((1-a)*c011+a*c111) );
}

//**************************************************************************************
// Random Functions. RandomFloat and RandomUInt
//**************************************************************************************
//  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//	are met:
//		1.	Redistributions of source code must retain the above copyright
//			notice, this list of conditions and the following disclaimer.
//		2.	Redistributions in binary form must reproduce the above copyright
//			notice, this list of conditions and the following disclaimer in the
//			documentation and/or other materials provided with the distribution.
//		3.	The names of its contributors may not be used to endorse or promote
//			products derived from this software without specific prior written
//			permission.
//	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
static unsigned long mt[N];		/* the array for the state vector  */
static int mti=N+1;				/* mti==N+1 means mt[N] is not initialized */

FORCEINLINE
void init_genrand(unsigned long seed)
{
    mt[0]= seed & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}
FORCEINLINE
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]= {0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* default initial seed */

        for (kk=0; kk<N-M; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk<N-1; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
    y = mt[mti++];
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}

FORCEINLINE
float RandomFloat()			/* generates a random number on [0,1)-Real-interval */
{
    return genrand_int32()*((float)1.0/(float)4294967296.0);	/* divided by 2^32 */
}

FORCEINLINE
float RandomFloat2()		/* generates a random number on [0,1]-Real-interval */
{
    return genrand_int32()*((float)1.0/(float)4294967295.0);	/* divided by 2^32-1 */
}

FORCEINLINE
unsigned long RandomUInt()
{
    return genrand_int32();
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK


//**************************************************************************************
// Root Finding Functions.
//**************************************************************************************
FORCEINLINE
float cbrt_5f(float f)
{
    unsigned int* p = (unsigned int *) &f;
    *p = *p/3 + 709921077;
    return f;
}

FORCEINLINE
float cbrta_newtonf(const float a, const float x)
{
    return a - (1.0f / 3.0f) * (a - x / (a*a));
}

FORCEINLINE
float my_cbrt(float d)
{
    float a = cbrt_5f(d);
    a = cbrta_newtonf(a, d);
    a = cbrta_newtonf(a, d);
    a = cbrta_newtonf(a, d);
    return cbrta_newtonf(a, d);
}

//Analytic Solution of x^3+a*x^2+b*x+c=0;
FORCEINLINE
int Cubic_Root(float a, float b, float c, float *root)
{
    float q=(3*b-a*a)*0.1111111111111f;
    float r=(9*a*b-27*c-2*a*a*a)*0.0185185185185f;
    float q3=q*q*q;
    float delta=q3+r*r;
    if(delta>0)
    {
        delta=sqrtf(delta);
        float s, t;
        if(r+delta<0)	s=-my_cbrt(-r-delta);
        else			s= my_cbrt( r+delta);
        if(r-delta<0)	t=-my_cbrt(-r+delta);
        else			t= my_cbrt( r-delta);
        root[0]=s+t-a*0.333333333f;
        return 1;
    }
    else
    {
        float rho=sqrtf(-q3);
        float length=my_cbrt(rho);
        float cos_theta=cosf(acosf(r/rho)*0.333333333f);
        float l_cos=length*cos_theta;
        float l_sin=length*sqrtf(1-cos_theta*cos_theta)*1.73205081f;
        float a_over_3=a*0.333333333f;
        root[0]=2*l_cos-a_over_3;
        root[1]= -l_cos-a_over_3-l_sin;
        root[2]= -l_cos-a_over_3+l_sin;
        return 3;
    }
    return 0;
}



//Analytic solution of ax^2+bx+c=0
template <class T> FORCEINLINE
int Quad_Root(T a, T b, T c, T *roots)
{
    T delta=b*b-4*a*c;
    if(delta>-1e-10f)
    {
        if(delta<0)	delta=0;
        delta=sqrtf(delta);
        roots[0]=(-b-delta)/(2*a);
        roots[1]=(-b+delta)/(2*a);
        return 2;
    }
    return 0;
}


//**************************************************************************************
// Misc Functions.
//**************************************************************************************
FORCEINLINE
void Default() {
    printf("Calling default function.\n");
}

FORCEINLINE
void Replace_Extension(char *filename, char *ext)
{
    int i,j;
    for(i=(int)strlen(filename)-1; i>=0 ; i--)
        if(filename[i]=='.') break;
    for(j=0; j<(int)strlen(ext); j++)
        filename[i+j+1]=ext[j];
    filename[i+j+1]='\0';
}

FORCEINLINE
void Color_Convert(float i, float *v, float min=-0.1, float max=0.1)
{

    i=(i-min)/(max-min);
    if(i>=1)	{
        v[0]=1;
        v[1]=0;
        v[2]=0;
        return;
    }
    if(i>=0.75)	{
        v[0]=1;
        v[1]=4-4*i;
        v[2]=0;
        return;
    }
    if(i>=0.5)	{
        v[0]=4*i-2;
        v[1]=1;
        v[2]=0;
        return;
    }
    if(i>=0.25)	{
        v[0]=0;
        v[1]=1;
        v[2]=2-4*i;
        return;
    }
    if(i>=0)	{
        v[0]=0;
        v[1]=4*i;
        v[2]=1;
        return;
    }
    v[0]=0;
    v[1]=0;
    v[2]=1;
}


FORCEINLINE
int Quick_Sort_Partition( float * array, int * index, int p, int r)
{
    float x=array[r];
    int i=p-1;
    for(int j=p; j<=r-1; j++)
        if( array[j]<=x)
        {
            i++;
            Swap( array[i], array[j]);
            Swap( index[i], index[j]);
        }
    Swap(array[i+1], array[r]);
    Swap(index[i+1], index[r]);
    return i+1;
}

FORCEINLINE
void Quick_Sort(float * array, int * index, int p, int r)
{
    if(p<r)
    {
        int q=Quick_Sort_Partition(array, index, p, r);
        Quick_Sort (array, index, p, q-1);
        Quick_Sort (array, index, q+1, r);
    }
}



#endif
