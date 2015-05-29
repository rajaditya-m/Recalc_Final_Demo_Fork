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

#ifndef _MATHDEFS_H_
#define _MATHDEFS_H_

#define	OK   1
#define	ERR -1

#ifndef TRUE
	#define	TRUE 1
#endif 
#ifndef FALSE
#define	FALSE 0
#endif 

#define	EXIT 0

//#define	MAX_LINE 250
//#define	VECSIZE	4
//#define	MAX_ARGS 50

//typedef	int myBOOL;
//#define	STRLEN 100
//typedef	char STR[STRLEN];

#ifdef WIN32
	#include <windows.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <cstring>

#ifndef	M_PI
#define	M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif

#ifndef MAXFLOAT
#define  MAXFLOAT        ((float)3.40282346638528860e+38)  
#endif

#ifndef MINFLOAT
#define MINFLOAT         (-(float)3.40282346638528860e+38)
#endif 

// some commands for a cleaner compilation
#pragma warning( disable : 4251 4018 ) // 4251 STL doesn't compile cleanly on Visual C++ due to problem with DLL declarations
									   // 4018 signed/unsigned mismatch when looping through STL vectors with integer declarations instead of unsigned integer declarations


#endif // _DEFS_H_

