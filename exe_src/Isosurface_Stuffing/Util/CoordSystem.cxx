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

#include <cstdio>
#include "CoordSystem.h"

void CoordSystem::print()
{
    printf("Origin: ") ;
    printVector(origin)	;
    printf("\tx: ") ;
    printVector(x) ;
    printf("\ty: ") ;
    printVector(y) ;
    printf("\tz: ") ;
    printVector(z) ;
    return ;
}

