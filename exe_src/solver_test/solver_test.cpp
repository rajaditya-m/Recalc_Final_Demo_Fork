//**************************************************************************************
//  Copyright (C) 2002 - 2011, Huamin Wang
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
// main.cpp
//**************************************************************************************
#include "SOLVER_MATH.h"
#include "BLOCK_MATRIX_GRAPH.h"

int main(int argc, char *argv[])
{
  (void) argc;
  (void) argv;
  using namespace solver;
	BLOCK_MATRIX_GRAPH<double> graph;
	BLOCK_MATRIX_GRAPH<double> temp_graph;
	temp_graph.Get(graph);
	
	double x[128], x1[128], r[128], b[128];
		
	for(int i=0; i<128; i++)
	{
		b[i]=8;
		x[i]=0;
	}
	//The Cholesky solver will destroy the matrix
	//so we use a temp one to solve it
	temp_graph.Solve(b, x);	

	graph.Multiply(x, r);
	for(int i=0; i<13; i++)
		printf("solution %d; %f\n", i, r[i]);

	for(int i=0; i<128; i++)
	{
		b[i]=8;
		x1[i]=0;
	}
	graph.CG_Solve(b, x1);


	for(int i=0; i<13; i++)
		printf("x %d: %f, %f\n", i, x[i], x1[i]);

	getchar();


	return 0;
}
