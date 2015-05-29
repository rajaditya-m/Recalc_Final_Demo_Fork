#ifndef MATRIXOPS_H_
#define MATRIXOPS_H_

#include <iostream>
#include <cstdlib>
#include <Eigen\Dense>
#include <Eigen\Sparse>

//MACROS FOR IDX
#define COLMAJORIDX(r,c,nrows) ((c)*(nrows)+(r)) 

//Your run of the mill dot product 
// Nothing special just move along 
double vecDotPdk(double *A, double *B, int len); 

//Your run of the mill scalar product
// Evaluates A = alpha * A
void vecScalarPdk(double *A, double alpha, int len);

//Again your run of the mill subtraction A = A - B
// Nothing to see here move along 
// Remember the result is stored back in A
void vecMinus(double *A, double *B, int len);

//Simple standard euclidean norm 
double vec2Norm(double *A, int len);

// Divide by Scalar operation 
void vecScalarDivide(double* A, double alpha, int len);



//The input to this algorithm is a MxN matrix in column major order with column representing the vectors that need to be orthogonalized 
// The output to this is two matrices 
// A matrix called Q which represents the orthogonal basis of dimensions MxN 
// A matrix called R which represents the Product (Q' * Mat)
// Assume that Q is already allocated with sufficient space 
// And so is R 
void QRDecomposition(int rows, int cols, double *U, double *Q, double *R);

void dimensionalityReductionByPCASVD(Eigen::MatrixXd dataPoints, int targetDims,Eigen::MatrixXd &fullToSub,Eigen::MatrixXd &mean);

void massPCA(Eigen::MatrixXd dataPoints, Eigen::SparseMatrix<double> massMatrix, int targetDims,Eigen::MatrixXd &fullToSub);

void massOrthogonalization(Eigen::MatrixXd dataPoints, Eigen::SparseMatrix<double> massMatrix, Eigen::MatrixXd &massOrtho);

void computeRotationInvariantPCA(Eigen::MatrixXd &dataPoints, int targetDims, Eigen::MatrixXd &fullToSub);

#endif





