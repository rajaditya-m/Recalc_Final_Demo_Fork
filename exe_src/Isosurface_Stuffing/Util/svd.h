/* file:        svd.h
** author:      Andrea Vedaldi
** desctiption: Declaration of function svd.
**/

#ifndef _SVD_H_

/*#ifdef __cplusplus
extern "C" {
#endif*/

class Matrix;

int svd(const Matrix& A, int M, int N, Matrix& U, Matrix& S, Matrix& V) ;

/*#ifdef __cplusplus
}
#endif*/

// _SVD_H_
#endif
