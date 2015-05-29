/* file:        inverse.h
** author:      Andrea Vedaldi
** desctiption: Declaration of function inverse.
**/

#ifndef _INVERSE_H_
 
/*#ifdef __cplusplus
extern "C" {
#endif*/

class Matrix;

//int inverse(const float *A, int N, float* B) ;
int inverse(const Matrix& A, int N, Matrix& B) ;

/*#ifdef __cplusplus
}
#endif*/

// _INVERSE_H_
#endif 

