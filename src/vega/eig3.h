/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, copied from the public domain Java library JAMA. */

#ifndef _eig_h
#define _eig_h
#pragma once

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
template <int n>
void eigen_decomposition(double A[][n], double v[][n], double d[]);

#endif
