/* file:        Matrix.h
** author:      Matt Gong
** description: Declaration of class Matrix.
*/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <vector>

class Matrix
{
public:
  
  Matrix();
  Matrix( int M, int N );
  Matrix( double* data, int M, int N );
  virtual ~Matrix();
  
  double& 
  operator()(int row, int column);

  const double&
  operator()(int row, int column) const;

  Matrix
  operator*(const Matrix& A) const;
  
  Matrix
  operator*(double scalar) const;  

  Matrix
  operator/(double scalar) const;  

  Matrix&
  operator*=(double scalar);

  Matrix&
  operator/=(double scalar);

  Matrix
  operator-() const;

  Matrix
  operator+(const Matrix& A) const;

  Matrix
  operator-(const Matrix& A) const;
  
  const std::vector<double>&
  getData() const;
  
  double&
  getData(int index);

  void
  fillArray( double*& data ) const;
  
  void
  fill4x4Array( double matrix[4][4] );

  void
  fill4x4ArrayOpenGL( double matrix[4][4] );

  void
  initializeData( double*& data );

  Matrix
  transpose() const;
  
  Matrix
  diag() const;

  double
  determinant() const;
  
  double
  magnitude() const;

  Matrix 
  cross(const Matrix& p) const;

  double 
  dot(const Matrix& p) const;

public:

  int m_M;
  int m_N;

private:

  // Data stored in column-major order:
  //
  // i.e. 3x4 array:
  // row 1: m_data[0]  m_data[3]  m_data[6]  m_data[9]
  // row 2: m_data[1]  m_data[4]  m_data[7]  m_data[10]
  // row 3: m_data[2]  m_data[5]  m_data[8]  m_data[11]
  //
  std::vector<double> m_data;
  
};

std::ostream& operator<<(std::ostream& os, const Matrix& A);

#endif 
