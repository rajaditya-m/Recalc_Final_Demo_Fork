/* file:        Matrix.cpp
** author:      Matt Gong
** description: Implementation of matrix class
**/

#include <iomanip>
#include <vector>
#include "Matrix.h"
#include "math.h"

//-----------------------------------------------------
Matrix::Matrix()
: m_M (0),
  m_N (0)
{
}

//-----------------------------------------------------
Matrix::Matrix( int M, int N )
: m_M (M),
  m_N (N)
{
  m_data.resize( M * N, 0.0 );
}

//-----------------------------------------------------
Matrix::Matrix( double* data, int M, int N )
: m_M (M),
  m_N (N)
{
  m_data.resize( M * N );
  for ( int i = 0; i < M*N; ++i )
  {
    m_data[i] = data[i];
  }
}

//-----------------------------------------------------
Matrix::~Matrix()
{
}

//-----------------------------------------------------
double& 
Matrix::operator()(int row, int column)
{ 
  return m_data[row + column * m_M];
}

//-----------------------------------------------------
const double&
Matrix::operator()(int row, int column) const
{ 
  return m_data[row + column * m_M];
}

//-----------------------------------------------------
Matrix
Matrix::operator*(const Matrix& A) const
{
  Matrix C( m_M, A.m_N );
  int i,j,k ;
  for(j = 0 ; j < A.m_N ; ++j) {
    for(i = 0 ; i < m_M ; ++i) {
      C(i,j) = 0 ;
      for(k = 0 ; k < m_N ; ++k) {
        C(i,j) += (*this)(i,k)*A(k,j) ;
      }      
    }
  }
  return C;
}

//-----------------------------------------------------
Matrix
Matrix::operator*(double scalar) const
{
  Matrix C( m_M, m_N );
  C = *this;
  C *= scalar;
  
  return C;
}

//-----------------------------------------------------
Matrix
Matrix::operator/(double scalar) const
{
  Matrix C( m_M, m_N );
  C = *this;
  C /= scalar;
  
  return C;
}

//-----------------------------------------------------
Matrix&
Matrix::operator*=(double scalar)
{
  for ( int i = 0; i < (int) m_data.size(); ++i )
  {
    m_data[i] *= scalar;
  }
  
  return *this;
}

//-----------------------------------------------------
Matrix&
Matrix::operator/=(double scalar)
{
  for ( int i = 0; i < (int) m_data.size(); ++i )
  {
    m_data[i] /= scalar;
  }
  
  return *this;
}

//-----------------------------------------------------
Matrix
Matrix::operator-() const
{
  Matrix C( m_M, m_N );
  C = *this;
  C *= -1;
  
  return C;
}

//-----------------------------------------------------
Matrix
Matrix::operator+(const Matrix& A) const
{
  Matrix C( m_M, m_N );
  
  if ( m_M == A.m_M && m_N == A.m_N )
  {
    for ( int i = 0; i < (int) m_data.size(); ++i ) 
    {
      C.m_data[i] = m_data[i] + A.m_data[i];
    }
  }
  
  return C;
}

//-----------------------------------------------------
Matrix
Matrix::operator-(const Matrix& A) const
{
  Matrix C( m_M, m_N );
  
  if ( m_M == A.m_M && m_N == A.m_N )
  {
    for ( int i = 0; i < (int) m_data.size(); ++i ) 
    {
      C.m_data[i] = m_data[i] - A.m_data[i];
    }
  }
  
  return C;
}

//-----------------------------------------------------
const std::vector<double>&
Matrix::getData() const
{
  return m_data;
}

//-----------------------------------------------------
double&
Matrix::getData(int index)
{
  return m_data[index];
}

//-----------------------------------------------------
void
Matrix::fillArray( double*& data ) const
{
  // Fills the given 1D array with the data
  // Note: C++ stores its multidimensional arrays in row major form.
  // The data vector is stored in column-major form.
  
  for ( int i = 0; i < m_M*m_N; ++i )
  {
    data[i] = m_data[i];
  }  
}

//-----------------------------------------------------
void
Matrix::fill4x4Array( double matrix[4][4] )
{
  // Fills the given 4x4 matrix in row-major form
  
  for ( int i = 0; i < 4; ++i ) 
  {
    for ( int j = 0; j < 4; ++j )
    {
      matrix[i][j] = (*this)(i, j);
    }
  }
}

//-----------------------------------------------------
void
Matrix::fill4x4ArrayOpenGL( double matrix[4][4] )
{
  // Fills the given 4x4 matrix using openGL format in column-major form
  
  for ( int i = 0; i < 4; ++i ) 
  {
    for ( int j = 0; j < 4; ++j )
    {
      matrix[j][i] = (*this)(i, j);
    }
  }
}

//-----------------------------------------------------
void
Matrix::initializeData( double*& data )
{
  // Initializes this matrix with the given 1D array of data.
  // Note: C++ stores its multidimensional arrays in row major form.
  // The data vector is stored in column-major form.
  
  for ( int i = 0; i < m_M*m_N; ++i )
  {
    m_data[i] = data[i];
  }  
}
  
//-----------------------------------------------------
Matrix
Matrix::transpose() const 
{
  Matrix A( m_N, m_M );
  for(int j = 0 ; j < m_M ; ++j) 
    for(int i = 0 ; i < m_N ; ++i) 
      A(i,j) = (*this)(j,i) ;
  return A ;
}
  
//-----------------------------------------------------
Matrix
Matrix::diag() const 
{
  Matrix A( (m_N>m_M) ? m_N : m_M, (m_N>m_M) ? m_N : m_M );
  for(int j = 0 ; j < m_N ; ++j)
    for(int i = 0 ; i < m_N ; ++i)        
      A(i,j) = 0.0 ;
  for(int i = 0 ; i < m_N ; ++i)
    A(i,i) = m_data[i];
  return A ;
}

//-----------------------------------------------------
double
Matrix::determinant() const
{
  //Return the determinant of the upper 3 by 3 portion of the matrix.
         
  return ((*this)(0,0) * (*this)(1,1) * (*this)(2,2)) + 
         ((*this)(0,2) * (*this)(1,0) * (*this)(2,1)) +
         ((*this)(0,1) * (*this)(1,2) * (*this)(2,0)) -
         ((*this)(0,2) * (*this)(1,1) * (*this)(2,0)) -
         ((*this)(0,0) * (*this)(1,2) * (*this)(2,1)) -
         ((*this)(0,1) * (*this)(1,0) * (*this)(2,2));
}

//-----------------------------------------------------
double
Matrix::magnitude() const
{
  // If this is a vector, find it's magnitude
  
  return sqrt(
    ((*this)(0,0) * (*this)(0,0)) + 
    ((*this)(1,0) * (*this)(1,0)) + 
    ((*this)(2,0) * (*this)(2,0)));
}

//-----------------------------------------------------
Matrix 
Matrix::cross(const Matrix& p) const
{
  // If this is a vector, cross it with the given vector
  
  Matrix C(3,1);
  
  C(0,0) = (*this)(1,0) * p(2,0) - (*this)(2,0) * p(1,0);
  C(1,0) = (*this)(2,0) * p(0,0) - (*this)(0,0) * p(2,0);
  C(2,0) = (*this)(0,0) * p(1,0) - (*this)(1,0) * p(0,0);
    
  return C;
}

//-----------------------------------------------------
double 
Matrix::dot(const Matrix& p) const
{
  return (*this)(0,0) * p(0,0) + (*this)(1,0) * p(1,0) + (*this)(2,0) * p(2,0);
}

//-----------------------------------------------------
std::ostream& 
operator<<(std::ostream& os, const Matrix& A)
{
  os<<std::endl<<std::setprecision(6)<<std::fixed ;
  for(int i = 0 ; i < A.m_M  ; ++i) {
    for(int j = 0 ; j <  A.m_N  ; ++j) {
      os<<std::setw(8)<<A(i,j)<<" " ;
    }
    os<<std::endl ;
  }
  return os ;
}