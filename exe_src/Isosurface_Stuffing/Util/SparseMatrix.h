/* file:        SparseMatrix.h
** author:      Matt Gong
** description: Sparse version of a Matrix - memory and performance efficient.
*/
#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <map>
#include <vector>

#define CLOSE_TO_ZERO 1.0E-10

// Typedefs
typedef std::map<int, double> Row;
typedef std::pair<int, double> RowPair;
typedef Row::iterator RowIterator;
typedef Row::const_iterator RowConstIterator;

typedef std::map<int, double> Column;
typedef std::pair<int, double> ColumnPair;
typedef Column::iterator ColumnIterator;
typedef Column::const_iterator ColumnConstIterator;

typedef std::vector<Row> Rows;
typedef Rows::iterator RowsIterator;
typedef Rows::const_iterator RowsConstIterator;

typedef std::vector<Column> Columns;
typedef Columns::iterator ColumnsIterator;
typedef Columns::const_iterator ColumnsConstIterator;

class SparseMatrix
{
public:
  
   SparseMatrix();
   SparseMatrix( int M, int N );
   virtual ~SparseMatrix();

   const Rows&
   getRows() const;
   const Columns&
   getColumns() const;

   double 
   getEntry(int row, int column) const;
   void
   setEntry(int row, int column, double value);
   void
   incrementEntry(int row, int column, double increment_value);

   SparseMatrix
   operator*(const SparseMatrix& A) const;

   SparseMatrix
   operator*(double scalar) const;  

   SparseMatrix
   operator/(double scalar) const;  

   SparseMatrix&
   operator*=(double scalar);

   SparseMatrix&
   operator/=(double scalar);

   SparseMatrix
   operator-() const;

   SparseMatrix
   operator+(const SparseMatrix& A) const;

   SparseMatrix
   operator-(const SparseMatrix& A) const;

   void
   fillArray( double*& data ) const;

   void
   initializeData( double*& data );
   
   void
   clear( int M, int N );

public:

   int m_M;
   int m_N;

private:

   // Data stored as a vector of rows and columns
   // Each vector stores a map from indices to value.
   // (e.g. 0->5.92, 4->-67.30, 26->45.17, etc.)
   Rows m_rows;
   Columns m_columns;
 
};

#endif 
