/* file:        SparseMatrix.cpp
** author:      Matt Gong
** description: Sparse version of a matrix - memory and performance efficient.
*/

#include <string.h>
#include "SparseMatrix.h"

//-----------------------------------------------------
SparseMatrix::SparseMatrix()
: m_M (0),
  m_N (0)
{
}

//-----------------------------------------------------
SparseMatrix::SparseMatrix( int M, int N )
: m_M (M),
  m_N (N)
{
   m_rows.resize( M );
   m_columns.resize( N );
}

//-----------------------------------------------------
SparseMatrix::~SparseMatrix()
{
}

//-----------------------------------------------------
const Rows&
SparseMatrix::getRows() const
{
   return m_rows;
}

//-----------------------------------------------------
const Columns&
SparseMatrix::getColumns() const
{
   return m_columns;
}

//-----------------------------------------------------
double
SparseMatrix::getEntry( int row, int column ) const
{ 
   // Get the row
   const Row& row_ref = m_rows[row];

   // Attempt to find the entry
   RowConstIterator iter = row_ref.find( column );
   
   if ( iter != row_ref.end() ) 
   {
      return iter->second;
   }
   else 
   {
      return 0.0;
   }
}

//-----------------------------------------------------
void
SparseMatrix::setEntry( int row, int column, double value )
{
   // Get row
   Row& row_ref = m_rows[row];
 
   // Attempt to find the entry
   RowIterator row_iter = row_ref.find( column );
  
   // If the value is close to zero 
   if ( value < CLOSE_TO_ZERO && value > -CLOSE_TO_ZERO ) 
   {
      // Erase entry if the old value exists
      if ( row_iter != row_ref.end() )
      {
         row_ref.erase( row_iter );
         m_columns[column].erase( row );
      }
   }
   else
   {
      // Value is non-zero
   
      // If the previous value exists, update values in row/column
      if ( row_iter != row_ref.end() )
      {
         row_iter->second = value;
            
         ColumnIterator column_iter = m_columns[column].find( row );
         column_iter->second = value;
      }
      else
      {      
         // Previous value did not exist - insert new value
         row_ref.insert( RowPair(column, value) );
         m_columns[column].insert( ColumnPair(row, value) );
      }
   } 
}

//-----------------------------------------------------
void
SparseMatrix::incrementEntry( int row, int column, double increment_value )
{ 
   // Get row
   Row& row_ref = m_rows[row];
 
   // Attempt to find the entry
   RowIterator row_iter = row_ref.find( column );
   
   // Determine new value
   double new_value = 0.0;
   if ( row_iter != row_ref.end() ) {
      new_value = row_iter->second + increment_value;
   }
   else {
      new_value = increment_value;
   }
 
   // If the new value is close to zero
   if ( new_value < CLOSE_TO_ZERO && new_value > -CLOSE_TO_ZERO ) 
   {
      // Erase entry if the old value exists
      if ( row_iter != row_ref.end() )
      {
         row_ref.erase( row_iter );
         m_columns[column].erase( row );
      }
   }
   else
   {
      // New value is non-zero
   
      // If the previous value exists, update values in row/column
      if ( row_iter != row_ref.end() )
      {
         row_iter->second = new_value;
            
         ColumnIterator column_iter = m_columns[column].find( row );
         column_iter->second = new_value;
      }
      else
      {      
         // Previous value did not exist - insert new value
         row_ref.insert( RowPair(column, new_value) );
         m_columns[column].insert( ColumnPair(row, new_value) );
      }
   }  
}

//-----------------------------------------------------
SparseMatrix
SparseMatrix::operator*(const SparseMatrix& A) const
{
   SparseMatrix C( m_M, A.m_N );
  
   // Loop over all rows
   int index = 0;
   int index2 = 0;
   double value = 0.0;
   double value2 = 0.0;
   RowConstIterator row_iter;
   RowConstIterator row_iter2;
   for ( int r = 0; r < m_M; ++r )
   {
      // Get the row
      const Row& row = m_rows[r];
   
      // Loop over all entries in the row
      row_iter = row.begin();
      for ( ; row_iter != row.end(); ++row_iter )
      {
         index = row_iter->first;
         value = row_iter->second;
         
         // Get the row in the other matrix for this index
         const Row& other_row = A.m_rows[index];
         
         // Loop over all entries in other row
         row_iter2 = other_row.begin();
         for ( ; row_iter2 != other_row.end(); ++row_iter2 )
         {
            index2 = row_iter2->first;
            value2 = row_iter2->second;
            
            // Increment sum in our product C matrix
            C.incrementEntry( r, index2, value * value2 );
         }
      }
   } 
  
   return C;   
}

//-----------------------------------------------------
SparseMatrix
SparseMatrix::operator*(double scalar) const
{
   SparseMatrix C( *this );
   C *= scalar;

   return C;
}

//-----------------------------------------------------
SparseMatrix
SparseMatrix::operator/(double scalar) const
{
   SparseMatrix C( *this );
   C /= scalar;

   return C;
}

//-----------------------------------------------------
SparseMatrix&
SparseMatrix::operator*=(double scalar)
{
   // Multiply value in all row entries
   RowsIterator row_iter = m_rows.begin();
   for ( ; row_iter != m_rows.end(); ++row_iter )
   {
      Row& row_map = *row_iter;
      
      RowIterator entry_iter = row_map.begin();
      for ( ; entry_iter != row_map.end(); ++entry_iter )
      {
         entry_iter->second *= scalar;
      }
   }

   // Multiply value in all column entries
   ColumnsIterator column_iter = m_columns.begin();
   for ( ; column_iter != m_columns.end(); ++column_iter )
   {
      Column& column_map = *column_iter;
      
      ColumnIterator entry_iter = column_map.begin();
      for ( ; entry_iter != column_map.end(); ++entry_iter )
      {
         entry_iter->second *= scalar;
      }
   }
  
   return *this;
}

//-----------------------------------------------------
SparseMatrix&
SparseMatrix::operator/=(double scalar)
{
   // Divide value in all row entries
   RowsIterator row_iter = m_rows.begin();
   for ( ; row_iter != m_rows.end(); ++row_iter )
   {
      Row& row_map = *row_iter;
      
      RowIterator entry_iter = row_map.begin();
      for ( ; entry_iter != row_map.end(); ++entry_iter )
      {
         entry_iter->second /= scalar;
      }
   }

   // Divide value in all column entries
   ColumnsIterator column_iter = m_columns.begin();
   for ( ; column_iter != m_columns.end(); ++column_iter )
   {
      Column& column_map = *column_iter;
      
      ColumnIterator entry_iter = column_map.begin();
      for ( ; entry_iter != column_map.end(); ++entry_iter )
      {
         entry_iter->second /= scalar;
      }
   }
  
   return *this;
}

//-----------------------------------------------------
SparseMatrix
SparseMatrix::operator-() const
{
   SparseMatrix C( *this );
   C *= -1;

   return C;
}

//-----------------------------------------------------
SparseMatrix
SparseMatrix::operator+(const SparseMatrix& A) const
{
   SparseMatrix C( m_M, m_N );

   // Loop over all rows in C
   int index = 0;
   double value = 0.0;
   RowConstIterator entry_iter;
   RowConstIterator entry_iter2;
   for ( int r = 0; r < m_M; ++r )
   {
      // Get the row in this matrix
      const Row& row = m_rows[r];
      
      // Get the row in A
      const Row& other_row = A.m_rows[r];
      
      // Loop over all entries in this row
      entry_iter = row.begin();
      for ( ; entry_iter != row.end(); ++entry_iter )
      {
         index = entry_iter->first;
         value = entry_iter->second;
         
         // Attempt to find the entry in the other row
         entry_iter2 = other_row.find( index );
         
         // If it is there, add them together
         if ( entry_iter2 != other_row.end() ) {
            value += entry_iter2->second;
         }
      
         C.setEntry( r, index, value );
      } 

      // Loop over all entries in other row
      entry_iter = other_row.begin();
      for ( ; entry_iter != other_row.end(); ++entry_iter )
      {
         index = entry_iter->first;
         value = entry_iter->second;
         
         // Attempt to find the entry in this matrix
         entry_iter2 = row.find( index );
         
         // If it is there, skip because we added it in the other loop
         if ( entry_iter2 != row.end() ) {
            continue;
         }
      
         C.setEntry( r, index, value );
      }
   }

   return C;  
}

//-----------------------------------------------------
SparseMatrix
SparseMatrix::operator-(const SparseMatrix& A) const
{
   SparseMatrix C( m_M, m_N );

   // Loop over all rows in C
   int index = 0;
   double value = 0.0;
   RowConstIterator entry_iter;
   RowConstIterator entry_iter2;
   for ( int r = 0; r < m_M; ++r )
   {
      // Get the row in this matrix
      const Row& row = m_rows[r];
      
      // Get the row in A
      const Row& other_row = A.m_rows[r];
      
      // Loop over all entries in this row
      entry_iter = row.begin();
      for ( ; entry_iter != row.end(); ++entry_iter )
      {
         index = entry_iter->first;
         value = entry_iter->second;
         
         // Attempt to find the entry in the other row
         entry_iter2 = other_row.find( index );
         
         // If it is there, subtract them
         if ( entry_iter2 != other_row.end() ) {
            value -= entry_iter2->second;
         }
      
         C.setEntry( r, index, value );
      } 

      // Loop over all entries in other row
      entry_iter = other_row.begin();
      for ( ; entry_iter != other_row.end(); ++entry_iter )
      {
         index = entry_iter->first;
         value = entry_iter->second;
         
         // Attempt to find the entry in this matrix
         entry_iter2 = row.find( index );
         
         // If it is there, skip because we added it in the other loop
         if ( entry_iter2 != row.end() ) {
            continue;
         }
      
         C.setEntry( r, index, -value );
      }
   }

   return C;   
}

//-----------------------------------------------------
void
SparseMatrix::fillArray( double*& data ) const
{
   // Fills the given 1D array with the data
   // The data vector is assumed to be stored in column-major form.

   // First clear all data
   memset( data, 0, m_M * m_N * sizeof(double) );
  
   // Loop over all rows
   RowConstIterator entry_iter;
   for ( int r = 0; r < m_M; ++r )
   {   
      const Row& row = m_rows[r];
   
      // Loop over all row entries
      entry_iter = row.begin();
      for ( ; entry_iter != row.end(); ++entry_iter )
      {
         // Set the data
         data[r + entry_iter->first * m_M] = entry_iter->second;
      }
   }   
}

//-----------------------------------------------------
void
SparseMatrix::initializeData( double*& data )
{
   // Initializes this matrix with the given 1D array of data.
   // The data vector is assumed to be in column-major form.
  
   // Loop over all rows
   for ( int r = 0; r < m_M; ++r )
   {
      // Loop over all columns
      for ( int c = 0; c < m_N; ++c )
      {
         setEntry( r, c, data[r + c * m_M] );
      }
   } 
}

//-----------------------------------------------------
void
SparseMatrix::clear( int M, int N )
{
   // Clear the matrix and reset with the given dimensions
   
   m_M = M;
   m_N = N;
   m_rows.clear();
   m_columns.clear();
   m_rows.resize( M );
   m_columns.resize( N );
}

