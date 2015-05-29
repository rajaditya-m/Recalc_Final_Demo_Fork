#ifndef MIXEDSPARSEMATRIX_H
#define MIXEDSPARSEMATRIX_H
#include <Eigen/Dense>
#include "multi_domain_tet.h"
#include "sparseMatrix.h"
//#include "print_macro.h"

class MixedMultiDomainTet;
class PardisoSolver;
class MixedSparseMatrix : public SparseMatrix
{
public:
  typedef Eigen::Matrix<int, 4, 4> Mat4i;
  typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Mati;
  MixedSparseMatrix(MixedMultiDomainTet* multi_tet);
  virtual ~MixedSparseMatrix() {}
  void Solve(double* rhs, double* x);
  void UpdateMatrix(const double dt);
  void InitColumnIndex();
  inline void AddSubMatrix(int row_start, int col_idx, int row_num, int col_num, double* mat) {
    int offset = 0;
    for (int r = 0, row = row_start; r < row_num; ++r, ++row) {
      for (int c = 0, col = col_idx; c < col_num; ++c, ++col, ++offset) {
        columnEntries[row][col] += mat[offset];
      }
    }
  }

  inline void AddSubMatrix(int row_start, int col_idx, int row_num, int col_num, int mat_col_size, double* mat) {
    for (int r = 0, row = row_start; r < row_num; ++r, ++row) {
      for (int c = 0, col = col_idx; c < col_num; ++c, ++col) {
        columnEntries[row][col] += mat[r * mat_col_size + c];
      }
    }
  }

  inline void AddSubMatrixTransposed(int row_start, int col_idx, int row_num, int col_num, double* mat) {
    for (int r = 0, row = row_start; r < row_num; ++r, ++row) {
      for (int c = 0, col = col_idx, idx = r; c < col_num; ++col, ++c, idx += col_num) {
        columnEntries[row][col] += mat[idx];
      }
    }
  }

  inline void AddSubMatrixTransposed(int row_start, int col_idx, int row_num, int col_num, int mat_col_size, double* mat) {
    for (int r = 0, row = row_start; r < row_num; ++r, ++row) {
      for (int c = 0, col = col_idx; c < col_num; ++col, ++c) {
        columnEntries[row][col] += mat[c * mat_col_size + r];
      }
    }
  }

  inline void Add3x3Matrix(int row_start, int col_idx, double* mat) {
    int offset = 0;
    for (int r = 0, row = row_start; r < 3; ++r, ++row) {
      for (int c = 0, col = col_idx; c < 3; ++c, ++col, ++offset) {
        columnEntries[row][col] += mat[offset];
      }
    }
  }

  template <class EigenMat3>
  inline void Add3x3Matrix(int row_start, int col_idx, const EigenMat3& mat) {
    for (int r = 0, row = row_start; r < 3; ++r, ++row) {
      for (int c = 0, col = col_idx; c < 3; ++c, ++col) {
        columnEntries[row][col] += mat(r, c);
      }
    }
  }

  MixedMultiDomainTet* multi_tet_;
  PardisoSolver* pardiso_solver_;
  Mati domain_column_idx_;
  std::vector<int> diag_column_idx_;
  std::vector<Mat4i> full_tet_column_idx_;
  std::vector<Mat4i> interface_tet_column_idx_;
};

#endif // MIXEDSPARSEMATRIX_H
