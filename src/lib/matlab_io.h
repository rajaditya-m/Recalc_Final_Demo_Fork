#pragma once
#include <fstream>
#include <functional>
#include "print_macro.h"

namespace dj {

#ifndef VECTOR_IO_FORMAT_CONSTANT
#define VECTOR_IO_FORMAT_CONSTANT
enum {
  kRow,
  kColumn
};
const char kSpace = '\t';
#endif

template <class T>
bool WriteVectorToMatlab(int size, T* vector, std::ostream& out, int format = kColumn)
{
  char space_char = (format == kColumn) ? '\n' : kSpace;
  for (int i = 0; i < size; i++) {
    out << vector[i] << space_char;
  }
  if (space_char == kSpace) {
    out << std::endl;
  }
  return true;
}

template <class T>
bool WriteVectorToMatlab(int size, T* vector, const char* file_name, int format = kColumn)
{
  std::ofstream output_file(file_name);
  ASSERT(output_file.is_open(), P(file_name));
  WriteVectorToMatlab<T>(size, vector, output_file, format);
  output_file.close();
  return true;
}

template <class EigenMatrix>
bool WriteEigenMatrixToMatlab(EigenMatrix& mat, std::ostream& out) {
  for (int r = 0; r < mat.rows(); ++r) {
    for (int c = 0; c < mat.cols(); ++c) {
      out << mat(r, c) << " ";
    }
    out << "\n";
  }
  return true;
}

template <class EigenMatrix>
bool WriteEigenMatrixToMatlab(EigenMatrix& mat, const char* file_name) {
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  WriteEigenMatrixToMatlab(mat, out);
  out.close();
  return true;
}

//------------------------------------------------------------------------------
// Read vector from file
//------------------------------------------------------------------------------
template <class T>
void ReadVectorFromMatlab(std::istream& in, std::vector<T>& vector, int expected_size = -1)
{
  while (!in.eof()) {
    T tmp;
    in >> tmp;
    if (in.eof()) break;
    vector.push_back(tmp);
  }
  ASSERT(expected_size == -1 || int(vector.size()) == expected_size);
}

template <class T>
inline void ReadVectorFromMatlab(const char* file_name, std::vector<T>& vector, int expected_size = -1)
{
  std::ifstream input_file(file_name);
  ASSERT(input_file.is_open(), P(file_name));
  ReadVectorFromMatlab<T>(input_file, vector, expected_size);
}


/*
   Matlab sparse matrix text file format :
   row_num col_num matrix_element_value
    .
    .
    .
*/
template <class T>
void ExportImplicitMatrixToMatlab(const char* file,
                                  std::function<void (T*, T*)> StiffnessMatrix,
                                  int size,
                                  T zero_threshold = T(1e-15))
{
  std::vector<T> vec(size, T(0));
  std::vector<T> col(size, T(0));
  std::ofstream out(file);
  if (!out.is_open()) {
    std::cerr << CURRENT_LINE << " => failed to open file " << file << std::endl;
    exit(0);
  }
  for (int c = 0; c < size; ++c) {
    vec[c] = 1;
    StiffnessMatrix(&vec[0], &col[0]);
    for (int r = 0; r < size; ++r) {
      if (col[r] > zero_threshold || col[r] < -zero_threshold) {
        // Matlab matrix is 1-indexed
        out << r + 1 << " " << c + 1 << " " << col[r] << "\n";
      }
    }
    vec[c] = 0;
  }
  out.close();
}

// Print implicit matrix
template <class T>
void WriteImplicitMatrixToMatlab(std::ostream& out, std::function<void (T*, T*)> StiffnessMatrix, int size) {
  using std::cout;
//  out <<  "[" ;
  for (int i = 0; i < size; ++i) {
    std::vector<T> vec(size, 0);
    vec[i] = 1;
    std::vector<T> row(size);
    StiffnessMatrix(&vec[0], &row[0]);
    for (int j = 0; j < size; ++j) {
      out << row[j] << " ";
    }
//    out << ";";
    out << "\n";
  }
//  out << "]\n" ;
}

template <class T>
void WriteImplicitMatrixToMatlab(const char* file_name, std::function<void (T*, T*)> StiffnessMatrix, int size) {
  std::ofstream out(file_name);
  if (!out.is_open()) {
    std::cerr << "WriteImplicitMatrixToMatlab() => failed to open file " << file_name << std::endl;
    exit(0);
  }
  for (int i = 0; i < size; ++i) {
    std::vector<T> vec(size, 0);
    vec[i] = 1;
    std::vector<T> row(size);
    StiffnessMatrix(&vec[0], &row[0]);
    for (int j = 0; j < size; ++j) {
      out << row[j] << " ";
    }
//    out << ";";
    out << "\n";
  }
  out.close();
}



} // namespace dj


