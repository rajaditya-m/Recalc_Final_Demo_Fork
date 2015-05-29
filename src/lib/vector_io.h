#pragma once
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <vector>
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
inline void ReadBuffer(T* dest, char*& buf, int count) {
  memcpy(dest, buf, sizeof(T) * count);
  buf += sizeof(T) * count;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// 1D vector

//------------------------------------------------------------------------------
// Text format
template <class T>
void Write1DVectorText(std::ostream& out, const T* vector, int size, int write_size = true, int format = kColumn) {
  char space_char = (format == kColumn) ? '\n' : kSpace;
  if (write_size) out << size << std::endl;
  for (int i = 0; i < size; i++) {
    out << vector[i] << space_char;
  }
  if (space_char == kSpace) {
    out << std::endl;
  }
}


template <class T>
inline void Write1DVectorText(const char* file_name,  const T* vector, int size, int write_size = true, int format = kColumn) {
  std::ofstream output_file(file_name);
  ASSERT(output_file.is_open(), P(file_name));
  Write1DVectorText(output_file, vector, size, write_size, format);
  output_file.close();
}

template <class T>
inline void Write1DVectorText(std::ostream& out, const std::vector<T>& vector, int write_size = true, int format = kColumn) {
  Write1DVectorText<T>(out, &vector[0], int(vector.size()), write_size, format);
}

template <class T>
inline void Write1DVectorText(const char* file_name,  const std::vector<T>& vector, int write_size = true, int format = kColumn) {
  Write1DVectorText(file_name, &vector[0], int(vector.size()), write_size, format);
}


//------------------------------------------------------------------------------
// Read
template <class T>
void Read1DVectorText(std::ifstream& in, T* vector, int expected_size, bool read_size = true) {
  int size;
  if (read_size) {
    in >> size;
    ASSERT(size == expected_size);
  } else {
    size = expected_size;
  }
  for (int i = 0; i < expected_size; ++i) {
    in >> vector[i];
  }
}

template <class T>
inline void Read1DVectorText(const char* file_name, T* vector, int expected_size, bool read_size = true) {
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  Read1DVectorText(in, vector, expected_size, read_size);
  in.close();
}

template <class T>
void Read1DVectorText(std::ifstream& in, std::vector<T>& vector, int expected_size = -1, bool read_size = true) {
  int size;
  if (read_size) {
    in >> size;
    ASSERT(expected_size == -1 || size == expected_size);
  } else {
    size = expected_size;
  }
  vector.resize(size);
  for (int i = 0; i < size; ++i) {
    in >> vector[i];
  }
}

template <class T>
void Read1DVectorText(const char* file_name, std::vector<T>& vector, int expected_size = -1, bool read_size = true) {
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  Read1DVectorText(in, vector, expected_size, read_size);
  in.close();
}


// Write1D binary
template <class T>
inline void Write1DVectorBinary(std::ostream& out, const T* vector, int size, bool write_size = true) {
  if (write_size) out.write((char*) &size, sizeof(int));
  out.write((char*) vector, sizeof(T) * size);
}

template <class T>
inline void Write1DVectorBinary(const char* file_name, const T* vector, int size, bool write_size = true) {
  std::ofstream out(file_name,std::ios::binary);
  ASSERT(out.is_open(), P(file_name));
  Write1DVectorBinary(out, vector, size, write_size);
  out.close();
}

template <class T>
inline void Write1DVectorBinary(std::ostream& out, const std::vector<T>& vector, bool write_size = true) {
  Write1DVectorBinary(out, &vector[0], int(vector.size()), write_size);
}

template <class T>
inline void Write1DVectorBinary(const char* file_name, const std::vector<T>& vector, bool write_size = true) {
  Write1DVectorBinary(file_name, &vector[0], int(vector.size()), write_size);
}

//------------------------------------------------------------------------------
// Read1D binary
template <class T>
inline void Read1DVectorBinary(std::istream& in, T* vector, int expected_size, bool read_size = true) {
  int size = -1;
  if (read_size) {
    in.read((char*) &size, sizeof(T));
    ASSERT(expected_size == size);
  } else {
    size = expected_size;
  }
  in.read((char*) &vector[0], sizeof(T) * size);
}

template <class T>
inline void Read1DVectorBinary(const char* file_name, T* vector, int expected_size, bool read_size = true) {
  std::ifstream in(file_name, std::ios::binary);
  ASSERT(in.is_open());
  Read1DVectorBinary(in, vector, expected_size, read_size);
  in.close();
}

template <class T>
inline void Read1DVectorBinary(std::istream& in, std::vector<T>& vector, int expected_size = -1, bool read_size = true) {
  int size = -1;
  if (read_size) {
    in.read((char*) &size, sizeof(T));
    ASSERT(expected_size == -1 || expected_size == size);
  } else {
    size = expected_size;
  }
  vector.resize(size);
  in.read((char*) &vector[0], sizeof(T) * size);
}

template <class T>
inline void Read1DVectorBinary(const char* file_name, std::vector<T>& vector, int expected_size = -1, bool read_size = true) {
  std::ifstream in(file_name, std::ios::binary);
  ASSERT(in.is_open());
  Read1DVectorBinary(in, vector, expected_size, read_size);
  in.close();
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// 2D Matrix
template <class T>
inline void Write2DVector(std::ostream& out, std::vector<std::vector<T> >& vector) {
  int size1 = (int) vector.size();
  out.write((char*) &size1, sizeof(int));
  for (int i = 0; i < int(vector.size()); ++i) {
    int size2 = (int) vector[i].size();
    out.write((char*) &size2, sizeof(int));
    if (size2 > 0) {
      out.write((char*) &vector[i][0], sizeof(T) * size2);
    }
  }
}

template <class T>
inline void Read2DVector(std::istream& in, std::vector<std::vector<T> >& vector) {
  int size1 = -1;
  in.read((char*) &size1, sizeof(int));
  vector.resize(size1);
  for (int i = 0; i < int(vector.size()); ++i) {
    int size2 = -1;
    in.read((char*) &size2, sizeof(int));
    vector[i].resize(size2);
    if (size2 > 0) {
      in.read((char*) &vector[i][0], sizeof(T) * size2);
    }
  }
}


template <class T>
inline bool WriteMatrix(int nx, int ny, T* matrix, std::ostream& out) {
  out << "[";
  for (int row = 0, offset = 0; row < nx; ++row) {
    for (int col = 0; col < ny; ++col, ++offset) {
      out << matrix[offset] << " ";
    }
    out << "; ";
  }
  out << "]" << std::endl;
  return true;
} //#WriteVector#



} // namespace dj
