#pragma once
#include <iostream>
#include <fstream>
#include <vector>

class BinaryFileReader
{
public:
  BinaryFileReader(const char* file_name);

  template <class T>
  inline void Read(T* data, int element_count)
  {
    if (!Valid()) {
      std::cerr << "BinaryFileReader::Read() => no more data to read from file" << std::endl;
      exit(0);
    }
    int size_in_bytes = sizeof(T) * element_count;
    memcpy(data, ptr_, size_in_bytes);
    ptr_ += size_in_bytes;
  }


  template <class T>
  inline void Read2DVector(std::vector<std::vector<T> >& vector)
  {
    int size1 = *((int*) ptr_);
    ptr_ += sizeof(int);
    vector.resize(size1);
    for (int i = 0; i < int(vector.size()); ++i) {
      int size2 = *((int*) ptr_);
      ptr_ += sizeof(int);
      vector[i].resize(size2);
      Read<T>(&vector[i][0], size2);
    }
  }

  template <class T>
  inline void Read1DVector(std::vector<T>& vector)
  {
    int size = *((int*) ptr_);
    ptr_ += sizeof(int);
    vector.resize(size);
    Read<T>(&vector[0], size);
  }

  bool Valid()
  {
    return (ptr_ - &buf_[0]) < (int) buf_.size();
  }

private:
  char* ptr_;
  std::vector<char> buf_;
};
