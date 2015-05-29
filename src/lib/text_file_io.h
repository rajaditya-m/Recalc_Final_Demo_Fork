#pragma once
#include <fstream>

class TextFileWriter
{
public:
  TextFileWriter(const char* file_name);

  template <class T>
  TextFileWriter& operator<<(T& v)
  {
    text_file_ << v;
    return *this;
  }

  ~TextFileWriter();

  void Close()
  {
    text_file_.close();
  }

private:
  std::ofstream text_file_;
};

class TextFileReader
{
public:
  TextFileReader(const char* file_name);

  template <class T>
  TextFileReader& operator>>(T& v)
  {
    text_file_ >> v;
    return *this;
  }

  ~TextFileReader();

  void Close();

private:
  std::ifstream text_file_;
};
