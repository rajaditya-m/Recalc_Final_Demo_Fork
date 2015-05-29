#include "binary_file_io.h"
#include "print_macro.h"

BinaryFileReader::BinaryFileReader(const char *file_name)
{
  std::ifstream in(file_name, std::ios::binary | std::ios::in);
  ASSERT(in.is_open(), P(file_name));
  // Get file size
  in.seekg(0, in.end);
  int size = int(in.tellg());
  // Allocate buffer
  buf_.resize(size);
  ptr_ = &buf_[0];
  // Read file data to buffer
  in.seekg(0, in.beg);
  in.read(ptr_, size);
  in.close();
}
