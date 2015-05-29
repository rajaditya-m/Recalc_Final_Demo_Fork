#include "text_file_io.h"
#include "print_macro.h"

TextFileWriter::TextFileWriter(const char *file_name) : text_file_(file_name)
{
  ASSERT(text_file_.is_open(), P(file_name));
}

TextFileWriter::~TextFileWriter()
{
  text_file_.close();
}


TextFileReader::TextFileReader(const char *file_name) : text_file_(file_name)
{
  ASSERT(text_file_.is_open(), P(file_name));
}

TextFileReader::~TextFileReader()
{
  text_file_.close();
}

void TextFileReader::Close()
{
  text_file_.close();
}
