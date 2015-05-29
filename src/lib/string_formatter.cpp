#include "string_formatter.h"
#include <cstdio>

namespace dj {

std::string GetStr(const char *fmt, ...) {
  const int kMaxBufferSize = 1024;
  char textString[kMaxBufferSize] = {'\0'};
  // -- Empty the buffer properly to ensure no leaks.
  //  memset(textString, '\0', sizeof(textString));
  va_list args;
  va_start(args, fmt);
  vsnprintf(textString, kMaxBufferSize, fmt, args);
  va_end (args);
  return std::string(textString);
}

const char* GetChar(const char *fmt, ...) {
  const int kMaxBufferSize = 1024;
  static char textString[kMaxBufferSize] = {'\0'};
  // -- Empty the buffer properly to ensure no leaks.
  //  memset(textString, '\0', sizeof(textString));
  va_list args;
  va_start(args, fmt);
  vsnprintf(textString, kMaxBufferSize, fmt, args);
  va_end (args);
  return &textString[0];
}

std::string Format(const char *fmt) {
  std::stringstream ss;
  while (*fmt) {
    if (*fmt == '%') {
      if (*(fmt + 1) == '%') {
        ++fmt;
      } else {
        std::cerr << CURRENT_LINE << " => invalid format string: \"" << fmt << "\"" << std::endl;
        std::abort();
      }
    }
    ss << *fmt++;
  }
  return ss.str();
}

}
