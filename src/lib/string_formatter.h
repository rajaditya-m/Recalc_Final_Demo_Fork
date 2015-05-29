#ifndef STRING_FORMATTER_H
#define STRING_FORMATTER_H
#pragma once
#include <string>
#include <iostream>
#include <stdarg.h>
#include <sstream>
#include <string.h>
#include <type_traits>
#include "macro_constant.h"

namespace dj {

//std::string GetStr(const char *fmt, ...);
const char* GetChar(const char *fmt, ...);
std::string GetStr(const char *fmt, ...);
std::string Format(const char *fmt);

inline bool IsSpecifier(char c) {
  const char specifier[] = "iduoxXfFeEgGaAcspn";
  const int length = sizeof(specifier) / sizeof(specifier[0]);
  for (int j = 0; j < length; ++j) {
    if (specifier[j] == c) {
      return true;
    }
  }
  return false;
}

template <class Type, bool IsClassType>
struct FormatSelector {
  static std::string Print(const char *fmt, ...) {
    const int kMaxBufferSize = 1024;
    char textString[kMaxBufferSize] = {'\0'};
    va_list args;
    va_start(args, fmt);
    vsnprintf(textString, kMaxBufferSize, fmt, args);
    va_end (args);
    return std::string(textString);
  }
};

template <class Type>
struct FormatSelector<Type, true> {
  template <class T>
  static std::string Print(const char *fmt, T val) {
    (void) fmt; (void) val;
    ASSERT(false, std::cerr << ("something wrong with formatter parameter") << std::endl;);
    return "";
  }
};

template <typename T, typename... Args>
std::string Format(const char *fmt, T value, Args... args) {
  std::stringstream ss;
  while (*fmt) {
    if (*fmt == '%') {
      if (*(fmt + 1) == '\0') {
        break; // error
      } else if (*(fmt + 1) == '%') {
        ++fmt;
      } else {
        if (*(fmt + 1) == 'z') {
          ss << value;
          fmt += 2;
        } else {
          int i = 1;
          for (; fmt[i] != '\0' && !IsSpecifier(fmt[i]); ++i) {
          }
          if (fmt[i] == '\0') {
            break;
          }
          ++i;
          std::string format_str(fmt, i);
          ss << FormatSelector<T, std::is_class<T>::value >::Print(format_str.c_str(), value);
          fmt += i;
        }
        ss << Format(fmt, args...); // call even when *s == 0 to detect extra arguments
        return ss.str();
      }
    }
    ss << *fmt++;
  }
  std::cerr << CURRENT_LINE << " => invalid parameter. \"" << fmt << "\"" << std::endl;
  std::abort();
  return "";
}

template <typename T, typename... Args>
void Format2Char(char* buf, const char* fmt, Args... args) {
  std::string str = Format(fmt, args...);
  memcpy(buf, &str[0], sizeof(char) * str.size());
}

} // namespace string_formatter

#endif // STRING_FORMATTER_H

