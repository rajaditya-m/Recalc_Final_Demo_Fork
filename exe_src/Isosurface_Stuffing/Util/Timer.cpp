
/* file:        Timer.cpp
** author:      Matt Gong
** description: Simple Portable Timer
*/

// System includes
#ifndef _WIN32
#include <sys/time.h>
#else
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#define _ATL_CSTRING_EXPLICIT_CONSTRUCTORS      // some CString constructors will be explicit
#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN            // Exclude rarely-used stuff from Windows headers
#endif
#include <windows.h>
#endif
#include <stdint.h>
#if defined(__linux)
#  define HAVE_POSIX_TIMER
#  include <time.h>
#  ifdef CLOCK_MONOTONIC
#     define CLOCKID CLOCK_MONOTONIC
#  else
#     define CLOCKID CLOCK_REALTIME
#  endif
#elif defined(__APPLE__)
#  define HAVE_MACH_TIMER
#  include <mach/mach_time.h>
#elif defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif
static uint64_t ns() {
  static uint64_t is_init = 0;
#if defined(__APPLE__)
  static mach_timebase_info_data_t info;
  if (0 == is_init) {
    mach_timebase_info(&info);
    is_init = 1;
  }
  uint64_t now;
  now = mach_absolute_time();
  now *= info.numer;
  now /= info.denom;
  return now;
#elif defined(__linux)
  static struct timespec linux_rate;
  if (0 == is_init) {
    clock_getres(CLOCKID, &linux_rate);
    is_init = 1;
  }
  uint64_t now;
  struct timespec spec;
  clock_gettime(CLOCKID, &spec);
  now = spec.tv_sec * 1.0e9 + spec.tv_nsec;
  return now;
#elif defined(_WIN32)
  static LARGE_INTEGER win_frequency;
  if (0 == is_init) {
    QueryPerformanceFrequency(&win_frequency);
    is_init = 1;
  }
  LARGE_INTEGER now;
  QueryPerformanceCounter(&now);
  return (uint64_t) ((1e9 * now.QuadPart)  / win_frequency.QuadPart);
#endif
}
#include "Timer.h"

//--------------------------------------------------------------------------
Timer::Timer()
  :  m_frequency_hertz(0),
     m_start_time(0)
{}

//--------------------------------------------------------------------------
Timer::~Timer()
{}

//--------------------------------------------------------------------------
void
Timer::start() {
  if (m_frequency_hertz == 0) {
#ifdef WIN32
    LARGE_INTEGER hertz_per_second;
    QueryPerformanceFrequency(&hertz_per_second);
    m_frequency_hertz = hertz_per_second.QuadPart;
#else
    // Solaris uses nanoseconds (1 billionth of a second)
    m_frequency_hertz = 1.0E9;
#endif
  }
  m_start_time = getNowInTicks();
}

//--------------------------------------------------------------------------
void
Timer::restart() {
  m_frequency_hertz = 0;
  m_start_time = 0;
  start();
}

//--------------------------------------------------------------------------
Timer::LONGLONG
Timer::getNowInTicks() {
  LONGLONG answer;
#ifdef WIN32
  LARGE_INTEGER temp;
  QueryPerformanceCounter(&temp);
  answer = (LONGLONG)temp.QuadPart;
#else
  answer = (LONGLONG) ns();
#endif
  return answer;
}

//--------------------------------------------------------------------------
long
Timer::getMillisSinceStart() {
  if (m_frequency_hertz == 0)		// force initialization if not done yet
    start();

  LONGLONG now_ticks = getNowInTicks();
  LONGLONG delta_ticks = now_ticks - m_start_time;

  double secs = ((double) (LONGLONG) delta_ticks)
                / ((double) (LONGLONG) m_frequency_hertz);
  double millis = secs * 1000.0;
  long answer = (long)millis;
  return answer;
}

//--------------------------------------------------------------------------
Timer::LONGLONG
Timer::getTicksSinceStart() {
  if (m_frequency_hertz == 0)		// force initialization if not done yet
    start();

  LONGLONG now_ticks = getNowInTicks();
  LONGLONG delta_ticks = now_ticks - m_start_time;

  return delta_ticks;
}

//--------------------------------------------------------------------------
Timer::LONGLONG
Timer::getTicksPerSecond() {
  return m_frequency_hertz;
}

