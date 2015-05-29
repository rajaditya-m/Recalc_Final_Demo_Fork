/* ----------------------------------------------------------------------- */
/*
  Easy embeddable cross-platform high resolution timer function. For each
  platform we select the high resolution timer. You can call the 'ns()'
  function in your file after embedding this.
*/
#ifndef TIMER__H_
#define TIMER__H_
#pragma once

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

class WorldTime {
public:
  static WorldTime* Instance() {
    static WorldTime timer;
    return &timer;
  }

  inline uint64_t NanoSecond() {
#if defined(__APPLE__)
    uint64_t now;
    now = mach_absolute_time();
    now *= info.numer;
    now /= info.denom;
    return now;
#elif defined(__linux)
    uint64_t now;
    struct timespec spec;
    clock_gettime(CLOCKID, &spec);
    now = spec.tv_sec * 1.0e9 + spec.tv_nsec;
    return now;
#elif defined(_WIN32)
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    return (uint64_t) ((1e9 * now.QuadPart)  / win_frequency.QuadPart);
#endif
  }

private:
  WorldTime() {
#if defined(__APPLE__)
    mach_timebase_info(&info);
#elif defined(__linux)
    clock_getres(CLOCKID, &linux_rate);
#elif defined(_WIN32)
    QueryPerformanceFrequency(&win_frequency);
#endif
  }

#if defined(__APPLE__)
  mach_timebase_info_data_t info;
#elif defined(__linux)
  struct timespec linux_rate;
#elif defined(_WIN32)
  LARGE_INTEGER win_frequency;
#endif
};


class Timer {
public:
  Timer() : time_(WorldTime::Instance()) {
    last_started_time_ = time_->NanoSecond();
  }

  void Start() {
    last_started_time_ = time_->NanoSecond();
  }

  float ElaspsedTime() {
    uint64_t current_time = time_->NanoSecond();
    return (current_time - last_started_time_) * 1e-9f;
  }
private:

  uint64_t last_started_time_;
  WorldTime* time_;
};

namespace unused {
#if 0
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
#endif


#if 0
//**************************************************************************************
// Copyright 2004 Huamin Wang.
//**************************************************************************************
// TIMER classes
//**************************************************************************************
#ifndef __TIMER_H__
#define __TIMER_H__
#include <iostream>
#include <sys/timeb.h>
#include <time.h>

class Timer {
public:
  struct timeb start_time;

  Timer() {Start();}

  ~Timer() {}

  void Start()
  {ftime( &start_time);}

  float Now() {
    struct timeb current_time;
    ftime(&current_time);
    return (float) current_time.time + 0.001f * current_time.millitm;
  }

  float GetTime() {
    struct timeb current_time;
    ftime( &current_time);
    return (float)(current_time.time - start_time.time) + 0.001f * (current_time.millitm - start_time.millitm);
  }
};

#ifdef __linux
bool operator<(const timespec& time1, const timespec& time2);
bool operator==(const timespec& time1, const timespec& time2);
bool operator>(const timespec& time1, const timespec& time2);
const timespec operator+(const timespec& time, const timespec& increment);
const timespec operator-(const timespec& end, const timespec& start);
float operator/(const timespec& dividend, const timespec& divisor);
unsigned long long GetNanoSeconds(const timespec& time);
const timespec GetCurrentTimeSpec(clockid_t clk_id = CLOCK_REALTIME);
std::ostream& operator<<(std::ostream& out, const timespec& time);
class AccurateTimer {
public:
  // CLOCK_REALTIME, a system-wide Realtime clock.
  // CLOCK_PROCESS_CPUTIME_ID, high-resolution timer provided by the CPU for each process.
  // CLOCK_THREAD_CPUTIME_ID, high-resolution timer provided by the CPU for each of the threads.
  AccurateTimer(clockid_t clk_id = CLOCK_REALTIME) : clk_id_ (clk_id) {
    Start(clk_id);
  }
  void Start(clockid_t clk_id = CLOCK_REALTIME) {
    clk_id_ = clk_id;
    clock_gettime(clk_id_, &start_time_);
  }

  float GetTime() {
    timespec end_time;
    clock_gettime(clk_id_, &end_time);
    end_time = end_time - start_time_;
    Start(clk_id_);
    return end_time.tv_sec + end_time.tv_nsec * 1e-9;
  }
private:
  timespec start_time_;
  clockid_t clk_id_;
};
#endif

#endif //__TIMER_H__
#endif
}


#endif
