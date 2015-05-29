
/* file:        Timer.h
** author:      Matt Gong
** description: Simple Timer
*/

#ifndef TIMER_H
#define TIMER_H

class Timer {
public:
  typedef unsigned long long LONGLONG;
   Timer();
   virtual ~Timer();

   void
   start();

   void
   restart();

   long
   getMillisSinceStart();

   LONGLONG
   getTicksSinceStart();

   LONGLONG
   getTicksPerSecond();

private:

   LONGLONG m_start_time;         // High frequency ticks
   LONGLONG m_frequency_hertz;    // Ticks per second

   LONGLONG 
   getNowInTicks();
};

#endif
