#ifndef FPS_CALCULATOR_H_
#define FPS_CALCULATOR_H_
#include "timer.h"

class FPSCalculator {
public:
  FPSCalculator(void)
    : frame_count_(0)
    , total_frame_(0) {
  }

  float CalculateFPS() {
    total_frame_++;
    frame_count_++;
    //  Calculate time passed
    float time_interval = timer_.ElaspsedTime();
    if (time_interval > 1e-2f) {
      // Calculate the number of frames per second
      fps_ = frame_count_ / (time_interval);
      // Reset frame count
      frame_count_ = 0;
      timer_.Start();
    }
    return fps_;
  }

  float fps(void) { return fps_; }
  unsigned int total_frame(void) { return total_frame_; }

private:
  Timer timer_;
  int frame_count_;
  unsigned int total_frame_;
  float fps_;
};
#endif // FPSCALCULATOR_H_
