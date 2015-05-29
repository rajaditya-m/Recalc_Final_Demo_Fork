#include "progress_bar.h"
#include <cstdio>
#include <cassert>

ProgressBar::ProgressBar(int total_step)
  : total_step_(total_step)
  , progress_(1)
  , max_char_num_(50)
{
  assert(total_step_ > 0);
  total_step_digits_ = 0;
  for (int i = total_step_; i > 0; i /= 10, total_step_digits_++) {
  }
  printf("\n");
}

void ProgressBar::Progress(int new_progress)
{
  progress_ = new_progress + 1;
  float completed_ratio = progress_ * 1.0f / total_step_;
  int percentage = (int)(completed_ratio * 100);
  printf("\r[%3d%%][", percentage);
  int char_num = int(completed_ratio * max_char_num_);
  for (int i = 0; i < char_num; ++i) {
    printf("=");
  }
  for (int i = 0; i < max_char_num_ - char_num; ++i) {
    printf(" ");
  }
  printf("][%*d/%*d]", total_step_digits_, progress_, total_step_digits_, total_step_);
  if (progress_ >= total_step_) {
    printf("\n");
  }
  fflush(stdout);
}
