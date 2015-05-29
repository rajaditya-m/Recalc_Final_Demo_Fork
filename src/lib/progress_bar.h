#pragma once

class ProgressBar
{
public:
  ProgressBar(int total_step);
  void Progress(int new_progress);
  ~ProgressBar() {}

private:
  int total_step_;
  int progress_;
  int max_char_num_;
  int total_step_digits_;
};

