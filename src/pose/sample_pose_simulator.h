#ifndef SAMPLE_POSE_SIMULATOR_H
#define SAMPLE_POSE_SIMULATOR_H
#include "input_handler.h"

class SamplePoseSimulator : public InputHandler
{
public:
  void Init(int argc, char *argv[]);

  virtual int Idle();
  virtual int HandleKeyPress(QKeyEvent *e);
  virtual int Render();
  virtual ~SamplePoseSimulator();
private:
  SamplePoseSimulator();
  DECLARE_SINGLETON_CLASS(SamplePoseSimulator)
};

#endif // SAMPLE_POSE_SIMULATOR_H
