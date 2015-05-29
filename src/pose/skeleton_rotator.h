#ifndef SKELETON_ROTATOR_H
#define SKELETON_ROTATOR_H
#include "input_handler.h"


class SkeletonRotator : public InputHandler
{
public:
  virtual int HandleKeyPress(QKeyEvent *e);
  void Init();
  virtual int Idle();

private:
  SkeletonRotator();
  DECLARE_SINGLETON_CLASS(SkeletonRotator)
};

#endif // SKELETON_ROTATOR_H
