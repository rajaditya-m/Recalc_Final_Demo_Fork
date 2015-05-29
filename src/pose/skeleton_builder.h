#ifndef SKELETON_BUILDER_H
#define SKELETON_BUILDER_H
#include "input_handler.h"

class SkeletonBuilder : public InputHandler
{
public:
  virtual int Render();
  virtual int HandleKeyPress(QKeyEvent *e);

private:
  SkeletonBuilder();
  DECLARE_SINGLETON_CLASS(SkeletonBuilder)
  float pos_[3];
};

#endif // SKELETON_BUILDER_H
