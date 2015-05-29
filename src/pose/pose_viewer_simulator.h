#pragma once
#include "input_handler.h"

class PoseViwerSimulator : public InputHandler
{
public:
  virtual int HandleKeyPress(QKeyEvent *e);
  virtual int Render();
private:
  PoseViwerSimulator();
  virtual ~PoseViwerSimulator();
  DECLARE_SINGLETON_CLASS(PoseViwerSimulator)
};
