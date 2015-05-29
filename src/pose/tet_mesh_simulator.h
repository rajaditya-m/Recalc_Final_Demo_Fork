#pragma once
#include "input_handler.h"

class TetMeshSimulator : public InputHandler
{
public:
  virtual int Idle();
  virtual int Render();
  virtual int HandleKeyPress(QKeyEvent *e);

private:
  TetMeshSimulator();
  virtual ~TetMeshSimulator();
  DECLARE_SINGLETON_CLASS(TetMeshSimulator)
};
