#ifndef CUBICASIMULATOR_H
#define CUBICASIMULATOR_H
#pragma once
#include "input_handler.h"

class CubicaSimulator : public InputHandler
{
public:
  // Mouse handlers
  virtual int HandleMouseMove(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMousePress(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMouseRelease(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  // Keyboard handlers
  virtual int HandleKeyPress(QKeyEvent* e);
  virtual int HandleKeyRelease(QKeyEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int Render();
  virtual int Idle();

private:
  CubicaSimulator();
  virtual ~CubicaSimulator();
  DECLARE_SINGLETON_CLASS(CubicaSimulator);
};

#endif // CUBICASIMULATOR_H
