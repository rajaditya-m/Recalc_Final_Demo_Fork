#ifndef MIXEDMULTIDOMAINSIMULATOR_H
#define MIXEDMULTIDOMAINSIMULATOR_H
#pragma once
#include "input_handler.h"
#include "singleton.h"

class MixedMultiDomainSimulator : public InputHandler
{
public:
  // Mouse handlers
  virtual int HandleMouseMove(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMousePress(QMouseEvent* e);
  virtual int HandleMouseRelease(QMouseEvent* e);
  // Keyboard handlers
  virtual int HandleKeyPress(QKeyEvent* e);
  virtual int HandleKeyRelease(QKeyEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int Render();
  virtual ~MixedMultiDomainSimulator();
  virtual int Idle();

private:
  int ui_mode_;
  double depth_;
  double clicked_pos_[3];
  MixedMultiDomainSimulator();
  DECLARE_SINGLETON_CLASS(MixedMultiDomainSimulator)
};

#endif // MIXEDMULTIDOMAINSIMULATOR_H
