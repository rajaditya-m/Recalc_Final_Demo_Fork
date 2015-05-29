#ifndef MULTI_DOMAIN_SIMULATOR_H
#define MULTI_DOMAIN_SIMULATOR_H
#pragma once
#include "input_handler.h"

class MultiDomainSimulator : public InputHandler
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
  MultiDomainSimulator();
  virtual ~MultiDomainSimulator();
  DECLARE_SINGLETON_CLASS(MultiDomainSimulator);
};

#endif // MULTI_DOMAIN_SIMULATOR_H
