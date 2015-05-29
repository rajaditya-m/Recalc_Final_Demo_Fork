#pragma once
#include "input_handler.h"

#include "single_domain_cubature.h"
#include "single_domain_basis_generator.h"

class SubspaceSimulator : public InputHandler
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
  SubspaceSimulator();
  virtual ~SubspaceSimulator();
  DECLARE_SINGLETON_CLASS(SubspaceSimulator)
};

