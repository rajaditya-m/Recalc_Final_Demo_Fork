#ifndef SUBSPACEPOSESAMPLESIMULATOR_H
#define SUBSPACEPOSESAMPLESIMULATOR_H

#pragma once
#include "input_handler.h"
class SubspacePoseSampleSimulator : public InputHandler
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
  SubspacePoseSampleSimulator();
  virtual ~SubspacePoseSampleSimulator();
  DECLARE_SINGLETON_CLASS(SubspacePoseSampleSimulator);
};

#endif // SUBSPACEPOSESAMPLESIMULATOR_H
