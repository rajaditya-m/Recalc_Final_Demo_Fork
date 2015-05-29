#ifndef RANDOMSUBSPACEPOSEVIEWER_H
#define RANDOMSUBSPACEPOSEVIEWER_H

#pragma once
#include "input_handler.h"
class RandomSubspacePoseViewer : public InputHandler
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
  virtual int Idle() { return kNotHandled; }
  virtual ~RandomSubspacePoseViewer();

private:
  RandomSubspacePoseViewer();
  DECLARE_SINGLETON_CLASS(RandomSubspacePoseViewer);
};

#endif // RANDOMSUBSPACEPOSEVIEWER_H
