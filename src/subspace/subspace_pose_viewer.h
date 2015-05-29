#ifndef SUBSPACEPOSEVIEWER_H
#define SUBSPACEPOSEVIEWER_H
#pragma once
#include "input_handler.h"
class SubspacePoseViewer : public InputHandler
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
  virtual ~SubspacePoseViewer();

private:
  SubspacePoseViewer();
  DECLARE_SINGLETON_CLASS(SubspacePoseViewer);
};

#endif // SUBSPACEPOSEVIEWER_H
