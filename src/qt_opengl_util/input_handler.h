#ifndef INPUT_HANDLER_H
#define INPUT_HANDLER_H

#include <QMouseEvent>
#include <QKeyEvent>
#include "singleton.h"
//#pragma once

class InputHandler
{
public:
  enum {
    kNotHandled = -1
  };
  static unsigned int mode;
  // Mouse handlers
  virtual int HandleMouseMove(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMousePress(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMouseRelease(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  // Keyboard handlers
  virtual int HandleKeyPress(QKeyEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleKeyRelease(QKeyEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int Render() { return kNotHandled; }
  virtual int Idle() { return kNotHandled; }

  bool CanHandle(int priority) const { return kPriority_ >= priority && (kType_ & mode); }
  unsigned int Type() const { return kType_; }
  int Priority() const { return kPriority_; }
  bool operator<(const InputHandler& other_handler) const { return this->kPriority_ < other_handler.kPriority_; }
  bool operator<=(const InputHandler& other_handler) const { return this->kPriority_ <= other_handler.kPriority_; }
  bool operator>(const InputHandler& other_handler) const { return this->kPriority_ > other_handler.kPriority_; }
  bool operator>=(const InputHandler& other_handler) const { return this->kPriority_ >= other_handler.kPriority_; }
protected:
  const int kPriority_;
  const unsigned int kType_;
  InputHandler(int priority = 0, int type = 1) : kPriority_(priority), kType_(type) {}
  virtual ~InputHandler() {}
  // Put this line in derived class
  // DECLARE_SINGLETON_CLASS(InputHandler)
};

#endif // INPUT_HANDLER_H


