#ifndef QT_OBJECT_SELECTOR_H
#define QT_OBJECT_SELECTOR_H
#include "input_handler.h"

template <class Float> class SelectableObject;

template <class Float>
class QtObjectSelector : public InputHandler
{
public:
  QtObjectSelector(SelectableObject<Float>* obj,
                   int mouse_button = Qt::LeftButton,
                   int modifier = Qt::ControlModifier,
                   int priority = 10);

  virtual ~QtObjectSelector() {}
  virtual int HandleMousePress(QMouseEvent *e);
  virtual int HandleMouseMove(QMouseEvent *e);
  virtual int Render();
  virtual int HandleMouseRelease(QMouseEvent* e);

private:
  SelectableObject<Float>* obj_;
  Float selected_pos_[3];
  Float current_pos_[3];
  Float depth_;
  unsigned int mouse_button_;
  unsigned int modifier_;
};

//extern template class QtObjectSelector<float>;
//extern template class QtObjectSelector<double>;
#endif // QT_OBJECT_SELECTOR_H
