#include "qt_object_selector.h"
#include <QMouseEvent>
#include "opengl_helper.h"
#include "selectable_object.h"
#include "print_macro.h"
#include "vector_lib.h"

template <class Float>
QtObjectSelector<Float>::QtObjectSelector(SelectableObject<Float> *obj, int mouse_button, int modifier, int priority)
  : InputHandler(priority)
  , obj_(obj)
  , mouse_button_(mouse_button)
  , modifier_(modifier)
{
}

template <class Float>
int QtObjectSelector<Float>::HandleMousePress(QMouseEvent *e) {
  QPoint pos = e->pos();
  if (e->button() == mouse_button_ && e->modifiers() == modifier_) {
    Float start[3], end[3];
    Float clicked_world_pos[3];
    GetPixelWorldPosition(pos.x(), pos.y(), clicked_world_pos);
    GetSelectionRay(pos.x(), pos.y(), start, end);
    obj_->SelectPoint(start, end, clicked_world_pos, selected_pos_);
    if (obj_->IsSelected()) {
      P(obj_->SelectedVertex(), selected_pos_[0], selected_pos_[1], selected_pos_[2]);
      depth_ = GetPointDepth(selected_pos_);
      return kPriority_;
    }
  }
  return -1;
}

template <class Float>
int QtObjectSelector<Float>::HandleMouseMove(QMouseEvent *e) {
  if (obj_->IsSelected()) {
    QPoint pos = e->pos();
    GetPixelWorldPosition(pos.x(), pos.y(), depth_, current_pos_);
    obj_->Move(current_pos_);
    return Priority();
  }
  return -1;
}

template <class Float>
int QtObjectSelector<Float>::Render() {
  if (obj_->IsSelected()) {
    obj_->GetVertexPosition(selected_pos_);
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
//    glPushMatrix();
//    Traslate(selected_pos_[0], selected_pos_[1], selected_pos_[2]);
//    glutWireSphere(.010f, 10, 10);
//    glPopMatrix();

    if (obj_->MouseMoved()) {
      glColor3f(0, 1, 0);
      DrawArrow(selected_pos_, current_pos_);
    }
    glPopAttrib();
  }
  return -1;
}

template <class Float>
int QtObjectSelector<Float>::HandleMouseRelease(QMouseEvent *e) {
  if (obj_->IsSelected()) {
    QPoint pos = e->pos();
    GetPixelWorldPosition(pos.x(), pos.y(), depth_, current_pos_);
    obj_->Unselect(current_pos_);
    return Priority();
  } else {
    return -1;
  }
}

template class QtObjectSelector<float>;
template class QtObjectSelector<double>;
