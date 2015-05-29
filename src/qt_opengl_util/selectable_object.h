#ifndef SELECTABLE_OBJECT_H
#define SELECTABLE_OBJECT_H
#include "print_macro.h"

template <class Float>
class SelectableObject {
public:
  //  typedef double Float;
  SelectableObject(int null_vertex = -1) : kNullVertex(null_vertex) {
    Float pos[3] = {0, 0, 0};
    Unselect(pos);
    mouse_moved_ = false;
  }

  SelectableObject(const SelectableObject<Float>& other) : kNullVertex(other.kNullVertex) {
    (void) other;
  }
  SelectableObject<Float>& operator=(const SelectableObject<Float>& other) {
    (void) other;
    return *this;
  }

  virtual ~SelectableObject() {}
  inline int MouseMoved() { return mouse_moved_; }
  inline int SelectedVertex() { return ui_selected_vertex_; }
  inline int SelectPoint(Float* ray_start, Float* ray_end, Float* clicked_world_pos, Float* selected_pos) {
    ui_selected_vertex_ = Select(ray_start, ray_end, clicked_world_pos, selected_pos);
    return ui_selected_vertex_;
  }

  virtual int Select(Float* ray_start, Float* ray_end, Float* clicked_world_pos, Float* selected_pos) = 0;

  virtual bool GetVertexPosition(Float* pos) { (void) pos; return false; }
  virtual bool IsSelected() { return ui_selected_vertex_ != kNullVertex; }
  virtual bool Move(Float* to) {
    if (IsSelected()) {
      mouse_moved_ = true;
      last_mouse_pos_[0] = to[0];
      last_mouse_pos_[1] = to[1];
      last_mouse_pos_[2] = to[2];
    }
    return false;
  }

  virtual void Unselect(Float* to) {
    (void) to;
    ui_selected_vertex_ = kNullVertex;
    mouse_moved_ = false;
  }

protected:
  const int kNullVertex;
  int ui_selected_vertex_;
  bool mouse_moved_;
  Float last_mouse_pos_[3];
};

#endif // SELECTABLE_OBJECT_H
