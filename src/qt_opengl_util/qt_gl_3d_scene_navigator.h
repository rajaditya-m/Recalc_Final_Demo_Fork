#ifndef QT_GL_3D_INPUT_HANLDER_H
#define QT_GL_3D_INPUT_HANLDER_H
#include "input_handler.h"
#include "singleton.h"

class OpenGLQt;
namespace dj {
class Camera;
}
class QtGL3DSceneNavigator : public InputHandler {
public:
  // Mouse handlers
  virtual int HandleMouseMove(QMouseEvent* event);
  virtual int HandleMousePress(QMouseEvent *event);
  virtual int HandleMouseRelease(QMouseEvent* event);

  virtual int HandleKeyPress(QKeyEvent* event);
  virtual int HandleKeyRelease(QKeyEvent* event);

  virtual int Render();
  void set_gl_widget(OpenGLQt* widget);
  void set_camera(dj::Camera* camera);
private:
  bool draw_axis_;
  OpenGLQt* gl_widget_;
//  dj::Scene* scene_;
  dj::Camera* camera_;
  bool render_focus_;
  unsigned long long last_click_time_;
  QPoint last_mouse_pos_;
  QtGL3DSceneNavigator();
  int mouse_motion_mode_;// = kDefault;
  DECLARE_SINGLETON_CLASS(QtGL3DSceneNavigator)
};

#endif // QT_GL_3D_INPUT_HANLDER_H
