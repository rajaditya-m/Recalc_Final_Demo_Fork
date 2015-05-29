#include <cmath>
#include "qt_gl_3d_scene_navigator.h"
#include <QPoint>
#include "scene.h"
#include "opengl_helper.h"
#include "open_gl_qt.h"
#include "print_macro.h"
#include "camera.h"
#include "global.h"

enum MouseState {
  kDefault,
  kShiftAndLeftMouseButton,
  kShiftAndRightMouseButton,
  kLeftMouseButton,
  kShift,
  kTranslate
};

int QtGL3DSceneNavigator::HandleMouseMove(QMouseEvent *event) {
  QPoint pos = event->pos();
  if (mouse_motion_mode_ == kLeftMouseButton ||
      mouse_motion_mode_ == kShiftAndLeftMouseButton ||
      mouse_motion_mode_ == kShiftAndRightMouseButton) {
    if (event->modifiers() != Qt::ShiftModifier && (event->buttons() & Qt::LeftButton)) {
      float angle_y = (float)(pos.x() - last_mouse_pos_.x()) * 360.0f / gl_widget_->width();
      float angle_x = (float)(pos.y() - last_mouse_pos_.y()) * 180.0f / gl_widget_->height();
      camera_->RotateX(-angle_x);
      camera_->RotateY(-angle_y);
    } else if (event->modifiers() == Qt::ShiftModifier && (event->buttons() & Qt::RightButton)) {
      float dx = (last_mouse_pos_.x() - pos.x()) * +0.003f;
      float dy = (last_mouse_pos_.y() - pos.y()) * -0.003f;
      camera_->MoveAlongX(dx);
      camera_->MoveAlongY(dy);
    } else if (event->modifiers() == Qt::ShiftModifier && (event->buttons() & Qt::LeftButton)) {
      float dx = (pos.x() - last_mouse_pos_.x());
      float dy = (pos.y() - last_mouse_pos_.y());
      if (std::abs(dx) > std::abs(dy)) {
        float angle_z = dx * 180.f / gl_widget_->height();
        camera_->RotateZ(angle_z);
      } else {
        //        camera_->Zoom(dy * 0.1);
      }
    }
  }

  last_mouse_pos_ = pos;
  return Priority();
}

#include "vector_lib.h"

int QtGL3DSceneNavigator::HandleMousePress(QMouseEvent* event) {
  QPoint pos = event->pos();
  bool handled = false;
  if (event->button() == Qt::LeftButton && event->modifiers() == Qt::ShiftModifier) {
    mouse_motion_mode_ = kShiftAndLeftMouseButton;
    handled = true;
  } else if (event->button() == Qt::LeftButton && event->modifiers() == Qt::NoModifier) {
    unsigned long long current_time = WorldTime::Instance()->NanoSecond();
    auto time_diff = current_time - last_click_time_;
    if (time_diff < 250000000) {
      float depth = GetPixelDepth(pos.x(), pos.y());
      if (depth < 1.0f - 1e-4f) {
        float world_pos[3];
        GetPixelWorldPosition(pos.x(), pos.y(), depth, world_pos);
        const float* focus = camera_->focus();
        const float* eye_pos = camera_->eye_pos();
        float new_eye_pos[3] = {
          world_pos[0] - focus[0] + eye_pos[0],
          world_pos[1] - focus[1] + eye_pos[1],
          world_pos[2] - focus[2] + eye_pos[2],
        };
        camera_->set_focus(world_pos);
        camera_->set_eye_pos(new_eye_pos);
        handled = true;
      }
    }
    if (!handled) {
      mouse_motion_mode_ = kLeftMouseButton;
    }
    last_click_time_ = current_time;
    handled = true;
  } else if (event->button() == Qt::RightButton && event->modifiers() == Qt::ShiftModifier) {
    mouse_motion_mode_ = kShiftAndRightMouseButton;
    handled = true;
  }
  //    else if (event->modifiers() & Qt::ControlModifier) {
  //      mouse_motion_mode_ = kTranslate;
  //    }
  last_mouse_pos_ = pos;
  return (handled) ? Priority() : -1;
}

int QtGL3DSceneNavigator::HandleMouseRelease(QMouseEvent *event) {
  QPoint pos = event->pos();
  mouse_motion_mode_ = kDefault;
  last_mouse_pos_ = pos;
  return -1;
}

int QtGL3DSceneNavigator::HandleKeyPress(QKeyEvent *event) {
  const float kTranslateStep = 0.01f;
  switch (event->key()) {
    case Qt::Key_Q:
      exit(0);
      break;
    case Qt::Key_P: {
      static int c = 0;
      std::stringstream in;
      in << c << ".jpeg";
      gl_widget_->ScreenShot(in.str().c_str(), "jpeg");
      c++;
      break;
    }
    case Qt::Key_W:
      camera_->Save();
      emit gl_widget_->updateStatusMessage("Scene file saved");
      break;
    case Qt::Key_A:
      if (event->modifiers() == Qt::ShiftModifier) {
        draw_axis_ = !draw_axis_;
      } else {
        camera_->Zoom(-camera_->zoom_step());
      }
      break;
    case Qt::Key_Z:
      camera_->Zoom(camera_->zoom_step());
      break;
    case Qt::Key_X:
      break;
    case Qt::Key_Down:
      camera_->MoveAlongY(-kTranslateStep);
      break;
    case Qt::Key_Up:
      camera_->MoveAlongY(+kTranslateStep);
      break;
    case Qt::Key_Left:
      camera_->MoveAlongX(-kTranslateStep);
      break;
    case Qt::Key_Right:
      camera_->MoveAlongX(+kTranslateStep);
      break;
    case Qt::Key_Shift:
      render_focus_ = true;
      break;
    case Qt::Key_Space:
      if (event->modifiers() & Qt::ShiftModifier) {
        global::pause_per_frame = !global::pause_per_frame;
      } else {
        global::simulate = !global::simulate;
      }
      break;
    default:
      break;
  }
  return -1;
}

int QtGL3DSceneNavigator::HandleKeyRelease(QKeyEvent *event) {
  Q_UNUSED(event);
  if (event->key() == Qt::Key_Shift) {
    render_focus_ = false;
  }
  return -1;
}

int QtGL3DSceneNavigator::Render() {
  if (mouse_motion_mode_ != kDefault || render_focus_) {
    glPushMatrix();
    glTranslatef(camera_->focus()[0], camera_->focus()[1], camera_->focus()[2]);
    glScalef(0.2f, 0.2f, 0.2f);
    //    DrawSphere(0.01, 10, 10);
    DrawAxis();
    glPopMatrix();
  }
  if (draw_axis_) {
    DrawAxis();
  }
  return -1;
}

void QtGL3DSceneNavigator::set_gl_widget(OpenGLQt *widget) {
  gl_widget_ = widget;
  camera_ = widget->camera();
}

void QtGL3DSceneNavigator::set_camera(dj::Camera* camera) {
  camera_ = camera;
}

QtGL3DSceneNavigator::QtGL3DSceneNavigator() : InputHandler(0), draw_axis_(false), render_focus_(false) {
  last_click_time_ = 0;
}
