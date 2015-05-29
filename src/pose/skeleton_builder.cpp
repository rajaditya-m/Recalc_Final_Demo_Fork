#include <QKeyEvent>
#include "skeleton_builder.h"
#include "opengl_helper.h"
#include "rainbow_color.h"
#include "vector_lib.h"

int SkeletonBuilder::Render()
{
  glPointSize(8);
  glBegin(GL_POINTS);
  glColor3fv(kYellow());
  glVertex3fv(pos_);
  glEnd();
  return -1;
}

int SkeletonBuilder::HandleKeyPress(QKeyEvent *e)
{
  const float kStepSize = 0.02f;
  switch (e->key()) {
    case Qt::Key_T: pos_[0] -= kStepSize; break;
    case Qt::Key_Y: pos_[0] += kStepSize; break;
    case Qt::Key_U: pos_[1] -= kStepSize; break;
    case Qt::Key_I: pos_[1] += kStepSize; break;
    case Qt::Key_O: pos_[2] -= kStepSize; break;
    case Qt::Key_P: pos_[2] += kStepSize; break;
    default: return -1;
  }
  dj::Vec3f pos(pos_);
  P(pos);
  return 0;
}

SkeletonBuilder::SkeletonBuilder()
  : InputHandler(20)
{
  pos_[0] = pos_[1] = pos_[2] = 0.0f;
}
