#include "user_interaction.h"
#include "opengl_helper.h"
#include "global.h"
#include "tet.h"
#include "skeleton.h"
//#include <GL/glut.h>
//#include <GL/glu.h>
#include "skeleton_node.h"
//#include "simulator.h"
#include "open_gl_qt.h"
// cuda
#include <QPoint>
#include <QMouseEvent>

#if 1
namespace body {
using namespace body;
int moving_bone = -1;
double offset[3];
double moving_target[3];
bool selected = false;
double depth_value;
using namespace global;
int SelectJoint(const double* start, const double* end) {
  int selected = skeleton->SelectJoint(start, end);
  if (selected >= 0) {
    skeleton->SetNewRoot(selected);
    skeleton->InitializeBoneRotation();
#ifdef CUDA_ENABLED
    if (cuda::skeleton) {
      cuda::skeleton->Update();
    }
#endif
    global::current_body->ReinitializeOffsetToSkeleton();
//    arch::current_cloth->ReinitializeOffsetToSkeleton();
  }
  return selected;
}

void InverseKinematics(const int max_iteration = 500, double step_size = 2.5) {
  if (!skeleton || !selected) return;
  SkeletonNode* moving_joint = skeleton->GetChildJoint(moving_bone);
  dj::Vec3d tmp_joint_offset(moving_joint->offset());
  moving_joint->set_offset(offset);
  int moving_joint_idx = moving_joint - skeleton->GetNode(0);
  skeleton->InverseKinematics(moving_joint_idx, moving_target, max_iteration, step_size);
  moving_joint->set_offset(tmp_joint_offset());
  skeleton->UpdatePosition();
#ifdef CUDA_ENABLED
  if (cuda::skeleton) {
    cuda::skeleton->Update();
  }
#endif
}

bool HandleMousePress(QMouseEvent * event) {
  using namespace body;
  QPoint pos = event->pos();
  if (event->button() == Qt::LeftButton) {
    double start[3], end[3];
    GetSelectionRay(pos.x(), pos.y(), start, end);
    if (event->modifiers() == Qt::ShiftModifier) {
      // Set new root
      int joint = body::SelectJoint(start, end);
      return joint >= 0;
    } else if (event->modifiers() == Qt::ControlModifier) {
      // Select a body point
      //      moving_bone = current_body->SelectVertex(start, end, moving_target, offset);
      //      moving_bone = current_body->SelectVertex(start, end, moving_target, offset);
//      body::moving_bone = arch::current_body->SelectVertex(start, end, body::moving_target, body::offset);
      if (body::moving_bone >= 0) {
        body::selected = true;
        body::depth_value = GetPointDepth(moving_target);
        global::simulate = false;
//        simulator::SwitchLevel(0);
        return true;
      }
    }
  }
  return false;
}

bool HandleMouseMove(QMouseEvent * event) {
  using namespace body;
  QPoint pos = event->pos();
  if (selected) {
    GetPixelWorldPosition(pos.x(), pos.y(), depth_value, moving_target);
    skeleton->InitializeBoneRotation();
    arch::current_body->ReinitializeOffsetToSkeleton();
//    arch::current_cloth->ReinitializeOffsetToSkeleton();
    InverseKinematics(100);
    skeleton->AssembleSkeletonTransformation();
    arch::current_body->ApplySkeletonTransformationToVertex();
//    arch::current_body->StrainLimiting(0, 40);

//    arch::current_cloth->ApplySkeletonTransformation(true);
//    arch::current_obj_collider->SurfaceCollision(body_collision::kCollisionThreshold, time_step);
//    for (int i = 0; i < 1; ++i) {
//      arch::current_cloth->StrainLimiting(10, 10, 2);
//      arch::current_obj_collider->SurfaceCollision(body_collision::kCollisionThreshold, time_step);
//    }
//    arch::current_cloth->ClearVelocity();
    return true;
  }
  return false;
}

bool HandleMouseRelease(QMouseEvent * event) {
  using namespace body;
  Q_UNUSED(event);
  if (selected) {
    global::simulate = true;
  }
//  simulator::simulate = false;
  selected = false;
  moving_bone = -1;
  return false;
}


} // namespace body

namespace cloth {
using namespace global;
int selected_vertex = -1;
double depth_value;
bool HandleMousePress(QMouseEvent * event) {
  QPoint pos = event->pos();
  if (event->button() == Qt::LeftButton) {
    dj::Vec3f startf(0, 0, 0), endf(0, 0, 0);
    GetSelectionRay(pos.x(), pos.y(), startf(), endf());
    dj::Vec3f start(startf[0], startf[1], startf[2]);
    dj::Vec3f end(endf[0], endf[1], endf[2]);
    if (event->modifiers() == Qt::ShiftModifier) {
      return false;
//      int joint =  current_cloth->SelectVertex(start(), end());
//      selected_vertex = joint;
//      if (selected_vertex >= 0) {
//        depth_value = GetPointDepth(current_cloth->vertex_ + joint * 3);
//      }
//      return joint >= 0;
    }
  }
  return false;
}



}
#endif
