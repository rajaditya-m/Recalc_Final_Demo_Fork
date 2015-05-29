#include <string.h>
#include <vector>
#include "camera.h"
#include "config_file.h"
#include "macro_constant.h"
#include "vector_lib.h"

namespace dj {

Camera::Camera(Real* eye_pos, Real* focus, Real* up,
       Real field_of_view, Real near_plane, Real far_plane) {
  Construct(eye_pos, focus, up, field_of_view, near_plane, far_plane);
  camera_config_file_ = NULL;
}

Camera::Camera(const char *config_file) {
  strcpy(config_file_name_, config_file);
  camera_config_file_ = new ConfigFile(config_file);
  std::vector<Real>* eye_pos = camera_config_file_->Get<std::vector<Real>*>("eye pos");
  std::vector<Real>* focus = camera_config_file_->Get<std::vector<Real>*>("focus");
  std::vector<Real>* up = camera_config_file_->Get<std::vector<Real>*>("up vector");
  Real near_plane = camera_config_file_->Get<Real>("near plane");
  Real far_plane = camera_config_file_->Get<Real>("far plane");
  Real field_of_view = camera_config_file_->Get<Real>("field of view");
  Construct(&(*eye_pos)[0], &(*focus)[0], &(*up)[0], field_of_view, near_plane, far_plane);
}

Camera::Camera(const Camera &camera) {
  Construct(camera.eye_pos_, camera.focus_, camera.up_, camera.field_of_view_, camera.near_plane_, camera.far_plane_);
  camera_config_file_ = NULL;
}

Camera::~Camera() {
  if (camera_config_file_) delete camera_config_file_;
}

void Camera::LookAt(const Camera::Real *eye_pos, const Camera::Real *focus, const Camera::Real *up) {
  for (int i = 0; i < 3; ++i) {
    eye_pos_[i] = eye_pos[i];
    focus_[i] = focus[i];
    up_[i] = up[i];
  }
  UpdateLocalFrame();
}

bool Camera::Save() {
  if (config_file_name_ == NULL) return false;
  std::vector<Real>* eye_pos = camera_config_file_->Get<std::vector<Real>*>("eye pos");
  std::vector<Real>* focus = camera_config_file_->Get<std::vector<Real>*>("focus");
  std::vector<Real>* up = camera_config_file_->Get<std::vector<Real>*>("up vector");
  for (int i = 0; i < 3; ++i) {
    (*eye_pos)[i] = eye_pos_[i];
    (*up)[i] = up_[i];
    (*focus)[i] = focus_[i];
  }
  camera_config_file_->Save(config_file_name_);
  return true;
}

void Camera::Reset() {
  LookAt(initial_eye_pos_, initial_focus_, initial_up_);
}

void Camera::RotateX(Real angle_in_degree) {
  Rotate(x_, angle_in_degree);
}

void Camera::RotateY(Real angle_in_degree) {
  Rotate(y_, angle_in_degree);
}

void Camera::RotateZ(Real angle_in_degree) {
  Rotate(z_, angle_in_degree);
}

void Camera::Zoom(Real zoom) {
  Real camer_to_focus_distance = dj::Distance3(&eye_pos_[0], &focus_[0]);
  const Real kMinZoom = Real(0.1);
  Real new_distance = camer_to_focus_distance + zoom;
  if (new_distance < kMinZoom) new_distance = kMinZoom;
  for (int i = 0; i < 3; ++i) {
    eye_pos_[i] = focus_[i] + new_distance * z_[i];
  }
}

void Camera::MoveAlongX(float dx) {
  for (int i = 0; i < 3; ++i) {
    eye_pos_[i] += x_[i] * dx;
    focus_[i] += x_[i] * dx;
  }
}

void Camera::MoveAlongY(float dy) {
  for (int i = 0; i < 3; ++i) {
    eye_pos_[i] += y_[i] * dy;
    focus_[i] += y_[i] * dy;
  }
}

void Camera::MoveAlongZ(float dz) {
  for (int i = 0; i < 3; ++i) {
    eye_pos_[i] += z_[i] * dz;
    focus_[i] += z_[i] * dz;
  }
}

void Camera::Construct(const Real *eye_pos, const Real *focus, const Real *up,
                       Real field_of_view, Real near_plane, Real far_plane) {
  zoom_step_  = 0.1f;
  LookAt(eye_pos, focus, up);
  field_of_view_ = field_of_view;
  near_plane_ = near_plane;
  far_plane_ = far_plane;
  for (int i = 0; i < 3; ++i) {
    initial_eye_pos_[i] = eye_pos[i];
    initial_focus_[i] = focus[i];
    initial_up_[i] = up[i];
  }
}


void Camera::UpdateLocalFrame() {
  for (int i = 0; i < 3; ++i) {
    z_[i] = eye_pos_[i] - focus_[i];
    y_[i] = up_[i];
  }
  dj::Cross3(y_, z_, x_);
  dj::Normalize3(x_);
  dj::Normalize3(y_);
  dj::Normalize3(z_);
}

void Camera::Rotate(Camera::Real *rotation_axis, Camera::Real angle_in_degree) {
  Real normalized_axis[3] = {rotation_axis[0], rotation_axis[1], rotation_axis[2]};
  dj::Normalize3(normalized_axis);
  Real rotation_matrix[3][3];
  Real angle_in_radian = dj::Degree2Radian(angle_in_degree);
  Real cosine = cos(angle_in_radian);
  Real sine = sin(angle_in_radian);
  rotation_matrix[0][0] = cosine + normalized_axis[0] * normalized_axis[0] * (1 - cosine);
  rotation_matrix[0][1] = normalized_axis[0] * normalized_axis[1] * (1 - cosine) - normalized_axis[2] * sine;
  rotation_matrix[0][2] = normalized_axis[0] * normalized_axis[2] * (1 - cosine) + normalized_axis[1] * sine;

  rotation_matrix[1][0] = normalized_axis[0] * normalized_axis[1] * (1 - cosine) + normalized_axis[2] * sine;
  rotation_matrix[1][1] = cosine + normalized_axis[1] * normalized_axis[1] * (1 - cosine);
  rotation_matrix[1][2] = normalized_axis[1] * normalized_axis[2] * (1 - cosine) - normalized_axis[0] * sine;

  rotation_matrix[2][0] = normalized_axis[0] * normalized_axis[2] * (1 - cosine) - normalized_axis[1] * sine;
  rotation_matrix[2][1] = normalized_axis[1] * normalized_axis[2] * (1 - cosine) + normalized_axis[0] * sine;
  rotation_matrix[2][2] = cosine + normalized_axis[2] * normalized_axis[2] * (1 - cosine);


  Real tmp[3];
  Real focus2eye[3];
  dj::SubVec3(eye_pos_, focus_, focus2eye);
  dj::MulMatrix3x3Vec<Real>((Real (*)[3]) rotation_matrix, focus2eye, tmp);
  dj::AddVec3(focus_, tmp, eye_pos_);
  dj::MulMatrix3x3Vec<Real>((Real (*)[3]) rotation_matrix, up_, tmp);
  set_up(tmp);
  UpdateLocalFrame();
}

void Camera::set_up(Camera::Real *up) {
  up_[0] = up[0];
  up_[1] = up[1];
  up_[2] = up[2];
}

void Camera::set_eye_pos(Camera::Real *eye_pos) {
  eye_pos_[0] = eye_pos[0];
  eye_pos_[1] = eye_pos[1];
  eye_pos_[2] = eye_pos[2];
}

} // namespace dj

