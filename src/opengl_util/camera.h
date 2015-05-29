#ifndef CAMERA_H
#define CAMERA_H
#pragma once

class ConfigFile;
namespace dj {

class Camera {
public:
  typedef float Real;
  Camera(Real* eye_pos, Real* focus, Real* up,
         Real field_of_view = Real(30),
         Real near_plane = Real(0.2),
         Real far_plane = Real(10000));
  Camera(const char* config_file);
  Camera(const Camera& camera);
  ~Camera();
  void LookAt(const Real *eye_pos, const Real *focus, const Real *up);

  // getter and setter methods
  const Real *eye_pos() const {
    return &eye_pos_[0];
  }
  const Real* up() const  {
    return &up_[0];
  }
  const Real* focus() const  {
    return &focus_[0];
  }
  Real zoom_step() const
  {
    return zoom_step_;
  }
  void set_zoom_step(Real zoom_step) { zoom_step_ = zoom_step; }
  void focus(Real* focus) const {
    focus[0] = focus_[0];
    focus[1] = focus_[1];
    focus[2] = focus_[2];
  }
  void set_focus(Real* focus) {
    focus_[0] = focus[0];
    focus_[1] = focus[1];
    focus_[2] = focus[2];
    UpdateLocalFrame();
  }
  void set_focus(Real x, Real y, Real z) {
    focus_[0] = x;
    focus_[1] = y;
    focus_[2] = z;
    UpdateLocalFrame();
  }
  void set_near_plane(Real near_plane) { near_plane_ = near_plane; }
  void set_far_plane(Real far_plane) { far_plane_ = far_plane; }
  void set_field_of_view(Real field_of_view) { field_of_view_ = field_of_view; }
  Real near_plane() { return near_plane_; }
  Real far_plane() { return far_plane_; }
  Real field_of_view() { return field_of_view_; }

  // Roll, pitch, and yaw in camera's local frame
  // camerate a origin with +y axis as the up vector
  void RotateX(Real angle_in_degree);
  void RotateY(Real angle_in_degree);
  void RotateZ(Real angle_in_degree);

  // Move camera position with wrt. camera's local frame
  void MoveAlongX(Real dx);
  void MoveAlongY(Real dy);
  void MoveAlongZ(Real dz);
  void Zoom(Real zoom);

  bool Save();
  void Reset();
  void set_eye_pos(Real* eye_pos);

private:
  void Construct(const Real *eye_pos, const Real *focus, const Real *up,
                 Real field_of_view, Real near_plane, Real far_plane);
  void UpdateLocalFrame();
  void Rotate(Real* rotation_axis, Real angle_in_degree);
  void set_up(Real* up);

  Real eye_pos_[3];
  Real focus_[3];
  Real up_[3];

  Real field_of_view_;
  Real near_plane_;
  Real far_plane_;

  Real initial_eye_pos_[3];
  Real initial_focus_[3];
  Real initial_up_[3];

  // camera's local frame vector expressed in world space
  Real x_[3];
  Real y_[3];
  Real z_[3];

  Real zoom_step_;
  // config file
  char config_file_name_[1024];
  ConfigFile* camera_config_file_;
};

} // namespace dj

#endif // CAMERA_H
