#include "rigid_body.h"
#include "affine_transformer.h"
#include "opengl_helper.h"
#include "rainbow_color.h"

RigidBody::RigidBody(TetrahedralMeshIO *mesh_io, const char *file_name, AffineTransformer<Real> *transformer)
  : Super(mesh_io, file_name, transformer)
{
  Construct();
}

RigidBody::RigidBody(TetrahedralMesh *mesh)
  : Super(mesh)
{
  Construct();
}

void RigidBody::Unselect(int vert, Real *pos)
{
  Force f;
  f.first = Vec3(pos[0] - vert_[vert][0], pos[1] - vert_[vert][1], pos[2] - vert_[vert][2]);
  f.first *= 30;
  f.second = vert;
  external_force_.push_back(f);
}

void RigidBody::Construct()
{
  gravity_ = Vec3(global::gravity);
  quaternion_ = Quaternion<Real>(1, 0, 0, 0);
  tet_density_ = std::vector<Real>(tet_num_, 10);
  translational_vel_ = Vec3(0, 0, 0);
  translational_damping_ = 0.02;
  rotational_damping_ = 0.02;
  rest_pos_ = vert_;
  ComputeTetrahedraVolume();
  // TODO uncomment this after making subspace mass diagonal
  if (0) {
    //  ComputeLumpedMass();
    L("using lumped mass");
  } else {
    mass_ = std::vector<Real>(v_num_, 1);
    L("Using uniform mass of 1");
  }
  total_mass_ = v_num_;
  ComputeCenterOfMass();
  CentralizeMassOfCenterAtOrigin();
  ComputeInertialTensor();
  rest_pos_.swap(vert_);
  ComputeCenterOfMass();
}

void RigidBody::Simulate(Real dt)
{
  //  Vec3 gravity_force(global::gravity);
  //  gravity_force *= total_mass_;
  Vec3 net_force = gravity_ * total_mass_;
  Vec3 net_torque(Real(0), Real(0), Real(0));
  ApplyGroundCollisionForce(dt, net_force, net_torque);
  ApplyExternalForce(net_force, net_torque);
  // Linear momentum
  acceleration_ = net_force / total_mass_;
  translational_vel_ += acceleration_ * dt;
  translational_vel_ *= (1 - translational_damping_);
  center_of_mass_ += translational_vel_ * dt;
  // Angular momentum
  Mat3 rotation_matrix;
  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  Mat3 cur_inv_inertial = rotation_matrix.transpose() * inv_inertia_tensor_ * rotation_matrix;
  net_torque *= dt;
  Vec3 change_in_angular_vel_;
  Eigen::Map<Vector3> domega_wraper(&change_in_angular_vel_[0]);
  Eigen::Map<Vector3> torque_wraper(&net_torque[0]);
  domega_wraper = cur_inv_inertial * torque_wraper;
  angular_vel_ += change_in_angular_vel_;
  angular_vel_ *= (1 - rotational_damping_);
  quaternion_ = quaternion_ + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[0], angular_vel_[1], angular_vel_[2]) * quaternion_;
  quaternion_.Normalize();
  UpdatePosition();
}

void RigidBody::ApplyExternalForce(Vec3& net_force, Vec3& net_torque)
{
  for (int i = 0; i < (int) external_force_.size(); ++i) {
    int v = external_force_[i].second;
    Vec3& force = external_force_[i].first;
    net_force += force;
    Vec3 r = vert_[v] - center_of_mass_;
    Vec3 tmp_torque;
    dj::Cross3(&r[0], &force[0], &tmp_torque[0]);
    net_torque += tmp_torque;
  }
  external_force_.clear();
}

void RigidBody::ApplyGroundCollisionForce(Real dt, TetrahedralMesh::Vec3 & net_force, TetrahedralMesh::Vec3 &net_torque)
{
  (void) dt;
  Real kFloor = 0.02;
  Real kStiffness = 900;
  Vec3 tmp_torque;
  for (int v = 0; v < v_num_; ++v) {
    if (vert_[v][1] < kFloor) {
      Vec3 force(Real(0), Real(0), Real(0));
      force[1] = (kFloor - vert_[v][1]) * mass_[v] * kStiffness;
      net_force += force;
      Vec3 r = vert_[v] - center_of_mass_;
      dj::Cross3(&r[0], &force[0], &tmp_torque[0]);
      net_torque += tmp_torque;
    }
  }
}

void RigidBody::UpdatePosition()
{
  Mat3 rotation_matrix;
  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  rotation_matrix.transposeInPlace();
  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
#if 1
    Vec3 offset;
    Eigen::Map<Vector3> offset_wrapper(&offset[0]);
    Eigen::Map<Vector3> rest_pos_wrapper(&rest_pos_[v][0]);
    offset_wrapper = rotation_matrix * rest_pos_wrapper;
    vert_[v] = offset + center_of_mass_;
#else
    vert_[v] = rest_pos_[v] + center_of_mass_;
#endif
  }
}

void RigidBody::ComputeInertialTensor()
{
  inertia_tensor_ = Mat3::Zero();
  for (int v = 0; v < v_num_; ++v) {
    Mat3 tmp = Mat3::Zero();
    Vec3 r = vert_[v] - center_of_mass_;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        tmp(i, j) = -r[i] * r[j];
      }
    }
    Real dot = r * r;
    tmp(0, 0) += dot;
    tmp(1, 1) += dot;
    tmp(2, 2) += dot;
    inertia_tensor_ += mass_[v] * tmp;
  }
  inv_inertia_tensor_ = inertia_tensor_.inverse();
}

void RigidBody::CentralizeMassOfCenterAtOrigin()
{
  AffineTransformer<Real>::Translate(&vert_[0][0], v_num_,  -center_of_mass_[0], -center_of_mass_[1], -center_of_mass_[2]);
  center_of_mass_ = Vec3(0, 0, 0);
}

void RigidBody::ComputeCenterOfMass()
{
  center_of_mass_ = Vec3(Real(0), Real(0), Real(0));
  for (int i = 0; i < v_num_; ++i) {
    center_of_mass_ += mass_[i] * vert_[i];
  }
  center_of_mass_ *= Real(1) / total_mass_;
}

void RigidBody::Render(int render_mode)
{
  //  glPushMatrix();
  //  glTranslated(center_of_mass_[0], center_of_mass_[1], center_of_mass_[1]);
  Super::Render(render_mode);
  glPointSize(6);
  glColor3fv(kRed());
  glBegin(GL_POINTS);
  Vertex3v(&center_of_mass_[0]);
  glEnd();
  //  glPopMatrix();
}
