#ifndef RIGID_BODY_H
#define RIGID_BODY_H
#include <Eigen/Dense>
#include "quaternion.h"
#include "simulation_tetrahedral_mesh.h"

template <class T = float>
struct EigenVectorChooser {
  typedef Eigen::Vector3f Vector3;
};
template <>
struct EigenVectorChooser<double> {
  typedef Eigen::Vector3d Vector3;
};

class TetrahedralMeshIO;
class RigidBody : public SimulationTetrahedralMesh
{
public:
  typedef SimulationTetrahedralMesh Super;
  typedef Eigen::Matrix<Real, 3, 3> Mat3;
  typedef std::pair<Vec3, int> Force;
  typedef EigenVectorChooser<Real>::Vector3 Vector3;

  RigidBody(TetrahedralMeshIO* mesh_io, const char* file_name, AffineTransformer<Real>* transformer = NULL);
  explicit RigidBody(TetrahedralMesh* mesh);

  virtual void Unselect(int vert, Real* pos);
  void Construct();
//  void ApplyGravity(Real dt);
  void Simulate(Real dt);
  void ApplyExternalForce(Vec3 &net_force, Vec3 &net_torque);
  void ApplyGroundCollisionForce(Real dt, Vec3 &net_force, Vec3 &net_torque);
  void UpdatePosition();
  void ComputeInertialTensor();
  void CentralizeMassOfCenterAtOrigin();
  void ComputeCenterOfMass();

  void Render(int render_mode);

  Quaternion<Real> quaternion_;
  std::vector<Vec3> rest_pos_;
  std::vector<Force> external_force_;
  Vec3 gravity_;
  Vec3 center_of_mass_;
  Mat3 inertia_tensor_;
  Mat3 inv_inertia_tensor_;
  Vec3 acceleration_;
  Vec3 translational_vel_;
  Vec3 angular_vel_;
  Real translational_damping_;
  Real rotational_damping_;
};

#endif // RIGID_BODY_H
