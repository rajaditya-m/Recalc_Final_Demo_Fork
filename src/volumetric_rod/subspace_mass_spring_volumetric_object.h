#ifndef SUPSPACE_MASS_SPRING_VOLUMETRIC_OBJECT_H
#define SUPSPACE_MASS_SPRING_VOLUMETRIC_OBJECT_H
//#include "tetrahedral_mesh.h"
#include "rigid_body.h"
#include "global.h"
#include <Eigen/Dense>

class TetCollider;
class TetrahedralMeshIO;
template <class Float> class ConjugateGradientSolver;

template <class T = float>
struct EigenTypeSelector {
  typedef Eigen::MatrixXf MyMatrixX;
  typedef Eigen::VectorXf MyVectorX;
  typedef Eigen::Matrix3f MyMatrix3;
  typedef Eigen::Vector3f MyVector3;
};

template <>
struct EigenTypeSelector<double> {
  typedef Eigen::MatrixXd MyMatrixX;
  typedef Eigen::VectorXd MyVectorX;
  typedef Eigen::Matrix3d MyMatrix3;
  typedef Eigen::Vector3d MyVector3;
};

template <class T> class AffineTransformer;

class ConfigFile;
class SubspaceMassSpringVolumetricObject : public RigidBody
{
public:
  static const int kCrossSectionVerextNum = 400;
  static const int kSegmentCrossSectionNum = 8;
  //  static const int kCrossSectionVerextNum = 120;
  //  static const int kSegmentCrossSectionNum = 2;
  static const int kVertexNumPerSegment = kCrossSectionVerextNum * kSegmentCrossSectionNum;
  typedef RigidBody Super;
  typedef EigenTypeSelector<Real>::MyMatrixX MatrixX;
  typedef EigenTypeSelector<Real>::MyVectorX VectorX;
  typedef EigenTypeSelector<Real>::MyVector3 Vector3;
  typedef EigenTypeSelector<Real>::MyMatrix3 Matrix3;
  typedef Eigen::Map<MatrixX> MapMatrixX;
  typedef Eigen::Map<VectorX> MapVectorX;
  typedef Eigen::Map<Matrix3> MapMatrixX3;
  typedef Eigen::Map<Vector3> MapVector3;

  SubspaceMassSpringVolumetricObject(TetrahedralMeshIO* mesh_io, const char* file_name, AffineTransformer<Real>* transformer = NULL);
  //  SubspaceMassSpringVolumetricObject(TetrahedralMeshIO *mesh_io,
  //                                     std::vector<const char*>  file_names,
  //                                     std::vector<AffineTransformer<Real>*> transformers);
  explicit SubspaceMassSpringVolumetricObject(TetrahedralMesh* mesh);

  void LoadSamplePose(const char* file_name, std::vector<VectorX> &poses);
  void ViewSamlePose(const char* pos_file);
  void GenerateSamplePose(const char* pose_file, int pos_num);
  void Construct();
  void LoadMassMatrix(const char* mass_file);

  virtual ~SubspaceMassSpringVolumetricObject();
  bool VerifySubspaceDiagonality(SubspaceMassSpringVolumetricObject::MatrixX& subspace);
  void LoadSubspace(const char *file_name);
  void BuildMomentumMatrix();

  void SubspaceSimulation(Real dt);
  void ImplicitStepOneSegment(Real dt);
  void ComputeForceGradient(std::vector<Real>& force_gradient);
  void ComputeForce(std::vector<Vec3>& tmp_vert_buffer_, Real dt);
  void SimulateInLocalFrame(Real dt);
  void Simulate(Real dt);
  void ReconstructObjectFromSubspace();
  void ConstructFullSubspace();
  void ComputeRestEdgeLength();
  void InitializeSubspaceData();
  void UpdatePosition();
  void Render(int render_mode);

  void FullSimulation(Real dt);


  int subspace_coordinate_per_segment_;
  int segment_num_;
  int global_subspace_coord_num_;
  MatrixX reduced_linear_momentum_matrix_;
  MatrixX reduced_angular_momentum_matrix_[3];
  MatrixX subspace_;
  MatrixX subspace_transpose_;
  MatrixX full_subspace_;
  MatrixX full_subspace_transpose_;
  VectorX q_;
  VectorX vel_q_;
  std::vector<Vector3> net_ext_force_;
  std::vector<Vec3> local_u_;
  std::vector<Real> segment_mass_;
  std::vector<Vec3> segment_center_;
  std::vector<Quaternion<Real> > segment_rotation_;
  std::vector<Vec3> displacement_;
  std::vector<Vec3> segment_translational_vel_;
  std::vector<Vec3> prev_vert_;
  std::vector<Vec3> vel_;
  std::vector<Vec3> tmp_vert_buffer_;
  std::vector<Real> rest_edge_length_;
  ConjugateGradientSolver<Real>* cg_solver_;
  TetCollider* collider_;
  Real stiffness_;
  Real damping_;
  ConfigFile* conf_;
};

#endif // SUPSPACE_MASS_SPRING_VOLUMETRIC_OBJECT_H
