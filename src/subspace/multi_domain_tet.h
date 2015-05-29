#ifndef MULTI_DOMAIN_TET_H
#define MULTI_DOMAIN_TET_H
#include <set>
#include <vector>
#include <functional>
#include <utility>
#include <unordered_map>
#include <Eigen/Dense>

#include "tet.h"
#include "quaternion.h"


class ObjMesh;
class ObjMeshRender;
template <class Float> class ConjugateGradientSolver;
template <class Float> class BiconjugateGradientStablized;
namespace solver {
template <class T> class BLOCK_MATRIX_GRAPH;
}
class MultiDomainTet : public Tet {
public:
  enum BasisFileFormat {
    kText,
    kBinary
  };
  struct VVCubature {
    int cubature_tet;
    double cubature_weight;
    int i, j;
    VVCubature(int tet = -1, double w = 0, int i_ = -1, int j_ = -1) {
      cubature_tet = tet;
      cubature_weight = w;
      i = i_;
      j = j_;
    }
  };


  typedef Tet Super;
  typedef ConjugateGradientSolver<double> Solver;
  //  typedef BiconjugateGradientStablized<double> Solver;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
  typedef Eigen::Matrix<int, 5, 1> Vec5i;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatRow;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatCol;
  typedef Eigen::VectorXd Vec;
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Mat3; // don't change to column major
  typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> Mat6;
  typedef Eigen::Vector3d Vec3;

  typedef Eigen::Map<Vec3> MapVec3;
  typedef Eigen::Map<Mat3> MapMat3;
  typedef Eigen::Map<Vec> MapVec;
  typedef Eigen::Map<Mat> MapMat;
  typedef Eigen::Map<MatCol> MapMatCol;
  typedef std::pair<int, double> CubaturePoint;
  typedef std::pair<int, Vec3> ExtForce;
  MultiDomainTet(const char *filename, double _limit_threshold = 0, AffineTransformer<double>* affine_transformer = NULL);
  virtual ~MultiDomainTet();

  template <class Vector0, class Vector1>
  void ProjectInterfaceSubspaceForce(int e,
                                     int p0, int p1,
                                     Vector0& force0,
                                     Vector1& force1);
  void UpdateCrossProductMatrix();
  void EnableFullRigidSimulation();
  void GenerateCubicaStylePartition(const char * folder);
  void ComputeLagradianMatrix();
  void ComputeVertexBasis();
  void ProjectPositionToSubspace();
  void LoadLocalMass(const char * prefix);
  void LoadSubspace(std::function<const char* (int)> GetFileName, int basis_file_format = kText);
  void LoadPartitionInfo(const char * partition_folder);
  void LoadPartitionInfo(std::string partition_folder);
  void VerifySubspaceDiagonality();
  void RemoveConstrainedDofFromSubspace();
  void ComputeInterfaceTorqueMatrix();
  void ComputeMomentumAndTorqueMatrix();
  void ComputePartAffineTransformamtion();
  void ComputeMassCenterPos();
  void InitializePartAffineTransformation();
  void InitializeMassCenterOffset();
  void ComputeGlobalPositionFromSubspace();
  void AssembleGlobalBasis();
  void InitInterfaceTetrahedra();
  void RotateBasis();
  void LoadMass(const char * file_name);
  void UpdatePosition();

  void PrecomputeFastSandwichTransform();
  void LoadCubature(const char * path_prefix);
  void ComputeInvInertiaTensor();
  void MultiBodyFullSimulation(double dt);
  void MultiBodyFullSimulationExplicit(double dt);
  void FullSimulation(double dt);
  void AddFictitiousForceAndGravity(const std::vector<Vec3>& acceleration,
                                    const std::vector<Vec3>& angular_acceleration,
                                    Vec &subspace_rhs);
  void GetExternalForce(double dt, std::vector<MultiDomainTet::ExtForce> *ext_force_ptr = nullptr);
  Vec ComputeSubspaceRhs(const double dt, std::vector<Vec3> &acceleration, std::vector<Vec3>& angular_acceleration, double internal_force_scaling_factor = 1.0);
  void ComputeSubspaceStiffnessMatrix(const double dt, double internal_force_scaling_factor = 1.0);
  MatCol ComputeFullRigidStiffnessMatrixWithCubature(const double dt);
  MatCol ComputeFullRigidStiffnessMatrix(const double dt);
  std::vector<Vec3> ComputeFullRigidMotionRhs(const double dt, std::vector<ExtForce> *ext_force_ptr);
  std::vector<Vec3> ComputeFullRigidMotionRhsWithCubature(const double dt, std::vector<ExtForce> *ext_force_ptr);
  void SimulateRigidMotionFull(const double dt, Vec& new_rigid_velocity, std::vector<ExtForce>* ext_force_ptr);
  void SimulateRigidMotion(const double dt, Vec& new_rigid_velocity, std::vector<ExtForce>* exe_force);
  void SimulateWithGlobalSubspace(double dt);
  void SubspaceMultiBodySimulation(double dt);
  void SubspaceMultiBodySimulationOneStepWithCubature(double dt);
  void SubspaceMultiBodySimulationOneStep(double dt);
  void SubspaceMultiBodySimulationMultiIteration(double dt, int iteration);
  void SubspaceMultiBodySimulationWithCubature(double dt);
  MatCol ComputeGlobalStiffnessMatrix(double dt, MatCol &RigidStiffnessMatrix);
  void Simulate(double dt);
  void BuildTopology();
  void Render(int render_mode, double * pos);
  // For each partition find the vertex that is closet to the center of mass and set it as fixed vertex
  void SetConstrainedVertex();
  void BuildParallelAccelerationStructure();
  void set_fixed_domain(int fixed_domain);
  void UpdateTexturedSurfaceMeshVBO();
  virtual void EnalbeTexturedSurfaceRendering(const char *shader_source);
  virtual void SetFixedDomains(std::set<int> fixed_domains);

  int part_num_;
  int total_basis_num_;
  Vec3 gravity_;
  Vec q_;
  Vec vel_q_;
  MatCol global_basis_;
  MatCol global_basis_transpose_;

  std::vector<Vec3> current_offset_from_mass_center_; // in global frame
  std::vector<Mat3> cross_product_matrix_;
  std::vector<Vec3> acceleration_;
  std::vector<Vec3> angular_acceleration_;

  std::vector<Mat> vert_rigid_basis_;
  std::vector<Mat> vert_rigid_basis_transpose_;
  std::vector<MatCol> rigid_basis_;
  std::vector<MatCol> rigid_basis_transpose_;

  std::vector<int> part_basis_size_; // num of basis in each part
  std::vector<int> basis_offset_; // part basis starting posiiton in global basis
  std::vector<MatCol> basis_;
  std::vector<MatCol> basis_transpose_;
  std::vector<Mat> vert_basis_;
  std::vector<Mat> vert_basis_transpose_;

  std::vector<MatCol> initial_basis_;
  std::vector<MatCol> initial_basis_transpose_;
  std::vector<Vec3> mass_center_initial_offset_;
  std::unordered_map<int, int> bone_id2domain_id_;
  std::vector<int> domain_attached_bone_; // list of part id
  std::vector<int> vert_num_per_part_;
  std::vector<int> vert_bone_id_;
  std::vector<int> vert_part_id_;
  std::vector<int> local_vert_idx_;
  std::vector<std::vector<int> > vert_local_id2global_id_;

  std::vector<Mat6> lagrangian_matrix_;
  std::vector<Mat> interface_torque_matrix_;
  std::vector<Mat> interface_torque_matrix_tranpose_;
  std::vector<Mat> momentum_matrix_;
  std::vector<Mat> torque_matrix_;
  std::vector<Mat> momentum_matrix_transpose_;
  std::vector<Mat> torque_matrix_transpose_;

  // cubature related data
  std::vector<int> all_cubature_tet_;
  std::vector<int> all_cubature_vert_;
  std::vector<int> interface_cubature_tet_;
  std::vector<int> interface_cubature_vert_;
  std::vector<CubaturePoint > all_cubature_;
  std::vector<std::vector<CubaturePoint> > domain_cubature_;
  std::vector<std::vector<CubaturePoint> > interface_cubature_;

  std::vector<std::vector<CubaturePoint> > interface_rigid_cubature_;

  inline void FullCholSolver2ConstrainedCholSolver(solver::BLOCK_MATRIX_GRAPH<double> *full_chol_solver,
                                                   solver::BLOCK_MATRIX_GRAPH<double> *constrained_chol_solver);

  // fast sandwich transform
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // inertial force: F_v = -m_v * a,
  //                F(q) = \sum_{all v} -(U_v^T * m_v * a) = (\sum_{all v} -(m_v * U_v)) * a
  //            sandwich = (\sum_{all v} -(m_v * U_v))
  //                F(q) = intertial_force_sandwich_^{rx3) * acceleration^{3x1}
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Coriolis force: F = -2m w\times v
  //              F(q) = \sum_{all v} -2(U_v^T * m_v * [w] * U_v * dq/dt)
  //          sandwich = \sum_{all v} -(U_v^T * m_v * [w] * U_v)
  //              F(q) = 2 * (\sum_{i\in [0, 8]} coriolis_force_sandwich_[i] * [w](i / 3, i % 3)) * dq/dt
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Euler force: F_v = -m_v * dw/dt \times r_v
  //             F(q) = \sum_{all v} -(U_v^T * m_v * [dw/dt] * (u + r0)
  //                  = \sum_{all v} {-(U_v^T * m_v * [dw/dt] * U * q) + (U_v^T * m_v * [r0] * dw/dt)}
  //  first term computed from Coriolis force sandwich
  //             sandwich = \sum_{all v} (U_v^T * m_v * [r0])
  //            F(q) = (\sum_{i\in [0,8]} coriolis_force_sandwich_[i] * [dw/dt](i / 3, i % 3)) * q
  //                   + euler_force_sandwich_^{r*3} * dw/dt^{3x1}
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // centerifugal force: F = -m * w \times (w \times r)
  //                  F(q) = \sum_{all v} -(U_v^T * m_v * [w] * [w] * (u + r0)
  //                       = \sum_{all v} {(-U_v^T * m_v * [w][w] * U_v * q) + (-U_v^T * m_v * [w] * [w] * r0)
  //  first term computed from Coriolis force sandwich
  //              sandwich = \sum_{all v} (-U_v^T * m_v * [w]^2 * r0)
  //                  F(q) = (\sum_{i\in [0,8]} coriolis_force_sandwich_[i] * ([w]*[w])(i / 3, i % 3)) * q
  //                         + (\sum_{i\in [0,8]} centrifugal_force_sandwich_[i] * ([w]*[w])(i / 3, i % 3))
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  std::vector<Mat> inertial_force_sandwich_;
  std::vector<Mat> coriolis_force_sandwich_;
  std::vector<Vec> centrifugal_force_sandwich_;
  std::vector<Mat> euler_force_sandwich_;

  // Rotation of a part is \sum mi*(ri*Ui*q)*ri^T = \sum mi*ri*ri^T + \sum mi*Ui*q*ri^T
  std::vector<Mat3> rotation_sandwich0_;
  std::vector<Mat3> rotation_sandwich1_;

  std::vector<Vec3> vert_offset_from_mass_center_;
  std::vector<double> mass_per_part_;
  std::vector<Vec3> angular_vel_;
  std::vector<Vec3> translational_vel_;
  std::vector<Quaternion<double> > quaternion_;
  std::vector<Mat3> current_inertia_tensor_;
  std::vector<Mat3> current_inv_inertia_tensor_;
  std::vector<Mat3> inv_inertia_tensor_;
  std::vector<Mat3> inertia_tensor_;
  std::vector<Vec3> part_translation_;
  std::vector<Mat3> part_rotation_;
  std::vector<Mat3> part_rotation_transpose_;
  std::vector<Vec3> center_of_mass_;
  std::vector<Mat3> initial_rotation_transpose_; // bone rotation transpose at rest pose
  std::vector<Vec3> initial_center_of_mass_;

  std::vector<Vec3> local_r_;
  std::vector<Vec3> local_u_;
  std::vector<std::vector<int> > topology_; // topology information about different domains
  std::vector<Mat> tmp_v_matrix_;
  std::vector<Mat> tmp_e_matrix_;
  solver::BLOCK_MATRIX_GRAPH<double>* full_chol_solver_;
  solver::BLOCK_MATRIX_GRAPH<double>* chol_solver_;
  std::vector<std::pair<int, Vec3> > ext_force_;
  Solver* cg_solver_;

  std::vector<Vec5i> interface_domains_; // (num of domains, domain0, domain1, domain2, domain3)
  int interface_num_;

  // parallel acceleration structure
  std::vector<std::vector<std::pair<int, int> > > vv_list_; // list of vert-vert pair in each domain
  std::vector<std::vector<std::vector<VVCubature> > > vv_cubature_;

  std::vector<std::vector<int> > v_list_;
  std::vector<std::vector<std::vector<VVCubature> > > v_cubature_;

  std::vector<std::vector<std::pair<int, int> > > interface_vv_list_; // list of vert-vert pair in each domain
  std::vector<std::vector<std::vector<VVCubature> > > interface_vv_cubature_;
  std::vector<std::vector<int> > interface_v_list_;
  std::vector<std::vector<std::vector<VVCubature> > > interface_v_cubature_;

  std::vector<std::vector<std::pair<int, int> > > domain_incident_interface_;
  std::vector<Mat3> tmp_interface_rigid_k_;

  solver::BLOCK_MATRIX_GRAPH<double>* constrained_chol_solver_;
  std::vector<int> fixed_domains_;
  std::vector<int> is_fixed_domain_;
  std::vector<int> constrained_dof_2_non_constrained_dof_;
  std::vector<int> non_constrained_dof_2_constrained_dof_;
  int fixed_domain_;

  bool full_rigid_simulation_;
  std::vector<int> interface_tet_;


  // embeded obj mesh
  std::vector<int> group_triangle_offset_;
  std::vector<Vec3> vert_tangent_;
  std::vector<Vec3> tri_tangent_;
  unsigned int tangent_vbo_;
  unsigned int shader_;
  unsigned int texture_id_;
  unsigned int bump_texture_id_;
  std::vector<std::vector<int> > groups_; // number of rendering group
  std::vector<int> tex_handle_;  // the texture id for each group
  std::vector<double> tex_coord_; // texture coordinate for 3 vertices of each triangle
  std::vector<int> tex_idx_;
  std::vector<int> vert_obj2tet_;
  std::vector<int> vert_tet2obj_;
  std::vector<int> tri_obj2tet_;
  ObjMesh *obj_;
  ObjMeshRender *obj_render_;
  void ComputeTangentVector();
  void EmbededObjMesh(ObjMesh * obj, ObjMeshRender * renderer, bool import_mapping, const char* folder = nullptr);
  void BuildRenderGroups();
  void BuildVertTriMapping(const char* export_folder = nullptr);
  void UpdateObjMeshPosition();

  template <class Vec1, class Vec2>
  inline static void Expand(const std::vector<int>& constrained_dof_2_non_constrained_dof,
                            Vec1& constrained_vec,
                            Vec2& non_constrained_vec) {
    memset(&non_constrained_vec[0], 0, sizeof(double) * int(non_constrained_vec.size()));
    for (int i = 0; i < int(constrained_vec.size()); ++i) {
      non_constrained_vec[constrained_dof_2_non_constrained_dof[i]] = constrained_vec[i];
    }
  }

  template <class Vec1, class Vec2>
  inline static void Contract(const std::vector<int>& constrained_dof_2_non_constrained_dof,
                              Vec2& non_constrained_vec,
                              Vec1& constrained_vec) {
    for (int i = 0; i < int(constrained_vec.size()); ++i) {
      //    P(i);
      //    P(constrained_dof_2_non_constrained_dof[i]);
      constrained_vec[i] = non_constrained_vec[constrained_dof_2_non_constrained_dof[i]];
    }
  }
};

inline MultiDomainTet::Mat3 GetSkewSymmetrixMatrix(const MultiDomainTet::Vec3& vec) {
  MultiDomainTet::Mat3 mat = MultiDomainTet::Mat3::Zero();
  mat(0, 1) = -vec[2];
  mat(0, 2) = +vec[1];
  mat(1, 0) = +vec[2];
  mat(1, 2) = -vec[0];
  mat(2, 0) = -vec[1];
  mat(2, 1) = +vec[0];
  return mat;
}


#endif // MULTI_DOMAIN_TET_H
