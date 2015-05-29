#ifndef SUBSPACETET_H
#define SUBSPACETET_H
#include <Eigen/Dense>
#include <set>
#include "quaternion.h"
#include "tet.h"


class SubspaceTet : public Tet {
public: public:
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
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Mat;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatCol;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatRow;
  typedef Eigen::VectorXd Vec;
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Mat3;
  typedef Eigen::Vector3d Vec3;

  typedef Eigen::Map<Vec3> MapVec3;
  typedef Eigen::Map<Mat3> MapMat3;
  typedef Eigen::Map<Vec> MapVec;
  typedef Eigen::Map<Mat> MapMat;
  typedef Eigen::Map<MatRow> MapMatRow;
  typedef Eigen::Map<MatCol> MapMatCol;
  typedef std::pair<int, double> CubaturePoint;

  SubspaceTet(const char *filename, double _limit_threshold, AffineTransformer<double>* affine_transformer,bool initialize_fem_module = true,bool load_interface = true, int split=0);
  virtual ~SubspaceTet() {}

  void UpdateRigidMotionAndLocalDeformation(double dt);
  double ComputeCubatureError(Vec &subspace_force, Mat& subspace_k);
  void SimulateWithRigidMotion(double dt);
  void SimulateWithCubature(double dt);
  void LoadCubature(const char* file_name);
  void LoadCubature(std::vector<int> &p,std::vector<double> &w);
  void BuildParallelAccelerationStructures();
  void UpdateLocalDeformation(double* new_vel_q, double dt);
  void GetReducedTangentStiffnessMatrix(Mat& reduced_k);
  void AddFictitiousForceAndGravity(double* force);
  template <class EigenMatrix>
  void LoadSubspaceFromEigenMatrix(EigenMatrix& basis) { this->basis_ = basis; OnSubspaceLoaded(); }
  void LoadSubspace(const char *basis_file);
  void LoadBinarySubspace(const char *basis_file);
  void LoadBinarySubspace(std::vector<double> &vec, int numVerts, int numBasis);
  void OnSubspaceLoaded();
  void OnSubspaceLoaded2();
  void VerifySubspaceDiagonality();
  void LoadMass(const char* file_name);
  void ComputeInertialTensor();
  void ProjectPositionToSubspace();
  void AnimateBasis();
  void Simulate(double dt);
  void SimulateWithReduceMassMatrix(double dt);
  void PrecomputeFastSanwichTransform();
  void ComputeMomentumAndTorqueMatrix();
  void ComputeCenterOfMass();
  void ComputeInterfaceTetData();
  void ComputeInterfaceCOMData();
  void ComputeInterfaceRotationMatrix();
  void ComputeCurrentInterfaceCOMData();
  void saveOldBasis();
  void LoadBasis(int basis_id);
  void NextBasis();
  void PrevBasis();
  void UpdatePosition();
  int getSimMode() { return sim_mode_;}
  void setSimMode(int x);// { sim_mode_ = x; inv_fem_->setSimMode(x);}
  inline static SubspaceTet::Mat3 GetSkewSymmetrixMatrix(const SubspaceTet::Vec3& vec) {
    SubspaceTet::Mat3 mat = SubspaceTet::Mat3::Zero();
    mat(0, 1) = -vec[2];
    mat(0, 2) = +vec[1];
    mat(1, 0) = +vec[2];
    mat(1, 2) = -vec[0];
    mat(2, 0) = -vec[1];
    mat(2, 1) = +vec[0];
    return mat;
  }

  Mat momentum_matrix_;
  Mat momentum_matrix_transpose_;
  Mat torque_matrix_;
  Mat torque_matrix_transpose_;
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
  Mat inertial_force_sandwich_;
  Mat coriolis_force_sandwich_[9];
  Vec centrifugal_force_sandwich_[9];
  Mat euler_force_sandwich_;

  Vec3 center_of_mass_;
  double magnitude_;
  int current_basis_;
  int basis_num_;
  std::vector<Vec3> vert_offset_from_mass_center_;
  std::vector<CubaturePoint> cubature_;
  Vec subspace_gravity_;
  std::vector<int> cubature_tets_;
  std::vector<Mat> vert_basis_;
  std::vector<Mat> vert_basis_transpose_;
  Vec q_;
  Vec vel_q_;
  Mat reduce_mass_matrix_;
  Mat basis_;
  Mat basis_transpose_;

  Mat old_basis_;
  Mat old_basis_transpose_;
  std::vector<Mat> old_vert_basis_;
  std::vector<Mat> old_vert_basis_transpose_;

  //Acceleration strucutres for cubatures
  std::vector<std::pair<int, int> >  vv_list_1_; // list of vert-vert pair in each domain
  std::vector<std::vector<VVCubature> >  vv_cubature_1_;
  std::vector<int>  v_list_1_;
  std::vector<std::vector<VVCubature> >  v_cubature_1_;

  std::vector<std::pair<int, int> >  vv_list_2_; // list of vert-vert pair in each domain
  std::vector<std::vector<VVCubature> >  vv_cubature_2_;
  std::vector<int>  v_list_2_;
  std::vector<std::vector<VVCubature> >  v_cubature_2_;



  Mat precompPdk_;
  bool hasPreCompPDK_;

  Vec3 translation_vel_;
  Vec3 translation_acc_;
  Vec3 angular_vel_;
  Vec3 angular_acc_;
  Vec3 gravity_;
  Mat3 rotation_;
  Mat3 inertial_tensor_;
  Mat3 inverse_intertial_tensor_;
  Quaternion<double> quaternion_;
  double total_mass_;

  std::vector<int> interface_1_tets;
  std::set<int> interface_1_verts_;
  Vec3 interface_1_restCOM_;
  Vec3 interface_1_currentCOM_;
  Vec3 domain_2_restCOM_;
  Vec3 domain_2_currentCOM_;
  Mat3 Apq;

  int sim_mode_;
  int stop ;

};

#endif // SUBSPACETET_H
