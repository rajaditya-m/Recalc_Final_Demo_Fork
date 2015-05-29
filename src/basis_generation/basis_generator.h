#ifndef BASIS_GENERATOR_H
#define BASIS_GENERATOR_H
#include <cstdio>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Sparse>

double ComputeDiff(int n3, int r, double* basis0, double* basis1);

double DiffFromFile(const char* file_name0, const char* file_name1);
class TetMesh;
class VolumetricMesh;
class SparseMatrix;
class StVKHessianTensor;

class BasisGenerator {
public:
  enum DataOrigin {
    kModalAnalysis,
    kSimulationData
  };

  BasisGenerator(int vert_num, double* vert,
                 int tet_num, int* tets,
                 double youngs_modulus = 1.5e5,
                 double poisson_ratio = 0.45,
                 double density = 1000.0,
                 double* mass = NULL,
                 std::set<int> *fixed_vertex_set = NULL);
  BasisGenerator(TetMesh* mesh, double* mass = NULL, std::set<int>* fixed_vertex_set = NULL);

  void Construct(double* mass, std::set<int>* fixed_vertex_set);
 ~BasisGenerator();
  void SetFixedVertices(std::set<int> fixed_verts);
  void RemoveSixRigidModes(int numVectors, double * x);
  void SetFixedVertexBasisZero(int basis_num, double *basis);
  void preLoad(const char* prefix);
  void SetInterfaceVertices(std::vector<std::pair<int,int> > &iv, int split );
  void setPrefix(const char* pfx) { prefix_ = pfx;}
  void setStiffyMultiplier(int sm) { stiffyMultiplier_ = sm;    }
  void setStitchStiffnessMatrix(SparseMatrix* spm);// { stitched_stiffness_matrix_ = new SparseMatrix(*spm); }

  bool ComputeLinearModes(int num_of_desired_modes);
  void ComputeModalDerivatives();
  bool ComputeNonLinearModes(int numNonLinearModes);
  bool MassPCA(int input_basis_num, int output_basis_num, std::vector<double>& basis, int data_origin = kModalAnalysis);
  bool MassPCA2(int output_basis_num, Eigen::MatrixXd dataMatrix,std::vector<double> &basis);
  void ComputeAllModes(int linear_mode_num, int non_linear_mode_num);

  void massOrthogonalization(Eigen::MatrixXd linBasis, Eigen::SparseMatrix<double> massMatrix, Eigen::MatrixXd &massOrtho);
  void massOrthogonalizationFast(Eigen::MatrixXd &linBasis);

  void SaveLinearModes(const char* file_name_prefix);

  void RegenerateAllModes(int modeId);
  inline void generateColsOfIMUUT(int colIdx, double *x, double *eigs, int numRows, int numCols);
  inline void generateColsOfIMUUT2(int colIdx1, int colIdx2, double *x, double *eigs, int numRows, int numCols, double multiplier);

  // linear modes
  int rigid_mode_num_;
  int linear_mode_num_; // total number of linear modes including rigid modes
  std::vector<double> frequencies_;
  std::vector<double> linear_modes_; // in column-major order
  std::vector<double> pure_eigen_vectors_;
  std::vector<double> pure_eigen_values_;

  std::vector<double> stitched_eigen_vectors_1_;
  std::vector<double> stitched_eigen_vectors_2_;
  std::vector<double> stitched_linear_modes_;
  std::vector<double> stitched_eigen_values_1_;
  std::vector<double> stitched_eigen_values_2_;
  std::vector<double> stitched_frequencies_1_;
  std::vector<double> stitched_frequencies_2_;
  std::vector<double> stitched_non_linear_modes_;
  std::vector<double> stitched_cubature_weights_;
  int stitched_non_linear_mode_num_;
  int stitched_linear_mode_num_;

  // modal derivatives
  int modal_derivative_num_;
  std::vector<double> modal_derivatives_; // in column-major order

  // non linear modes
  int non_linear_mode_num_;
  std::vector<double> non_linear_modes_; // in column-major order
  std::vector<double> eigen_values_;

  std::set<int> fixed_vertices_;
  int numVertsToRemove_;
  std::vector<int> vertsToRemove_;
  VolumetricMesh* mesh_;
  SparseMatrix* mass_matrix_;
  Eigen::SparseMatrix<double> mass_matrix_eig_;

  std::vector<std::pair<int,int> > interface_vertex_;
  std::vector<int> nonZeroColumns_1_;
  std::vector<int> columnInformation_1_;
  std::vector<int> nonZeroColumns_2_;
  std::vector<int> columnInformation_2_;
  const char* prefix_;

  Eigen::MatrixXd rhsOriginal_;
  int numColsOriginalRHS_;

  SparseMatrix* stitched_stiffness_matrix_;
  StVKHessianTensor *stVKStiffnessHessian;

 std::vector<double> massInv_;
 std::vector<double> massFlat_;
 std::vector<double> massSqrt_;

 int stiffyMultiplier_;

};
#endif // BASIS_GENERATOR_H
