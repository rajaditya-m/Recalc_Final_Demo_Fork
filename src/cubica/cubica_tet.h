#ifndef CUBICATET_H
#define CUBICATET_H
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include "tet.h"

namespace solver {
template <class T> class BLOCK_MATRIX_GRAPH;
}
class BufferManager;
class SubspaceTet;
class CubicaTet : public Tet
{
public:
  typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mati;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatRow;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatCol;
  typedef Eigen::VectorXd Vec;
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Mat3;
  typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> Mat6;
  typedef Eigen::Vector3d Vec3;

  typedef Eigen::Map<Vec3> MapVec3;
  typedef Eigen::Map<Mat3> MapMat3;
  typedef Eigen::Map<Vec> MapVec;
  typedef Eigen::Map<Mat> MapMat;
  typedef Eigen::Map<MatCol> MapMatCol;
  typedef Tet Super;

  CubicaTet(const char* full_mesh_file, const char* folder);
  void ComputeSubspaceForceProjectionMatrix();
  void ComputeInterfaceTorqueMatrix();
  void LoadSubspace(std::function<const char* (int)> GetSubspaceFileName);
  void LoadCubature(std::function<const char* (int)> GetCubatureFileName);
  void PreComputeInterDomainCouplingMatrix();
//  virtual int Select(double* ray_start, double* ray_end, double* clicked_world_pos, double*selected_pos);
//  virtual bool GetVertexPosition(double* pos);
  template <class Vector0, class Vector1>
  void ProjectSubspaceForce(int e, Vector0 &force0, Vector1 &force1);

  ~CubicaTet();
  void Simulate(double dt);
  void NextRenderMode();
  void Render();
  void ComputeInterDomainSubspaceForce(int e, Vec& force0, Vec& force1);
  void ComputeVertexInterfacialArea();
  std::vector<int> parts_;
  int part_num_;
  std::vector<std::vector<int> > part_connectivity_;
  int block_edge_num_;

  int total_basis_num_;
  std::vector<int> basis_offset_;
  Vec vel_q_;
  Vec q_;

  Vec3 gravity_;
  std::vector<Mat> mat0_;
  std::vector<Mat> mat1_;
  std::vector<Mat> mat2_;
  std::vector<Vec> vec0_;
  std::vector<Vec> vec1_;

  std::vector<Mat> subspace_force_projection_matrix_;
  std::vector<Vec3> interface_center_of_mass_;
  std::vector<Mat> interface_torque_matrix_;
  std::vector<std::vector<int> > interface_vert_;
  std::vector<int> interface_vert_num_;
  std::vector<int> edge_list_;
  Mati block_edge_;
  std::vector<std::vector<int> > vert_local_id2global_id_;
  std::vector<double> vert_interfacial_area_;
  std::vector<int> vert_partition_info_;
  std::vector<std::vector<int> > is_interface_vert_;
  solver::BLOCK_MATRIX_GRAPH<double>* solver_;
  std::vector<SubspaceTet*> domain_container_;
  SubspaceTet** domain_;
  int ui_selected_domain_;
};

#endif // CUBICATET_H
