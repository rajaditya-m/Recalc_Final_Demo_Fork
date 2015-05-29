#include "cubica_tet.h"
#include <cstdio>
#include <tuple>
#include <fstream>
#include <fixed_vector.h>
#include "opengl_helper.h"
#include "vector_io.h"
#include "tet_mesh_simulator_bridge.h"
#include "BLOCK_MATRIX_GRAPH.h"
#include "config_file.h"
#include "matlab_io.h"
#include "subspace_tet.h"
#include "vega_tet_mesh_io.h"
#include "global.h"
#include "vector_lib.h"
#include "print_macro.h"
#include "basis_io.h"
#ifdef HAS_BASIS_GENERATION_MODULE
#include "basis_generator.h"
#endif

CubicaTet::CubicaTet(const char* mesh_file, const char *folder)
  : Super(mesh_file, 0, NULL) {
  ui_selected_domain_ = -1;
  gravity_ = Vec3(global::gravity[0], global::gravity[1], global::gravity[2]);
  char file_name[512];
  if (0) {
    sprintf(file_name, "%s/global_mesh.veg", folder);
    std::vector<double> verts(rest_pos_, rest_pos_ + vertex_num_ *  3);
    std::vector<int> tets(tet_, tet_ + tet_number * 4);
    VegaTetMeshIO::Instance()->Write(file_name, verts, tets);
    exit(0);
  }
  // read the list of partitions and their partition ids
  {
    sprintf(file_name, "%s/partition_id.txt", folder);
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    in >> part_num_;
    parts_.resize(part_num_);
    for (int p = 0; p < part_num_; ++p) {
      in >> parts_[p];
    }
    in.close();
  }
  domain_container_.resize(part_num_);
  domain_ = &domain_container_[0];
  // read the tet mesh for each domain
  for (int p = 0; p < part_num_;  ++p) {
    //    KK;P(p);
    sprintf(file_name, "%s/partition_%d", folder, parts_[p]);
    P(p);
    domain_[p] = new SubspaceTet(file_name, 0, NULL);
    KK;
    //    sprintf(file_name, "%s/partition_%d.veg", folder, parts_[p]);
    //    std::vector<double> verts(domain_[p]->rest_pos_, domain_[p]->rest_pos_ + domain_[p]->vertex_num_ * 3);
    //    std::vector<int> tets(domain_[p]->tet_, domain_[p]->tet_ + domain_[p]->tet_number * 4);
    //    VegaTetMeshIO::Instance()->Write(file_name, verts, tets);
  }
  //  exit(0);
  vert_partition_info_ = std::vector<int>(vertex_num_ * 4, -1);
  block_edge_ = Mati::Zero(part_num_, part_num_);
  block_edge_.fill(-1);
  {
    sprintf(file_name, "%s/vert_partition.txt", folder);
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    int v_num;
    in >> v_num;
    ASSERT(v_num == vertex_num_);
    for (int i = 0; i < vertex_num_ * 4; ++i) {
      in >> vert_partition_info_[i];
    }
    in.close();
  }

  vert_local_id2global_id_.resize(part_num_);
  is_interface_vert_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    is_interface_vert_[p] = std::vector<int>(domain_[p]->vertex_num_, 0);
    vert_local_id2global_id_[p].resize(domain_[p]->vertex_num_);
  }

  part_connectivity_ = std::vector<std::vector<int> >(part_num_, std::vector<int>(part_num_, 0));
  for (int v = 0; v < vertex_num_; ++v) {
    // a vertex belongs to two domains
    if (vert_partition_info_[v * 4 + 2] != -1) {
      ASSERT(vert_partition_info_[v * 4 + 0] != -1);
      int p0 = vert_partition_info_[v * 4 + 0];
      int v0 = vert_partition_info_[v * 4 + 1];
      int p1 = vert_partition_info_[v * 4 + 2];
      int v1 = vert_partition_info_[v * 4 + 3];
      part_connectivity_[p0][p1] = 1;
      part_connectivity_[p1][p0] = 1;
      is_interface_vert_[p0][v0] = 1;
      is_interface_vert_[p1][v1] = 1;
      vert_local_id2global_id_[p0][v0] = v;
      vert_local_id2global_id_[p1][v1] = v;
    } else {
      int p0 = vert_partition_info_[v * 4 + 0];
      int v0 = vert_partition_info_[v * 4 + 1];
      vert_local_id2global_id_[p0][v0] = v;
    }
  }

  std::vector<int> dummy_block_size(part_num_, 1);
  auto tmp_connectivity_ = part_connectivity_;
  solver::BLOCK_MATRIX_GRAPH<double> solve(tmp_connectivity_, dummy_block_size);
  block_edge_num_ = solve.e_number0;
  edge_list_ = std::vector<int>(solve.E, solve.E + solve.e_number0 * 2);
  for (int i = 0; i < part_num_; ++i) {
    for (int j = 0; j < part_num_; ++j) {
      block_edge_(i, j) = tmp_connectivity_[i][j];
    }
  }

  interface_vert_num_ = std::vector<int>(block_edge_num_, 0);
  interface_vert_.resize(block_edge_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    if (vert_partition_info_[v * 4 + 2] != -1) {
      ASSERT(vert_partition_info_[v * 4 + 0] != -1);
      int p0 = vert_partition_info_[v * 4 + 0];
      int p1 = vert_partition_info_[v * 4 + 2];
      int e = block_edge_(p0, p1);
      interface_vert_num_[e]++;
      interface_vert_[e].push_back(v);
      interface_vert_[e].push_back(vert_partition_info_[v * 4 + 1]);
      interface_vert_[e].push_back(vert_partition_info_[v * 4 + 3]);
    }
  }
  subspace_force_projection_matrix_.resize(block_edge_num_ * 4);
  ComputeVertexInterfacialArea();
}

void CubicaTet::ComputeSubspaceForceProjectionMatrix() {
  typedef Eigen::Matrix<double, 6, 1> Vec6;
  typedef Eigen::Matrix<double, 6, 6> Mat6;
  for (int e = 0; e < block_edge_num_; ++e) {
    int p0 = edge_list_[e * 2 + 0];
    int p1 = edge_list_[e * 2 + 1];
    Mat momentum_matrix0 = domain_[p0]->rotation_ * domain_[p0]->momentum_matrix_;
    Mat momentum_matrix1 = domain_[p1]->rotation_ * domain_[p1]->momentum_matrix_;
    Mat interface_torque_matrix0 = domain_[p0]->rotation_ * interface_torque_matrix_[e * 2 + 0];
    Mat interface_torque_matrix1 = domain_[p1]->rotation_ * interface_torque_matrix_[e * 2 + 1];

    Mat momentum_matrix0_transpose = momentum_matrix0.transpose();
    Mat momentum_matrix1_transpose = momentum_matrix1.transpose();
    Mat interface_torque_matrix0_transpose = interface_torque_matrix0.transpose();
    Mat interface_torque_matrix1_transpose = interface_torque_matrix1.transpose();

    Mat6 lagrangian_matrix;
    lagrangian_matrix.block<3, 3>(0, 0) = momentum_matrix0 * momentum_matrix0_transpose +
                                          momentum_matrix1 * momentum_matrix1_transpose;

    lagrangian_matrix.block<3, 3>(0, 3) = momentum_matrix0 * interface_torque_matrix0_transpose +
                                          momentum_matrix1 * interface_torque_matrix1_transpose;

    lagrangian_matrix.block<3, 3>(3, 0) = interface_torque_matrix0 * momentum_matrix0_transpose +
                                          interface_torque_matrix1 * momentum_matrix1_transpose;

    lagrangian_matrix.block<3, 3>(3, 3) = interface_torque_matrix0 * interface_torque_matrix0_transpose +
                                          interface_torque_matrix1 * interface_torque_matrix1_transpose;

    Mat6 inv_lagrangian_matrix = (-1) * lagrangian_matrix.inverse();
    //    Mat ldlt = lagrangian_matrix.ldlt();
    //    for (int i = 0; i < 6; ++i) {
    //      Vec6 vec6_zero = Vec6::Zero();
    //      vec6_zero(i) = 1;
    //      inv_lagrangian_matrix.col(i) = lagrangian_matrix.ldlt().solve(vec6_zero);
    //      vec6_zero(i) = 0;
    //    }

    //    {
    //      Mat6 zero = lagrangian_matrix.inverse() * lagrangian_matrix - Mat6::Identity();
    //      ASSERT(zero.norm() < 1e-3, P(zero.norm(), zero.maxCoeff(), zero.minCoeff()));
    //    }

    Mat tmp00 = inv_lagrangian_matrix.block<3, 3>(0, 0) * momentum_matrix0 +
                inv_lagrangian_matrix.block<3, 3>(0, 3) * interface_torque_matrix0;

    Mat tmp01 = inv_lagrangian_matrix.block<3, 3>(0, 0) * momentum_matrix1 +
                inv_lagrangian_matrix.block<3, 3>(0, 3) * interface_torque_matrix1;

    Mat tmp10 = inv_lagrangian_matrix.block<3, 3>(3, 0) * momentum_matrix0 +
                inv_lagrangian_matrix.block<3, 3>(3, 3) * interface_torque_matrix0;

    Mat tmp11 = inv_lagrangian_matrix.block<3, 3>(3, 0) * momentum_matrix1 +
                inv_lagrangian_matrix.block<3, 3>(3, 3) * interface_torque_matrix1;

    if (0) {
      using namespace dj;
      WriteEigenMatrixToMatlab(lagrangian_matrix, "/tmp/lag");
      WriteEigenMatrixToMatlab(inv_lagrangian_matrix, "/tmp/invlag");
      WriteEigenMatrixToMatlab(momentum_matrix0, "/tmp/m0");
      WriteEigenMatrixToMatlab(momentum_matrix1, "/tmp/m1");
      WriteEigenMatrixToMatlab(interface_torque_matrix0, "/tmp/t0");
      WriteEigenMatrixToMatlab(interface_torque_matrix1, "/tmp/t1");
      WriteEigenMatrixToMatlab(subspace_force_projection_matrix_[e * 4 + 0], "/tmp/r0");
      WriteEigenMatrixToMatlab(subspace_force_projection_matrix_[e * 4 + 1], "/tmp/r1");
      WriteEigenMatrixToMatlab(subspace_force_projection_matrix_[e * 4 + 2], "/tmp/r2");
      WriteEigenMatrixToMatlab(subspace_force_projection_matrix_[e * 4 + 3], "/tmp/r3");
    }

    subspace_force_projection_matrix_[e * 4 + 0] = momentum_matrix0_transpose * tmp00 + interface_torque_matrix0_transpose * tmp10;
    subspace_force_projection_matrix_[e * 4 + 1] = momentum_matrix0_transpose * tmp01 + interface_torque_matrix0_transpose * tmp11;

    subspace_force_projection_matrix_[e * 4 + 2] = momentum_matrix1_transpose * tmp00 + interface_torque_matrix1_transpose * tmp10;
    subspace_force_projection_matrix_[e * 4 + 3] = momentum_matrix1_transpose * tmp01 + interface_torque_matrix1_transpose * tmp11;

    for (int i = 0; i < subspace_force_projection_matrix_[e * 4 + 0].rows(); ++i) {
      subspace_force_projection_matrix_[e * 4 + 0](i, i) += 1;
    }

    for (int i = 0; i < subspace_force_projection_matrix_[e * 4 + 3].rows(); ++i) {
      subspace_force_projection_matrix_[e * 4 + 3](i, i) += 1;
    }
  }
}

void CubicaTet::ComputeInterfaceTorqueMatrix() {
  auto GetInterfaceTorqueMatrix = [&](int p, Vec3 & center_of_mass, Mat & matrix) {
    matrix = Mat::Zero(3, domain_[p]->basis_num_);
    for (int v = 0; v < domain_[p]->vertex_num_; ++v) {
      double mass = domain_[p]->mass_[v];
      for (int c = 0; c < domain_[p]->basis_num_; ++c) {
        Vec3 col(domain_[p]->basis_(v * 3 + 0, c),
                 domain_[p]->basis_(v * 3 + 1, c),
                 domain_[p]->basis_(v * 3 + 2, c));
        col *= mass;
        MapVec3 pos(domain_[p]->rest_pos_ + v * 3);
        Vec3 r = pos - center_of_mass;
        Vec3 cross = r.cross(col);
        matrix(0, c) += cross[0];
        matrix(1, c) += cross[1];
        matrix(2, c) += cross[2];
      }
    }
  };

  interface_center_of_mass_.clear();
  interface_torque_matrix_.resize(block_edge_num_ * 2);
  for (int e = 0; e < block_edge_num_; ++e) {
    int p0 = edge_list_[e * 2 + 0];
    int p1 = edge_list_[e * 2 + 1];
    int vert_num = int(interface_vert_[e].size()) / 3;
    Vec3 center_of_mass(0, 0, 0);
    double total_mass = 0;
    for (int i = 0; i < vert_num; ++i) {
      int v = interface_vert_[e][i * 3 + 1];
      double mass = domain_[p0]->mass_[v];
      center_of_mass += (mass * MapVec3(domain_[p0]->rest_pos_ + v * 3));
      total_mass += mass;
    }
    center_of_mass /= total_mass;
    interface_center_of_mass_.push_back(center_of_mass);

    GetInterfaceTorqueMatrix(p0, center_of_mass, interface_torque_matrix_[e * 2 + 0]);
    GetInterfaceTorqueMatrix(p1, center_of_mass, interface_torque_matrix_[e * 2 + 1]);
  }
}

void CubicaTet::LoadSubspace(std::function<const char *(int)> GetSubspaceFileName) {
  basis_offset_ = std::vector<int>(part_num_ + 1, 0);
  total_basis_num_ = 0;
  for (int p = 0; p < part_num_; ++p) {
    const char* file_name = GetSubspaceFileName(parts_[p]);
    std::vector<double> basis;
    int v_num, basis_num;
    ReadBasisInBinary(file_name, v_num, basis_num, basis);
    total_basis_num_ += basis_num;
    basis_offset_[p + 1] = basis_offset_[p] + basis_num;
    //    ReadBasisInText(file_name, v_num, basis_num, basis);
    ASSERT(v_num == domain_[p]->vertex_num_);
    MapMatCol mat_basis(&basis[0], v_num * 3, basis_num);
    domain_[p]->LoadSubspaceFromEigenMatrix(mat_basis);
  }
  ComputeInterfaceTorqueMatrix();
  PreComputeInterDomainCouplingMatrix();
}

void CubicaTet::LoadCubature(std::function<const char *(int)> GetCubatureFileName) {
  (void) GetCubatureFileName;
  for (int p = 0; p < part_num_; ++p) {
    domain_[p]->LoadCubature(NULL);
  }
}

void CubicaTet::PreComputeInterDomainCouplingMatrix() {
  // inter-domain coupling force on domain q, ('ed quanity is from neighbor domain
  // F_vv'(q, q') = -\sum_{all v} k_v*(R*U_v*q + R*r0 + MassCenter - R'U'_v*q'- R'*r0'- R'MassCenter')
  // f(q, q') = -\sum_{all v} k_v*(U_v^T*U_v*q + U_v^T*r0 + U^T*R^T*MassCenter - U_v^T*R^T*R'U'_v*q'-U_v^T*R^T*R'*r0'-U_v^T*R^TMassCenter')
  // = \sum_{all v} k_v*(-U_v^T*U_v*q + U_v^T*R^T*R'*U'_v*q' + U_v^T*R^T*(-MassCenter+MassCenter')  - U_v^T*r0 + U_v^T*R^T*R'*r0')
  // = \sum_{all v} mat0*q            + mat1*{R^T*R'}*q'     + mat2*{R^T}*(-MassCenter+MassCenter') + vec0 + vec1*{R^T*R'} * r0'
  // mat0 = \sum_{all v} -k_v*U_v^T*U_v
  // mat1 = \sum_{all v} +k_v*U_v^T*R^T*R'*U'_v
  // mat2 = \sum_{all v} +k_v*U_v^T*R^T
  // vec0 = \sum_{all v} -k_v*U_v^T*r0
  // vec1 = \sum_{all v} +k_v*U_v^T*R^T*R'*r0'
  mat0_.resize(block_edge_num_ * 2);
  mat1_.resize(block_edge_num_ * 2 * 9);
  mat2_.resize(block_edge_num_ * 2 * 9);
  vec0_.resize(block_edge_num_ * 2);
  vec1_.resize(block_edge_num_ * 2 * 9);
  for (int e = 0; e < block_edge_num_; ++e) {
    int p0 = edge_list_[e * 2 + 0];
    int p1 = edge_list_[e * 2 + 1];
    mat0_[e * 2 + 0] = Mat::Zero(domain_[p0]->basis_num_, domain_[p0]->basis_num_);
    mat0_[e * 2 + 1] = Mat::Zero(domain_[p1]->basis_num_, domain_[p1]->basis_num_);

    vec0_[e * 2 + 0] = Vec::Zero(domain_[p0]->basis_num_);
    vec0_[e * 2 + 1] = Vec::Zero(domain_[p1]->basis_num_);

    for (int i = 0; i < 9; ++i) {
      mat1_[e * 18 + i] = Mat::Zero(domain_[p0]->basis_num_, domain_[p1]->basis_num_);
      mat1_[e * 18 + 9 + i] = Mat::Zero(domain_[p1]->basis_num_, domain_[p0]->basis_num_);

      mat2_[e * 18 + i] = Mat::Zero(domain_[p0]->basis_num_, 3);
      mat2_[e * 18 + 9 + i] = Mat::Zero(domain_[p1]->basis_num_, 3);

      vec1_[e * 18 + i] = Vec::Zero(domain_[p0]->basis_num_);
      vec1_[e * 18 + 9 + i] = Vec::Zero(domain_[p1]->basis_num_);
    }

    int vert_num = int(interface_vert_[e].size()) / 3;
    for (int v = 0; v < vert_num; ++v) {
      int global_v = interface_vert_[e][v * 3 + 0];
      int v0 = interface_vert_[e][v * 3 + 1];
      int v1 = interface_vert_[e][v * 3 + 2];
      mat0_[e * 2 + 0] -= vert_interfacial_area_[global_v] * (domain_[p0]->vert_basis_transpose_[v0] * domain_[p0]->vert_basis_[v0]);
      mat0_[e * 2 + 1] -= vert_interfacial_area_[global_v] * (domain_[p1]->vert_basis_transpose_[v1] * domain_[p1]->vert_basis_[v1]);

      Vec3 r0 = MapVec3(domain_[p0]->rest_pos_ + v0 * 3) - domain_[p0]->center_of_mass_;
      Vec3 r1 = MapVec3(domain_[p1]->rest_pos_ + v1 * 3) - domain_[p1]->center_of_mass_;

      vec0_[e * 2 + 0] -= vert_interfacial_area_[global_v] * (domain_[p0]->vert_basis_transpose_[v0] * r0);
      vec0_[e * 2 + 1] -= vert_interfacial_area_[global_v] * (domain_[p1]->vert_basis_transpose_[v1] * r1);
      for (int row = 0, idx = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col, ++idx) {
          Mat3 R = Mat3::Zero();
          R(row, col) = 1;
          mat1_[e * 18 + idx] += vert_interfacial_area_[global_v] * (domain_[p0]->vert_basis_transpose_[v0] * (R * domain_[p1]->vert_basis_[v1]));
          mat1_[e * 18 + 9 + idx] += vert_interfacial_area_[global_v] * (domain_[p1]->vert_basis_transpose_[v1] * (R * domain_[p0]->vert_basis_[v0]));

          mat2_[e * 18 + idx] += vert_interfacial_area_[global_v] * (domain_[p0]->vert_basis_transpose_[v0] * R);
          mat2_[e * 18 + 9 + idx] += vert_interfacial_area_[global_v] * (domain_[p1]->vert_basis_transpose_[v1] * R);

          vec1_[e * 18 + idx] += vert_interfacial_area_[global_v] * (domain_[p0]->vert_basis_transpose_[v0] * R * r1);
          vec1_[e * 18 + 9 + idx] += vert_interfacial_area_[global_v] * (domain_[p1]->vert_basis_transpose_[v1] * R * r0);
        }
      }
    }
  }
}

//int CubicaTet::Select(double *ray_start, double *ray_end, double *clicked_world_pos, double *selected_pos)
//{
//  for (int p = 0; p < vertex_num_; ++p) {
//    int vert = domain_[p]->Select(ray_start, ray_end, clicked_world_pos, selected_pos);
//    if (vert >= 0) {
//      ui_selected_domain_ = p;
//      return vert;
//    }
//  }
//  return -1;
//}

//bool CubicaTet::GetVertexPosition(double *pos)
//{
//  if (ui_selected_domain_ >= 0) {
//    SubspaceTet* mesh = domain_[ui_selected_domain_];
//    int vert = mesh->SelectedVertex();
////    ASSERT(vert >= 0);
//    pos[0] = mesh->X[vert * 3 + 0];
//    pos[1] = mesh->X[vert * 3 + 1];
//    pos[2] = mesh->X[vert * 3 + 2];
//    return true;
//  } else {
//    return false;
//  }
//}


template <class Vector0, class Vector1>
void CubicaTet::ProjectSubspaceForce(int e, Vector0 &force0, Vector1 &force1) {
#if 0
  force0 = subspace_force_projection_matrix_[e * 4 + 0] * force0 + subspace_force_projection_matrix_[e * 4 + 1] * force1;
  force1 = subspace_force_projection_matrix_[e * 4 + 2] * force0 + subspace_force_projection_matrix_[e * 4 + 3] * force1;
#else
  typedef Eigen::Matrix<double, 6, 1> Vec6;
  int p0 = edge_list_[e * 2 + 0];
  int p1 = edge_list_[e * 2 + 1];
  //  P(domain_container_.size());
  Mat3 rot0 = domain_[p0]->rotation_;
  Mat3 rot1 = domain_[p1]->rotation_;
  Mat3 rot0_tran = rot0.transpose();
  Mat3 rot1_tran = rot1.transpose();
  Mat6 lagrangian_matrix;
  Mat mom0 = domain_[p0]->momentum_matrix_;
  Mat mom1 = domain_[p1]->momentum_matrix_;
  Mat mom0_tran = mom0.transpose();
  Mat mom1_tran = mom1.transpose();
  Mat tor0 = interface_torque_matrix_[e * 2 + 0];
  Mat tor1 = interface_torque_matrix_[e * 2 + 1];
  Mat tor0_tran = tor0.transpose();
  Mat tor1_tran = tor1.transpose();
  lagrangian_matrix.block<3, 3>(0, 0) = rot0 * (mom0 * mom0_tran) * rot0_tran + rot1 * (mom1 * mom1_tran) * rot1_tran;
  lagrangian_matrix.block<3, 3>(0, 3) = rot0 * (mom0 * tor0_tran) * rot0_tran + rot1 * (mom1 * tor1_tran) * rot1_tran;
  lagrangian_matrix.block<3, 3>(3, 0) = rot0 * (tor0 * mom0_tran) * rot0_tran + rot1 * (tor1 * mom1_tran) * rot1_tran;
  lagrangian_matrix.block<3, 3>(3, 3) = rot0 * (tor0 * tor0_tran) * rot0_tran + rot1 * (tor1 * tor1_tran) * rot1_tran;

  Vec6 residual;
  residual.block<3, 1>(0, 0) = rot0 * (mom0 * force0) + rot1 * (mom1 * force1);
  residual.block<3, 1>(3, 0) = rot0 * (tor0 * force0) + rot1 * (tor1 * force1);
  Vec6 lagranian_multiplier = lagrangian_matrix.ldlt().solve(residual);
  MapVec3 lambda0(&lagranian_multiplier[0]);
  if (0) {
    using namespace dj;
    WriteEigenMatrixToMatlab(lagrangian_matrix, "/tmp/lag");
    Mat m0 = rot0 * mom0;
    Mat m1 = rot1 * mom1;
    Mat t0 = rot0 * tor0;
    Mat t1 = rot1 * tor1;
    //      WriteEigenMatrixToMatlab(inv_lagrangian_matrix, "/tmp/invlag");
    WriteEigenMatrixToMatlab(m0, "/tmp/m0");
    WriteEigenMatrixToMatlab(m1, "/tmp/m1");
    WriteEigenMatrixToMatlab(t0, "/tmp/t0");
    WriteEigenMatrixToMatlab(t1, "/tmp/t1");
    WriteVectorToMatlab(6, &lagranian_multiplier[0], "/tmp/lambda");
    WriteVectorToMatlab(force0.size(), &force0[0], "/tmp/f0");
    WriteVectorToMatlab(force1.size(), &force1[0], "/tmp/f1");
  }
  MapVec3 lambda1(&lagranian_multiplier[3]);
  force0 -= (mom0_tran * (rot0_tran * lambda0)) + (tor0_tran * (rot0_tran * lambda1));
  force1 -= (mom1_tran * (rot1_tran * lambda0)) + (tor1_tran * (rot1_tran * lambda1));
  if (0) {
    using namespace dj;
    WriteVectorToMatlab(force0.size(), &force0[0], "/tmp/truth0");
    WriteVectorToMatlab(force1.size(), &force1[0], "/tmp/truth1");
  }
  //  force0 =  -(mom0_tran * (rot0_tran * lambda0)) + (tor0_tran * (rot0_tran * lambda1));
  //  force1 =  -(mom1_tran * (rot1_tran * lambda0)) + (tor1_tran * (rot1_tran * lambda1));
#endif
}

CubicaTet::~CubicaTet() {
  for (SubspaceTet * tet : domain_container_) {
    delete tet;
  }
}

void CubicaTet::Simulate(double dt) {
  const float kInternalForceScaling = conf.Get<double>("internal force scaling");
  //  dt = 0.033; P(dt);
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));
  // gravity
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += domain_[p]->total_mass_ * gravity_;
  }

  // collision
  std::vector<std::tuple<int, int, Vec3> > collision_force;
  collision_force.reserve(vertex_num_ / 100);
  const double kFloor = 0.00;
  const double kFloorStiffness = 3000;
  for (int p = 0; p < part_num_; ++p) {
    for (int v = 0; v < domain_[p]->vertex_num_; ++v) {
      //    if (0)
      if (domain_[p]->X[v * 3 + 1] < kFloor) {
        Vec3 force(0, 0, 0);
        force[1] = (kFloor - domain_[p]->X[v * 3 + 1]) * domain_[p]->mass_[v] * kFloorStiffness;
        collision_force.emplace_back(p, v, force);
        net_force[p][1] += force[1];

        Vec3 r = MapVec3(domain_[p]->X + v * 3) - domain_[p]->center_of_mass_;
        net_torque[p] += r.cross(force);
      }
    }
  }
  Vec3 force;
  int ui_force_vert = GetUIForce(&force[0]);
  if (ui_force_vert >= 0) {
    int p = vert_partition_info_[ui_force_vert * 4 + 0];
    int local_v = vert_partition_info_[ui_force_vert * 4 + 1];
    collision_force.emplace_back(p, local_v, force);
  }

  #define ASSUME_RIGID
  //    ComputeSubspaceForceProjectionMatrix();
  // inter-domain elastic force
  //  if (0)
  for (int e = 0; e < block_edge_num_; ++e) {
    int p0 = edge_list_[e * 2 + 0];
    int p1 = edge_list_[e * 2 + 1];
    // inter-domain coupling force on domain q, ('ed quanity is from neighbor domain
    // F_vv'(q, q') = -\sum_{all v} k_v*(R*U_v*q + R*r0 + MassCenter - R'U'_v*q'- R'*r0'- R'MassCenter')
    // f(q, q') = -\sum_{all v} k_v*(U_v^T*U_v*q + U_v^T*r0 + U^T*R^T*MassCenter - U_v^T*R^T*R'U'_v*q'-U_v^T*R^T*R'*r0'-U_v^T*R^TMassCenter')
    // = \sum_{all v} k_v*(-U_v^T*U_v*q + U_v^T*R^T*R'*U'_v*q' + U_v^T*R^T*(-MassCenter+MassCenter')  - U_v^T*r0 + U_v^T*R^T*R'*r0')
    // = \sum_{all v} mat0*q            + mat1*{R^T*R'}*q'     + mat2*{R^T}*(-MassCenter+MassCenter') + vec0 + vec1*{R^T*R'}
    // mat0 = \sum_{all v} -k_v*U_v^T*U_v
    // mat1 = \sum_{all v} +k_v*U_v^T*R^T*R'*U'_v
    // mat2 = \sum_{all v} +k_v*U_v^T*R^T
    // vec0 = \sum_{all v} -k_v*U_v^T*r0
    // vec1 = \sum_{all v} +k_v*U_v^T*R^T*R'*r0'
#ifndef ASSUME_RIGID
    Vec subspace_force0 = mat0_[e * 2 + 0] * domain_[p0]->q_ + vec0_[e * 2 + 0];
    Vec subspace_force1 = mat0_[e * 2 + 1] * domain_[p1]->q_ + vec0_[e * 2 + 1];
    Mat mat10 = Mat::Zero(domain_[p0]->basis_num_, domain_[p0]->basis_num_);
    Mat mat11 = Mat::Zero(domain_[p1]->basis_num_, domain_[p1]->basis_num_);
    Mat mat20 = Mat::Zero(domain_[p0]->basis_num_, 3);
    Mat mat21 = Mat::Zero(domain_[p1]->basis_num_, 3);
    Vec vec10 = Vec::Zero(domain_[p0]->basis_num_);
    Vec vec11 = Vec::Zero(domain_[p1]->basis_num_);
    Mat3& rot0 = domain_[p0]->rotation_;
    Mat3 rot0_tran = rot0.transpose();
    Mat3& rot1 = domain_[p1]->rotation_;
    Mat3 rot1_tran = rot1.transpose();
    Mat3 rot0_tran_dot_rot1 = rot0_tran * rot1;
    Mat3 rot1_tran_dot_rot0 = rot1_tran * rot0;
    for (int row = 0, idx = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col, ++idx) {
        mat10 += rot0_tran_dot_rot1(row, col) * mat1_[e * 18 + idx];
        mat11 += rot1_tran_dot_rot0(row, col) * mat1_[e * 18 + 9 + idx];

        mat20 += rot0_tran(row, col) * mat2_[e * 18 + idx];
        mat21 += rot1_tran(row, col) * mat2_[e * 18 + 9 + idx];

        vec10 += rot0_tran_dot_rot1(row, col) * vec1_[e * 18 + idx];
        vec11 += rot1_tran_dot_rot0(row, col) * vec1_[e * 18 + 9 + idx];
      }
    }
    Vec diff_mass_center = domain_[p1]->center_of_mass_ - domain_[p0]->center_of_mass_;
    subspace_force0 += mat10 * domain_[p0]->q_ + mat20 * (diff_mass_center) + vec10;
    subspace_force1 += mat11 * domain_[p1]->q_ - mat21 * (diff_mass_center) + vec11;
#else
    Vec subspace_force0 = vec0_[e * 2 + 0];
    Vec subspace_force1 = vec0_[e * 2 + 1];
    Mat mat10 = Mat::Zero(domain_[p0]->basis_num_, domain_[p0]->basis_num_);
    Mat mat11 = Mat::Zero(domain_[p1]->basis_num_, domain_[p1]->basis_num_);
    Mat mat20 = Mat::Zero(domain_[p0]->basis_num_, 3);
    Mat mat21 = Mat::Zero(domain_[p1]->basis_num_, 3);
    Vec vec10 = Vec::Zero(domain_[p0]->basis_num_);
    Vec vec11 = Vec::Zero(domain_[p1]->basis_num_);
    Mat3& rot0 = domain_[p0]->rotation_;
    Mat3 rot0_tran = rot0.transpose();
    Mat3& rot1 = domain_[p1]->rotation_;
    Mat3 rot1_tran = rot1.transpose();
    Mat3 rot0_tran_dot_rot1 = rot0_tran * rot1;
    Mat3 rot1_tran_dot_rot0 = rot1_tran * rot0;
    for (int row = 0, idx = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col, ++idx) {
        mat10 += rot0_tran_dot_rot1(row, col) * mat1_[e * 18 + idx];
        mat11 += rot1_tran_dot_rot0(row, col) * mat1_[e * 18 + 9 + idx];

        mat20 += rot0_tran(row, col) * mat2_[e * 18 + idx];
        mat21 += rot1_tran(row, col) * mat2_[e * 18 + 9 + idx];

        vec10 += rot0_tran_dot_rot1(row, col) * vec1_[e * 18 + idx];
        vec11 += rot1_tran_dot_rot0(row, col) * vec1_[e * 18 + 9 + idx];
      }
    }
    Vec diff_mass_center = domain_[p1]->center_of_mass_ - domain_[p0]->center_of_mass_;
    subspace_force0 += vec10 + mat20 * (diff_mass_center);
    subspace_force1 += vec11 - mat21 * (diff_mass_center);
#endif
    // Verify fast sandwich transform
    if (0) {
      Vec f0, f1;
      ComputeInterDomainSubspaceForce(e, f0, f1);
      Vec diff0 = f0 - subspace_force0;
      Vec diff1 = f1 - subspace_force1;
      //      P(diff0.norm(), f0.norm(), subspace_force0.norm());
      //      P(diff1.norm(), f1.norm(), subspace_force1.norm());
      ASSERT(diff0.norm() < 1e-8, P(diff0.norm(), f0.norm(), subspace_force0.norm()));
      ASSERT(diff1.norm() < 1e-8, P(diff1.norm(), f1.norm(), subspace_force1.norm()));
    }
    //    {
    //    Vec f0 = subspace_force0, f1 = subspace_force1;
    //    PROJECT_SUBSPACE_FORCE(e, f0, f1);
    //    }
    ProjectSubspaceForce(e, subspace_force0, subspace_force1);
    if (0) {
      //      Vec diff0 = f0 - subspace_force0;
      //      Vec diff1 = f1 - subspace_force1;
      //      P(diff0.norm(), f0.norm(), subspace_force0.norm());
      //      P(diff1.norm(), f1.norm(), subspace_force1.norm());
      //      ASSERT(diff0.norm() < 1e-8);
      //      ASSERT(diff1.norm() < 1e-8);
    }
    // Verify projection matrix
    if (0) {
      {
        Vec f0 = rot0 * (domain_[p0]->momentum_matrix_ * subspace_force0);
        Vec f1 = rot1 * (domain_[p1]->momentum_matrix_ * subspace_force1);
        Vec diff = f0 + f1;
        ASSERT(diff.norm() < 1e-7, P(diff.norm(), f0.norm(), f1.norm()));
      }
      {
        Vec f0 = rot0 * (interface_torque_matrix_[e * 2 + 0] * subspace_force0);
        Vec f1 = rot1 * (interface_torque_matrix_[e * 2 + 1] * subspace_force1);
        Vec diff = f0 + f1;
        ASSERT(diff.norm() < 1e-7, P(diff.norm(), f0.norm(), f1.norm()));
        //      P(diff.norm());
      }
    }

    net_force[p0] += rot0 * (domain_[p0]->momentum_matrix_ * subspace_force0);
    net_torque[p0] += rot0 * (domain_[p0]->torque_matrix_ * subspace_force0);
    //    P(net_force[p0], net_torque[p0]);
    net_force[p1] += rot1 * (domain_[p1]->momentum_matrix_ * subspace_force1);
    net_torque[p1] += rot1 * (domain_[p1]->torque_matrix_ * subspace_force1);
  }

  Vec rigid_rhs = Vec::Zero(part_num_ * 6);
  Mat rigid_k = Mat::Zero(part_num_ * 6, part_num_ * 6);
  for (int p = 0; p < part_num_; ++p) {
    SubspaceTet* mesh = domain_[p];
    MapVec3 force_rhs(&rigid_rhs[p * 6 + 0]);
    MapVec3 torque_rhs(&rigid_rhs[p * 6 + 3]);
    //    net_force[p].setZero();
    force_rhs = net_force[p] * dt + mesh->total_mass_ * mesh->translation_vel_;
    Mat3 cur_inertial_tensor = mesh->rotation_ * mesh->inertial_tensor_ * mesh->rotation_.transpose();
    //    net_torque[p].setZero();
    torque_rhs = net_torque[p] * dt +  cur_inertial_tensor * mesh->angular_vel_;
  }

  auto Dr_Dw = [](Mat3 & rotation, Mat3 & dr_dx, Mat3 & dr_dy, Mat3 & dr_dz) {
    dr_dx(0, 0) = dr_dx(0, 1) = dr_dx(0, 2) = 0;
    dr_dx(1, 0) = -rotation(2, 0);
    dr_dx(1, 1) = -rotation(2, 1);
    dr_dx(1, 2) = -rotation(2, 2);
    dr_dx(2, 0) = +rotation(1, 0);
    dr_dx(2, 1) = +rotation(1, 1);
    dr_dx(2, 2) = +rotation(1, 2);

    dr_dy(1, 0) = dr_dy(1, 1) = dr_dy(1, 2) = 0;
    dr_dy(0, 0) = +rotation(2, 0);
    dr_dy(0, 1) = +rotation(2, 1);
    dr_dy(0, 2) = +rotation(2, 2);
    dr_dy(2, 0) = -rotation(0, 0);
    dr_dy(2, 1) = -rotation(0, 1);
    dr_dy(2, 2) = -rotation(0, 2);

    dr_dz(2, 0) = dr_dz(2, 1) = dr_dz(2, 2) = 0;
    dr_dz(0, 0) = -rotation(1, 0);
    dr_dz(0, 1) = -rotation(1, 1);
    dr_dz(0, 2) = -rotation(1, 2);
    dr_dz(1, 0) = +rotation(0, 0);
    dr_dz(1, 1) = +rotation(0, 1);
    dr_dz(1, 2) = +rotation(0, 2);
  };


  // TODO delete
  //  if (0)
  for (int e = 0; e < block_edge_num_; ++e) {
    int p0 = edge_list_[e * 2 + 0];
    int p1 = edge_list_[e * 2 + 1];

    Mat t_mat2_0 = Mat::Zero(domain_[p0]->basis_num_, 3);
    Mat t_mat2_1 = Mat::Zero(domain_[p1]->basis_num_, 3);

    Mat r_mat1_vec1_0_r0 = Mat::Zero(domain_[p0]->basis_num_, 3);
    Mat r_mat1_vec1_0_r1 = Mat::Zero(domain_[p0]->basis_num_, 3);
    Mat r_mat1_vec1_1_r0 = Mat::Zero(domain_[p1]->basis_num_, 3);
    Mat r_mat1_vec1_1_r1 = Mat::Zero(domain_[p1]->basis_num_, 3);


    Vec diff_mass_center = domain_[p1]->center_of_mass_ - domain_[p0]->center_of_mass_;

    Mat3& rot0 = domain_[p0]->rotation_;
    Mat3& rot1 = domain_[p1]->rotation_;
    Mat3 dr0_dx, dr0_dy, dr0_dz;
    Mat3 dr1_dx, dr1_dy, dr1_dz;
    Dr_Dw(rot0, dr0_dx, dr0_dy, dr0_dz);
    Dr_Dw(rot1, dr1_dx, dr1_dy, dr1_dz);

    Mat3 dr0_T_r1_dx0 = (dr0_dx.transpose() * domain_[p1]->rotation_);
    Mat3 dr0_T_r1_dy0 = (dr0_dy.transpose() * domain_[p1]->rotation_);
    Mat3 dr0_T_r1_dz0 = (dr0_dz.transpose() * domain_[p1]->rotation_);

    Mat3 dr0_T_r1_dx1 = (domain_[p0]->rotation_.transpose() * dr1_dx);
    Mat3 dr0_T_r1_dy1 = (domain_[p0]->rotation_.transpose() * dr1_dy);
    Mat3 dr0_T_r1_dz1 = (domain_[p0]->rotation_.transpose() * dr1_dz);

    Mat3 dr1_T_r0_dx0 = (domain_[p1]->rotation_.transpose() * dr0_dx);
    Mat3 dr1_T_r0_dy0 = (domain_[p1]->rotation_.transpose() * dr0_dy);
    Mat3 dr1_T_r0_dz0 = (domain_[p1]->rotation_.transpose() * dr0_dz);

    Mat3 dr1_T_r0_dx1 = (dr1_dx.transpose() * domain_[p0]->rotation_);
    Mat3 dr1_T_r0_dy1 = (dr1_dy.transpose() * domain_[p0]->rotation_);
    Mat3 dr1_T_r0_dz1 = (dr1_dz.transpose() * domain_[p0]->rotation_);

    // = \sum_{all v} R * (mat0*q + mat1*{R^T*R'}*q' + mat2*{R^T}*(-MassCenter+MassCenter') + vec0 + vec1*{R^T*R'})
    for (int row = 0, idx = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col, ++idx) {
        t_mat2_0 += mat2_[e * 18 + 0 + idx] * domain_[p0]->rotation_(col, row); // rotation tranpose
        t_mat2_1 += mat2_[e * 18 + 9 + idx] * domain_[p1]->rotation_(col, row); // rotation tranpose
#ifndef ASSUME_RIGID
        r_mat1_vec1_0_r0.col(0) += (mat1_[e * 18 + idx] * domain_[p0]->q_ + vec1_[e * 18 + idx]) * dr0_T_r1_dx0(row, col) +
                                   mat2_[e * 18 + idx] * (dr0_dx(col, row) * diff_mass_center); // tranpose

        r_mat1_vec1_0_r0.col(1) += (mat1_[e * 18 + idx] * domain_[p0]->q_ + vec1_[e * 18 + idx]) * dr0_T_r1_dy0(row, col) +
                                   mat2_[e * 18 + idx] * (dr0_dy(col, row) * diff_mass_center); // tranpose

        r_mat1_vec1_0_r0.col(2) += (mat1_[e * 18 + idx] * domain_[p0]->q_ + vec1_[e * 18 + idx]) * dr0_T_r1_dz0(row, col) +
                                   mat2_[e * 18 + idx] * (dr0_dz(col, row) * diff_mass_center); // tranpose

        r_mat1_vec1_0_r1.col(0) += (mat1_[e * 18 + idx] * domain_[p0]->q_ + vec1_[e * 18 + idx]) * dr0_T_r1_dx1(row, col);
        r_mat1_vec1_0_r1.col(1) += (mat1_[e * 18 + idx] * domain_[p0]->q_ + vec1_[e * 18 + idx]) * dr0_T_r1_dy1(row, col);
        r_mat1_vec1_0_r1.col(2) += (mat1_[e * 18 + idx] * domain_[p0]->q_ + vec1_[e * 18 + idx]) * dr0_T_r1_dz1(row, col);


        r_mat1_vec1_1_r0.col(0) += (mat1_[e * 18 + 9 + idx] * domain_[p1]->q_ + vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dx0(row, col);
        r_mat1_vec1_1_r0.col(1) += (mat1_[e * 18 + 9 + idx] * domain_[p1]->q_ + vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dy0(row, col);
        r_mat1_vec1_1_r0.col(2) += (mat1_[e * 18 + 9 + idx] * domain_[p1]->q_ + vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dz0(row, col);

        r_mat1_vec1_1_r1.col(0) += (mat1_[e * 18 + 9 + idx] * domain_[p1]->q_ + vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dx1(row, col) +
                                   mat2_[e * 18 + 9 + idx] * (-dr1_dx(col, row) * diff_mass_center);// transpose

        r_mat1_vec1_1_r1.col(1) += (mat1_[e * 18 + 9 + idx] * domain_[p1]->q_ + vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dy1(row, col) +
                                   mat2_[e * 18 + 9 + idx] * (-dr1_dy(col, row) * diff_mass_center);// transpose

        r_mat1_vec1_1_r1.col(2) += (mat1_[e * 18 + 9 + idx] * domain_[p1]->q_ + vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dz1(row, col) +
                                   mat2_[e * 18 + 9 + idx] * (-dr1_dz(col, row) * diff_mass_center);// transpose
#else
        r_mat1_vec1_0_r0.col(0) += (vec1_[e * 18 + idx]) * dr0_T_r1_dx0(row, col) +
                                   mat2_[e * 18 + idx] * (dr0_dx(col, row) * diff_mass_center); // tranpose

        r_mat1_vec1_0_r0.col(1) += (vec1_[e * 18 + idx]) * dr0_T_r1_dy0(row, col) +
                                   mat2_[e * 18 + idx] * (dr0_dy(col, row) * diff_mass_center); // tranpose

        r_mat1_vec1_0_r0.col(2) += (vec1_[e * 18 + idx]) * dr0_T_r1_dz0(row, col) +
                                   mat2_[e * 18 + idx] * (dr0_dz(col, row) * diff_mass_center); // tranpose

        r_mat1_vec1_0_r1.col(0) += (vec1_[e * 18 + idx]) * dr0_T_r1_dx1(row, col);
        r_mat1_vec1_0_r1.col(1) += (vec1_[e * 18 + idx]) * dr0_T_r1_dy1(row, col);
        r_mat1_vec1_0_r1.col(2) += (vec1_[e * 18 + idx]) * dr0_T_r1_dz1(row, col);


        r_mat1_vec1_1_r0.col(0) += (vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dx0(row, col);
        r_mat1_vec1_1_r0.col(1) += (vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dy0(row, col);
        r_mat1_vec1_1_r0.col(2) += (vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dz0(row, col);

        r_mat1_vec1_1_r1.col(0) += (vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dx1(row, col) +
                                   mat2_[e * 18 + 9 + idx] * (-dr1_dx(col, row) * diff_mass_center);// transpose

        r_mat1_vec1_1_r1.col(1) += (vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dy1(row, col) +
                                   mat2_[e * 18 + 9 + idx] * (-dr1_dy(col, row) * diff_mass_center);// transpose

        r_mat1_vec1_1_r1.col(2) += (vec1_[e * 18 + 9 + idx]) * dr1_T_r0_dz1(row, col) +
                                   mat2_[e * 18 + 9 + idx] * (-dr1_dz(col, row) * diff_mass_center);// transpose
#endif
      }
    }
    t_mat2_0 *= -1;
    for (int col = 0; col < 3; ++col) {
      {
        auto f0 = t_mat2_0.col(col);
        auto f1 = t_mat2_1.col(col);
        ProjectSubspaceForce(e, f0, f1);
      }
      {
        auto f0 = r_mat1_vec1_0_r0.col(col);
        auto f1 = r_mat1_vec1_1_r0.col(col);
        ProjectSubspaceForce(e, f0, f1);
      }
      {
        auto f0 = r_mat1_vec1_0_r1.col(col);
        auto f1 = r_mat1_vec1_1_r1.col(col);
        ProjectSubspaceForce(e, f0, f1);
      }
    }
    // rotation
    rigid_k.block<3, 3>(p0 * 6 + 0, p0 * 6 + 3) += rot0 * domain_[p0]->momentum_matrix_ * r_mat1_vec1_0_r0;
    rigid_k.block<3, 3>(p0 * 6 + 3, p0 * 6 + 3) += rot0 * domain_[p0]->torque_matrix_ * r_mat1_vec1_0_r0;

    rigid_k.block<3, 3>(p0 * 6 + 0, p1 * 6 + 3) += rot0 * domain_[p0]->momentum_matrix_ * r_mat1_vec1_0_r1;
    rigid_k.block<3, 3>(p0 * 6 + 3, p1 * 6 + 3) += rot0 * domain_[p0]->torque_matrix_ * r_mat1_vec1_0_r1;

    rigid_k.block<3, 3>(p1 * 6 + 0, p0 * 6 + 3) += rot1 * domain_[p1]->momentum_matrix_ * r_mat1_vec1_1_r0;
    rigid_k.block<3, 3>(p1 * 6 + 3, p0 * 6 + 3) += rot1 * domain_[p1]->torque_matrix_ * r_mat1_vec1_1_r0;

    rigid_k.block<3, 3>(p1 * 6 + 0, p1 * 6 + 3) += rot1 * domain_[p1]->momentum_matrix_ * r_mat1_vec1_1_r1;
    rigid_k.block<3, 3>(p1 * 6 + 3, p1 * 6 + 3) += rot1 * domain_[p1]->torque_matrix_ * r_mat1_vec1_1_r1;
    // translation
    rigid_k.block<3, 3>(p0 * 6 + 0, p0 * 6) += rot0 * domain_[p0]->momentum_matrix_ * t_mat2_0;
    rigid_k.block<3, 3>(p0 * 6 + 3, p0 * 6) += rot0 * domain_[p0]->torque_matrix_ * t_mat2_0;

    rigid_k.block<3, 3>(p0 * 6 + 0, p1 * 6) -= rot0 * domain_[p0]->momentum_matrix_ * t_mat2_0;
    rigid_k.block<3, 3>(p0 * 6 + 3, p1 * 6) -= rot0 * domain_[p0]->torque_matrix_ * t_mat2_0;

    rigid_k.block<3, 3>(p1 * 6 + 0, p0 * 6) += rot1 * domain_[p1]->momentum_matrix_ * t_mat2_1;
    rigid_k.block<3, 3>(p1 * 6 + 3, p0 * 6) += rot1 * domain_[p1]->torque_matrix_ * t_mat2_1;

    rigid_k.block<3, 3>(p1 * 6 + 0, p1 * 6) -= rot1 * domain_[p1]->momentum_matrix_ * t_mat2_1;
    rigid_k.block<3, 3>(p1 * 6 + 3, p1 * 6) -= rot1 * domain_[p1]->torque_matrix_ * t_mat2_1;
  }

  rigid_k *= -dt * dt;
  for (int p = 0; p < part_num_; ++p) {
    SubspaceTet* mesh = domain_[p];
    Mat3 cur_inertial_tensor = mesh->rotation_ * mesh->inertial_tensor_ * mesh->rotation_.transpose();
    rigid_k(p * 6 + 0, p * 6 + 0) += mesh->total_mass_;
    rigid_k(p * 6 + 1, p * 6 + 1) += mesh->total_mass_;
    rigid_k(p * 6 + 2, p * 6 + 2) += mesh->total_mass_;
    rigid_k.block<3, 3>(p * 6 + 3, p * 6 + 3) += cur_inertial_tensor;
  }
  //  PVEC(nt_force[0]);
  //  P(rigid_rhs.norm());
  //  P(rigid_k.norm());
  //  WriteEigenMatrixToMatlab(rigid_k, "/tmp/log/k");
  Vec new_rigid_vel = rigid_k.colPivHouseholderQr().solve(rigid_rhs);
  //    new_rigid_vel.setZero();
  if (0) {
    Vec diff = rigid_k * new_rigid_vel - rigid_rhs;
    P(diff.norm());
    P(new_rigid_vel.norm());
  }
  for (int p = 0; p < part_num_; ++p) {
    domain_[p]->translation_acc_ = (MapVec3(&new_rigid_vel[p * 6]) - domain_[p]->translation_vel_) / dt;
    domain_[p]->angular_acc_ = (MapVec3(&new_rigid_vel[p * 6 + 3]) - domain_[p]->angular_vel_) / dt;
    //    P(p, dj::Vec3d(domain_[p]->translation_acc_.data()), dj::Vec3d(domain_[p]->angular_acc_.data()));
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
#if 1
  Vec rhs = Vec::Zero(total_basis_num_);
  // collision force
  for (int i = 0; i < int(collision_force.size()); ++i) {
    int& p = std::get<0>(collision_force[i]);
    int& v = std::get<1>(collision_force[i]);
    Vec3 force = domain_[p]->rotation_.transpose() * std::get<2>(collision_force[i]);
    MapVec subspace_force(&rhs[basis_offset_[p]], domain_[p]->basis_num_);
    subspace_force += domain_[p]->vert_basis_transpose_[v] * force;
  }
  // gravity and fititious force
  for (int p = 0; p < part_num_; ++p) {
    domain_[p]->AddFictitiousForceAndGravity(&rhs[basis_offset_[p]]);
  }

  Mat k = Mat::Zero(total_basis_num_, total_basis_num_);
  // reduced internal force and diagonal terms of reduce k
  for (int p = 0; p < part_num_; ++p) {
    int basis_num = domain_[p]->basis_num_;
    domain_[p]->inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
    MatCol part_reduced_k;
    domain_[p]->GetReducedTangentStiffnessMatrix(part_reduced_k);
    k.block(basis_offset_[p], basis_offset_[p], basis_num, basis_num) += part_reduced_k;
    MapVec subspace_force(&rhs[basis_offset_[p]], domain_[p]->basis_num_);
    // vega internal force is in opposite direction of actual force
    subspace_force -= kInternalForceScaling *
                      (domain_[p]->basis_transpose_ * MapVec((double*) domain_[p]->inv_fem_->internal_force_, domain_[p]->vertex_num_ * 3));
    subspace_force *= dt;
    subspace_force += domain_[p]->vel_q_;
  }
  // off diagonal terms of reudced k
  // = \sum_{all v}  (mat0*q + mat1*{R^T*R'}*q' + mat2*{R^T}*(-MassCenter+MassCenter') + vec0 + vec1*{R^T*R'})
  for (int e = 0; e < block_edge_num_; ++e) {
    int p0 = edge_list_[e * 2 + 0];
    int p1 = edge_list_[e * 2 + 1];
    Mat mat1_0 = Mat::Zero(domain_[p0]->basis_num_, domain_[p1]->basis_num_);
    Mat mat1_1 = Mat::Zero(domain_[p1]->basis_num_, domain_[p0]->basis_num_);

    Mat3& rot0 = domain_[p0]->rotation_;
    Mat3& rot1 = domain_[p1]->rotation_;
    Mat3 r0_T_dot_r1 = rot0.transpose() * rot1;
    Mat3 r1_T_dot_r0 = r0_T_dot_r1.transpose();

    // F = (mat0*q + mat1*{R^T*R'}*q' + mat2*{R^T}*(-MassCenter+MassCenter') + vec0 + vec1*{R^T*R'})
    for (int row = 0, idx = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col, ++idx) {
        mat1_0 += mat1_[e * 18 + 0 + idx] * r0_T_dot_r1(row, col);
        mat1_1 += mat1_[e * 18 + 9 + idx] * r1_T_dot_r0(row, col);
      }
    }
    k.block(basis_offset_[p0], basis_offset_[p0], domain_[p0]->basis_num_, domain_[p0]->basis_num_) -= mat0_[e * 2 + 0];
    k.block(basis_offset_[p0], basis_offset_[p1], domain_[p0]->basis_num_, domain_[p1]->basis_num_) -= mat1_0;

    k.block(basis_offset_[p1], basis_offset_[p1], domain_[p1]->basis_num_, domain_[p1]->basis_num_) -= mat0_[e * 2 + 1];
    k.block(basis_offset_[p1], basis_offset_[p0], domain_[p1]->basis_num_, domain_[p0]->basis_num_) -= mat1_1;
  }
  k *= (kInternalForceScaling * dt * dt);
  for (int i = 0; i < total_basis_num_; ++i) {
    k(i, i) += 1;
  }
  Vec new_vel_q = k.colPivHouseholderQr().solve(rhs);
//  P(new_vel_q.norm());
  if (0) {
    Vec diff = k * new_vel_q - rhs;
    P(diff.norm());
  }
  //  new_vel_q.setZero();

  for (int p = 0; p < part_num_; ++p) {
    MapVec part_vel(&new_vel_q[basis_offset_[p]], domain_[p]->basis_num_);
    domain_[p]->vel_q_ = part_vel;
    domain_[p]->q_ += domain_[p]->vel_q_ * dt;
    domain_[p]->vel_q_ *= 0.96;
  }
#endif

  profiler.Start("update");
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    domain_[p]->UpdateRigidMotionAndLocalDeformation(dt);
  }
  profiler.End("update");
}

void CubicaTet::NextRenderMode() {
  for (SubspaceTet * tet : domain_container_) {
    tet->NextRenderMode();
  }
}

void CubicaTet::Render() {
  //  Super::Render(); return;
  for (SubspaceTet * tet : domain_container_) {
    tet->Render();
  }

  if (1) {
    glPointSize(8);
    glColor3f(1, 1, 1);
    glBegin(GL_POINTS);
    for (int p = 0; p < part_num_; ++p) {
      Vertex3v(&domain_[p]->center_of_mass_[0]);
    }
    glEnd();
  }
}

void CubicaTet::ComputeInterDomainSubspaceForce(int e, CubicaTet::Vec & force0, CubicaTet::Vec & force1) {
  int p0 = edge_list_[e * 2 + 0];
  int p1 = edge_list_[e * 2 + 1];
  int v_num = int(interface_vert_[e].size()) / 3;
  force0 = Vec::Zero(domain_[p0]->basis_num_);
  force1 = Vec::Zero(domain_[p1]->basis_num_);
  for (int i = 0; i < v_num; ++i) {
    const int global_v = interface_vert_[e][i * 3 + 0];
    const int v0 = interface_vert_[e][i * 3 + 1];
    const int v1 = interface_vert_[e][i * 3 + 2];
    MapVec3 pos0(domain_[p0]->X + v0 * 3);
    MapVec3 pos1(domain_[p1]->X + v1 * 3);
    Vec3 diff = vert_interfacial_area_[global_v] * (pos1 - pos0);
    force0 += domain_[p0]->vert_basis_transpose_[v0] * domain_[p0]->rotation_.transpose() * diff;
    //    Vec tmp = domain_[p0]->vert_basis_transpose_[v0] * domain_[p0]->rotation_.transpose() * diff;
    diff *= -1;
    force1 += domain_[p1]->vert_basis_transpose_[v1] * domain_[p1]->rotation_.transpose() * diff;
  }
  P(force0.norm(), force1.norm());
}

void CubicaTet::ComputeVertexInterfacialArea() {
  const double kStiffness = 1.0e7;
  vert_interfacial_area_ = std::vector<double>(vertex_num_, 0);
  for (int p = 0; p < part_num_; ++p) {
    for (int t = 0; t < domain_[p]->triangle_num_; ++t) {
      int* verts = domain_[p]->T + t * 3;
      bool is_surface_tri = true;
      for (int i = 0; i < 3; ++i) {
        if (!is_interface_vert_[p][verts[i]]) {
          is_surface_tri = false;
          break;
        }
      }
      if (is_surface_tri) {
        double area = dj::ComputeTriangleArea(domain_[p]->rest_pos_ + verts[0] * 3,
                                              domain_[p]->rest_pos_ + verts[1] * 3,
                                              domain_[p]->rest_pos_ + verts[2] * 3);
        area /= 3;
        for (int i = 0; i < 3; ++i) {
          int global_v = vert_local_id2global_id_[p][verts[i]];
          vert_interfacial_area_[global_v] += area;
        }
      }
    }
  }
  for (int v = 0; v < vertex_num_; ++v) {
    vert_interfacial_area_[v] /= 2;
    vert_interfacial_area_[v] *= kStiffness;
  }
}
