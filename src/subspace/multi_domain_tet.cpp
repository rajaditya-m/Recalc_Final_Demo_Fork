#include <Eigen/Geometry>
#include <fstream>
#include <set>
#include <unordered_set>
#include <cstdio>
#include "tet_gen_mesh_io.h"
#include "skeleton_node.h"
#include "skeleton.h"
#include "string_formatter.h"
#include "multi_domain_tet.h"
#include "conjugate_gradient_solver.h"
#include "opengl_helper.h"
#include "biconjugate_gradient_stablized.h"
#include "global.h"
#include "vector_io.h"
#include "print_macro.h"
#include "textured_triangle_mesh_renderer.h"
#include "sparseMatrix.h"
#include "tet_mesh_simulator_bridge.h"
#include "vector_io.h"
#include "config_file.h"
#include "matlab_io.h"
#include "BLOCK_MATRIX_GRAPH.h"
#include "affine_transformer.h"
#include "rainbow_color.h"
#include "open_gl_qt.h"

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// embeded mesh
#include "objMesh.h"
#include "objMeshRender.h"


MultiDomainTet::MultiDomainTet(const char *filename, double _limit_threshold, AffineTransformer<double> *affine_transformer)
  : Super(filename, _limit_threshold, affine_transformer, true) {
  ext_force_.reserve(vertex_num_ / 100);
  chol_solver_ = nullptr;
  full_chol_solver_ = nullptr;
  constrained_chol_solver_ = nullptr;
  gravity_[0] = global::gravity[0];
  gravity_[1] = global::gravity[1];
  gravity_[2] = global::gravity[2];
  P(gravity_);
  cross_product_matrix_.resize(vertex_num_);
  current_offset_from_mass_center_.resize(vertex_num_);
  local_u_ = std::vector<Vec3>(vertex_num_, Vec3(0, 0, 0));
  local_r_ = std::vector<Vec3>(vertex_num_, Vec3(0, 0, 0));
  cg_solver_ = new Solver(vertex_num_ * 3);
  P(surface_vertices_.size());
  obj_ = nullptr;
  obj_render_ = nullptr;
  fixed_domain_ = -1;
  full_rigid_simulation_ = false;
  fixed_domains_.clear();
  position_changed_ = true;
}

MultiDomainTet::~MultiDomainTet() {
  delete chol_solver_;
  delete full_chol_solver_;
}

void MultiDomainTet::UpdateCrossProductMatrix() {
  OMP_FOR
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    current_offset_from_mass_center_[v] = MapVec3(X + v * 3) - center_of_mass_[p];
    cross_product_matrix_[v] = GetSkewSymmetrixMatrix(current_offset_from_mass_center_[v]);
  }
}

void MultiDomainTet::EnableFullRigidSimulation() {
  std::unordered_set<int> tets(all_cubature_tet_.begin(), all_cubature_tet_.end());
  tets.insert(interface_tet_.begin(), interface_tet_.end());
  all_cubature_tet_.clear();
  all_cubature_tet_.insert(all_cubature_tet_.end(), tets.begin(), tets.end());
  full_rigid_simulation_ = true;
}

template <class Vector0, class Vector1>
void MultiDomainTet::ProjectInterfaceSubspaceForce(int e,
                                                   int p0, int p1,
                                                   Vector0& force0,
                                                   Vector1& force1) {
#if 1
  typedef Eigen::Matrix<double, 6, 1> Vec6;
  int e0 = e * 2 + 0;
  int e1 = e * 2 + 1;
  lagrangian_matrix_[e].block<3, 3>(0, 0) =
    part_rotation_[p0] * (momentum_matrix_[p0] * momentum_matrix_[p0].transpose()) * part_rotation_[p0].transpose() +
    part_rotation_[p1] * (momentum_matrix_[p1] * momentum_matrix_[p1].transpose()) * part_rotation_[p1].transpose();

  lagrangian_matrix_[e].block<3, 3>(0, 3) =
    part_rotation_[p0] * (momentum_matrix_[p0] * interface_torque_matrix_[e0].transpose()) * part_rotation_[p0].transpose() +
    part_rotation_[p1] * (momentum_matrix_[p1] * interface_torque_matrix_[e1].transpose()) * part_rotation_[p1].transpose();

  lagrangian_matrix_[e].block<3, 3>(3, 0) =
    part_rotation_[p0] * (interface_torque_matrix_[e0] * momentum_matrix_[p0].transpose()) * part_rotation_[p0].transpose() +
    part_rotation_[p1] * (interface_torque_matrix_[e1] * momentum_matrix_[p1].transpose()) * part_rotation_[p1].transpose();

  lagrangian_matrix_[e].block<3, 3>(3, 3) =
    part_rotation_[p0] * (interface_torque_matrix_[e0] * interface_torque_matrix_[e0].transpose()) * part_rotation_[p0].transpose() +
    part_rotation_[p1] * (interface_torque_matrix_[e1] * interface_torque_matrix_[e1].transpose()) * part_rotation_[p1].transpose();

  Vec6 residual;
  residual.block<3, 1>(0, 0) = part_rotation_[p0] * (momentum_matrix_[p0] * force0) +
                               part_rotation_[p1] * (momentum_matrix_[p1] * force1);
  residual.block<3, 1>(3, 0) = part_rotation_[p0] * (interface_torque_matrix_[e0] * force0) +
                               part_rotation_[p1] * (interface_torque_matrix_[e1] * force1);
  Vec6 lagranian_multiplier = lagrangian_matrix_[e].ldlt().solve(residual);
  MapVec3 lambda0(&lagranian_multiplier[0]);
  MapVec3 lambda1(&lagranian_multiplier[3]);
  force0 -=  (momentum_matrix_transpose_[p0] * (part_rotation_[p0].transpose() * lambda0)) +
             (interface_torque_matrix_tranpose_[e0] * (part_rotation_[p0].transpose() * lambda1));

  force1 -=  (momentum_matrix_transpose_[p1] * (part_rotation_[p1].transpose() * lambda0)) +
             (interface_torque_matrix_tranpose_[e1] * (part_rotation_[p1].transpose() * lambda1));
#else
  //  Mat3 lagrangian_matrix =
  //    part_rotation_[p0] * (momentum_matrix_[p0] * momentum_matrix_[p0].transpose()) * part_rotation_[p0].transpose() +
  //    part_rotation_[p1] * (momentum_matrix_[p1] * momentum_matrix_[p1].transpose()) * part_rotation_[p1].transpose();
  //
  //  Vec3 residual =
  //    part_rotation_[p0] * (momentum_matrix_[p0] * force0) + part_rotation_[p1] * (momentum_matrix_[p1] * force1);
  //  Vec3 lagranian_multiplier = lagrangian_matrix.ldlt().solve(residual);
  //  Vec3& lambda0 = lagranian_multiplier;
  //  force0 -=  (momentum_matrix_transpose_[p0] * (part_rotation_[p0].transpose() * lambda0));
  //  force1 -=  (momentum_matrix_transpose_[p1] * (part_rotation_[p1].transpose() * lambda0));
#endif
}

void MultiDomainTet::GenerateCubicaStylePartition(const char* folder) {
  std::vector<int> tet_partition_map(tet_number, -1); // tell which partition each tet goes to
  std::vector<std::vector<int> > tet_partition(part_num_); // list of global tet id for each partition
  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    std::set<int> partitions;
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      int p = vert_part_id_[v];
      partitions.insert(p);
    }
    ASSERT(partitions.size() <= 2, L("Assumption that a tet can only be shared by two partitions is not true!!!"));
    if (partitions.size() == 1) {
      int p = *(partitions.begin());
      tet_partition[p].push_back(t);
      tet_partition_map[t] = p;
    } else {
      int p0 = *(partitions.begin());
      int p1 = *(++partitions.begin());
      if (p0 > p1) std::swap(p0, p1);
      tet_partition[p0].push_back(t);
      tet_partition_map[t] = p0;
    }
  }

  // for each vert, (part_id0, local_vert_id0) + (part_id1, local_verid1)
  std::vector<std::pair<int, int> > vert_partition_map(vertex_num_ * 2, std::make_pair(-1, -1));
  // for each partition, list the four local vertex index
  std::vector<std::vector<int> > part_tets(part_num_);
  // list of vertex position for each partition
  std::vector<std::vector<double> > part_vert_pos(part_num_);
  // for each partition, list of vertex's global id
  std::vector<std::vector<int> > part_vert_global_id(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    part_vert_global_id[p].clear();
    std::vector<int> vert_local_id(vertex_num_, -1);
    int vert_count = 0;
    for (int t : tet_partition[p]) {
      int* verts = tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int v = verts[i];
        if (vert_local_id[v] == -1) {
          // meet new vertex
          vert_local_id[v] = vert_count;
          part_tets[p].push_back(vert_local_id[v]);
          part_vert_global_id[p].push_back(v);
          part_vert_pos[p].push_back(rest_pos_[v * 3 + 0]);
          part_vert_pos[p].push_back(rest_pos_[v * 3 + 1]);
          part_vert_pos[p].push_back(rest_pos_[v * 3 + 2]);
          if (vert_partition_map[v * 2 + 0].first == -1) {
            vert_partition_map[v * 2 + 0].first = p;
            vert_partition_map[v * 2 + 0].second = vert_local_id[v];
          } else {
            ASSERT(vert_partition_map[v * 2 + 1].first == -1, L("a vertex is put into more than two domains!!!"));
            vert_partition_map[v * 2 + 1].first = p;
            vert_partition_map[v * 2 + 1].second = vert_local_id[v];
          }
          vert_count++;
        } else {
          // old vertex
          part_tets[p].push_back(vert_local_id[v]);
        }
      }
    }
  }

  char file_name[512];
  for (int p = 0; p < part_num_; ++p) {
    sprintf(file_name, "%s/partition_%d", folder, domain_attached_bone_[p]);
    TetGenMeshIO::Instance()->Write(file_name, part_vert_pos[p], part_tets[p]);
    {
      sprintf(file_name, "%s/partition_%d.vert_info.txt", folder, domain_attached_bone_[p]);
      std::ofstream out(file_name);
      ASSERT(out.is_open(), P(file_name));
      out << part_vert_global_id[p].size() << "\n";
      bool has_constraied_vertex = false;
      for (int i = 0; i < int(part_vert_global_id[p].size()); ++i) {
        int global_v = part_vert_global_id[p][i];
        out << global_v << " ";
        int constrained = (vert_basis_[global_v].norm() == 0) ? 1 : 0;
        if (constrained) has_constraied_vertex = true;
        out << constrained << " ";
        if (vert_partition_map[global_v * 2 + 1].first == -1) {
          // vertex only belong to one partition
          out << "-1 -1\n";
        } else {
          if (vert_partition_map[global_v * 2 + 0].first == p) {
            out << vert_partition_map[global_v * 2 + 1].first << " ";
            out << vert_partition_map[global_v * 2 + 1].second << "\n";
          } else {
            ASSERT(vert_partition_map[global_v * 2 + 1].first == p, L("can't find vertex in other partition"));
            out << vert_partition_map[global_v * 2 + 0].first << " ";
            out << vert_partition_map[global_v * 2 + 0].second << "\n";
          }
        }
      }
      out.close();
      ASSERT(has_constraied_vertex, P(p));
    }
    {
      sprintf(file_name, "%s/partition_%d.tet_info.txt", folder, domain_attached_bone_[p]);
      std::ofstream out(file_name);
      ASSERT(out.is_open(), P(file_name));
      int part_tet_num = int(tet_partition[p].size());
      out << part_tet_num << "\n";
      for (int local_t = 0; local_t < part_tet_num; ++local_t) {
        out << tet_partition[p][local_t] << "\n";
      }
      out.close();
    }
  }

  sprintf(file_name, "%s/tet_partition.txt", folder);
  {
    std::ofstream out(file_name);
    ASSERT(out.is_open(), P(file_name));
    out << tet_number << "\n";
    for (int t = 0; t < tet_number; ++t) {
      out << tet_partition_map[t] << "\n";
    }
    out.close();
  }

  {
    sprintf(file_name, "%s/vert_partition.txt", folder);
    std::ofstream out(file_name);
    ASSERT(out.is_open(), P(file_name));
    out << vertex_num_ << "\n";
    for (int v = 0; v < vertex_num_; ++v) {
      out << vert_partition_map[v * 2 + 0].first << " ";
      out << vert_partition_map[v * 2 + 0].second << " ";
      out << vert_partition_map[v * 2 + 1].first << " ";
      out << vert_partition_map[v * 2 + 1].second << "\n";
    }
    out.close();
  }
  {
    sprintf(file_name, "%s/partition_id.txt", folder);
    std::ofstream out(file_name);
    ASSERT(out.is_open(), P(file_name));
    out << part_num_ << "\n";
    for (int p = 0; p < part_num_; ++p) {
      out << domain_attached_bone_[p] << std::endl;
    }
    out.close();
  }
}

void MultiDomainTet::ComputeLagradianMatrix() {
  lagrangian_matrix_.resize(chol_solver_->e_number0);
  for (int e = 0; e < chol_solver_->e_number0; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    lagrangian_matrix_[e].block<3, 3>(0, 0) = momentum_matrix_[p0] * momentum_matrix_[p0].transpose() +
                                              momentum_matrix_[p1] * momentum_matrix_[p1].transpose();

    lagrangian_matrix_[e].block<3, 3>(0, 3) = momentum_matrix_[p0] * torque_matrix_[p0].transpose() +
                                              momentum_matrix_[p1] * torque_matrix_[p1].transpose();

    lagrangian_matrix_[e].block<3, 3>(3, 0) = torque_matrix_[p0] * momentum_matrix_[p0].transpose() +
                                              torque_matrix_[p1] * momentum_matrix_[p1].transpose();

    lagrangian_matrix_[e].block<3, 3>(3, 3) = torque_matrix_[p0] * torque_matrix_[p0].transpose() +
                                              torque_matrix_[p1] * torque_matrix_[p1].transpose();
  }
}

void MultiDomainTet::ComputeVertexBasis() {
  vert_basis_.resize(vertex_num_);
  vert_basis_transpose_.resize(vertex_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    int local_v = local_vert_idx_[v];
    vert_basis_[v] = basis_[p].middleRows<3>(local_v * 3);
    vert_basis_transpose_[v] = vert_basis_[v].transpose();
  }

  {
#if 0
    // TODO: comment out for production
    // verify
    for (int v = 0; v < vertex_num_; ++v) {
      int p = vert_compact_part_id_[v];
      int local_v = local_vert_idx_[v];
      for (int i = 0; i < part_basis_size_[p]; ++i) {
        ASSERT(vert_basis_[v](0, i) == basis_[p](local_v * 3 + 0, i));
        ASSERT(vert_basis_[v](1, i) == basis_[p](local_v * 3 + 1, i));
        ASSERT(vert_basis_[v](2, i) == basis_[p](local_v * 3 + 2, i));
      }
    }
#endif
  }
}

/// q = U^T * M * (x - X)
void MultiDomainTet::ProjectPositionToSubspace() {
  MapVec map_x(X, vertex_num_ * 3);
  dj::Swap(X, tmp_X);
  skeleton_->AssembleSkeletonTransformation();
  ApplySkeletonTransformationToVertex();
  dj::Swap(X, tmp_X);
  MapVec map_X(tmp_X, vertex_num_ * 3);
  Vec u = map_x - map_X;
  for (int v = 0; v < vertex_num_; ++v) {
    u[v * 3 + 0] *= mass_[v];
    u[v * 3 + 1] *= mass_[v];
    u[v * 3 + 2] *= mass_[v];
  }
  q_ = global_basis_transpose_ * u;
  MapVec dx((double*)inv_fem_->u_, vertex_num_ * 3);
  MapVec rest_pos(rest_pos_, vertex_num_ * 3);
  dx = map_x - rest_pos;
}

void MultiDomainTet::LoadLocalMass(const char *prefix) {
  for (int p = 0; p < part_num_; ++p) {
    char file_name[512];
    sprintf(file_name, "%s_%d.M", prefix, domain_attached_bone_[p]);
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    for (int v = 0; v < vert_num_per_part_[p]; ++v) {
      int global_v = vert_local_id2global_id_[p][v];
      in >> mass_[global_v];
    }
    in.close();
  }
}

void MultiDomainTet::LoadSubspace(std::function<const char* (int)> GetFileName, int basis_file_format) {
  basis_.resize(part_num_);
  part_basis_size_.resize(part_num_);
  total_basis_num_ = 0;
  basis_offset_.resize(part_num_ + 1);
  basis_offset_[0] = 0;
  for (int p = 0; p < part_num_; ++p) {
    //        sprintf(file_name, "%s/%d/asciiBasis_%d.txt", prefix, parts_[p], parts_[p]);
    //        sprintf(file_name, "%s/%d/ascii_Basis_%d.txt", prefix, parts_[p], parts_[p]);
    //        sprintf(file_name, "%s/%d/FromGlobal/asciiMat_%d_2.txt", prefix, parts_[p], parts_[p]);
    //    sprintf(file_name, "%s/%d/asciiMat_%d.txt", prefix, parts_[p], parts_[p]);
    const char* file_name = GetFileName(p);
    int row_num, basis_num;
    if (basis_file_format == kText) {
      std::ifstream in(file_name);
      ASSERT(in.is_open(), P(file_name));
      in >> row_num >> basis_num;
      ASSERT(row_num == vert_num_per_part_[p] * 3);
      basis_[p] = Mat(row_num, basis_num);
      for (int i = 0; i < row_num; ++i) {
        for (int j = 0; j < basis_num; ++j) {
          in >> basis_[p](i, j);
        }
      }
      in.close();
    } else {
      std::ifstream in(file_name, std::ios::binary);
      ASSERT(in.is_open(), P(file_name));
      in.read((char*) &row_num, sizeof(int));
      ASSERT(row_num == vert_num_per_part_[p] * 3);
      in.read((char*) &basis_num, sizeof(int));
      MatCol tmp_basis(row_num, basis_num);
      in.read((char*) tmp_basis.data(), sizeof(double) * row_num * basis_num);
      basis_[p] = tmp_basis;
      in.close();
    }
    L(file_name);
    P(p, row_num, basis_num);
    part_basis_size_[p] = basis_num;
    total_basis_num_ += basis_num;
    basis_offset_[p + 1] = total_basis_num_;
  }

  rigid_basis_.resize(part_num_);
  rigid_basis_transpose_.resize(part_num_);
  vert_rigid_basis_.resize(vertex_num_);
  vert_rigid_basis_transpose_.resize(vertex_num_);
  for (int p = 0; p < part_num_; ++p) {
    std::string file_name = dj::Format("%z/modal_basis/partition_%d.rigid_basis.bin", GetDataFolder(), p);
    int row_num, basis_num;
    std::ifstream in(file_name, std::ios::binary);
    ASSERT(in.is_open(), P(file_name));
    in.read((char*) &row_num, sizeof(int));
    ASSERT(row_num == vert_num_per_part_[p] * 3);
    in.read((char*) &basis_num, sizeof(int));
    MatCol tmp_basis(row_num, basis_num);
    in.read((char*) tmp_basis.data(), sizeof(double) * row_num * basis_num);
    rigid_basis_[p] = tmp_basis;
    rigid_basis_transpose_[p] = rigid_basis_[p].transpose();
    in.close();

    for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
      int v = vert_local_id2global_id_[p][local_v];
      vert_rigid_basis_[v] = rigid_basis_[p].middleRows<3>(local_v * 3);
      vert_rigid_basis_transpose_[v] = vert_rigid_basis_[v].transpose();
    }
  }

  q_ = Vec::Zero(total_basis_num_);
  vel_q_ = Vec::Zero(total_basis_num_);

  BuildTopology();

  //  for (int p0 = 0; p0 < 3; ++p0) {
  //    P(topology_[p0].size(), p0);
  //    for (int p1 = 0; p1 < 3; ++p1) {
  //      P(p0, p1, topology_[p0][p1]);
  //    }
  //  }
  auto tmp_topology = topology_;
  chol_solver_ = new solver::BLOCK_MATRIX_GRAPH<double>(topology_, part_basis_size_);
  full_chol_solver_ = new solver::BLOCK_MATRIX_GRAPH<double>(tmp_topology, part_basis_size_, 6);
  //  for (int p0 = 0; p0 < 3; ++p0) {
  //    P(topology_[p0].size(), p0);
  //    for (int p1 = 0; p1 < 3; ++p1) {
  //      P(p0, p1, topology_[p0][p1]);
  //    }
  //  }
  //  P(solver_->e_number);
  //  exit(0);

  tmp_v_matrix_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    tmp_v_matrix_[p] = Mat::Zero(basis_[p].rows(), basis_[p].cols());
  }

  tmp_e_matrix_.resize(chol_solver_->e_number);
  for (int e = 0; e < chol_solver_->e_number; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    tmp_e_matrix_[e] = Mat::Zero(basis_[p0].rows(), basis_[p1].cols());
  }

  basis_transpose_ = basis_;
  for (int i = 0; i < (int)basis_transpose_.size(); ++i) {
    basis_transpose_[i].transposeInPlace();
  }
  VerifySubspaceDiagonality();
  //  RemoveConstrainedDofFromSubspace();
  basis_transpose_ = basis_;
  //  basis_transpose_.resize(part_num_);
  for (int i = 0; i < (int)basis_transpose_.size(); ++i) {
    //    basis_transpose_[i] = basis_[i].transpose();
    basis_transpose_[i].transposeInPlace();
  }
  initial_basis_ = basis_;
  initial_basis_transpose_ = basis_transpose_;
  ComputeMomentumAndTorqueMatrix();
  ComputeVertexBasis();
  ComputeLagradianMatrix();
  if (0) {
    int max_idx = -1;
    double max_diff = 0;
    Mat max_mat;
    int count = 0;
    Mat avg = Mat::Zero(3, 3);
    for (int v = 0; v < vertex_num_; ++v) {
      if (vert_basis_transpose_[v].norm() < 1e-10) continue;
      int p = vert_part_id_[v];
      Mat r = momentum_matrix_[p] * vert_basis_transpose_[v];
      //      PMAT(r);
      r(0, 0) -= 1;
      r(1, 1) -= 1;
      r(2, 2) -= 1;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          r(i, j) = dj::Abs(r(i, j));
        }
      }
      avg += r;
      count++;
      double max = dj::Max(dj::Abs(r.maxCoeff()), dj::Abs(r.minCoeff()));
      if (max > max_diff) {
        max_diff = max;
        max_idx = v;
        max_mat = r;
      }
      //      break;
    }
    P(max_diff, max_idx);
    avg /= count;
    PMAT(max_mat);
    PMAT(avg);
    exit(0);
  }
}

void MultiDomainTet::LoadPartitionInfo(std::string partition_folder) {
  LoadPartitionInfo(partition_folder.c_str());
}

/// Read information about how mesh vertex is decomposed into different parts
void MultiDomainTet::LoadPartitionInfo(const char* partition_folder) {
  // list of partition ids
  {
    std::string vert_part_info_name = dj::Format("%z/vertex_partition_info.txt", partition_folder);
    std::ifstream in(vert_part_info_name);
    ASSERT(in.is_open(), P(vert_part_info_name));
    int v_num;
    in >> v_num;
    ASSERT(v_num == vertex_num_);
    local_vert_idx_.resize(vertex_num_);
    vert_bone_id_.resize(vertex_num_);
    vert_part_id_.resize(vertex_num_);
    part_num_ = 0;
    // partition info for each vertex
    for (int v = 0; v < vertex_num_; ++v) {
      int local_id, compact_partition_id;
      in >> compact_partition_id >> local_id;
      part_num_ = std::max(part_num_, compact_partition_id);
      vert_part_id_[v] = compact_partition_id;
      local_vert_idx_[v] = local_id;
      compact_partition_id = vert_part_id_[v];
      ASSERT(local_id >= 0);
    }
    in.close();
  }

  part_num_++;
  part_translation_.resize(part_num_);
  part_rotation_ = std::vector<Mat3>(part_num_, Mat3::Identity());
  part_rotation_transpose_ = std::vector<Mat3>(part_num_, Mat3::Identity());
  vert_num_per_part_ = std::vector<int>(part_num_, 0);
  mass_per_part_ = std::vector<double>(part_num_, 0);
  for (int v = 0; v < vertex_num_; ++v) {
    int compact_partition_id = vert_part_id_[v];
    ++vert_num_per_part_[compact_partition_id];
    mass_per_part_[compact_partition_id] += mass_[v];
  }

  std::string domain_bone_id_file = dj::Format("%z/domain_bone_id.txt", partition_folder) ;
  std::ifstream in(domain_bone_id_file);
  if (in.is_open()) {
    int part_num;
    in >> part_num;
    ASSERT(part_num == part_num_);
    domain_attached_bone_.resize(part_num_);
    for (int i = 0; i < part_num_; ++i) {
      in >> domain_attached_bone_[i];
      bone_id2domain_id_[domain_attached_bone_[i]] = i;
    }
    for (int v = 0; v < vertex_num_; ++v) {
      int bone_id = domain_attached_bone_[vert_part_id_[v]];
      vert_bone_id_[v] = bone_id;
      if (skeleton_) {
        closest_bone_[v] = bone_id;
        skeleton_->AttachedToBone(X + v * 3, closest_bone_[v], initial_offset_ + v * 3);
      }
    }
    in.close();
  }
  //  for (int i = 0; i < part_num_; ++i) {
  //    P(i, vert_num_per_part_[i]);
  //  }
  //  exit(0);

  vert_local_id2global_id_.resize(part_num_);
  for (int i = 0; i < (int)vert_local_id2global_id_.size(); ++i) {
    vert_local_id2global_id_[i].resize(vert_num_per_part_[i]);
  }

  for (int v = 0; v < vertex_num_; ++v) {
    int compact_part_id = vert_part_id_[v];
    int local_id = local_vert_idx_[v];
    ASSERT(local_id >= 0 && local_id < vert_num_per_part_[compact_part_id],
           P(local_id, vert_num_per_part_[compact_part_id]));
    vert_local_id2global_id_[compact_part_id][local_id] = v;
  }

  // read interface domains
  {
    std::string file = dj::Format("%z/interface_domains.txt", partition_folder);
    std::ifstream in(file.c_str());
    ASSERT(in.is_open(), P(file));
    int p_num;
    in >> p_num;
    ASSERT(p_num == part_num_);
    in >> interface_num_;
    interface_domains_.resize(interface_num_);
    for (int i = 0; i < int(interface_domains_.size()); ++i) {
      interface_domains_[i].fill(-1);
      interface_domains_[i][0] = 0;
      int dummy = 0;
      in >> dummy;
      for (int j = 0; j < 4; ++j) {
        in >> interface_domains_[i][j + 1];
        if (interface_domains_[i][j + 1] != -1) {
          interface_domains_[i][0]++;
        }
      }
    }
    in.close();
  }
  InitializePartAffineTransformation();
  InitializeMassCenterOffset();
  InitInterfaceTetrahedra();
  for (int p = 0; p < part_num_; ++p) {
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      local_r_[v] = MapVec3(rest_pos_ + v * 3) - center_of_mass_[p];
    }
  }
  //  ComputePartAffineTransformamtion();
  delete cg_solver_;
  // 3 dof from translation and 3 dof from rotation for each partition
  cg_solver_ = new Solver(vertex_num_ * 3);// + part_num_ * 6);

  if (0)
    for (int p = 0; p < part_num_; ++p) {
      std::ofstream out(dj::Format("%z/partition_%d.mass.txt", partition_folder, p));
      out << vert_num_per_part_[p] << std::endl;
      for (int i = 0; i < vert_num_per_part_[p]; ++i) {
        int global_v = vert_local_id2global_id_[p][i];
        out << mass_[global_v] << "\n";
      }
      out.close();
    }
  acceleration_ = std::vector<Vec3>(part_num_, Vec3(0, 0, 0));
  angular_acceleration_ = std::vector<Vec3>(part_num_, Vec3(0, 0, 0));
  UpdateCrossProductMatrix();
}

void MultiDomainTet::VerifySubspaceDiagonality() {
  for (int p = 0; p < part_num_; ++p) {
    Mat m_dot_u = basis_[p];
    ASSERT(vert_num_per_part_[p] * 3 == basis_[p].rows());
    for (int v = 0; v < (int)vert_num_per_part_[p]; ++v) {
      int global_v = vert_local_id2global_id_[p][v];
      ASSERT(local_vert_idx_[global_v] == v);
      for (int i = 0; i < basis_[p].cols(); ++i) {
        ASSERT(basis_transpose_[p](i, v) == basis_[p](v, i));
        m_dot_u(v * 3 + 0, i) *= mass_[global_v];
        m_dot_u(v * 3 + 1, i) *= mass_[global_v];
        m_dot_u(v * 3 + 2, i) *= mass_[global_v];
      }
      //      m_dot_u.row(v * 3 + 0) *= mass_[global_v];
      //      m_dot_u.row(v * 3 + 1) *= mass_[global_v];
      //      m_dot_u.row(v * 3 + 2) *= mass_[global_v];
    }

    if (0 && p == 9) {
      std::cout << basis_[p].row(0) << std::endl;
      std::cout << m_dot_u.row(0) << std::endl;
      P(vert_local_id2global_id_[p][0]);
      P(mass_[vert_local_id2global_id_[p][0]], domain_attached_bone_[p], basis_[p].rows());
      exit(0);
    }
    Mat zero = basis_transpose_[p] * m_dot_u - Mat::Identity(part_basis_size_[p], part_basis_size_[p]);
    double max = 0;
    for (int i = 0; i < zero.rows(); ++i) {
      for (int j = 0; j < zero.cols(); ++j) {
        max = dj::Max(max, dj::Abs(zero(i, j)));
      }
    }
    //    auto max = std::max(dj::Abs(zero.minCoeff()), dj::Abs(zero.maxCoeff()));
    //    P(domain_attached_bone_[p], max);
    //    if (p == 0) {
    //      std::ofstream out("/tmp/mass.txt");
    //      for (int v : vert_local_id2global_id_[p]) {
    //        out << mass_[v] << std::endl;
    //      }
    //      out.close();
    //      L("finish");
    //      exit(0);
    //    }
    ASSERT(max < 1e-5, P(p, max));
  }
  L("basis mass diagonality verified");
}

void MultiDomainTet::RemoveConstrainedDofFromSubspace() {
  //  return;
  for (int v = 0; v < vertex_num_; ++v) {
    if (is_constrainted_[v]) {
      int p = vert_part_id_[v];
      int local_v = local_vert_idx_[v];
      //      P(p, local_v, basis_[p].rows(), basis_[p].cols());
      for (int i = 0; i < basis_[p].cols(); ++i) {
        basis_[p](local_v * 3 + 0, i) = 0;
        basis_[p](local_v * 3 + 1, i) = 0;
        basis_[p](local_v * 3 + 2, i) = 0;
      }
    }
  }
}

void MultiDomainTet::ComputeInterfaceTorqueMatrix() {
  std::vector<MatCol>& basis = rigid_basis_;
  auto GetInterfaceTorqueMatrix = [&](int p, Vec3 & center_of_mass, Mat & matrix) {
    matrix = Mat::Zero(3, basis[p].cols());
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      double mass = mass_[v];
      for (int c = 0; c < basis[p].cols(); ++c) {
        Vec3 col(basis[p](i * 3 + 0, c), basis[p](i * 3 + 1, c), basis[p](i * 3 + 2, c));
        col *= mass;
        MapVec3 pos(rest_pos_ + v * 3);
        Vec3 r = pos - center_of_mass;
        Vec3 cross = r.cross(col);
        matrix(0, c) += cross[0];
        matrix(1, c) += cross[1];
        matrix(2, c) += cross[2];
      }
    }
  };

  interface_torque_matrix_.resize(interface_num_ * 2);
  interface_torque_matrix_tranpose_.resize(interface_num_ * 2);
  for (int e = 0; e < interface_num_; ++e) {
    std::unordered_set<int> verts;
    for (CubaturePoint & cubature : interface_cubature_[e]) {
      int t = cubature.first;
      verts.insert(tet_[t * 4 + 0]);
      verts.insert(tet_[t * 4 + 1]);
      verts.insert(tet_[t * 4 + 2]);
      verts.insert(tet_[t * 4 + 3]);
    }

    Vec3 interface_center_of_mass(0, 0, 0);
    double total_mass = 0;
    for (int v : verts) {
      MapVec3 pos(rest_pos_ + v * 3);
      interface_center_of_mass += (mass_[v] * pos);
      total_mass += mass_[v];
    }
    interface_center_of_mass /= total_mass;

    int p0 = interface_domains_[e][1];
    int p1 = interface_domains_[e][2];
    GetInterfaceTorqueMatrix(p0, interface_center_of_mass, interface_torque_matrix_[e * 2 + 0]);
    GetInterfaceTorqueMatrix(p1, interface_center_of_mass, interface_torque_matrix_[e * 2 + 1]);
    interface_torque_matrix_tranpose_[e * 2 + 0] = interface_torque_matrix_[e * 2 + 0].transpose();
    interface_torque_matrix_tranpose_[e * 2 + 1] = interface_torque_matrix_[e * 2 + 1].transpose();
  }
}

void MultiDomainTet::ComputeMomentumAndTorqueMatrix() {
  // momemtum matirx = \sum_{all vertex v} m_v * U_v, (U_v are the 3 rows in basis that correspond to the vertex)
  // torque matirx = \sum_{all vertex v} r.cross(m_v * U_v), r = vertex_rest_position - mass_center
  //    std::vector<Mat>& basis = basis_;
  std::vector<MatCol>& basis = rigid_basis_;
  momentum_matrix_.resize(part_num_);
  torque_matrix_.resize(part_num_);
  momentum_matrix_transpose_.resize(part_num_);
  torque_matrix_transpose_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    momentum_matrix_[p] = Mat::Zero(3, basis[p].cols());
    torque_matrix_[p] = Mat::Zero(3, basis[p].cols());
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      double mass = mass_[v];
      for (int c = 0; c < basis[p].cols(); ++c) {
        for (int r = 0; r < 3; ++r) {
          momentum_matrix_[p](r, c) += mass * basis[p](i * 3 + r, c);
        }
        Vec3 col(basis[p](i * 3 + 0, c), basis[p](i * 3 + 1, c), basis[p](i * 3 + 2, c));
        col *= mass;
        Vec3& r = vert_offset_from_mass_center_[v];
        Vec3 cross = r.cross(col);
        torque_matrix_[p](0, c) += cross[0];
        torque_matrix_[p](1, c) += cross[1];
        torque_matrix_[p](2, c) += cross[2];
      }
    }
    momentum_matrix_transpose_[p] = momentum_matrix_[p].transpose();
    torque_matrix_transpose_[p] = torque_matrix_[p].transpose();
    P(momentum_matrix_[p].norm(), torque_matrix_[p].norm(), p);
    //    PMAT(momentum_matrix_[p]);
    //    PMAT(torque_matrix_[p]);
  }
}

void MultiDomainTet::ComputePartAffineTransformamtion() {
  skeleton_->AssembleSkeletonTransformation();
  double* skeleton_transformation = &skeleton_->skeleton_transformation_[0];
  for (int p = 0; p < part_num_; ++p) {
    part_translation_[p] = initial_center_of_mass_[p];
    double* pos = &part_translation_[p][0];
    int attached_bone = domain_attached_bone_[p];
    dj::MulMatrix3x3Vec<double>((double(*)[3]) (skeleton_transformation + attached_bone * 12),
                                &mass_center_initial_offset_[p][0], pos);
    pos[0] += skeleton_transformation[attached_bone * 12 + 9];
    pos[1] += skeleton_transformation[attached_bone * 12 + 10];
    pos[2] += skeleton_transformation[attached_bone * 12 + 11];
    center_of_mass_[p] = part_translation_[p];
    part_translation_[p] -= initial_center_of_mass_[p];
    // rotation
    // R^T
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> current_rotation(skeleton_transformation + attached_bone * 12);
    //    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> init_rotation(&skeleton_->bone_rotation_[attached_bone * 9]);
    part_rotation_[p] = current_rotation * initial_rotation_transpose_[p];
    //    P(p);
    //    PMAT(init_rotation);
    //    PMAT(initial_rotation_[p]);
  }

  // verify
  {
#if 0
    std::vector<double> tmp_pos(vertex_num_ * 3);
    memcpy(&tmp_pos[0], X, sizeof(double) * 3 * vertex_num_);
    ApplySkeletonTransformationToVertex();

    for (int v = 0; v < vertex_num_; ++v) {
      Vec3 pos(rest_pos_ + v * 3);
      int part = vert_compact_part_id_[v];
      //      Vec3 cur_pos = part_rotation_[part] * (pos - initial_center_of_mass_[part]) + part_translation_[part] + initial_center_of_mass_[part];
      //      Vec3 cur_pos = part_rotation_[part] * (pos - initial_center_of_mass_[part]) + current_center_of_mass_[part];
      Vec3 cur_pos = part_rotation_[part] * (pos - initial_center_of_mass_[part]) + current_center_of_mass_[part];
      double diff = dj::Distance3(&cur_pos[0], X + v * 3);
      ASSERT(diff < 1e-5, P(diff, v) P(part_rotation_[part], part_translation_[part]));
      //      X[v * 3 + 0] = cur_pos[0];
      //      X[v * 3 + 1] = cur_pos[1];
      //      X[v * 3 + 2] = cur_pos[2];
    }
    memcpy(X, &tmp_pos[0], sizeof(double) * 3 * vertex_num_);
    L("bone rotation and translation verfied");
#endif
  }
}

void MultiDomainTet::InitializePartAffineTransformation() {
  // translate (center of mass)
  initial_center_of_mass_.resize(part_num_);
  std::fill(initial_center_of_mass_.begin(), initial_center_of_mass_.end(), Vec3(0, 0, 0));
  for (int v = 0; v < vertex_num_; ++v) {
    int part = vert_part_id_[v];
    initial_center_of_mass_[part][0] += rest_pos_[v * 3 + 0] * mass_[v];
    initial_center_of_mass_[part][1] += rest_pos_[v * 3 + 1] * mass_[v];
    initial_center_of_mass_[part][2] += rest_pos_[v * 3 + 2] * mass_[v];
  }

  for (int p = 0; p < part_num_; ++p) {
    initial_center_of_mass_[p] /= mass_per_part_[p];
  }
  center_of_mass_ = initial_center_of_mass_;

  vert_offset_from_mass_center_.resize(vertex_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    MapVec3 init_pos(rest_pos_ + v * 3);
    vert_offset_from_mass_center_[v] = init_pos - initial_center_of_mass_[p];
  }

  // rotation
  if (skeleton_) {
    initial_rotation_transpose_.resize(part_num_);
    skeleton_->AssembleSkeletonTransformation();
    double* skeleton_transformation = &skeleton_->skeleton_transformation_[0];
    for (int p = 0; p < part_num_; ++p) {
      int attached_bone = domain_attached_bone_[p];
      initial_rotation_transpose_[p] = Mat3(skeleton_transformation + attached_bone * 12);
      initial_rotation_transpose_[p].transposeInPlace();
    }
  }

  for (int p = 0; p < part_num_; ++p) {
    part_rotation_[p] = Mat3::Identity();
  }

  angular_vel_ = std::vector<Vec3>(part_num_, Vec3(0, 0, 0));
  translational_vel_ = std::vector<Vec3>(part_num_, Vec3(0, 0, 0));
  quaternion_ = std::vector<Quaternion<double> >(part_num_, Quaternion<double>(1, 0, 0, 0));
  ComputeInvInertiaTensor();
}

void MultiDomainTet::ComputeMassCenterPos() {
  double* skeleton_transformation = &skeleton_->skeleton_transformation_[0];
  for (int p = 0; p < part_num_; ++p) {
    int attached_bone = domain_attached_bone_[p];
    MapMat3 rotation(skeleton_transformation + attached_bone * 12);
    MapVec3 offset(skeleton_transformation + attached_bone * 12 + 9);
    center_of_mass_[p] = rotation * mass_center_initial_offset_[p] + offset;
  }
}

void MultiDomainTet::InitializeMassCenterOffset() {
  if (skeleton_) {
    mass_center_initial_offset_.resize(part_num_);
    for (int i = 0; i < part_num_; ++i) {
      skeleton_->AttachedToBone(&initial_center_of_mass_[i][0], domain_attached_bone_[i], &mass_center_initial_offset_[i][0]);
    }
  }
}

void MultiDomainTet::ComputeGlobalPositionFromSubspace() {
  //  ComputePartAffineTransformamtion();
  OMP_FOR
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    int local_v = local_vert_idx_[v];
    Vec3 local_u(0, 0, 0);
    for (int r = 0; r < part_basis_size_[p]; ++r) {
      local_u[0] += q_[basis_offset_[p] + r] * basis_[p](local_v * 3 + 0, r);
      local_u[1] += q_[basis_offset_[p] + r] * basis_[p](local_v * 3 + 1, r);
      local_u[2] += q_[basis_offset_[p] + r] * basis_[p](local_v * 3 + 2, r);
    }
    Mat3& rot = part_rotation_[p];
    MapVec3 map_x(X + v * 3);
    map_x = local_u + rot * vert_offset_from_mass_center_[v] + center_of_mass_[p];
  }
}

void MultiDomainTet::AssembleGlobalBasis() {
  global_basis_ = MatCol::Zero(vertex_num_ * 3, total_basis_num_);
  for (int p = 0; p < part_num_; ++p) {
    for (int i = basis_offset_[p]; i < basis_offset_[p + 1]; ++i) {
      int basis_idx = i - basis_offset_[p];
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int global_v = vert_local_id2global_id_[p][local_v];
        global_basis_(global_v * 3 + 0, i) = basis_[p](local_v * 3 + 0, basis_idx);
        global_basis_(global_v * 3 + 1, i) = basis_[p](local_v * 3 + 1, basis_idx);
        global_basis_(global_v * 3 + 2, i) = basis_[p](local_v * 3 + 2, basis_idx);
      }
    }
  }
  global_basis_transpose_ = global_basis_.transpose();
  {
#if 0
    MatCol m_dot_u = global_basis_;
    for (int i = 0; i < vertex_num_; ++i) {
      m_dot_u.row(i * 3 + 0) *= mass_[i];
      m_dot_u.row(i * 3 + 1) *= mass_[i];
      m_dot_u.row(i * 3 + 2) *= mass_[i];
    }
    MatCol zero = global_basis_transpose_ * m_dot_u - MatCol::Identity(total_basis_num_, total_basis_num_);
    if (0) {
      std::ofstream out("/tmp/armadillo_150k.full_basis.txt");
      out << global_basis_.rows() << "\n" << global_basis_.cols() << std::endl;
      for (int i = 0; i < global_basis_.rows(); ++i) {
        for (int j = 0; j < global_basis_.cols(); ++j) {
          out << global_basis_(i, j) << " ";
        }
        out << "\n";
      }
      out.close();
      exit(0);
    }
    ASSERT(dj::Abs(zero.maxCoeff()) < 1e-5, P(dj::Abs(zero.maxCoeff())));
    ASSERT(dj::Abs(zero.minCoeff()) < 1e-5, P(dj::Abs(zero.minCoeff())));
    L("global basis mass diagonality verified");
#endif
  }
}

void MultiDomainTet::InitInterfaceTetrahedra() {
  interface_tet_.clear();
  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    int vert_part[4] = {
      vert_part_id_[verts[0]],
      vert_part_id_[verts[1]],
      vert_part_id_[verts[2]],
      vert_part_id_[verts[3]],
    };
    bool is_boundary_tet = false;
    for (int i = 1; i < 4; ++i) {
      if (vert_part[i] != vert_part[0]) {
        is_boundary_tet = true;
        break;
      }
    }
    if (is_boundary_tet) {
      interface_tet_.push_back(t);
    }
  }
}

void MultiDomainTet::RotateBasis() {
  ComputePartAffineTransformamtion();
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Mat3& rot = part_rotation_[p];
    for (int i = 0; i < part_basis_size_[p]; ++i) {
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        MapVec3 u(&initial_basis_[p](local_v * 3, i));
        MapVec3 rotated_u(&basis_[p](local_v * 3, i));
        rotated_u = rot * u;
      }
    }
    basis_transpose_[p] = basis_[p].transpose();
  }
  AssembleGlobalBasis();
}

void MultiDomainTet::LoadMass(const char *file_name) {
  std::ifstream in(file_name);
  std::vector<double> mm(vertex_num_);
  ASSERT(in.is_open(), P(file_name));
  int v_num;
  in >> v_num;
  ASSERT(v_num == vertex_num_);
  double total0 = 0, total1 = 0;
  for (int v = 0; v < vertex_num_; ++v) {
    double mass;
    in >> mass;
    total0 += mass;
    total1 += mass_[v];
    mm[v] = mass;
    ASSERT(dj::Abs(mass - mass_[v]) < 1e-8, P(mass, mass_[v]));
    mass_[v] = mass;
  }
  P(total0, total1);
  in.close();

  if (0) {
    int p = 0;
    std::ofstream out("/tmp/mm");
    for (int v : vert_local_id2global_id_[p]) {
      out << mm[v] << std::endl;
    }
    out.close();
    L("finish");
    exit(0);
  }
}

void MultiDomainTet::UpdatePosition() {
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Mat3& rotation = part_rotation_[p];
    MapVec sub_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
      int v = vert_local_id2global_id_[p][local_v];
      Vec3 local_u = vert_basis_[v] * sub_q;
      MapVec3 map_x(X + v * 3);
      map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
    }
  }
}


void MultiDomainTet::PrecomputeFastSandwichTransform() {
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
  inertial_force_sandwich_.resize(part_num_);
  euler_force_sandwich_.resize(part_num_);
  coriolis_force_sandwich_.resize(part_num_ * 9);
  centrifugal_force_sandwich_.resize(part_num_ * 9);
  rotation_sandwich0_.resize(part_num_);
  rotation_sandwich1_.resize(total_basis_num_);
  for (int p = 0; p < part_num_; ++p) {
    int basis_num = part_basis_size_[p];
    inertial_force_sandwich_[p] = Mat::Zero(basis_num, 3);
    euler_force_sandwich_[p] = Mat::Zero(basis_num, 3);
    for (int idx = 0; idx < 9; ++idx) {
      coriolis_force_sandwich_[p * 9 + idx] = Mat::Zero(basis_num, basis_num);
      centrifugal_force_sandwich_[p * 9 + idx] = Vec::Zero(basis_num);
    }
    rotation_sandwich0_[p].setZero();
    for (int i = 0; i < part_basis_size_[p]; ++i) {
      rotation_sandwich1_[basis_offset_[p] + i].setZero();
    }
    for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
      int v = vert_local_id2global_id_[p][local_v];
      inertial_force_sandwich_[p] += -mass_[v] * vert_basis_transpose_[v]; // (\sum_{all v} -(m_v * U_v))
      Mat3 skew_symmetrix_matrix = GetSkewSymmetrixMatrix(vert_offset_from_mass_center_[v]);
      // sandwich = \sum_{all v} (U_v^T * m_v * [r0])
      euler_force_sandwich_[p] += vert_basis_transpose_[v] * (mass_[v] * skew_symmetrix_matrix);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          rotation_sandwich0_[p](r, c) += mass_[v] * vert_offset_from_mass_center_[v][r] * vert_offset_from_mass_center_[v][c];
          for (int i = 0; i < part_basis_size_[p]; ++i) {
            rotation_sandwich1_[basis_offset_[p] + i](r, c) += mass_[v] * vert_basis_[v](r, i) * vert_offset_from_mass_center_[v][c];
          }
          int idx = r * 3 + c;
          Mat3 rotation = Mat3::Zero();
          rotation(r, c) = mass_[v];
          // sandwich = \sum_{all v} (-U_v^T * m_v * [w]^2 * r0)
          centrifugal_force_sandwich_[p * 9 + idx] -= vert_basis_transpose_[v] * (rotation * vert_offset_from_mass_center_[v]);
          // sandwich = \sum_{all v} -(U_v^T * m_v * [w] * U_v)
          coriolis_force_sandwich_[p * 9 + idx] -= vert_basis_transpose_[v] * rotation * vert_basis_[v];
        }
      }
    }
  }
}

void MultiDomainTet::LoadCubature(const char *directory) {
  const double kMinCubatureWeight = 1e-10;
  P(kMinCubatureWeight);

  std::unordered_set<int> cubature_vert_set;
  all_cubature_tet_.clear();
  domain_cubature_.resize(part_num_);
  interface_cubature_tet_.clear();
  std::unordered_set<int> interface_cubature_vert_set;
  interface_cubature_.resize(interface_num_);
  interface_rigid_cubature_.resize(interface_num_);
  if (conf.Get<int>("use all cubature points")) {
    L("using all tets as cubature point for testing");
    all_cubature_tet_.clear();
    all_cubature_vert_.clear();
    interface_cubature_tet_.clear();
    interface_cubature_vert_.clear();
    all_cubature_.clear();
    domain_cubature_.clear();
    interface_cubature_.clear();

    interface_cubature_.resize(interface_num_);
    domain_cubature_.resize(part_num_);
    all_cubature_vert_.clear();
    for (int v = 0; v < vertex_num_; ++v) {
      all_cubature_vert_.push_back(v);
    }
    interface_cubature_vert_set.clear();
    for (int t = 0; t < tet_number; ++t) {
      int* verts = tet_ + t * 4;
      std::set<int> domain_set;
      for (int i = 0; i < 4; ++i) {
        int v = verts[i];
        int p = vert_part_id_[v];
        domain_set.insert(p);
      }
      // interface tet
      if (domain_set.size() > 1) {
        for (int i = 0; i < 4; ++i) {
          interface_cubature_vert_set.insert(verts[i]);
        }
        interface_cubature_tet_.push_back(t);
        std::vector<int> domains(domain_set.begin(), domain_set.end());
        std::sort(domains.begin(), domains.end());
        int e = 0;
        for (; e < interface_num_; ++e) {
          bool found = true;
          for (int i = 0; i < int(domains.size()); ++i) {
            if (interface_domains_[e][i + 1] != domains[i]) {
              found = false;
              break;
            }
          }
          if (found) break;
        }
        ASSERT(e < interface_num_);
        //        if (t == 378) {
        //          P(e);
        //          PVEC(domains);
        //          exit(0);
        //        }
        interface_cubature_[e].emplace_back(t, 1.0);
      } else {
        int p = vert_part_id_[verts[0]];
        domain_cubature_[p].emplace_back(t, 1.0);
      }
      all_cubature_.push_back(CubaturePoint(t, 1.0));
      all_cubature_tet_.push_back(t);
    }
    interface_cubature_vert_.clear();
    interface_cubature_vert_.insert(interface_cubature_vert_.end(), interface_cubature_vert_set.begin(), interface_cubature_vert_set.end());
    interface_rigid_cubature_ = interface_cubature_;
  } else {
    for (int p = 0; p < part_num_; ++p) {
      // read the mapping from local tet_id to global tet id
      std::string cub_file_name = dj::Format("%z/partition_%d.cubature.global_tet.txt", directory, p);
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // read cubature tets and weights
      std::ifstream cubature_file(cub_file_name);
      ASSERT(cubature_file.is_open(), P(cub_file_name));
      int cubature_num;
      cubature_file >> cubature_num;
      domain_cubature_[p].reserve(cubature_num);
      int num_of_discarded_cubature = 0;
      for (int i = 0; i < cubature_num; ++i) {
        int global_tet_id;
        double weight;
        cubature_file >> global_tet_id;
        cubature_file >> weight;
        if (weight >= kMinCubatureWeight) {
          int tet_id = global_tet_id;
          ASSERT(tet_id < tet_number);
          all_cubature_tet_.push_back(tet_id);
          cubature_vert_set.insert(tet_[tet_id * 4 + 0]);
          cubature_vert_set.insert(tet_[tet_id * 4 + 1]);
          cubature_vert_set.insert(tet_[tet_id * 4 + 2]);
          cubature_vert_set.insert(tet_[tet_id * 4 + 3]);
          domain_cubature_[p].push_back(std::make_pair(tet_id, weight));
        } else {
          num_of_discarded_cubature++;
        }
      }

      std::string info =  dj::Format("inner domain %d: %d selected cubatures points "
                                     "and %d discarded with weight less than %e",
                                     p,
                                     domain_cubature_[p].size(),
                                     num_of_discarded_cubature,
                                     kMinCubatureWeight);
      L(info);
    }

    for (int e = 0; e < interface_num_; ++e) {
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // read cubature tets and weights
      std::string cubature_file_name = dj::Format("%z/interface_%d_%d_%d_%d.cubature.global_tet.txt",
                                                  directory,
                                                  interface_domains_[e][1], interface_domains_[e][2],
                                                  interface_domains_[e][3], interface_domains_[e][4]);
      std::ifstream cubature_file(cubature_file_name);
      ASSERT(cubature_file.is_open(), P(cubature_file_name));
      int cubature_num;
      cubature_file >> cubature_num;
      interface_cubature_[e].reserve(cubature_num);
      int num_of_discarded_cubature = 0;
      for (int i = 0; i < cubature_num; ++i) {
        int global_tet_id;
        double weight;
        cubature_file >> global_tet_id;
        cubature_file >> weight;
        if (weight >= kMinCubatureWeight) {
          int tet_id = global_tet_id;
          ASSERT(tet_id < tet_number);
          cubature_vert_set.insert(tet_[tet_id * 4 + 0]);
          cubature_vert_set.insert(tet_[tet_id * 4 + 1]);
          cubature_vert_set.insert(tet_[tet_id * 4 + 2]);
          cubature_vert_set.insert(tet_[tet_id * 4 + 3]);
          interface_cubature_vert_set.insert(tet_[tet_id * 4 + 0]);
          interface_cubature_vert_set.insert(tet_[tet_id * 4 + 1]);
          interface_cubature_vert_set.insert(tet_[tet_id * 4 + 2]);
          interface_cubature_vert_set.insert(tet_[tet_id * 4 + 3]);
          all_cubature_tet_.push_back(tet_id);
          interface_cubature_tet_.push_back(tet_id);
          interface_cubature_[e].push_back(std::make_pair(tet_id, weight));
        } else {
          num_of_discarded_cubature++;
        }
      }
      std::string info =  dj::Format("interface cubature %d, %d, %d, %d: %d selected cubatures points "
                                     "and %d discarded with weight less than %e",
                                     interface_domains_[e][1], interface_domains_[e][2],
                                     interface_domains_[e][3], interface_domains_[e][4],
                                     interface_cubature_[e].size(),
                                     num_of_discarded_cubature,
                                     kMinCubatureWeight);
      L(info);
      if (0) {
        std::string cubature_file_name = dj::Format("%z/rigid_interface_%d_%d_%d_%d.cubature.global_tet.txt",
                                                    directory,
                                                    interface_domains_[e][1], interface_domains_[e][2],
                                                    interface_domains_[e][3], interface_domains_[e][4]);
        std::ifstream cubature_file(cubature_file_name);
        ASSERT(cubature_file.is_open(), P(cubature_file_name));
        int cubature_num;
        cubature_file >> cubature_num;
        interface_rigid_cubature_[e].reserve(cubature_num);
        int num_of_discarded_cubature = 0;
        for (int i = 0; i < cubature_num; ++i) {
          int global_tet_id;
          double weight;
          cubature_file >> global_tet_id;
          cubature_file >> weight;
          if (weight >= kMinCubatureWeight) {
            int tet_id = global_tet_id;
            ASSERT(tet_id < tet_number);
            all_cubature_tet_.push_back(tet_id);
            interface_rigid_cubature_[e].push_back(std::make_pair(tet_id, weight));
          } else {
            num_of_discarded_cubature++;
          }
        }
        std::string info = dj::Format("rigid interface cubature %d, %d, %d, %d: %d selected cubatures points "
                                      "and %d discarded with weight less than %e",
                                      interface_domains_[e][1], interface_domains_[e][2],
                                      interface_domains_[e][3], interface_domains_[e][4],
                                      interface_rigid_cubature_[e].size(),
                                      num_of_discarded_cubature,
                                      kMinCubatureWeight);
        L(info);
      }
    }
    std::unordered_set<int> all_cubature_tet_set(all_cubature_tet_.begin(), all_cubature_tet_.end());
    if (full_rigid_simulation_) {
      all_cubature_tet_set.insert(interface_tet_.begin(), interface_tet_.end());
    }
    all_cubature_tet_ = std::vector<int>(all_cubature_tet_set.begin(), all_cubature_tet_set.end());
    all_cubature_vert_.clear();
    all_cubature_vert_.insert(all_cubature_vert_.end(), cubature_vert_set.begin(), cubature_vert_set.end());
    interface_cubature_vert_.clear();
    interface_cubature_vert_.insert(interface_cubature_vert_.end(), interface_cubature_vert_set.begin(), interface_cubature_vert_set.end());
    all_cubature_.clear();
    for (int p = 0; p < part_num_; ++p) {
      all_cubature_.insert(all_cubature_.end(), domain_cubature_[p].begin(), domain_cubature_[p].end());
    }
    for (int e = 0; e < interface_num_; ++e) {
      all_cubature_.insert(all_cubature_.end(), interface_cubature_[e].begin(), interface_cubature_[e].end());
    }
    interface_rigid_cubature_ = interface_cubature_;
  }
  ComputeInterfaceTorqueMatrix();
  BuildParallelAccelerationStructure();
  P(all_cubature_.size(), tet_number);
  P(all_cubature_vert_.size(), vertex_num_);
}

void MultiDomainTet::ComputeInvInertiaTensor() {
  inv_inertia_tensor_.resize(part_num_);
  inertia_tensor_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    Mat3 inertia_tensor = Mat3::Zero();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      Mat3 tmp = Mat3::Zero();
      Vec3 r(rest_pos_[v * 3 + 0] - initial_center_of_mass_[p][0],
             rest_pos_[v * 3 + 1] - initial_center_of_mass_[p][1],
             rest_pos_[v * 3 + 2] - initial_center_of_mass_[p][2]);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          tmp(i, j) = -r[i] * r[j];
        }
      }
      double dot = r.dot(r);
      tmp(0, 0) += dot;
      tmp(1, 1) += dot;
      tmp(2, 2) += dot;
      inertia_tensor += mass_[v] * tmp;
    }
    inertia_tensor_[p] = inertia_tensor;
    inv_inertia_tensor_[p] = inertia_tensor.inverse();
  }
  current_inertia_tensor_ = inertia_tensor_;
  current_inv_inertia_tensor_ = inv_inertia_tensor_;
}

void MultiDomainTet::MultiBodyFullSimulation(double dt) {
  static unsigned int steps = 0;
  if (1 && steps == 0) {
    constrainted_vertex_.clear();
    if (1)
      for (int v = 0 ; v < vertex_num_; ++v) {
        if (X[v * 3 + 2] < 1e-3) constrainted_vertex_.push_back(v);
      }
    if (0)
      for (int p = 1; p < part_num_; ++p) {
        int min_v = -9999999;
        double min_dist = 1e10;
        for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
          int v = vert_local_id2global_id_[p][local_v];
          double dist = dj::Distance3(X + v * 3, &center_of_mass_[p][0]);
          if (dist < min_dist) {
            min_v = v;
            min_dist = dist;
          }
        }
        ASSERT(min_v > 0);
        P(p, min_v);
        constrainted_vertex_.push_back(min_v);
        constrainted_vertex_.push_back(min_v - 32);
        constrainted_vertex_.push_back(min_v + 1);
      }
    //    constrainted_vertex_.push_back(177);
    //    constrainted_vertex_.push_back(179);
    //  constrainted_vertex_.clear();
    memset(&is_constrainted_[0], 0, sizeof(int) * vertex_num_);
    for (int v : constrainted_vertex_) {
      is_constrainted_[v] = true;
    }
  }
  steps++;

  //  gravity_[0] = 0; gravity_[1] = 0; gravity_[2] = 0;
  //  PVEC(gravity_);
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  Real dt_2 = dt * dt;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // net force on each part
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  // net torque on each part
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));

  // gravity on each part
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += mass_per_part_[p] * gravity_;
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external force in world space to each vertex
  //      FullSimulation(dt);
  //      return;
  std::vector<Vec3> ext_force(vertex_num_, Vec3(0, 0, 0));
  if (dt * steps < 1.0) {
    int applied_vertex = conf.Get<int>("selected_vertex");
    if (applied_vertex >= 0 && applied_vertex < vertex_num_) {
      double* force = conf.Get<double*>("constant force");
      ext_force[applied_vertex][0] += force[0];
      ext_force[applied_vertex][1] += force[1];
      ext_force[applied_vertex][2] += force[2];
    }
  } else {
    //    FullSimulation(dt);
    //    return;
  }
  const double kFloor = -0.20;
  const double kFloorStiffness = 3000;
  // Ground collision
  if (conf.Get<int>("ground collision"))
    for (int v = 0; v < vertex_num_; ++v) {
      if (X[v * 3 + 1] < kFloor) {
        Vec3 collision_force(0, 0, 0);
        //      collision_force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
        collision_force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
        ext_force[v][1] += collision_force[1];
        int p = vert_part_id_[v];
        net_force[p][1] += collision_force[1];

        MapVec3 map_x(X + v * 3);
        Vec3 r = map_x - center_of_mass_[p];
        net_torque[p] += r.cross(collision_force);
      }
    }

  bool enable_interface_force = true;
  {
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      ext_force[v] += force;
      //      has_ui_force = true;
      int p = vert_part_id_[v];
      net_force[p] += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[p] += r.cross(force);
    }
  }
  // elastic force from boundary tets
  if (enable_interface_force)
    for (int t : interface_tet_) {
      // the force is computed from vega, the direction of the force is opposite to the actual direction
      double* force = inv_fem_->element_force_ + t * 12;
      for (int i = 0; i < 4; ++i) {
        int v = tet_[t * 4 + i];
        int p = vert_part_id_[v];
        MapVec3 map_force(force + i * 3);
        net_force[p] -= map_force;
        MapVec3 map_x(X + v  * 3);
        Vec3 r = map_x - center_of_mass_[p];
        net_torque[p] -= r.cross(map_force);
      }
    }

  std::vector<Vec3> rhs(part_num_ * 2);
  std::vector<Mat3> current_inertia_tensor_ = inertia_tensor_;
  if (0) {
    Vec3 acc0 = net_force[0] / mass_per_part_[0];
    Vec3 acc1 = net_force[1] / mass_per_part_[1];
    Vec3 dv0 = acc0 * dt;
    Vec3 dv1 = acc1 * dt;
    P(acc0, mass_per_part_[0], dv0);
    P(acc1, mass_per_part_[1], dv1);
  }
  for (int p = 0; p < part_num_; ++p) {
    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_[p].transpose();
    rhs[p * 2 + 0] = net_force[p] * dt + translational_vel_[p] * mass_per_part_[p];
    rhs[p * 2 + 1] = net_torque[p] * dt + current_inertia_tensor_[p] * angular_vel_[p];
  }
  //  P(net_force[1], translational_vel_[1]);

  std::vector<Vec3> implicit_vel(vertex_num_);
  auto K = [&](double * x, double * result) {
    for (int v = 0; v < vertex_num_; ++v) {
      int p = vert_part_id_[v];
      MapVec3 map_pos(X + v * 3);
      Vec3 r = map_pos - center_of_mass_[p];
      double angular_vel[3];//= angular_vel_[p].cross(r);
      dj::Cross3(x + p * 6 + 3, &r[0], &angular_vel[0]);
      implicit_vel[v][0] = x[p * 6 + 0] + angular_vel[0];
      implicit_vel[v][1] = x[p * 6 + 1] + angular_vel[1];
      implicit_vel[v][2] = x[p * 6 + 2] + angular_vel[2];
    }

    for (int p = 0; p < part_num_; ++p) {
      result[p * 6 + 0] = x[p * 6 + 0] * mass_per_part_[p];
      result[p * 6 + 1] = x[p * 6 + 1] * mass_per_part_[p];
      result[p * 6 + 2] = x[p * 6 + 2] * mass_per_part_[p];

      MapVec3 map_result(result + p * 6 + 3);
      MapVec3 map_x(x + p * 6 + 3);
      map_result = current_inertia_tensor_[p] * map_x;
    }
    //        return;
    // elastic force from boundary tets
    if (enable_interface_force)
      for (int t : interface_tet_) {
        // the force is computed from vega, the direction of the force is opposite to the actual direction
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          Vec3 force(0, 0, 0);
          int v = tet_[t * 4 + i];
          // compute single tet element force on vertex v
          for (int j = 0; j < 4; ++j) {
            int vj = tet_[t * 4 + j];
            force += element_k.block<3, 3>(i * 3, j * 3) * implicit_vel[vj];
          }
          force *= dt_2;
          int p = vert_part_id_[v];
          MapVec3(result + p * 6 + 0) += force;

          MapVec3 map_x(X + v  * 3);
          Vec3 r = map_x - center_of_mass_[p];
          MapVec3(result + p * 6 + 3) += r.cross(force);
        }
      }
  };

  MatCol RigidStiffnessMatrix = MatCol::Zero(part_num_ * 6, part_num_ * 6);
  Vec identity = Vec::Zero(part_num_ * 6);
  for (int i = 0; i < part_num_ * 6; ++i) {
    identity[i] = 1;
    K(&identity[0], RigidStiffnessMatrix.col(i).data());
    identity[i] = 0;
  }
  Vec new_rigid_velocity = Vec::Zero(part_num_ * 6);
  if (fixed_domain_ >= 0 && fixed_domain_ < part_num_) {
    int size0 = fixed_domain_ * 6;
    int size1 = (part_num_ - fixed_domain_ - 1) * 6;
    Vec fixed_rhs = Vec::Zero((part_num_ - 1) * 6);
    MatCol fixed_k = MatCol::Zero((part_num_ - 1) * 6, (part_num_ - 1) * 6);
    if (fixed_domain_ > 0) {
      memcpy(&fixed_rhs[0], &rhs[0][0], sizeof(double) * size0);
      fixed_k.block(0, 0, size0, size0) = RigidStiffnessMatrix.block(0, 0, size0, size0);
    }
    if (fixed_domain_ < part_num_ - 1) {
      memcpy(&fixed_rhs[size0], &rhs[0][0] + size0 + 6, sizeof(double) * size1);
      fixed_k.block(size0, size0, size1, size1) = RigidStiffnessMatrix.block(size0 + 6, size0 + 6, size1, size1);
    }
    Vec fixed_vel = fixed_k.colPivHouseholderQr().solve(fixed_rhs);
    if (fixed_domain_ > 0) {
      memcpy(&new_rigid_velocity[0], &fixed_vel[0], sizeof(double) * size0);
    }
    if (fixed_domain_ < part_num_ - 1) {
      memcpy(&new_rigid_velocity[size0 + 6], &fixed_vel[size0], sizeof(double) * size1);
    }
  } else {
    MapVec map_rhs(&rhs[0][0], part_num_ * 6);
    new_rigid_velocity = RigidStiffnessMatrix.colPivHouseholderQr().solve(map_rhs);
  }


  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    MapVec3 new_translational_vel(&new_rigid_velocity[p * 6]);
    MapVec3 new_angular_vel(&new_rigid_velocity[p * 6 + 3]);
    if (p == 0) {
      new_translational_vel.setZero();
      new_angular_vel.setZero();
    }
    //    new_translational_vel.setZero();
    //    new_angular_vel.setZero();
    //    new_translational_vel.setZero();
    //    if (!has_ui_force) {
    //      new_translational_vel.setZero();
    //    }
    Vec3 dv = new_translational_vel - translational_vel_[p];
    Vec3 dw = new_angular_vel - angular_vel_[p];
    acceleration[p] = (dv) / dt;
    angular_acceleration[p] = (dw) / dt;
    //    P(p, acceleration[p], angular_acceleration[p]);

    //    translational_vel_[p] = new_translational_vel;
    //    angular_vel_[p] = new_angular_vel;
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
#if 1
  // internal force
  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    if (!enable_interface_force) {
      bool is_boundary = false;
      for (int i = 1; i < 4; ++i) {
        if (vert_part_id_[verts[i]] != vert_part_id_[verts[0]]) {
          is_boundary = true;
        }
      }
      if (is_boundary) continue;
    }
    double* element_force = inv_fem_->element_force_ + t * 12;
    for (int i = 0; i < 4; ++i) {
      ext_force[verts[i]][0] -= element_force[i * 3 + 0];
      ext_force[verts[i]][1] -= element_force[i * 3 + 1];
      ext_force[verts[i]][2] -= element_force[i * 3 + 2];
    }
  }
  if (0)
    for (int t : interface_tet_) {
      int* verts = tet_ + t * 4;
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = vert_part_id_[vi];
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int pj = vert_part_id_[vj];
          Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
          ext_force[verts[i]] -= part_rotation_[pi].transpose() *
                                 element_k.block<3, 3>(i * 3, j * 3)
                                 * (translational_vel_[pj] + angular_vel_[pj].cross(rj)) * dt;
        }
      }
    }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // apply fictitious force
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Vec3& acc = acceleration[p];
    Vec3& angular_acc = angular_acceleration[p];
    Mat3 rotation_transpose = part_rotation_[p].transpose();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      // inertial force: F= -m * a
      Vec3 fictitious_force = -mass_[v] * acc;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[p];
      // Euler force: F = -m * dw/dt \times r
      if (1) {
        Vec3 euler_force = -mass_[v] * angular_acc.cross(r);
        fictitious_force += euler_force;
      }
      // centerifugal force: F = -m * w \times (w \times r)
      if (1)  {
        Vec3 centerifugal_force = -mass_[v] * angular_vel_[p].cross(angular_vel_[p].cross(r));
        fictitious_force += centerifugal_force;
      }
      // Coriolis force: F = -2m w\times v
      if (1) {
        MapVec3 map_vel(velocity_ + v * 3);
        Vec3 croilis_force = -2 * mass_[v] * angular_vel_[p].cross(map_vel);
        fictitious_force += croilis_force;
      }
      ext_force[v] += fictitious_force;
      ext_force[v] += gravity_ * mass_[v];
      // transpose external force to part's local frame
      ext_force[v] = rotation_transpose * ext_force[v];
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // simulate internal dynamics
  // rhs = dt * R^T * F + M * vel
  {
    std::vector<Vec3> rhs(vertex_num_);
    for (int v = 0; v < vertex_num_; ++v) {
      MapVec3 map_vel(velocity_ + v * 3);
      rhs[v] = dt * ext_force[v] + mass_[v] * map_vel;
    }

    std::vector<Vec3> tmp_rhs(vertex_num_);
    // A = (M + dt^2 * R^T * K * R)
    static_assert(sizeof(Vec3) == sizeof(double) * 3, "invalid assumption");
#if 1
    auto StiffnessMatrix = [&](Real * x, Real * result) {
      // result = R * x
      for (int v = 0; v < vertex_num_; ++v) {
        Mat3& rotation = part_rotation_[vert_part_id_[v]];
        MapVec3 map_x(x + v * 3);
        MapVec3 map_result(result + v * 3);
        map_result = rotation * map_x;
      }

      memset(&tmp_rhs[0][0], 0, sizeof(double) * vertex_num_ * 3);
      for (int t = 0; t < tet_number; ++t) {
        int* verts = tet_ + t * 4;
        if (!enable_interface_force) {
          bool is_boundary = false;
          for (int i = 1; i < 4; ++i) {
            if (vert_part_id_[verts[i]] != vert_part_id_[verts[0]]) {
              is_boundary = true;
            }
          }
          if (is_boundary) continue;
        }
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            tmp_rhs[verts[i]] += element_k.block<3, 3>(i * 3, j * 3) * MapVec3(result + verts[j] * 3);
          }
        }
      }


      // tmp_rhs = K * R * x;
      //      inv_fem_->tangent_stiffness_matrix_->MultiplyVector(result, &tmp_rhs[0][0]);
      //    OMP_FOR
      for (int v = 0; v < vertex_num_; ++v) {
        Mat3 rotation_transpose = part_rotation_[vert_part_id_[v]].transpose();
        tmp_rhs[v] = rotation_transpose * tmp_rhs[v];
        tmp_rhs[v] *= dt_2;
        result[v * 3 + 0] = mass_[v] * x[v * 3 + 0] + tmp_rhs[v][0];
        result[v * 3 + 1] = mass_[v] * x[v * 3 + 1] + tmp_rhs[v][1];
        result[v * 3 + 2] = mass_[v] * x[v * 3 + 2] + tmp_rhs[v][2];

        if (0)
          if (vert_part_id_[v] == 0) {
            result[v * 3 + 0] = 0;
            result[v * 3 + 1] = 0;
            result[v * 3 + 2] = 0;
          }
      }
      // constrained vertex
      for (int v : constrainted_vertex_) {
        result[v * 3 + 0] = 0;
        result[v * 3 + 1] = 0;
        result[v * 3 + 2] = 0;
      }
    };
#else
    auto StiffnessMatrix = [&](Real * x, Real * result) {
      // tmp_rhs = K * R * x;
      inv_fem_->tangent_stiffness_matrix_->MultiplyVector(x, &tmp_rhs[0][0]);
      //    OMP_FOR
      for (int v = 0; v < vertex_num_; ++v) {
        Mat3 rotation_transpose = part_rotation_[vert_part_id_[v]].transpose();
        tmp_rhs[v] = rotation_transpose * tmp_rhs[v];
        tmp_rhs[v] *= dt_2;
        result[v * 3 + 0] = mass_[v] * x[v * 3 + 0] + tmp_rhs[v][0];
        result[v * 3 + 1] = mass_[v] * x[v * 3 + 1] + tmp_rhs[v][1];
        result[v * 3 + 2] = mass_[v] * x[v * 3 + 2] + tmp_rhs[v][2];
        if (0)
          if (vert_part_id_[v] == 0) {
            result[v * 3 + 0] = 0;
            result[v * 3 + 1] = 0;
            result[v * 3 + 2] = 0;
          }
      }
      // constrained vertex
      for (int v : constrainted_vertex_) {
        result[v * 3 + 0] = 0;
        result[v * 3 + 1] = 0;
        result[v * 3 + 2] = 0;
      }
    };
#endif
    for (int v : constrainted_vertex_) {
      rhs[v][0] = 0;
      rhs[v][1] = 0;
      rhs[v][2] = 0;
      velocity_[v * 3 + 0] = 0;
      velocity_[v * 3 + 1] = 0;
      velocity_[v * 3 + 2] = 0;
    }

    if (0)
      for (int v = 0; v < vertex_num_; ++v) {
        if (vert_part_id_[v] == 0) {
          rhs[v][0] = 0;
          rhs[v][1] = 0;
          rhs[v][2] = 0;
          velocity_[v * 3 + 0] = 0;
          velocity_[v * 3 + 1] = 0;
          velocity_[v * 3 + 2] = 0;
        }
      }

    cg_solver_->Resize(vertex_num_ * 3);
    auto info = cg_solver_->Solve((double*) (&rhs[0][0]), velocity_, StiffnessMatrix, 2000, 1e-15);
    (void) info;
    //    P(info.first, info.second);
    //    ASSERT(info.second == info.second);
  }
#endif
#if 0
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    center_of_mass_[p] += translational_vel_[p] * dt;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();
    //    const double kDamping = 0.97;
    const double kDamping = 1.000;
    // damping
    translational_vel_[p] *= kDamping;
    angular_vel_[p] *= kDamping;

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
    quaternion_[p].Quaternion2Matrix(rotation_matrix.data());
    part_rotation_[p] = rotation_matrix;
    //    part_rotation_[p].setIdentity();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      local_u_[v][0] += velocity_[v * 3 + 0] * dt;
      local_u_[v][1] += velocity_[v * 3 + 1] * dt;
      local_u_[v][2] += velocity_[v * 3 + 2] * dt;
#if 0
      local_u_[v][0] = 0;
      local_u_[v][1] = 0;
      local_u_[v][2] = 0;
#endif

      // damping
      velocity_[v * 3 + 0] *= kDamping;
      velocity_[v * 3 + 1] *= kDamping;
      velocity_[v * 3 + 2] *= kDamping;

      MapVec3 map_x(X + v * 3);
      map_x = rotation_matrix * (local_u_[v] + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
      //      map_x = rotation_matrix * (vert_offset_from_mass_center_[v]) + center_of_mass_[p];
      inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
      inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
      inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
    }
  }
#else
  for (int v = 0; v < vertex_num_; ++v) {
    local_u_[v] += MapVec3(velocity_ + v * 3) * dt;
    local_r_[v] = local_u_[v] + vert_offset_from_mass_center_[v];
  }
  std::vector<Vec3> cm_offset(part_num_, Vec3(0, 0, 0));
  std::vector<Mat3> rotation(part_num_, Mat3::Zero());
  if (1) {
    for (int p = 0; p < part_num_; ++p) {
      if (p == fixed_domain_) continue;
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int v = vert_local_id2global_id_[p][local_v];
        cm_offset[p] += mass_[v] * local_u_[v];
        for (int i = 0; i < 3; ++i) {
          rotation[p](i, 0) += mass_[v] * local_r_[v][i] * vert_offset_from_mass_center_[v][0];
          rotation[p](i, 1) += mass_[v] * local_r_[v][i] * vert_offset_from_mass_center_[v][1];
          rotation[p](i, 2) += mass_[v] * local_r_[v][i] * vert_offset_from_mass_center_[v][2];
        }
      }
      cm_offset[p] /= mass_per_part_[p];
      Mat3 rot;
      Mat3 scaling;
      dj::PolarDecompose<double>(rotation[p].data(), &rot(0, 0), &scaling(0, 0), 1e-8);
      rotation[p] = rot;
      Mat3 I = rot.transpose() * rot - Mat3::Identity();
      //      PMAT(rot);
      ASSERT(I.norm() < 1e-8);
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int v = vert_local_id2global_id_[p][local_v];
        Vec3 old_pos = local_u_[v];
        //        local_u_[v] = local_r_[v] - vert_offset_from_mass_center_[v];
        local_u_[v] = rot.transpose() * local_r_[v] - vert_offset_from_mass_center_[v];
        local_u_[v] -= cm_offset[p];
        local_r_[v] = local_u_[v] + vert_offset_from_mass_center_[v];
        MapVec3(velocity_ + v * 3) += (local_u_[v] - old_pos) / dt;
      }
    }
  }

  //    memcpy(velocity_, &rhs[0][0], sizeof(double) * vertex_num_ * 3);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //  double sum_vel = 0;
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    translational_vel_[p] += acceleration[p] * dt;
    translational_vel_[p] += part_rotation_[p] * cm_offset[p] / dt;
    center_of_mass_[p] += translational_vel_[p] * dt;
    //    const double kDamping = 0.995;
    const double kDamping = 1.0;
    translational_vel_[p] *= kDamping;
    angular_vel_[p] += angular_acceleration[p] * dt;
    if (1)
      if (p != fixed_domain_) {
        double angle_in_radian;
        Vec3 axis;
        auto rot = Quaternion<double>::Matrix2Quaternion(&rotation[p](0, 0));
        rot.GetRotation(&angle_in_radian, &axis[0]);
        //    dj::RotationMatrix2AngleAxis(rotation[p].data(), angle_in_radian, &axis[0]);
        //        PMAT(rotation[p]);
        //        P(axis, angle_in_radian);
        axis *= angle_in_radian;
        axis = part_rotation_[p] * axis;
        angular_vel_[p] += axis / dt;
      }
    angular_vel_[p] *= kDamping;
    quaternion_[p] = quaternion_[p] +
                     (0.5 * dt) *
                     Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) *
                     quaternion_[p];
    quaternion_[p].Normalize();

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
    quaternion_[p].Quaternion2Matrix(rotation_matrix.data());
    part_rotation_[p] = rotation_matrix;
    part_rotation_transpose_[p] = part_rotation_[p].transpose();
    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
    current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();
    //    part_rotation_[p].setIdentity();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      MapVec3 map_x(X + v * 3);
      map_x = rotation_matrix * (local_r_[v]) + center_of_mass_[p];
      inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
      inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
      inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
    }
  }

#endif
  //  inv_fem_->UpdateOffset();
  if (0) {
    P(dj::Vec3d(&translational_vel_[0][0]), dj::Vec3d(&angular_vel_[0][0]), dj::Vec3d(&center_of_mass_[0][0]));
    P(dj::Vec3d(&translational_vel_[1][0]), dj::Vec3d(&angular_vel_[1][0]), dj::Vec3d(&center_of_mass_[1][0]));
    //    PMAT(part_rotation_[0]);
    //    PMAT(part_rotation_[1]);
  }
}

void MultiDomainTet::MultiBodyFullSimulationExplicit(double dt) {
  //#define EXPLICIT_NO_ROTATION
  static unsigned int steps = 0;
  if (1 && steps == 0) {
    constrainted_vertex_.clear();
    if (1)
      for (int v = 0 ; v < vertex_num_; ++v) {
        if (X[v * 3 + 2] < 1e-3) constrainted_vertex_.push_back(v);
      }
    if (0)
      for (int p = 1; p < part_num_; ++p) {
        int min_v = -9999999;
        double min_dist = 1e10;
        for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
          int v = vert_local_id2global_id_[p][local_v];
          double dist = dj::Distance3(X + v * 3, &center_of_mass_[p][0]);
          if (dist < min_dist) {
            min_v = v;
            min_dist = dist;
          }
        }
        ASSERT(min_v > 0);
        P(p, min_v);
        constrainted_vertex_.push_back(min_v);
        constrainted_vertex_.push_back(min_v - 32);
        constrainted_vertex_.push_back(min_v + 1);
      }
    //    constrainted_vertex_.push_back(177);
    //    constrainted_vertex_.push_back(179);
    //  constrainted_vertex_.clear();
    memset(&is_constrainted_[0], 0, sizeof(int) * vertex_num_);
    for (int v : constrainted_vertex_) {
      is_constrainted_[v] = true;
    }
  }
  steps++;


  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  SparseMatrix* tangent_stiffness = inv_fem_->tangent_stiffness_matrix_;
  Real dt_2 = dt * dt;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // net force on each part
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  // net torque on each part
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));
  // gravity on each part
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += mass_per_part_[p] * gravity_;
  }
  // elastic force from boundary tets
  for (int t : interface_tet_) {
    // the force is computed from vega, the direction of the force is opposite to the actual direction
    double* force = inv_fem_->element_force_ + t * 12;
    for (int i = 0; i < 4; ++i) {
      int v = tet_[t * 4 + i];
      int p = vert_part_id_[v];
      MapVec3 map_force(force + i * 3);
      net_force[p] -= map_force;
      MapVec3 map_x(X + v  * 3);
      Vec3 r = map_x - center_of_mass_[p];
      net_torque[p] -= r.cross(map_force);
    }
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external force in world space to each vertex
  std::vector<Vec3> ext_force(vertex_num_, Vec3(0, 0, 0));

  if (dt * steps < 2.0) {
    ext_force[621] += Vec3(0, -500, 0);
  } else {
    //        FullSimulation(dt);
    //        return;
  }
  const double kFloor = -0.02;
  const double kFloorStiffness = 38000;
  for (int v = 0; v < vertex_num_; ++v) {
    // gravity and internal elastic force
    ext_force[v][0] += inv_fem_->rhs_[v][0];
    ext_force[v][1] += inv_fem_->rhs_[v][1];
    ext_force[v][2] += inv_fem_->rhs_[v][2];
    // Ground collision
    if (0)
      if (X[v * 3 + 1] < kFloor) {
        Vec3 collision_force(0, 0, 0);
        collision_force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
        //      collision_force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
        ext_force[v][1] += collision_force[1];
        net_force[vert_part_id_[v]][1] += collision_force[1];

        MapVec3 map_x(X + v * 3);
        Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
        net_torque[vert_part_id_[v]] += r.cross(collision_force);
      }
  }

  bool has_ui_force = true;
  (void) has_ui_force;
  {
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      has_ui_force = true;
      ext_force[v] += force;
      int p = vert_part_id_[v];
      net_force[p] += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[p] += r.cross(force);
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // compute rigid translational and rotational acceleration
  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    // TODO delete
    if (p == 0) {
      acceleration[p].setZero();
      angular_acceleration[p].setZero();
      continue;
    }
    if (0) {
      Vec3 acc0 = net_force[0] / mass_per_part_[0];
      Vec3 acc1 = net_force[1] / mass_per_part_[1];
      Vec3 dv0 = acc0 * dt;
      Vec3 dv1 = acc1 * dt;
      P(acc0, mass_per_part_[0], dv0);
      P(acc1, mass_per_part_[1], dv1);
    }
    acceleration[p] = net_force[p] / mass_per_part_[p];
    Mat3 cur_inv_inertial = part_rotation_[p] * inv_inertia_tensor_[p] * part_rotation_[p].transpose();
    angular_acceleration[p] = cur_inv_inertial * net_torque[p];
#ifdef EXPLICIT_NO_ROTATION
    angular_acceleration[p].setZero();
#endif
    //    angular_acceleration[p].setZero();
    //    acceleration[p].setZero();
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // apply fictitious force
  //  L("no fictitious_force force");
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Vec3& acc = acceleration[p];
    Vec3& angular_acc = angular_acceleration[p];
    Mat3 rotation_transpose = part_rotation_[p].transpose();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      //       inertial force: F= -m * a
      Vec3 fictitious_force = -mass_[v] * acc;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[p];
      // Euler force: F = -m * dw/dt \times r
      if (1) {
        Vec3 euler_force = -mass_[v] * angular_acc.cross(r);
        fictitious_force += euler_force;
      }
      // centerifugal force: F = -m * w \times (w \times r)
      if (1)  {
        Vec3 centerifugal_force = -mass_[v] * angular_vel_[p].cross(angular_vel_[p].cross(r));
        fictitious_force += centerifugal_force;
      }
      // Coriolis force: F = -2m w\times v
      if (1) {
        MapVec3 map_vel(velocity_ + v * 3);
        Vec3 croilis_force = -2 * mass_[v] * angular_vel_[p].cross(map_vel);
        fictitious_force += croilis_force;
      }

      ext_force[v] += fictitious_force;
      // transpose external force to local part local frame
      ext_force[v] = rotation_transpose * ext_force[v];
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // simulate internal dynamics
  // rhs = dt * R^T * F + M * vel
  std::vector<Vec3> rhs(vertex_num_);
  //  if (0)
  for (int t : interface_tet_) {
    int* verts = tet_ + t * 4;
    Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
    for (int i = 0; i < 4; ++i) {
      int vi = verts[i];
      int pi = vert_part_id_[vi];
      for (int j = 0; j < 4; ++j) {
        int vj = verts[j];
        int pj = vert_part_id_[vj];
        Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
        ext_force[verts[i]] -= part_rotation_[pi].transpose() *
                               element_k.block<3, 3>(i * 3, j * 3)
                               * (translational_vel_[pj] + angular_vel_[pj].cross(rj)) * dt;
      }
    }
  }
  for (int v = 0; v < vertex_num_; ++v) {
    MapVec3 map_vel(velocity_ + v * 3);
    rhs[v] = dt * ext_force[v] + mass_[v] * map_vel;
  }

  // A = (M + dt^2 * R^T * K * R)
  std::vector<Vec3> tmp_rhs(vertex_num_);
  static_assert(sizeof(Vec3) == sizeof(double) * 3, "invalid assumption");
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    //    memcpy(result, x, sizeof(double) * vertex_num_ * 3); return;
    // result = R * x
    //    OMP_FOR
    for (int v = 0; v < vertex_num_; ++v) {
      Mat3& rotation = part_rotation_[vert_part_id_[v]];
      MapVec3 map_x(x + v * 3);
      MapVec3 map_result(result + v * 3);
      map_result = rotation * map_x;
    }
    // tmp_rhs = K * R * x;
    tangent_stiffness->MultiplyVector(result, &tmp_rhs[0][0]);
    //    memcpy(result, &tmp_rhs[0][0], sizeof(double) * vertex_num_ * 3); return;
    //    OMP_FOR
    for (int v = 0; v < vertex_num_; ++v) {
      Mat3 rotation_transpose = part_rotation_[vert_part_id_[v]].transpose();
      tmp_rhs[v] = rotation_transpose * tmp_rhs[v];
      tmp_rhs[v] *= dt_2;
      result[v * 3 + 0] = mass_[v] * x[v * 3 + 0] + tmp_rhs[v][0];
      result[v * 3 + 1] = mass_[v] * x[v * 3 + 1] + tmp_rhs[v][1];
      result[v * 3 + 2] = mass_[v] * x[v * 3 + 2] + tmp_rhs[v][2];
    }
    // constrained vertex
    for (int v : constrainted_vertex_) {
      result[v * 3 + 0] = 0;
      result[v * 3 + 1] = 0;
      result[v * 3 + 2] = 0;
    }
  };
  memcpy(velocity_, &rhs[0][0], sizeof(double) * vertex_num_ * 3);
  for (int v : constrainted_vertex_) {
    rhs[v][0] = 0;
    rhs[v][1] = 0;
    rhs[v][2] = 0;
    velocity_[v * 3 + 0] = 0;
    velocity_[v * 3 + 1] = 0;
    velocity_[v * 3 + 2] = 0;
  }

  {
#if 0
    // verify
    MapVec b(&rhs[0][0], (vertex_num_) * 3);
    //    Vec b((vertex_num_ + part_num_ * 2) * 3);
    //    b.setZero();
    //    b[0] = 1;
    Vec tmp = b;
    StiffnessMatrix(&b[0], &tmp[0]);
    //    tangent_stiffness->MultiplyVector(&b[0], &tmp[0]);
    P(tmp);
    //    P(b);
    //    double dot = 0;
    //    for (int i = 0; i < vertex_num_ * 3; ++i) {
    //      dot += tmp[i] * tmp[i];
    //    }
    //    P(dot);
    P(b.dot(b), tmp.dot(tmp));
    exit(0);
#endif
  }
  cg_solver_->Resize(vertex_num_ * 3);
  auto info = cg_solver_->Solve((double*) (&rhs[0][0]), velocity_, StiffnessMatrix, 2000, 1e-20);
  P(info.first, info.second);


  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // compensate rigid motion after local deformation
#if 0
  std::vector<Vec3> momentum(part_num_, Vec3(0, 0, 0));
  std::vector<Vec3> angular_momentum(part_num_, Vec3(0, 0, 0));
  if (1)
    for (int p = 0; p < part_num_; ++p) {
      if (p == 0) continue;
      int fixed_vert_num = 0;
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int v = vert_local_id2global_id_[p][local_v];
        if (is_constrainted_[v]) fixed_vert_num++;
        MapVec3 map_vel(velocity_ + v * 3);
        momentum[p] += mass_[v] * map_vel;
        angular_momentum[p] += mass_[v] * (vert_offset_from_mass_center_[v].cross(map_vel));
      }
      Vec3 momentum_correction = momentum[p] * (vert_num_per_part_[p] / (vert_num_per_part_[p] - fixed_vert_num) / (mass_per_part_[p]));
      Vec3 angular_correction = inv_inertia_tensor_[p] * angular_momentum[p] *
                                (vert_num_per_part_[p] / (vert_num_per_part_[p] - fixed_vert_num));
      //    PMAT(inertia_tensor_[p]);
      //    PMAT(inv_inertia_tensor_[p]);
      //    P(angular_momentum[p]);
      //    P(angular_correction, momentum_correction, mass_per_part_[p]);
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int v = vert_local_id2global_id_[p][local_v];
        if (is_constrainted_[v]) continue;
        MapVec3 map_vel(velocity_ + v * 3);
        map_vel -= momentum_correction;
#ifndef EXPLICIT_NO_ROTATION
        //        map_vel += mass_[v] * vert_offset_from_mass_center_[v].cross(angular_correction);
#endif
      }
    }
  //    memcpy(velocity_, &rhs[0][0], sizeof(double) * vertex_num_ * 3);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //  double sum_vel = 0;
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    translational_vel_[p] += acceleration[p] * dt + part_rotation_[p] * momentum[p] / mass_per_part_[p];
    center_of_mass_[p] += translational_vel_[p] * dt;
    const double kDamping = 0.995;
    //    const double kDamping = 1.0;
    translational_vel_[p] *= kDamping;
    angular_vel_[p] += angular_acceleration[p] * dt;
#ifndef EXPLICIT_NO_ROTATION
    if (p == 1) {
      //      angular_vel_[p] += part_rotation_[p] * inv_inertia_tensor_[p] * angular_momentum[p];
    }
#endif
    angular_vel_[p] *= kDamping;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
    quaternion_[p].Quaternion2Matrix(rotation_matrix.data());
    part_rotation_[p] = rotation_matrix;
    part_rotation_transpose_[p] = part_rotation_[p].transpose();
    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
    current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();
    //    part_rotation_[p].setIdentity();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      //      if (has_ui_force) {
      local_u_[v][0] += velocity_[v * 3 + 0] * dt;
      local_u_[v][1] += velocity_[v * 3 + 1] * dt;
      local_u_[v][2] += velocity_[v * 3 + 2] * dt;
      local_r_[v] = local_u_[v] + vert_offset_from_mass_center_[v];

      //      local_u_[v][0] = 0;
      //      local_u_[v][1] = 0;
      //      local_u_[v][2] = 0;
      //      }
      //      sum_vel += local_u_[v][1];
      velocity_[v * 3 + 0] *= kDamping;
      velocity_[v * 3 + 1] *= kDamping;
      velocity_[v * 3 + 2] *= kDamping;
      //      MapVec3 map_rest_pos(rest_pos_ + v * 3);
      MapVec3 map_x(X + v * 3);
      map_x = rotation_matrix * (local_r_[v]) + center_of_mass_[p];
      inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
      inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
      inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
    }
  }

#else
  for (int v = 0; v < vertex_num_; ++v) {
    local_u_[v] += MapVec3(velocity_ + v * 3) * dt;
    local_r_[v] = local_u_[v] + vert_offset_from_mass_center_[v];
  }
  std::vector<Vec3> cm_offset(part_num_, Vec3(0, 0, 0));
  std::vector<Mat3> rotation(part_num_, Mat3::Zero());
  if (1) {
    for (int p = 0; p < part_num_; ++p) {
      if (p == 0) continue;
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int v = vert_local_id2global_id_[p][local_v];
        cm_offset[p] += mass_[v] * local_u_[v];
        for (int i = 0; i < 3; ++i) {
          rotation[p](i, 0) += mass_[v] * local_r_[v][i] * vert_offset_from_mass_center_[v][0];
          rotation[p](i, 1) += mass_[v] * local_r_[v][i] * vert_offset_from_mass_center_[v][1];
          rotation[p](i, 2) += mass_[v] * local_r_[v][i] * vert_offset_from_mass_center_[v][2];
        }
      }
      cm_offset[p] /= mass_per_part_[p];
      Mat3 rot;
      Mat3 scaling;
      dj::PolarDecompose<double>(rotation[p].data(), &rot(0, 0), &scaling(0, 0), 1e-8);
      rotation[p] = rot;
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int v = vert_local_id2global_id_[p][local_v];
        Vec3 old_pos = local_u_[v];
        local_u_[v] = rot.transpose() * local_r_[v] - vert_offset_from_mass_center_[v];
        //        local_u_[v] = local_r_[v] - vert_offset_from_mass_center_[v];
        local_u_[v] -= cm_offset[p];
        local_r_[v] = local_u_[v] + vert_offset_from_mass_center_[v];
        MapVec3(velocity_ + v * 3) += (local_u_[v] - old_pos) / dt;
      }
    }
  }

  //    memcpy(velocity_, &rhs[0][0], sizeof(double) * vertex_num_ * 3);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //  double sum_vel = 0;
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    translational_vel_[p] += acceleration[p] * dt + part_rotation_[p] * cm_offset[p] / dt;
    center_of_mass_[p] += translational_vel_[p] * dt;
    const double kDamping = 0.995;
    //    const double kDamping = 1.0;
    translational_vel_[p] *= kDamping;
#ifndef EXPLICIT_NO_ROTATION
    if (p == 1) {
      double angle_in_radian;
      Vec3 axis;
      auto rot = Quaternion<double>::Matrix2Quaternion(&rotation[p](0, 0));
      rot.GetRotation(&angle_in_radian, &axis[0]);
      //    dj::RotationMatrix2AngleAxis(rotation[p].data(), angle_in_radian, &axis[0]);
      PMAT(rotation[p]);
      P(axis, angle_in_radian);
      axis *= angle_in_radian;
      axis = part_rotation_[p] * axis;
      angular_vel_[p] += angular_acceleration[p] * dt;
      angular_vel_[p] += axis / dt;
    }
#endif
    angular_vel_[p] *= kDamping;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
    quaternion_[p].Quaternion2Matrix(rotation_matrix.data());
    part_rotation_[p] = rotation_matrix;
    part_rotation_transpose_[p] = part_rotation_[p].transpose();
    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
    current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();
    //    part_rotation_[p].setIdentity();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      MapVec3 map_x(X + v * 3);
      map_x = rotation_matrix * (local_r_[v]) + center_of_mass_[p];
      inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
      inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
      inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
    }
  }

#endif

  if (0) {
    P(dj::Vec3d(&translational_vel_[0][0]), dj::Vec3d(&angular_vel_[0][0]), dj::Vec3d(&center_of_mass_[0][0]));
    P(dj::Vec3d(&translational_vel_[1][0]), dj::Vec3d(&angular_vel_[1][0]), dj::Vec3d(&center_of_mass_[1][0]));
  }
}


void MultiDomainTet::FullSimulation(double dt) {
  static unsigned int steps = 0;
  steps++;
  //  std::fill(is_constrainted_.begin(), is_constrainted_.end(), false);
  //  for (int v = 0; v < vertex_num_; ++v) {
  //    //    if (X[v * 3 + 1] < 0.01) {
  //    //    }
  //    if (vert_part_id_[v] == 0) {
  //      is_constrainted_[v] = true;
  //    }
  //  }
  //  is_constrainted_[461] = false;
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  SparseMatrix* tangent_stiffness = inv_fem_->tangent_stiffness_matrix_;
  Real dt_2 = dt * dt;
  auto K = [&](Real * x, Real * result) {
    //    memcpy((Real*) result, (Real*) x, sizeof(Real) * 3 * v_num_);
    //    return;
    tangent_stiffness->MultiplyVector(x, result);
    OMP_FOR
    for (int v = 0; v < vertex_num_; ++v) {
      int v3 = v * 3;
      if (is_constrainted_[v]) {
        result[v3 + 0] = 0;
        result[v3 + 1] = 0;
        result[v3 + 2] = 0;
      } else {
        result[v3 + 0] = x[v3 + 0] * mass_[v] + dt_2 * result[v3 + 0];
        result[v3 + 1] = x[v3 + 1] * mass_[v] + dt_2 * result[v3 + 1];
        result[v3 + 2] = x[v3 + 2] * mass_[v] + dt_2 * result[v3 + 2];
      }
    }
  };

  double(*rhs)[3] = inv_fem_->rhs_;
  //  memset((double*) rhs, 0, sizeof(double) * vertex_num_ * 3);
  //  rhs[621][1] -= 500;
  //  if (dt * steps < 1.0) {
  //    rhs[621][1] -= 500;
  //  }
  for (int v = 0; v < vertex_num_; ++v) {
    if (is_constrainted_[v]) {
      rhs[v][0] = 0;
      rhs[v][1] = 0;
      rhs[v][2] = 0;
    } else {
      rhs[v][0] = dt * (rhs[v][0]) + inv_fem_->vel_[v][0] * mass_[v];
      rhs[v][1] = dt * (rhs[v][1]) + inv_fem_->vel_[v][1] * mass_[v];
      rhs[v][2] = dt * (rhs[v][2]) + inv_fem_->vel_[v][2] * mass_[v];
    }
  }

  double force[3];
  int applied_v = GetUIForce(force);
  if (applied_v > 0) {
    rhs[applied_v][0] += force[0] * dt;
    rhs[applied_v][1] += force[1] * dt;
    rhs[applied_v][2] += force[2] * dt;
  }
  if (0) {
    std::vector<int> visited(vertex_num_, false);
    double f = 200.8;
    //    double f[] = {30, 28, 25, 30};
    rhs[selected_vertex_][0] -= f * mass_[selected_vertex_] * dt;
    rhs[selected_vertex_][1] += 100 * mass_[selected_vertex_] * dt;
    std::queue<int> q;
    q.push(selected_vertex_);
    visited[selected_vertex_] = true;
    int node_at_next_level = 0;
    int count = 1;
    int depth = 1;
    while (depth <= 5 && !q.empty()) {
      int v = q.front();
      q.pop();
      for (int i : adjacent_vertex_[v]) {
        if (visited[i] == false) {
          q.push(i);
          rhs[i][0] -= f * mass_[i] * dt;
          rhs[i][1] += 100 * mass_[i] * dt;
          node_at_next_level++;
          visited[i] = true;
        }
      }
      count--;
      if (count == 0) {
        depth++;
        count = node_at_next_level;
        node_at_next_level = 0;
      }
    }
  }

  //  memcpy((double*)inv_fem_->vel_, (double*)rhs, sizeof(double) * 3 * vertex_num_);
  //  P(MapVec((double*) rhs, vertex_num_ * 3).norm());
  //  KK;
  auto info = cg_solver_->Solve((double*)rhs, (double*)inv_fem_->vel_, K, 1000, 1e-15);
  (void) info;
  P(info.first, info.second);

  double(*vel)[3] = inv_fem_->vel_;
  double(*u)[3] = inv_fem_->u_;
  //  OMP_FOR
  for (int v = 0; v < vertex_num_; ++v) {
    if (is_constrainted_[v]) {
      vel[v][0] = 0;
      vel[v][1] = 0;
      vel[v][2] = 0;
    } else {
      u[v][0] += vel[v][0] * dt;
      u[v][1] += vel[v][1] * dt;
      u[v][2] += vel[v][2] * dt;

      //      const double kDamping = 0.96;
      //      vel[v][0] *= kDamping;
      //      vel[v][1] *= kDamping;
      //      vel[v][2] *= kDamping;

      //    if (tet_mesh_->attached_joint_of_vertex_[v] >= 0)
      //    P(dj::Vec3d(u_[v]));
      X[v * 3 + 0] = u[v][0] + rest_pos_[v * 3 + 0];
      X[v * 3 + 1] = u[v][1] + rest_pos_[v * 3 + 1];
      X[v * 3 + 2] = u[v][2] + rest_pos_[v * 3 + 2];
    }
  }
}

void MultiDomainTet::AddFictitiousForceAndGravity(const std::vector<Vec3>& acceleration,
                                                  const std::vector<Vec3>& angular_acceleration,
                                                  MultiDomainTet::Vec & subspace_rhs) {
  for (int p = 0; p < part_num_; ++p) {
    Mat3 rotation_transpose = part_rotation_[p].transpose();
    Vec3 acc = rotation_transpose * (acceleration[p] - gravity_);
    Vec3 angular_acc = rotation_transpose * angular_acceleration[p];
    MapVec subspace_force(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
    MapVec map_vel_q(&vel_q_[basis_offset_[p]], part_basis_size_[p]);
    MapVec map_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // inertial force: F_v = -m_v * a,
    //                F(q) = \sum_{all v} -(U_v^T * m_v * a) = (\sum_{all v} -(m_v * U_v)) * a
    //            sandwich = (\sum_{all v} -(m_v * U_v))
    //                F(q) = intertial_force_sandwich_^{rx3) * acceleration^{3x1}
    subspace_force += inertial_force_sandwich_[p] * (acc); // gravity also added here
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Coriolis force: F = -2m w\times v
    //              F(q) = \sum_{all v} -2(U_v^T * m_v * [w] * U_v * dq/dt)
    //          sandwich = \sum_{all v} -(U_v^T * m_v * [w] * U_v)
    //              F(q) = 2 * (\sum_{i\in [0, 8]} coriolis_force_sandwich_[i] * [w](i / 3, i % 3)) * dq/dt
    Mat3 angular_vel_mat = GetSkewSymmetrixMatrix(rotation_transpose * angular_vel_[p]);
    double* data = angular_vel_mat.data();
    Mat coriolis_subspace_force = Mat::Zero(part_basis_size_[p], part_basis_size_[p]);
    for (int i = 0; i < 9; ++i) {
      coriolis_subspace_force += (coriolis_force_sandwich_[p * 9 + i] * data[i]);
    }
    subspace_force += 2 * (coriolis_subspace_force) * map_vel_q;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Euler force: F_v = -m_v * dw/dt \times r_v
    //             F(q) = \sum_{all v} -(U_v^T * m_v * [dw/dt] * (u + r0)
    //                  = \sum_{all v} {-(U_v^T * m_v * [dw/dt] * U * q) + (U_v^T * m_v * [r0] * dw/dt)}
    //  first term computed from Coriolis force sandwich
    //             sandwich = \sum_{all v} (U_v^T * m_v * [r0])
    //            F(q) = (\sum_{i\in [0,8]} coriolis_force_sandwich_[i] * [dw/dt](i / 3, i % 3)) * q
    //                   + euler_force_sandwich_^{r*3} * dw/dt^{3x1}
    Mat3 angular_acc_mat = GetSkewSymmetrixMatrix(angular_acc);
    data = angular_acc_mat.data();
    Mat euler_force = Mat::Zero(part_basis_size_[p], part_basis_size_[p]);
    for (int i = 0; i < 9; ++i) {
      euler_force += coriolis_force_sandwich_[p * 9 + i] * data[i];
    }
    subspace_force += euler_force * map_q + euler_force_sandwich_[p] * angular_acc;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // centerifugal force: F = -m * w \times (w \times r)
    //                  F(q) = \sum_{all v} -(U_v^T * m_v * [w] * [w] * (u + r0)
    //                       = \sum_{all v} {(-U_v^T * m_v * [w][w] * U_v * q) + (-U_v^T * m_v * [w] * [w] * r0)
    //  first term computed from Coriolis force sandwich
    //              sandwich = \sum_{all v} (-U_v^T * m_v * [w]^2 * r0)
    //                  F(q) = (\sum_{i\in [0,8]} coriolis_force_sandwich_[i] * ([w]*[w])(i / 3, i % 3)) * q
    //                         + (\sum_{i\in [0,8]} centrifugal_force_sandwich_[i] * ([w]*[w])(i / 3, i % 3))
    angular_vel_mat = angular_vel_mat * angular_vel_mat;
    data = angular_vel_mat.data();
    Mat centrifugal_force0 = Mat::Zero(part_basis_size_[p], part_basis_size_[p]);
    Vec centrifugal_force1 = Vec::Zero(part_basis_size_[p]);
    for (int i = 0; i < 9; ++i) {
      centrifugal_force0 += coriolis_force_sandwich_[p * 9 + i] * data[i];
      centrifugal_force1 += centrifugal_force_sandwich_[p * 9 + i] * data[i];
    }
    subspace_force += centrifugal_force0 * map_q + centrifugal_force1;
  }
}

void MultiDomainTet::GetExternalForce(double dt, std::vector<MultiDomainTet::ExtForce> *ext_force_ptr) {
  if (ext_force_ptr == nullptr) ext_force_ptr = &ext_force_;
  ext_force_ptr->clear();
  profiler.Start("collision force");
  // external collision force in world space to each vertex
  if (conf.Get<int>("ground collision")) {
    const double kFloor = -0.20;
    const double kFloorStiffness = 25000;
    for (int v = 0; v < vertex_num_; ++v) {
      // Ground collision
      if (X[v * 3 + 1] < kFloor) {
        Vec3 force(0, 0, 0);
        force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
        //      force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
        ext_force_ptr->emplace_back(v, force);
      }
    }
  }
  profiler.End("collision force");
  static unsigned int steps = 0;
  if (steps * dt < 1.0) {
    int applied_vertex = conf.Get<int>("selected_vertex");
    if (applied_vertex >= 0 && applied_vertex < vertex_num_) {
      double* force = conf.Get<double*>("constant force");
      ext_force_ptr->emplace_back(applied_vertex, Vec3(force[0], force[1], force[2]));
    }
  }
  steps++;
  // user interaction force
  {
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      ext_force_ptr->emplace_back(v, force);
    }
  }
}

MultiDomainTet::Vec MultiDomainTet::ComputeSubspaceRhs(const double dt,
                                                       std::vector<Vec3>& acceleration,
                                                       std::vector<Vec3>& angular_acceleration,
                                                       double internal_force_scaling_factor) {
  Vec subspace_rhs = Vec::Zero(total_basis_num_);
  for (std::pair<int, Vec3>& collision : ext_force_) {
    int& v = collision.first;
    Vec3& force = collision.second;
    int p = vert_part_id_[v];
    MapVec map_rhs(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
    map_rhs += vert_basis_transpose_[v] * (part_rotation_transpose_[p] * force);
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // apply fictitious force
  profiler.Start("fictitious_force");
  AddFictitiousForceAndGravity(acceleration, angular_acceleration, subspace_rhs);
  profiler.End("fictitious_force");

  profiler.Start("force");
  // internal force
  for (CubaturePoint & cubature : all_cubature_) {
    int& t = cubature.first;
    double& weight = cubature.second;
    // reduced internal force
    double* force = inv_fem_->element_force_ + t * 12;
    for (int local_v = 0; local_v < 4; ++local_v) {
      int v = tet_[t * 4 + local_v];
      int p = vert_part_id_[v];
      MapVec3 map_force(force + local_v * 3);
      MapVec subspace_force(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
      // negative sign before weight is because internal force from vega is in oppisite direction of actual force
      subspace_force -= vert_basis_transpose_[v] * (internal_force_scaling_factor * weight * (part_rotation_[p].transpose() * map_force));
    }
  }
  profiler.End("force");
  // implicit terms
  if (0)
    for (int e = 0; e < interface_num_; ++e) {
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            MapVec(&subspace_rhs[basis_offset_[pi]], part_basis_size_[pi]) -=
              vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
              (element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] + angular_vel_[pj].cross(rj)) * weight * dt);
          }
        }
      }
    }
  subspace_rhs *= dt;
  subspace_rhs += vel_q_;
  return std::move(subspace_rhs);
}

void MultiDomainTet::ComputeSubspaceStiffnessMatrix(const double dt, double internal_force_scaling_factor) {
  const double dt_2 = dt * dt;
  chol_solver_->SetMatrixZero();
  profiler.Start("mul p");
  //  if (0)
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    double* mat = chol_solver_->Av[p].A;
    Eigen::Map<MatRow> reduced_k(mat, part_basis_size_[p], part_basis_size_[p]);
    // all cubature vertex in this domain
    for (int j = 0; j < int(v_list_[p].size()); ++j) {
      int v = v_list_[p][j];
      ASSERT(vert_part_id_[v] == p);
      int p = vert_part_id_[v];
      Mat3 sub_k = Mat3::Zero();
      // all cubature tet incident on this v
      for (int k = 0; k < int(v_cubature_[p][j].size()); ++k) {
        VVCubature& cubature = v_cubature_[p][j][k];
        const int& t = cubature.cubature_tet;
        const double& weight = cubature.cubature_weight;
        const int& idx = cubature.i;
        //        if (p == 0) {
        //          printf("%d\t%d\t%lf\t%d\n", v, t, weight, idx);
        //          std::cout << v << "\t" << t << "\t" << weight << "\t" << idx << "\n";
        //          P(v, t, weight, idx);
        //        }
        double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            sub_k(r, c) += weight * tet_k[idx * 3 + r][idx * 3 + c];
          }
        }
      }
      reduced_k += vert_basis_transpose_[v]
                   * (part_rotation_transpose_[p] * sub_k * part_rotation_[p])
                   * vert_basis_[v];
    }
  }
  profiler.End("mul p");

  profiler.Start("mul e");
  //  if (0)
  OMP_FOR
  for (int e = 0; e < part_num_ + chol_solver_->e_number0; ++e) {
    int p0, p1;
    double* mat = NULL;
    if (e < part_num_) {
      p0 = p1 = e;
      mat = chol_solver_->Av[e].A;
    } else {
      p0 = chol_solver_->E[(e - part_num_) * 2 + 0];
      p1 = chol_solver_->E[(e - part_num_) * 2 + 1];
      mat = chol_solver_->Ae[e - part_num_].A;
      mat = chol_solver_->Ae[e - part_num_].A;
    }
    Eigen::Map<MatRow> reduced_k(mat, part_basis_size_[p0], part_basis_size_[p1]);
    // list of vert-vert pair in this domain
    for (int i = 0; i < int(vv_list_[e].size()); ++i) {
      int vi = vv_list_[e][i].first;
      int vj = vv_list_[e][i].second;
      // list of  cubature tets that contains this vert-vert pair
      Mat3 sub_k = Mat3::Zero();
      for (int k = 0; k < int(vv_cubature_[e][i].size()); ++k) {
        VVCubature& cubature = vv_cubature_[e][i][k];
        const int& t = cubature.cubature_tet;
        const double& weight = cubature.cubature_weight;
        const int& idx0 = cubature.i;
        const int& idx1 = cubature.j;
        double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            sub_k(r, c) += weight * tet_k[idx0 * 3 + r][idx1 * 3 + c];
          }
        }
      }
      Mat k = vert_basis_transpose_[vi]
              * (part_rotation_transpose_[p0] * sub_k * part_rotation_[p1])
              * vert_basis_[vj];
      reduced_k += k;
      if (e < part_num_) {
        reduced_k += k.transpose();
      }
    }
  }
  profiler.End("mul e");

  profiler.Start("subspace diag & scale");
  //  if (0)
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapMat result(chol_solver_->Av[p].A, part_basis_size_[p], part_basis_size_[p]);
    result *= dt_2 * internal_force_scaling_factor;
    for (int i = 0; i < part_basis_size_[p]; ++i) {
      result(i, i) += 1;
    }
  }

  //  if (0)
  OMP_FOR
  for (int e = 0; e < chol_solver_->e_number; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    MapMat result(chol_solver_->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1]);
    result *= (dt_2 * internal_force_scaling_factor);
  }
  profiler.End("subspace diag & scale");
}

MultiDomainTet::MatCol MultiDomainTet::ComputeFullRigidStiffnessMatrixWithCubature(const double dt) {
  MatCol stiffness_matrix = MatCol::Zero(part_num_ * 6, part_num_ * 6);
  for (int e = 0; e < interface_num_; ++e) {
    for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
      int& t = cubature.first;
      double& weight = cubature.second;
      // TODO delete
      //      weight = 0.4;
      // the force is computed from vega, the direction of the force is opposite to the actual direction
      Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
      int* verts = tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = vert_part_id_[vi];
        Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
        auto skew_symmetric_mat_i = GetSkewSymmetrixMatrix(ri);
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int pj = vert_part_id_[vj];
          Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
          auto skew_symmetric_mat_j = GetSkewSymmetrixMatrix(rj);
          stiffness_matrix.block<3, 3>(pi * 6 + 0, pj * 6 + 0) += element_k.block<3, 3>(i * 3, j * 3) * weight;
          stiffness_matrix.block<3, 3>(pi * 6 + 0, pj * 6 + 3) -= element_k.block<3, 3>(i * 3, j * 3) * skew_symmetric_mat_j * weight;
          stiffness_matrix.block<3, 3>(pi * 6 + 3, pj * 6 + 0) += skew_symmetric_mat_i * element_k.block<3, 3>(i * 3, j * 3) * weight;
          stiffness_matrix.block<3, 3>(pi * 6 + 3, pj * 6 + 3) -= skew_symmetric_mat_i * element_k.block<3, 3>(i * 3, j * 3) * skew_symmetric_mat_j * weight;
        }
      }
    }
  }
  stiffness_matrix *= (dt * dt);
  // Diagonal terms
  for (int p = 0; p < part_num_; ++p) {
    stiffness_matrix(p * 6 + 0, p * 6 + 0) += mass_per_part_[p];
    stiffness_matrix(p * 6 + 1, p * 6 + 1) += mass_per_part_[p];
    stiffness_matrix(p * 6 + 2, p * 6 + 2) += mass_per_part_[p];
    stiffness_matrix.block<3, 3>(p * 6 + 3, p * 6 + 3) += current_inertia_tensor_[p];
  }
  return std::move(stiffness_matrix);
}

MultiDomainTet::MatCol MultiDomainTet::ComputeFullRigidStiffnessMatrix(const double dt) {
  MatCol stiffness_matrix = MatCol::Zero(part_num_ * 6, part_num_ * 6);
  for (int t : interface_tet_) {
    Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
    int* verts = tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      int vi = verts[i];
      int pi = vert_part_id_[vi];
      Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
      auto skew_symmetric_mat_i = GetSkewSymmetrixMatrix(ri);
      for (int j = 0; j < 4; ++j) {
        int vj = verts[j];
        int pj = vert_part_id_[vj];
        Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
        auto skew_symmetric_mat_j = GetSkewSymmetrixMatrix(rj);
        stiffness_matrix.block<3, 3>(pi * 6 + 0, pj * 6 + 0) += element_k.block<3, 3>(i * 3, j * 3);
        stiffness_matrix.block<3, 3>(pi * 6 + 0, pj * 6 + 3) -= element_k.block<3, 3>(i * 3, j * 3) * skew_symmetric_mat_j;
        stiffness_matrix.block<3, 3>(pi * 6 + 3, pj * 6 + 0) += skew_symmetric_mat_i * element_k.block<3, 3>(i * 3, j * 3);
        stiffness_matrix.block<3, 3>(pi * 6 + 3, pj * 6 + 3) -= skew_symmetric_mat_i * element_k.block<3, 3>(i * 3, j * 3) * skew_symmetric_mat_j;
      }
    }
  }
  stiffness_matrix *= (dt * dt);
  // explicit
  //  const double kImplicitness = 1.0;
  //  stiffness_matrix *= kImplicitness;
  // Diagonal terms
  for (int p = 0; p < part_num_; ++p) {
    stiffness_matrix(p * 6 + 0, p * 6 + 0) += mass_per_part_[p];
    stiffness_matrix(p * 6 + 1, p * 6 + 1) += mass_per_part_[p];
    stiffness_matrix(p * 6 + 2, p * 6 + 2) += mass_per_part_[p];
    stiffness_matrix.block<3, 3>(p * 6 + 3, p * 6 + 3) += current_inertia_tensor_[p];
  }
  return std::move(stiffness_matrix);
}

std::vector<MultiDomainTet::Vec3> MultiDomainTet::ComputeFullRigidMotionRhs(const double dt,
                                                                            std::vector<MultiDomainTet::ExtForce> *ext_force_ptr) {
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // net force on each part
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  // net torque on each part
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));
  // elastic force from boundary tets
  for (int t : interface_tet_) {
    // the force is computed from vega, the direction of the force is opposite to the actual direction
    double* force = inv_fem_->element_force_ + t * 12;
    for (int i = 0; i < 4; ++i) {
      int v = tet_[t * 4 + i];
      int p = vert_part_id_[v];
      MapVec3 map_force(force + i * 3);
      net_force[p] -= map_force;
      MapVec3 map_x(X + v  * 3);
      Vec3 r = map_x - center_of_mass_[p];
      net_torque[p] -= r.cross(map_force);
    }

#if 0
    {
      Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = tet_[t * 4 + i];
        int pi = vert_part_id_[vi];
        Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
        for (int j = 0; j < 4; ++j) {
          int vj = tet_[t * 4 + j];
          int pj = vert_part_id_[vj];
          Vec3 vel = vert_basis_[vj] * MapVec(&vel_q_[basis_offset_[pj]], part_basis_size_[pj]);
          vel = part_rotation_[pj] * (vel * dt);
          Vec force = element_k.block<3, 3>(i * 3, j * 3) * vel;
          net_force[pi] -= force;
          net_torque[pi] -= ri.cross(force);
        }
      }
    }
#endif
  }

  // gravity on each part
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += mass_per_part_[p] * gravity_;
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external collision force in world space to each vertex
  {
    if (ext_force_ptr == nullptr) {
      ext_force_ptr = &ext_force_;
    }
    for (const ExtForce & f : *ext_force_ptr) {
      const int& v = f.first;
      const Vec3& force = f.second;
      int p = vert_part_id_[v];
      net_force[p] += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[p] += r.cross(force);
    }
  }

  std::vector<Vec3> rhs(part_num_ * 2);
  //  std::vector<Mat3> current_inertia_tensor_ = inertia_tensor_;
  for (int p = 0; p < part_num_; ++p) {
    //    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_[p].transpose();
    rhs[p * 2 + 0] = net_force[p] * dt + translational_vel_[p] * mass_per_part_[p];
    rhs[p * 2 + 1] = net_torque[p] * dt + current_inertia_tensor_[p] * angular_vel_[p];
  }
  return std::move(rhs);
}

std::vector<MultiDomainTet::Vec3> MultiDomainTet::ComputeFullRigidMotionRhsWithCubature(const double dt, std::vector<MultiDomainTet::ExtForce> *ext_force_ptr) {
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // net force on each part
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  // net torque on each part
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));
  // elastic force from boundary tets
  for (int e = 0; e < interface_num_; ++e) {
    for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
      int& t = cubature.first;
      double& weight = cubature.second;
      // the force is computed from vega, the direction of the force is opposite to the actual direction
      double* force = inv_fem_->element_force_ + t * 12;
      for (int i = 0; i < 4; ++i) {
        int v = tet_[t * 4 + i];
        int p = vert_part_id_[v];
        Vec3 cubature_force = MapVec3(force + i * 3) * weight;
        net_force[p] -= cubature_force;
        MapVec3 map_x(X + v  * 3);
        Vec3 r = map_x - center_of_mass_[p];
        net_torque[p] -= r.cross(cubature_force);
      }
    }
  }

  // gravity on each part
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += mass_per_part_[p] * gravity_;
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external collision force in world space to each vertex
  {
    if (ext_force_ptr == nullptr) {
      ext_force_ptr = &ext_force_;
    }
    for (const ExtForce & f : *ext_force_ptr) {
      const int& v = f.first;
      const Vec3& force = f.second;
      int p = vert_part_id_[v];
      net_force[p] += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[p] += r.cross(force);
    }
  }

  std::vector<Vec3> rhs(part_num_ * 2);
  for (int p = 0; p < part_num_; ++p) {
    rhs[p * 2 + 0] = net_force[p] * dt + translational_vel_[p] * mass_per_part_[p];
    rhs[p * 2 + 1] = net_torque[p] * dt + current_inertia_tensor_[p] * angular_vel_[p];
  }
  return std::move(rhs);
}

void MultiDomainTet::SimulateRigidMotionFull(const double dt,
                                             MultiDomainTet::Vec & new_rigid_velocity,
                                             std::vector<MultiDomainTet::ExtForce> *ext_force_ptr) {
  std::vector<Vec3> rhs = ComputeFullRigidMotionRhs(dt, ext_force_ptr);
  MatCol RigidStiffnessMatrix = ComputeFullRigidStiffnessMatrix(dt);
  profiler.Start("rigid solve");
  // has fixed domain
  if (fixed_domain_ >= 0 && fixed_domain_ < part_num_) {
    int size0 = fixed_domain_ * 6;
    int size1 = (part_num_ - fixed_domain_ - 1) * 6;
    Vec fixed_rhs = Vec::Zero((part_num_ - 1) * 6);
    MatCol fixed_k = MatCol::Zero((part_num_ - 1) * 6, (part_num_ - 1) * 6);
    if (fixed_domain_ > 0) {
      memcpy(&fixed_rhs[0], &rhs[0][0], sizeof(double) * size0);
      fixed_k.block(0, 0, size0, size0) = RigidStiffnessMatrix.block(0, 0, size0, size0);
    }
    if (fixed_domain_ < part_num_ - 1) {
      memcpy(&fixed_rhs[size0], &rhs[0][0] + size0 + 6, sizeof(double) * size1);
      fixed_k.block(size0, size0, size1, size1) = RigidStiffnessMatrix.block(size0 + 6, size0 + 6, size1, size1);
    }
    Vec fixed_vel = fixed_k.colPivHouseholderQr().solve(fixed_rhs);
    if (fixed_domain_ > 0) {
      memcpy(&new_rigid_velocity[0], &fixed_vel[0], sizeof(double) * size0);
    }
    if (fixed_domain_ < part_num_ - 1) {
      memcpy(&new_rigid_velocity[size0 + 6], &fixed_vel[size0], sizeof(double) * size1);
    }
  } else {
    MapVec map_rhs(&rhs[0][0], part_num_ * 6);
    new_rigid_velocity = RigidStiffnessMatrix.colPivHouseholderQr().solve(map_rhs);
  }
  profiler.End("rigid solve");
}

void MultiDomainTet::SimulateRigidMotion(const double dt,
                                         MultiDomainTet::Vec & new_rigid_velocity,
                                         std::vector<MultiDomainTet::ExtForce> *ext_force_ptr) {
  double rigid_coupling_scaling = conf.Get<double>("rigid coupling scaling");
  //  double rigid_coupling_scaling = 1;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // net force on each part
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  // net torque on each part
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));
  // gravity on each part
  std::vector<Vec> subspace_interface_force(part_num_);
  std::vector<MatCol>& basis = rigid_basis_;
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += mass_per_part_[p] * gravity_;
    subspace_interface_force[p] = Vec::Zero(basis[p].cols());
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // elastic force from boundary tets
  profiler.Start("boundary force");
  for (int e = 0; e < interface_num_; ++e) {
    // FIXME: ignor interface forces that come from interface of more than 2 domains
    if (interface_domains_[e][0] > 2) continue;
    const int p0 = interface_domains_[e][1];
    const int p1 = interface_domains_[e][2];
    Vec force0 = Vec::Zero(basis[p0].cols());
    Vec force1 = Vec::Zero(basis[p1].cols());
    for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
      int& t = cubature.first;
      double& weight = cubature.second;

      double* element_force = inv_fem_->element_force_ + t * 3 * 4;
      for (int local_v = 0; local_v < 4; ++local_v) {
        MapVec3 force(element_force + local_v * 3);
        int v = tet_[t * 4 + local_v];
        int p = vert_part_id_[v];
        // negative sign before weight is because internal force from vega is in oppisite direction of actual force
        if (p == p0) {
          force0 += vert_rigid_basis_transpose_[v] * (-weight * (part_rotation_transpose_[p] * force));
        } else {
          ASSERT(p == p1);
          force1 += vert_rigid_basis_transpose_[v] * (-weight * (part_rotation_transpose_[p] * force));
        }
      }
    }
    if (0) {
      L("before");
      int e0 = e * 2 + 0; int e1 = e * 2 + 1;
      Vec3 r0 = part_rotation_[p0] * momentum_matrix_[p0] * force0 + part_rotation_[p1] * momentum_matrix_[p1] * force1;
      Vec3 r1 = part_rotation_[p0] * interface_torque_matrix_[e0] * force0 + part_rotation_[p1] * interface_torque_matrix_[e1] * force1;
      P(force0.norm(), force1.norm());
      P(r0.norm(), r1.norm(), p0, p1);
    }
    ProjectInterfaceSubspaceForce(e, p0, p1, force0, force1);
    // verify projection
    if (0) {
      int e0 = e * 2 + 0; int e1 = e * 2 + 1;
      Vec3 r0 = part_rotation_[p0] * momentum_matrix_[p0] * force0 + part_rotation_[p1] * momentum_matrix_[p1] * force1;
      Vec3 r1 = part_rotation_[p0] * interface_torque_matrix_[e0] * force0 + part_rotation_[p1] * interface_torque_matrix_[e1] * force1;
      P(r0.norm(), r1.norm(), p0, p1);
      P(force0.norm(), force1.norm());
      ASSERT(r0.norm() < 1e-6, P(r0.norm()));
      ASSERT(r1.norm() < 1e-6, P(r1.norm()));
    }
    subspace_interface_force[p0] += force0 * rigid_coupling_scaling;
    subspace_interface_force[p1] += force1 * rigid_coupling_scaling;
  }
  profiler.End("boundary force");

  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Mat3& rotation = part_rotation_[p];
    net_force[p] += rotation * momentum_matrix_[p] * subspace_interface_force[p];
    net_torque[p] += rotation * torque_matrix_[p] * subspace_interface_force[p];
  }


  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external collision force in world space to each vertex
  {
    if (ext_force_ptr == nullptr) {
      ext_force_ptr = &ext_force_;
    }
    for (const ExtForce & f : *ext_force_ptr) {
      const int& v = f.first;
      const Vec3& force = f.second;
      int p = vert_part_id_[v];
      net_force[p] += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[p] += r.cross(force);
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // compute rigid transformation
  std::vector<Vec3> rhs(part_num_ * 2);
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    rhs[p * 2 + 0] = net_force[p] * dt + translational_vel_[p] * mass_per_part_[p];
    rhs[p * 2 + 1] = net_torque[p] * dt + current_inertia_tensor_[p] * angular_vel_[p];
  }

  profiler.Start("rigid k");
#if 0
  //  Mat truth_rigid_k;
  //  {
  std::vector<Vec3> implicit_vel(vertex_num_);
  std::vector<Vec3> implicit_force(vertex_num_);
  auto RigidMotionMatrix = [&](double * x, double * result) {
    for (int v : interface_cubature_vert_) {
      int p = vert_compact_part_id_[v];
      MapVec3 map_pos(X + v * 3);
      //        Vec3 r = map_pos - center_of_mass_[p];
      Vec3 r = part_rotation_[p] * vert_offset_from_mass_center_[v];
      double angular_vel[3] = {0, 0, 0};//= angular_vel_[p].cross(r);
      dj::Cross3(x + p * 6 + 3, &r[0], &angular_vel[0]);
      implicit_vel[v][0] = x[p * 6 + 0] + angular_vel[0];
      implicit_vel[v][1] = x[p * 6 + 1] + angular_vel[1];
      implicit_vel[v][2] = x[p * 6 + 2] + angular_vel[2];
      implicit_force[v].setZero();
    }

    //      if (0)
    for (int p = 0; p < part_num_; ++p) {
      result[p * 6 + 0] = x[p * 6 + 0] * mass_per_part_[p];
      result[p * 6 + 1] = x[p * 6 + 1] * mass_per_part_[p];
      result[p * 6 + 2] = x[p * 6 + 2] * mass_per_part_[p];

      MapVec3 map_result(result + p * 6 + 3);
      MapVec3 map_x(x + p * 6 + 3);
      map_result = current_inertia_tensor_[p] * map_x;
    }
    // TODO delete
    //              return;
    for (int p = 0; p < part_num_; ++p) {
      subspace_interface_force[p].setZero();
    }
    for (int e = 0; e < solver_->e_number0; ++e) {
      int p0 = solver_->E[e * 2 + 0];
      int p1 = solver_->E[e * 2 + 1];
      Vec force0 = Vec::Zero(part_basis_size_[p0]);
      Vec force1 = Vec::Zero(part_basis_size_[p1]);
      for (CubaturePoint & cubature : inter_domain_cubature_[e]) {
        int& t = cubature.first;
        double& weight = cubature.second;
        // the force is computed from vega, the direction of the force is opposite to the actual direction
        double (*element_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          Vec3 force(0, 0, 0);
          int v = tet_[t * 4 + i];
          int p = vert_compact_part_id_[v];
          // compute single tet element force on vertex v
          for (int j = 0; j < 4; ++j) {
            int vj = tet_[t * 4 + j];
            force[0] += dj::Dot3(&element_k[i * 3 + 0][j * 3], &implicit_vel[vj][0]);
            force[1] += dj::Dot3(&element_k[i * 3 + 1][j * 3], &implicit_vel[vj][0]);
            force[2] += dj::Dot3(&element_k[i * 3 + 2][j * 3], &implicit_vel[vj][0]);
          }
          if (p == p0) {
            force0 += vert_basis_transpose_[v] * (weight * (part_rotation_[p].transpose() * force));
          } else {
            force1 += vert_basis_transpose_[v] * (weight * (part_rotation_[p].transpose() * force));
          }

        }
      }
      ProjectInterfaceSubspaceForce(e, p0, p1, force0, force1);
      subspace_interface_force[p0] += force0 * rigid_coupling_scaling;
      subspace_interface_force[p1] += force1 * rigid_coupling_scaling;
    }

    for (int p = 0; p < part_num_; ++p) {
      Vec3 momemtum = part_rotation_[p] * momentum_matrix_[p] * subspace_interface_force[p];
      momemtum *= dt_2;
      result[p * 6 + 0] += momemtum[0];
      result[p * 6 + 1] += momemtum[1];
      result[p * 6 + 2] += momemtum[2];

      Vec3 torque = part_rotation_[p] * torque_matrix_[p] * subspace_interface_force[p];
      torque *= dt_2;
      result[p * 6 + 3] += torque[0];
      result[p * 6 + 4] += torque[1];
      result[p * 6 + 5] += torque[2];
    }
  };

  MatCol Rigid = MatCol::Zero(part_num_ * 6, part_num_ * 6);
  Vec identity = Vec::Zero(part_num_ * 6);
  for (int i = 0; i < part_num_ * 6; ++i) {
    identity[i] = 1;
    RigidMotionMatrix(&identity[0], Rigid.col(i).data());
    identity[i] = 0;
  }
  //    truth_rigid_k = Rigid;
  //  }
#else
  MatCol Rigid = MatCol::Zero(part_num_ * 6, part_num_ * 6);
  // Off-diagonal terms
  OMP_FOR
  for (int e = 0; e < interface_num_; ++e) {
    if (interface_domains_[e][0] > 2) continue;
    const int p0 = interface_domains_[e][1];
    const int p1 = interface_domains_[e][2];

    Mat df_dt[2][2] = {
      {Mat::Zero(basis[p0].cols(), 3), Mat::Zero(basis[p1].cols(), 3)},
      {Mat::Zero(basis[p1].cols(), 3), Mat::Zero(basis[p1].cols(), 3)}
    }; // translation

    Mat df_dr[2][2] = {
      {Mat::Zero(basis[p0].cols(), 3), Mat::Zero(basis[p0].cols(), 3)},
      {Mat::Zero(basis[p1].cols(), 3), Mat::Zero(basis[p1].cols(), 3)}
    }; // rotation

    for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
      int& t = cubature.first;
      double& weight = cubature.second;
      // the force is computed from vega, the direction of the force is opposite to the actual direction
      typedef Eigen::Matrix<double, 12, 12, Eigen::RowMajor> Mat12;
      Eigen::Map<Mat12> element_k((double*) (inv_fem_->element_k_ + t * 144));
      //      double (*element_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = tet_[t * 4 + i];
        int pi = vert_part_id_[vi];
        int idx0 = (pi == p0) ? 0 : 1;
        // compute single tet element force on vertex v
        for (int j = 0; j < 4; ++j) {
          int vj = tet_[t * 4 + j];
          int pj = vert_part_id_[vj];
          //          Vec3 r = center_of_mass_[pj] - MapVec3(X + vj * 3); // -r;
          Vec3 r = -part_rotation_[pj] * vert_offset_from_mass_center_[vj];
          Mat3 skew_symmetrix_mat = GetSkewSymmetrixMatrix(r);
          int idx1 = (pj == p0) ? 0 : 1;
          df_dt[idx0][idx1] += vert_rigid_basis_transpose_[vi] * (weight * part_rotation_transpose_[pi] * element_k.block<3, 3>(i * 3, j * 3));
          df_dr[idx0][idx1] += vert_rigid_basis_transpose_[vi] * (weight * part_rotation_transpose_[pi] * (element_k.block<3, 3>(i * 3, j * 3) * skew_symmetrix_mat));
        }
      }
    }

    for (int i = 0; i < 2; ++i) {
      for (int col = 0; col < 3; ++col) {
        {
          auto f0 = df_dt[0][i].col(col);
          auto f1 = df_dt[1][i].col(col);
          ProjectInterfaceSubspaceForce(e, p0, p1, f0, f1);
        }
        {
          auto f0 = df_dr[0][i].col(col);
          auto f1 = df_dr[1][i].col(col);
          ProjectInterfaceSubspaceForce(e, p0, p1, f0, f1);
        }
      }
    }
    //    Rigid.block<3, 3>(p0 * 6 + 0, p0 * 6 + 0) += part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][0];
    //    Rigid.block<3, 3>(p0 * 6 + 0, p0 * 6 + 3) += part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][0];
    //    Rigid.block<3, 3>(p0 * 6 + 3, p0 * 6 + 0) += part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][0];
    //    Rigid.block<3, 3>(p0 * 6 + 3, p0 * 6 + 3) += part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][0];
    //    Rigid.block<3, 3>(p1 * 6 + 0, p1 * 6 + 0) += part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][1];
    //    Rigid.block<3, 3>(p1 * 6 + 0, p1 * 6 + 3) += part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][1];
    //    Rigid.block<3, 3>(p1 * 6 + 3, p1 * 6 + 0) += part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][1];
    //    Rigid.block<3, 3>(p1 * 6 + 3, p1 * 6 + 3) += part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][1];
    tmp_interface_rigid_k_[e * 8 + 0] = part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][0];
    tmp_interface_rigid_k_[e * 8 + 1] = part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][0];
    tmp_interface_rigid_k_[e * 8 + 2] = part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][0];
    tmp_interface_rigid_k_[e * 8 + 3] = part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][0];
    tmp_interface_rigid_k_[e * 8 + 4] = part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][1];
    tmp_interface_rigid_k_[e * 8 + 5] = part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][1];
    tmp_interface_rigid_k_[e * 8 + 6] = part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][1];
    tmp_interface_rigid_k_[e * 8 + 7] = part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][1];


    Rigid.block<3, 3>(p0 * 6 + 0, p1 * 6 + 0) += part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][1];
    Rigid.block<3, 3>(p0 * 6 + 0, p1 * 6 + 3) += part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][1];
    Rigid.block<3, 3>(p0 * 6 + 3, p1 * 6 + 0) += part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][1];
    Rigid.block<3, 3>(p0 * 6 + 3, p1 * 6 + 3) += part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][1];
    Rigid.block<3, 3>(p1 * 6 + 0, p0 * 6 + 0) += part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][0];
    Rigid.block<3, 3>(p1 * 6 + 0, p0 * 6 + 3) += part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][0];
    Rigid.block<3, 3>(p1 * 6 + 3, p0 * 6 + 0) += part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][0];
    Rigid.block<3, 3>(p1 * 6 + 3, p0 * 6 + 3) += part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][0];
  }


  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    for (std::pair<int, int>& edge_idx : domain_incident_interface_[p]) {
      int& e = edge_idx.first;
      int i = edge_idx.second * 4;
      Rigid.block<3, 3>(p * 6 + 0, p * 6 + 0) += tmp_interface_rigid_k_[e * 8 + i + 0];
      Rigid.block<3, 3>(p * 6 + 0, p * 6 + 3) += tmp_interface_rigid_k_[e * 8 + i + 1];
      Rigid.block<3, 3>(p * 6 + 3, p * 6 + 0) += tmp_interface_rigid_k_[e * 8 + i + 2];
      Rigid.block<3, 3>(p * 6 + 3, p * 6 + 3) += tmp_interface_rigid_k_[e * 8 + i + 3];
    }
  }
  Rigid *= (rigid_coupling_scaling * dt * dt);
  //    Rigid *= 0;
  // Diagonal terms
  //  if (0)
  for (int p = 0; p < part_num_; ++p) {
    Rigid(p * 6 + 0, p * 6 + 0) += mass_per_part_[p];
    Rigid(p * 6 + 1, p * 6 + 1) += mass_per_part_[p];
    Rigid(p * 6 + 2, p * 6 + 2) += mass_per_part_[p];
    Rigid.block<3, 3>(p * 6 + 3, p * 6 + 3) += current_inertia_tensor_[p];
  }

  //  Mat diff = truth_rigid_k - Rigid;
  //  WriteEigenMatrixToMatlab(diff, " / tmp / log / diff");
  //  WriteEigenMatrixToMatlab(truth_rigid_k, " / tmp / log / truth");
  //  WriteEigenMatrixToMatlab(Rigid, " / tmp / log / test");
  //  Rigid = truth_rigid_k;
  //    P(diff.maxCoeff());
  //    P(diff.minCoeff());
  //  P(diff.norm());
  //  ASSERT(diff.norm() < 1e-6, P(diff.norm()));
#endif
  profiler.End("rigid k");
  new_rigid_velocity = Vec::Zero(part_num_ * 6);
  profiler.Start("rigid solve");
  if (fixed_domain_ >= 0 && fixed_domain_ < part_num_) {
    int size0 = fixed_domain_ * 6;
    int size1 = (part_num_ - fixed_domain_ - 1) * 6;
    Vec fixed_rhs = Vec::Zero((part_num_ - 1) * 6);
    MatCol fixed_k = MatCol::Zero((part_num_ - 1) * 6, (part_num_ - 1) * 6);
    if (fixed_domain_ > 0) {
      memcpy(&fixed_rhs[0], &rhs[0][0], sizeof(double) * size0);
      fixed_k.block(0, 0, size0, size0) = Rigid.block(0, 0, size0, size0);
    }
    if (fixed_domain_ < part_num_ - 1) {
      memcpy(&fixed_rhs[size0], &rhs[0][0] + size0 + 6, sizeof(double) * size1);
      fixed_k.block(size0, size0, size1, size1) = Rigid.block(size0 + 6, size0 + 6, size1, size1);
    }
    //    PMATCOL(fixed_k);
    Vec fixed_vel = fixed_k.colPivHouseholderQr().solve(fixed_rhs);
    //    P(fixed_k.norm());
    //    P(fixed_vel.norm());
    if (fixed_domain_ > 0) {
      memcpy(&new_rigid_velocity[0], &fixed_vel[0], sizeof(double) * size0);
    }
    if (fixed_domain_ < part_num_ - 1) {
      memcpy(&new_rigid_velocity[size0 + 6], &fixed_vel[size0], sizeof(double) * size1);
    }
  } else {
    MapVec map_rhs(&rhs[0][0], part_num_ * 6);
    new_rigid_velocity = Rigid.colPivHouseholderQr().solve(map_rhs);
  }
  profiler.End("rigid solve");
}

void MultiDomainTet::SimulateWithGlobalSubspace(double dt) {
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  MatCol k_dot_u = MatCol::Zero(vertex_num_ * 3, total_basis_num_);
  for (int b = 0; b < total_basis_num_; ++b) {
    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(global_basis_.col(b).data(), k_dot_u.col(b).data());
  }
  MatCol reduced_k = global_basis_transpose_ * k_dot_u;
  reduced_k *= dt * dt;
  reduced_k += MatCol::Identity(total_basis_num_, total_basis_num_);
  //  std::ofstream out(" / tmp / k");
  //  out << reduced_k;
  //  out.close();
  MapVec full_force((double*)inv_fem_->rhs_, vertex_num_ * 3);
  //  full_force[selected_vertex_ * 3 + 2] += pose_conf.Get<double>("force");
  //  full_force[selected_vertex_ * 3 + 1] -= pose_conf.Get<double>("force");
  Vec reduced_force = global_basis_transpose_ * full_force;
  reduced_force *= dt;
  reduced_force += vel_q_;
  vel_q_ = reduced_k.colPivHouseholderQr().solve(reduced_force);
  P(vel_q_.dot(vel_q_));
  q_ += vel_q_ * dt;
  vel_q_ *= 0.96;
  ComputeGlobalPositionFromSubspace();
  //  PV(vel_q_, vel_q_.size());
  MapVec map_X(X, vertex_num_ * 3);
  MapVec map_rest_pos_(rest_pos_, vertex_num_ * 3);
  MapVec map_u((double*)inv_fem_->u_, vertex_num_ * 3);
  map_u = map_X - map_rest_pos_;
  //  map_u = global_basis_ * q_;
  //  map_X = map_u + map_rest_pos_;
}

void MultiDomainTet::SubspaceMultiBodySimulation(double dt) {
  const double dt_2 = dt * dt;
  // accumulate external force and torque for each part
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // net force on each part
  std::vector<Vec3> net_force(part_num_, Vec3(0, 0, 0));
  // net torque on each part
  std::vector<Vec3> net_torque(part_num_, Vec3(0, 0, 0));
  // gravity on each part
  for (int p = 0; p < part_num_; ++p) {
    net_force[p] += mass_per_part_[p] * gravity_;
  }

  // elastic force from boundary tets
  for (int t : interface_tet_) {
    // the force is computed from vega, the direction of the force is opposite to the actual direction
    double* force = inv_fem_->element_force_ + t * 12;
    for (int i = 0; i < 4; ++i) {
      int v = tet_[t * 4 + i];
      int p = vert_part_id_[v];
      MapVec3 map_force(force + i * 3);
      net_force[p] -= map_force;
      MapVec3 map_x(X + v  * 3);
      Vec3 r = map_x - center_of_mass_[p];
      net_torque[p] -= r.cross(map_force);
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external force in world space to each vertex
  std::vector<Vec3> ext_force(vertex_num_, Vec3(0, 0, 0));
  const double kFloor = -0.10;
  const double kFloorStiffness = 3000;
  for (int v = 0; v < vertex_num_; ++v) {
    // gravity and internal elastic force
    ext_force[v][0] += inv_fem_->rhs_[v][0];
    ext_force[v][1] += inv_fem_->rhs_[v][1];
    ext_force[v][2] += inv_fem_->rhs_[v][2];

    // Ground collision
    //    if (0)
    if (X[v * 3 + 1] < kFloor) {
      Vec3 collision_force(0, 0, 0);
      collision_force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
      //      collision_force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
      ext_force[v][1] += collision_force[1];
      net_force[vert_part_id_[v]][1] += collision_force[1];

      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[vert_part_id_[v]] += r.cross(collision_force);
    }
  }
  {
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      ext_force[v] += force;
      int p = vert_part_id_[v];
      net_force[p] += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[p] += r.cross(force);
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // compute rigid transformation
  std::vector<Vec3> rhs(part_num_ * 2);
  std::vector<Mat3> current_inertia_tensor_ = inertia_tensor_;
  for (int p = 0; p < part_num_; ++p) {
    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_[p].transpose();
    rhs[p * 2 + 0] = net_force[p] * dt + translational_vel_[p] * mass_per_part_[p];
    rhs[p * 2 + 1] = net_torque[p] * dt + current_inertia_tensor_[p] * angular_vel_[p];
  }
  std::vector<Vec3> implicit_vel(vertex_num_);
  auto RigidMotionMatrix = [&](double * x, double * result) {
    for (int v = 0; v < vertex_num_; ++v) {
      int p = vert_part_id_[v];
      MapVec3 map_pos(X + v * 3);
      Vec3 r = map_pos - center_of_mass_[p];
      double angular_vel[3];//= angular_vel_[p].cross(r);
      dj::Cross3(x + p * 6 + 3, &r[0], &angular_vel[0]);
      implicit_vel[v][0] = x[p * 6 + 0] + angular_vel[0];
      implicit_vel[v][1] = x[p * 6 + 1] + angular_vel[1];
      implicit_vel[v][2] = x[p * 6 + 2] + angular_vel[2];
    }

    for (int p = 0; p < part_num_; ++p) {
      result[p * 6 + 0] = x[p * 6 + 0] * mass_per_part_[p];
      result[p * 6 + 1] = x[p * 6 + 1] * mass_per_part_[p];
      result[p * 6 + 2] = x[p * 6 + 2] * mass_per_part_[p];

      MapVec3 map_result(result + p * 6 + 3);
      MapVec3 map_x(x + p * 6 + 3);
      map_result = current_inertia_tensor_[p] * map_x;
    }
    //        return;
    // elastic force from boundary tets
    for (int t : interface_tet_) {
      // the force is computed from vega, the direction of the force is opposite to the actual direction
      double (*element_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        Vec3 force(0, 0, 0);
        int v = tet_[t * 4 + i];
        // compute single tet element force on vertex v
        for (int j = 0; j < 4; ++j) {
          int vj = tet_[t * 4 + j];
          force[0] += dj::Dot3(&element_k[i * 3 + 0][j * 3], &implicit_vel[vj][0]);
          force[1] += dj::Dot3(&element_k[i * 3 + 1][j * 3], &implicit_vel[vj][0]);
          force[2] += dj::Dot3(&element_k[i * 3 + 2][j * 3], &implicit_vel[vj][0]);
        }
        force *= dt_2;
        int p = vert_part_id_[v];
        result[p * 6 + 0] += force[0];
        result[p * 6 + 1] += force[1];
        result[p * 6 + 2] += force[2];

        MapVec3 map_x(X + v  * 3);
        Vec3 r = map_x - center_of_mass_[p];
        Vec3 implicit_torque = r.cross(force);
        result[p * 6 + 3] += implicit_torque[0];
        result[p * 6 + 4] += implicit_torque[1];
        result[p * 6 + 5] += implicit_torque[2];
      }
    }
  };

  MatCol Rigid = MatCol::Zero(part_num_ * 6, part_num_ * 6);
  Vec identity = Vec::Zero(part_num_ * 6);
  for (int i = 0; i < part_num_ * 6; ++i) {
    identity[i] = 1;
    RigidMotionMatrix(&identity[0], Rigid.col(i).data());
    identity[i] = 0;
  }

  //  Mat3 mm = Rigid.block<3, 3>(0, 0);
  //  PMAT(mm);
  //  L("exit");
  //  exit(0);
  MapVec map_rhs(&rhs[0][0], part_num_ * 6);
  Vec new_vel = Rigid.colPivHouseholderQr().solve(map_rhs);

  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    MapVec3 new_translational_vel(&new_vel[p * 6]);
    MapVec3 new_angular_vel(&new_vel[p * 6 + 3]);
    acceleration[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration[p] = (new_angular_vel - angular_vel_[p]) / dt;

    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // local deformation
#define LOCAL_DEFORMATION
#ifdef LOCAL_DEFORMATION
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // apply fictitious force
  //  L("no fictitious_force force");
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Vec3& acc = acceleration[p];
    Vec3& angular_acc = angular_acceleration[p];
    //    Mat3 rotation_transpose = part_rotation_[p].transpose();
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      //       inertial force: F= -m * a
      Vec3 fcititious_force = -mass_[v] * acc;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[p];
      // Euler force: F = -m * dw/dt \times r
      if (1) {
        Vec3 euler_force = -mass_[v] * angular_acc.cross(r);
        fcititious_force += euler_force;
      }
      // centerifugal force: F = -m * w \times (w \times r)
      if (1)  {
        Vec3 centerifugal_force = -mass_[v] * angular_vel_[p].cross(angular_vel_[p].cross(r));
        fcititious_force += centerifugal_force;
      }
      // Coriolis force: F = -2m w\times v
      if (1) {
        MapVec3 map_vel(velocity_ + v * 3);
        Vec3 croilis_force = -2 * mass_[v] * angular_vel_[p].cross(map_vel);
        fcititious_force += croilis_force;
      }

      ext_force[v] += fcititious_force;
      // transpose external force to local part local frame
      // ext_force[v] = rotation_transpose * ext_force[v];
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // simulate internal dynamics
  // rhs = dt * R^T * F + M * vel

  // Compute reduce stiffness matrix
  SparseMatrix* K = inv_fem_->tangent_stiffness_matrix_;
  for (int p = 0; p < part_num_; ++p) {
    tmp_v_matrix_[p] = Mat::Zero(tmp_v_matrix_[p].rows(), tmp_v_matrix_[p].cols());
  }
  for (int e = 0; e < (int)tmp_e_matrix_.size(); ++e) {
    tmp_e_matrix_[e] = Mat::Zero(tmp_e_matrix_[e].rows(), tmp_e_matrix_[e].cols());
  }
  for (int v0 = 0; v0 < vertex_num_; ++v0) {
    for (int j = 0; j < K->rowLength[v0 * 3]; j += 3) {
      int col = K->columnIndices[v0 * 3][j];
      int v1 = col / 3;
      int p0 = vert_part_id_[v0];
      int p1 = vert_part_id_[v1];
      if (p0 > p1) continue;
      int local_v0 = local_vert_idx_[v0];
      int local_v1 = local_vert_idx_[v1];
      double block[3][3] = {
        K->columnEntries[v0 * 3 + 0][j + 0], K->columnEntries[v0 * 3 + 0][j + 1], K->columnEntries[v0 * 3 + 0][j + 2],
        K->columnEntries[v0 * 3 + 1][j + 0], K->columnEntries[v0 * 3 + 1][j + 1], K->columnEntries[v0 * 3 + 1][j + 2],
        K->columnEntries[v0 * 3 + 2][j + 0], K->columnEntries[v0 * 3 + 2][j + 1], K->columnEntries[v0 * 3 + 2][j + 2],
      };
      // v_matrix
      if (p0 == p1) {
        for (int i = 0; i < basis_[p0].cols(); ++i) {
          for (int k = 0; k < 3; ++k) {
            tmp_v_matrix_[p0](local_v0 * 3 + k, i) += block[k][0] * basis_[p0](local_v1 * 3 + 0, i) +
                                                      block[k][1] * basis_[p0](local_v1 * 3 + 1, i) +
                                                      block[k][2] * basis_[p0](local_v1 * 3 + 2, i);
          }
        }
      } else {
        // e_matrix
        int e = topology_[p0][p1];
        for (int i = 0; i < basis_[p1].cols(); ++i) {
          for (int k = 0; k < 3; ++k) {
            tmp_e_matrix_[e](local_v0 * 3 + k, i) += block[k][0] * basis_[p1](local_v1 * 3 + 0, i) +
                                                     block[k][1] * basis_[p1](local_v1 * 3 + 1, i) +
                                                     block[k][2] * basis_[p1](local_v1 * 3 + 2, i);
          }
        }
      }
    }
  }

  for (int p = 0; p < part_num_; ++p) {
    MapMat result(chol_solver_->Av[p].A, basis_[p].cols(), basis_[p].cols());
    result = basis_transpose_[p] * tmp_v_matrix_[p];
    result *= dt_2;
    result += Mat::Identity(basis_[p].cols(), basis_[p].cols());
  }

  for (int e = 0; e < chol_solver_->e_number; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    MapMat result(chol_solver_->Ae[e].A, basis_[p0].cols(), basis_[p1].cols());
    result = basis_transpose_[p0] * tmp_e_matrix_[e];
    result *= dt_2;
  }

  Vec subspace_rhs = Vec::Zero(total_basis_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    double* f = &(ext_force[v][0]);
    int local_v = local_vert_idx_[v];
    for (int i = basis_offset_[p]; i < basis_offset_[p + 1]; ++i) {
      int basis_idx = i - basis_offset_[p];
      subspace_rhs[i] += basis_[p](local_v * 3 + 0, basis_idx) * f[0] +
                         basis_[p](local_v * 3 + 1, basis_idx) * f[1] +
                         basis_[p](local_v * 3 + 2, basis_idx) * f[2];
    }
  }
  subspace_rhs *= dt;
  subspace_rhs += vel_q_;

  Vec subspace_force_residual(total_basis_num_);
  chol_solver_->Multiply(&vel_q_[0], &subspace_force_residual[0]);
  subspace_force_residual = (subspace_force_residual - vel_q_) / (-1 * dt);

#if 1
  auto A = [&](double * x, double * result) {
    chol_solver_->Multiply(x, result);
  };
  cg_solver_->Resize(total_basis_num_);
  auto info = cg_solver_->Solve(&subspace_rhs[0], &vel_q_[0], A, 2000, 1e-6);
  P(info.first, info.second);
#else
  solver_->Solve(&subspace_rhs[0], &vel_q_[0]);
#endif
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#endif

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    center_of_mass_[p] += translational_vel_[p] * dt;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();
    translational_vel_[p] *= 0.98;
    angular_vel_[p] *= 0.98;

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
    quaternion_[p].Quaternion2Matrix(rotation_matrix.data());
    part_rotation_[p] = rotation_matrix;
  }

  q_ += vel_q_ * dt;
  vel_q_ *= 0.96;
  P(vel_q_.dot(vel_q_));

  // Rotate basis
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    Mat3& rot = part_rotation_[p];
    for (int i = 0; i < part_basis_size_[p]; ++i) {
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        MapVec3 u(&initial_basis_[p](local_v * 3, i));
        MapVec3 rotated_u(&basis_[p](local_v * 3, i));
        rotated_u = rot * u;
      }
    }
    basis_transpose_[p] = basis_[p].transpose();
  }
  // compute global position and update vertex offset
  ComputeGlobalPositionFromSubspace();
  for (int v = 0; v < vertex_num_; ++v) {
    inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
    inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
    inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
  }
}

void MultiDomainTet::SubspaceMultiBodySimulationOneStepWithCubature(double dt) {
  const double dt_2 = dt * dt;
  profiler.Start("simulation");
  profiler.Start("vega tets");
  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(all_cubature_tet_);
  profiler.End("vega tets");
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  profiler.Start("external force");
  GetExternalForce(dt, &ext_force_);
  profiler.End("external force");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
  // rhs = dt * R^T * F + M * vel
  profiler.Start("subspace rhs");
  Vec subspace_rhs = ComputeSubspaceRhs(dt, acceleration_, angular_acceleration_, internal_force_scaling_factor);
  profiler.End("subspace rhs");


  profiler.Start("rigid rhs");
  std::vector<Vec3> rigid_rhs = ComputeFullRigidMotionRhsWithCubature(dt, &ext_force_);
  profiler.End("rigid rhs");

  profiler.Start("rigid K");
  MatCol RigidStiffnessMatrix = ComputeFullRigidStiffnessMatrixWithCubature(dt);
  profiler.End("rigid K");

  profiler.Start("subspace K");
  ComputeSubspaceStiffnessMatrix(dt, internal_force_scaling_factor);
  profiler.End("subspace K");
  profiler.Start("full chol asm");
  full_chol_solver_->ResetTopology();
  full_chol_solver_->SetMatrixZero();
#if 0
  if (1) {
    profiler.Start("int. rigid-local");
    OMP_FOR
    for (int e = 0; e < interface_num_; ++e) {
      //    for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        const double factor = weight;// * dt_2;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (pi >= pj) continue;
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            Mat3 skew_symmetrix_mat_j = GetSkewSymmetrixMatrix(rj);
            MapMat full_mat = (pi == pj) ? MapMat(full_chol_solver_->Av[pi].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6)
                              : MapMat(full_chol_solver_->Ae[topology_[pi][pj]].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6);
            Mat tmp_mat = vert_basis_transpose_[vi] * (part_rotation_transpose_[pi] * (element_k.block<3, 3>(i * 3, j * 3) * factor));
            full_mat.block(0, part_basis_size_[pj] + 0, part_basis_size_[pi], 3) += tmp_mat;
            full_mat.block(0, part_basis_size_[pj] + 3, part_basis_size_[pi], 3) -= tmp_mat * skew_symmetrix_mat_j;

            tmp_mat = (factor * element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj]) * vert_basis_[vj];
            full_mat.block(part_basis_size_[pi] + 0, 0, 3, part_basis_size_[pj]) += tmp_mat;
            full_mat.block(part_basis_size_[pi] + 3, 0, 3, part_basis_size_[pj]) += skew_symmetrix_mat_i * tmp_mat;
          }
        }
      }
    }
    profiler.End("int. rigid-local");
    profiler.Start("domain rigid-local");
    //  if (0)
    for (int e = 0; e < interface_num_; ++e) {
      //    for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        const double factor = weight;// * dt_2;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (pi != pj) continue;
            //            if (pi == 0) {
            //              printf("%06d\t%06d\t%06d\n", vi, t, vj);
            //            }
            MapMat full_mat = (pi == pj) ? MapMat(full_chol_solver_->Av[pi].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6)
                              : MapMat(full_chol_solver_->Ae[topology_[pi][pj]].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6);

            Mat tmp_mat = (factor * element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj]) * vert_basis_[vj];
            full_mat.block(part_basis_size_[pi] + 0, 0, 3, part_basis_size_[pj]) += tmp_mat;
            full_mat.block(part_basis_size_[pi] + 3, 0, 3, part_basis_size_[pj]) += skew_symmetrix_mat_i * tmp_mat;
          }
        }
      }
    }
    profiler.End("domain rigid-local");
    //    exit(0);

  }
#else
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapMat full_mat =  MapMat(full_chol_solver_->Av[p].A, part_basis_size_[p] + 6, part_basis_size_[p] + 6);
    // all cubature vertex in this domain
    for (int j = 0; j < int(interface_v_list_[p].size()); ++j) {
      int v = interface_v_list_[p][j];
      ASSERT(vert_part_id_[v] == p);
      Mat3 sub_k = Mat3::Zero();
      // all cubature tet incident on this v
      for (int k = 0; k < int(interface_v_cubature_[p][j].size()); ++k) {
        VVCubature& cubature = interface_v_cubature_[p][j][k];
        const int& t = cubature.cubature_tet;
        const double& weight = cubature.cubature_weight;
        const int& idx = cubature.i;
        double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            sub_k(r, c) += weight * tet_k[idx * 3 + r][idx * 3 + c];
          }
        }
      }
      //      sub_k *= dt_2;
      Vec3 r = MapVec3(X + v * 3) - center_of_mass_[p];
      Mat3 skew_symmetric_mat = GetSkewSymmetrixMatrix(r);
      Mat tmp_mat = (sub_k * part_rotation_[p]) * vert_basis_[v];
      full_mat.block(part_basis_size_[p] + 0, 0, 3, part_basis_size_[p]) += tmp_mat;
      full_mat.block(part_basis_size_[p] + 3, 0, 3, part_basis_size_[p]) += skew_symmetric_mat * tmp_mat;
    }
  }
  //if (0)
  OMP_FOR
  for (int e = 0; e < part_num_ + full_chol_solver_->e_number0; ++e) {
    int p0, p1;
    double* mat = NULL;
    if (e < part_num_) {
      p0 = p1 = e;
      mat = full_chol_solver_->Av[e].A;
    } else {
      p0 = full_chol_solver_->E[(e - part_num_) * 2 + 0];
      p1 = full_chol_solver_->E[(e - part_num_) * 2 + 1];
      mat = full_chol_solver_->Ae[e - part_num_].A;
    }
    MapMat full_mat = MapMat(mat, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
    // list of vert-vert pair in this domain
    for (int i = 0; i < int(interface_vv_list_[e].size()); ++i) {
      int vi = interface_vv_list_[e][i].first;
      int vj = interface_vv_list_[e][i].second;
      // list of  cubature tets that contains this vert-vert pair
      Mat3 sub_k = Mat3::Zero();
      for (int k = 0; k < int(interface_vv_cubature_[e][i].size()); ++k) {
        VVCubature& cubature = interface_vv_cubature_[e][i][k];
        const int& t = cubature.cubature_tet;
        const double& weight = cubature.cubature_weight;
        //        if (p0 == 0 && p1 == 0) {
        //          printf("%06d\t%06d\t%06d\n", vi, t, vj);
        //        }
        const int& idx0 = cubature.i;
        const int& idx1 = cubature.j;
        double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            sub_k(r, c) += weight * tet_k[idx0 * 3 + r][idx1 * 3 + c];
          }
        }
      }
      //      sub_k *= dt_2;
      Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[vert_part_id_[vi]];
      Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[vert_part_id_[vj]];
      Mat3 skew_symmetric_mat_i = GetSkewSymmetrixMatrix(ri);
      Mat3 skew_symmetric_mat_j = GetSkewSymmetrixMatrix(rj);
      if (p0 == p1) {
        const int p = p0;
        // dvi / dvj
        Mat tmp_mat = (sub_k * part_rotation_[p]) * vert_basis_[vj];
        full_mat.block(part_basis_size_[p] + 0, 0, 3, part_basis_size_[p]) += tmp_mat;
        full_mat.block(part_basis_size_[p] + 3, 0, 3, part_basis_size_[p]) += skew_symmetric_mat_i * tmp_mat;

        //        tmp_mat = vert_basis_transpose_[vi] * (part_rotation_transpose_[p0] * sub_k);
        //        full_mat.block(0, part_basis_size_[p1] + 0, part_basis_size_[p0], 3) += tmp_mat;
        //        full_mat.block(0, part_basis_size_[p1] + 3, part_basis_size_[p0], 3) -= tmp_mat * skew_symmetric_mat_j;
        //        // dvj / dvi
        //                tmp_mat = vert_basis_transpose_[vj] * (part_rotation_transpose_[p1] * sub_k.transpose());
        //                full_mat.block(0, part_basis_size_[p1] + 0, part_basis_size_[p0], 3) += tmp_mat;
        //                full_mat.block(0, part_basis_size_[p1] + 3, part_basis_size_[p0], 3) -= tmp_mat * skew_symmetric_mat_i;

        tmp_mat = (sub_k.transpose() * part_rotation_[p]) * vert_basis_[vi];
        full_mat.block(part_basis_size_[p] + 0, 0, 3, part_basis_size_[p]) += tmp_mat;
        full_mat.block(part_basis_size_[p] + 3, 0, 3, part_basis_size_[p]) += skew_symmetric_mat_j * tmp_mat;
      } else {
        Mat tmp_mat = vert_basis_transpose_[vi] * (part_rotation_transpose_[p0] * sub_k);
        full_mat.block(0, part_basis_size_[p1] + 0, part_basis_size_[p0], 3) += tmp_mat;
        full_mat.block(0, part_basis_size_[p1] + 3, part_basis_size_[p0], 3) -= tmp_mat * skew_symmetric_mat_j;

        tmp_mat = (sub_k * part_rotation_[p1]) * vert_basis_[vj];
        full_mat.block(part_basis_size_[p0] + 0, 0, 3, part_basis_size_[p1]) += tmp_mat;
        full_mat.block(part_basis_size_[p0] + 3, 0, 3, part_basis_size_[p1]) += skew_symmetric_mat_i * tmp_mat;
      }
    }
  }
#endif
  if (0) {
    auto A = [&](double * a, double * b) {
      chol_solver_->Multiply(a, b);
    };
    dj::WriteImplicitMatrixToMatlab<double>("/tmp/full", A, total_basis_num_);
    exit(0);
  }

  profiler.Start("domain k asm");
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapMat mat(chol_solver_->Av[p].A, part_basis_size_[p], part_basis_size_[p]);
    MapMat full_mat(full_chol_solver_->Av[p].A, part_basis_size_[p] + 6, part_basis_size_[p] + 6);
    full_mat.block(part_basis_size_[p], 0, 6, part_basis_size_[p]) *= dt_2;
    full_mat.block(0, part_basis_size_[p], part_basis_size_[p], 6) = full_mat.block(part_basis_size_[p], 0, 6, part_basis_size_[p]).transpose();

    full_mat.block(0, 0, part_basis_size_[p], part_basis_size_[p]) = mat;
    full_mat.block<6, 6>(part_basis_size_[p], part_basis_size_[p]) = RigidStiffnessMatrix.block<6, 6>(p * 6, p * 6);
  }
  profiler.End("domain k asm");
  profiler.Start("rigid k asm");
  OMP_FOR
  for (int e = 0; e < full_chol_solver_->e_number0; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    MapMat mat(chol_solver_->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1]);
    MapMat full_mat(full_chol_solver_->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
    full_mat.block(0, part_basis_size_[p1], part_basis_size_[p0], 6) *= dt_2;
    full_mat.block(part_basis_size_[p0], 0, 6, part_basis_size_[p1]) *= dt_2;
    full_mat.block(0, 0, part_basis_size_[p0], part_basis_size_[p1]) = mat;
    full_mat.block<6, 6>(part_basis_size_[p0], part_basis_size_[p1]) = RigidStiffnessMatrix.block<6, 6>(p0 * 6, p1 * 6);
  }
  profiler.End("rigid k asm");

  profiler.End("full chol asm");
  if (0) {
    for (int p = 0; p < part_num_; ++p) {
      MapMat full_mat =  MapMat(full_chol_solver_->Av[p].A, part_basis_size_[p] + 6, part_basis_size_[p] + 6);
      full_mat.block(0, part_basis_size_[p], part_basis_size_[p], 6) = full_mat.block(part_basis_size_[p], 0, 6, part_basis_size_[p]).transpose();
      P(p, full_mat.norm());
    }
    for (int e = 0; e < full_chol_solver_->e_number; ++e) {
      int p0 = full_chol_solver_->E[e * 2 + 0];
      int p1 = full_chol_solver_->E[e * 2 + 1];
      MapMat full_mat = MapMat(full_chol_solver_->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
      P(e, full_mat.norm());
    }
  }
  if (0) {
    auto A = [&](double * a, double * b) {
      full_chol_solver_->Multiply(a, b);
    };
    dj::WriteImplicitMatrixToMatlab<double>("/tmp/full", A, total_basis_num_ + part_num_ * 6);
    //    exit(0);
  }
  profiler.Start("solver");
  Vec rhs_all = Vec::Zero(total_basis_num_ + part_num_ * 6);
  for (int p = 0; p < part_num_; ++p) {
    MapVec(&rhs_all[basis_offset_[p] + p * 6], part_basis_size_[p]) = MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
    MapVec(&rhs_all[basis_offset_[p] + p * 6 + part_basis_size_[p]], 6) = MapVec(&rigid_rhs[p * 2][0], 6);
  }
  Vec tmp_new_vel = rhs_all;
  Vec new_vel = rhs_all;
  if (constrained_chol_solver_ == nullptr) {
    int status = full_chol_solver_->Solve(&rhs_all[0], &tmp_new_vel[0]);
    ASSERT(status == solver::BLOCK_MATRIX_GRAPH<double>::kSuccess);
  } else {
    int dof = total_basis_num_ + (part_num_ - int(fixed_domains_.size())) * 6;
    Vec constrained_rhs(dof);
    Contract(constrained_dof_2_non_constrained_dof_, rhs_all, constrained_rhs);
    constrained_chol_solver_->ResetTopology();
    FullCholSolver2ConstrainedCholSolver(full_chol_solver_, constrained_chol_solver_);
    if (0) {
      dj::WriteVectorToMatlab(total_basis_num_ + part_num_ * 6, &rhs_all[0], "/tmp/truth");
      dj::WriteVectorToMatlab(dof, &constrained_rhs[0], "/tmp/test");
      auto A = [&](double * a, double * b) {
        constrained_chol_solver_->Multiply(a, b);
      };
      dj::WriteImplicitMatrixToMatlab<double>("/tmp/constrained", A, dof);
      exit(0);
    }
    Vec constrain_x(dof);
    int status = constrained_chol_solver_->Solve(&constrained_rhs[0], &constrain_x[0]);
    ASSERT(status == solver::BLOCK_MATRIX_GRAPH<double>::kSuccess);
    Expand(constrained_dof_2_non_constrained_dof_, constrain_x, tmp_new_vel);
  }
  profiler.End("solver");

  for (int p = 0; p < part_num_; ++p) {
    MapVec(&new_vel[basis_offset_[p]], part_basis_size_[p]) = MapVec(&tmp_new_vel[basis_offset_[p] + p * 6], part_basis_size_[p]);
    MapVec(&new_vel[total_basis_num_ + p * 6], 6) = MapVec(&tmp_new_vel[basis_offset_[p] + p * 6 + part_basis_size_[p]], 6);
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  profiler.Start("rigid update");
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapVec3 new_translational_vel(&new_vel[total_basis_num_ + p * 6 + 0]);
    MapVec3 new_angular_vel(&new_vel[total_basis_num_ + p * 6 + 3]);
    acceleration_[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration_[p] = (new_angular_vel - angular_vel_[p]) / dt;

    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
  }
  profiler.End("rigid update");

  profiler.Start("update");
  vel_q_ = MapVec(&new_vel[0], total_basis_num_);
  q_ += vel_q_ * dt;
  vel_q_ *= 0.98;
  //  P(vel_q_.dot(vel_q_));

  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    if (fixed_domain_ != p) {
      //      translational_vel_[p] += part_rotation_[p] * cm_offset[p] / dt;
    }
    center_of_mass_[p] += translational_vel_[p] * dt;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();
    translational_vel_[p] *= 0.98;
    angular_vel_[p] *= 0.98;

    quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
    part_rotation_transpose_[p] = part_rotation_[p].transpose();

    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];

    Mat3& rotation = part_rotation_[p];
    MapVec sub_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
      int v = vert_local_id2global_id_[p][local_v];
      Vec3 local_u = vert_basis_[v] * sub_q;
      MapVec3 map_x(X + v * 3);
      map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
    }
  }

  // compute global position and update vertex offset
  inv_fem_->UpdateOffset();
  UpdateObjMeshPosition();
  profiler.End("update");
  profiler.End("simulation");
  position_changed_ = true;
}


void MultiDomainTet::SubspaceMultiBodySimulationOneStep(double dt) {
  profiler.Start("total");
  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(all_cubature_tet_);
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  GetExternalForce(dt, &ext_force_);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
  // rhs = dt * R^T * F + M * vel
  std::vector<Vec3> rigid_rhs = ComputeFullRigidMotionRhs(dt, &ext_force_);
  Vec subspace_rhs = ComputeSubspaceRhs(dt, acceleration_, angular_acceleration_, internal_force_scaling_factor);
  ComputeSubspaceStiffnessMatrix(dt, internal_force_scaling_factor);
  MatCol RigidStiffnessMatrix = ComputeFullRigidStiffnessMatrix(dt);
  // export full matrix
  //  {
#if 0
  MatCol final_matrix = ComputeGlobalStiffnessMatrix(RigidStiffnessMatrix);
  profiler.Start("solver");
  Vec rhs_all = Vec::Zero(total_basis_num_ + part_num_ * 6);
  MapVec(&rhs_all[0], total_basis_num_) = subspace_rhs;
  if (fixed_domain_ != -1) {
    rigid_rhs[fixed_domain_ * 2 + 0].setZero();
    rigid_rhs[fixed_domain_ * 2 + 1].setZero();
  }
  MapVec(&rhs_all[total_basis_num_], part_num_ * 6) = MapVec(&rigid_rhs[0][0], part_num_ * 6);
  Vec new_vel = rhs_all;
  if (1) {
    //      dj::WriteEigenMatrixToMatlab(final_matrix, "/Volumes/ram/truth0");
    dj::WriteEigenMatrixToMatlab(final_matrix, "/tmp/truth0");
    exit(0);
  }
  // 1: use conjugate gradient, 0: use cholesky solver
#if 1
  auto A = [&](double * x, double * result) {
    MapVec(result, total_basis_num_ + part_num_ * 6) = final_matrix * MapVec(x, total_basis_num_ + part_num_ * 6);
    if (fixed_domain_ != -1) {
      MapVec(&result[total_basis_num_ + fixed_domain_ * 6], 6).setZero();
    }
  };
  cg_solver_->Resize(total_basis_num_ + part_num_ * 6);
  auto info = cg_solver_->Solve(&rhs_all[0], &new_vel[0], A, 780, 1e-10);
  P(info.first, info.second);
#else
  new_vel  = final_matrix.colPivHouseholderQr().solve(rhs_all);
#endif // LOCAL_DEFORMATION
  profiler.End("solver");
#else
  full_chol_solver_->ResetTopology();
  full_chol_solver_->SetMatrixZero();
  for (int p = 0; p < part_num_; ++p) {
    MapMat mat(chol_solver_->Av[p].A, part_basis_size_[p], part_basis_size_[p]);
    MapMat full_mat(full_chol_solver_->Av[p].A, part_basis_size_[p] + 6, part_basis_size_[p] + 6);
    full_mat.block(0, 0, part_basis_size_[p], part_basis_size_[p]) = mat;
    full_mat.block<6, 6>(part_basis_size_[p], part_basis_size_[p]) = RigidStiffnessMatrix.block<6, 6>(p * 6, p * 6);
  }
  for (int e = 0; e < full_chol_solver_->e_number0; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    MapMat mat(chol_solver_->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1]);
    MapMat full_mat(full_chol_solver_->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
    full_mat.block(0, 0, part_basis_size_[p0], part_basis_size_[p1]) = mat;
    full_mat.block<6, 6>(part_basis_size_[p0], part_basis_size_[p1]) = RigidStiffnessMatrix.block<6, 6>(p0 * 6, p1 * 6);
  }

  for (int e = 0; e < interface_num_; ++e) {
    for (CubaturePoint & cubature : interface_cubature_[e]) {
      int t = cubature.first;
      double weight = cubature.second;
      int* verts = tet_ + t * 4;
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = vert_part_id_[vi];
        Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
        Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int pj = vert_part_id_[vj];
          if (pi > pj) continue;
          Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
          Mat3 skew_symmetrix_mat_j = GetSkewSymmetrixMatrix(rj);
          MapMat full_mat = (pi == pj) ? MapMat(full_chol_solver_->Av[pi].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6)
                            : MapMat(full_chol_solver_->Ae[topology_[pi][pj]].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6);

          full_mat.block(0, part_basis_size_[pj] + 0, part_basis_size_[pi], 3) +=
            vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
            (element_k.block<3, 3>(i * 3, j * 3) * weight * dt * dt);

          full_mat.block(0, part_basis_size_[pj] + 3, part_basis_size_[pi], 3) -=
            vert_basis_transpose_[vi] * part_rotation_transpose_[pi] * (element_k.block<3, 3>(i * 3, j * 3) * (skew_symmetrix_mat_j * weight * dt * dt));

          full_mat.block(part_basis_size_[pi] + 0, 0, 3, part_basis_size_[pj]) +=
            element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * (weight * dt * dt);

          full_mat.block(part_basis_size_[pi] + 3, 0, 3, part_basis_size_[pj]) +=
            skew_symmetrix_mat_i * element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * (weight * dt * dt);
        }
      }
    }
  }
  profiler.Start("solver");
  Vec rhs_all = Vec::Zero(total_basis_num_ + part_num_ * 6);
  for (int p = 0; p < part_num_; ++p) {
    MapVec(&rhs_all[basis_offset_[p] + p * 6], part_basis_size_[p]) = MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
    MapVec(&rhs_all[basis_offset_[p] + p * 6 + part_basis_size_[p]], 6) = MapVec(&rigid_rhs[p * 2][0], 6);
  }
  Vec tmp_new_vel = rhs_all;
  int status = full_chol_solver_->Solve(&rhs_all[0], &tmp_new_vel[0]);
  ASSERT(status == solver::BLOCK_MATRIX_GRAPH<double>::kSuccess);
  profiler.End("solver");
  Vec new_vel = rhs_all;
  for (int p = 0; p < part_num_; ++p) {
    MapVec(&new_vel[basis_offset_[p]], part_basis_size_[p]) = MapVec(&tmp_new_vel[basis_offset_[p] + p * 6], part_basis_size_[p]);
    MapVec(&new_vel[total_basis_num_ + p * 6], 6) = MapVec(&tmp_new_vel[basis_offset_[p] + p * 6 + part_basis_size_[p]], 6);
  }
#endif
  //    dj::WriteEigenMatrixToMatlab(final_matrix, "/tmp/full");
  //    L("written");
  //    exit(0);
  //  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  profiler.Start("rigid update");
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapVec3 new_translational_vel(&new_vel[total_basis_num_ + p * 6 + 0]);
    MapVec3 new_angular_vel(&new_vel[total_basis_num_ + p * 6 + 3]);
    acceleration_[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration_[p] = (new_angular_vel - angular_vel_[p]) / dt;

    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
  }
  profiler.End("rigid update");

  profiler.Start("update");
  vel_q_ = MapVec(&new_vel[0], total_basis_num_);
  q_ += vel_q_ * dt;
  vel_q_ *= 0.98;
  //  P(vel_q_.dot(vel_q_));

  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    if (fixed_domain_ != p) {
      //      translational_vel_[p] += part_rotation_[p] * cm_offset[p] / dt;
    }
    center_of_mass_[p] += translational_vel_[p] * dt;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();
    translational_vel_[p] *= 0.98;
    angular_vel_[p] *= 0.98;

    quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
    part_rotation_transpose_[p] = part_rotation_[p].transpose();

    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];

    Mat3& rotation = part_rotation_[p];
    MapVec sub_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
      int v = vert_local_id2global_id_[p][local_v];
      Vec3 local_u = vert_basis_[v] * sub_q;
      MapVec3 map_x(X + v * 3);
      map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
    }
  }

  // compute global position and update vertex offset
  inv_fem_->UpdateOffset();
  UpdateObjMeshPosition();
  profiler.End("update");
  profiler.End("total");
  Build_VN();
}

void MultiDomainTet::SubspaceMultiBodySimulationMultiIteration(double dt, int iteration) {
  profiler.Start("total");
  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(all_cubature_tet_);
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  GetExternalForce(dt, &ext_force_);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
  // rhs = dt * R^T * F + M * vel
  std::vector<Vec3> rigid_rhs0 = ComputeFullRigidMotionRhs(dt, &ext_force_);
  Vec subspace_rhs0 = ComputeSubspaceRhs(dt, acceleration_, angular_acceleration_, internal_force_scaling_factor);
  ComputeSubspaceStiffnessMatrix(dt, internal_force_scaling_factor);
  MatCol RigidStiffnessMatrix = ComputeFullRigidStiffnessMatrix(dt);
  MatCol local_k = MatCol::Zero(total_basis_num_, total_basis_num_);
  for (int i = 0; i < total_basis_num_; ++i) {
    Vec one = Vec::Zero(total_basis_num_);
    one[i] = 1.0;
    chol_solver_->Multiply(&one[0], &local_k(0, i));
  }
  MatCol global_matrix = ComputeGlobalStiffnessMatrix(dt, RigidStiffnessMatrix);
  chol_solver_->Cholesky_Decomposition();
  auto old_translation_vel = translational_vel_;
  auto old_angular_vel = angular_vel_;
  for (int iter = 0; iter < iteration; ++iter) {
    std::vector<Vec3> rigid_rhs = rigid_rhs0;
    // elastic force from boundary tets
    for (int t : interface_tet_) {
      Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = tet_[t * 4 + i];
        int pi = vert_part_id_[vi];
        Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
        for (int j = 0; j < 4; ++j) {
          int vj = tet_[t * 4 + j];
          int pj = vert_part_id_[vj];
          Vec3 vel = vert_basis_[vj] * MapVec(&vel_q_[basis_offset_[pj]], part_basis_size_[pj]);
          vel = part_rotation_[pj] * (vel);
          Vec3 force = element_k.block<3, 3>(i * 3, j * 3) * vel * dt * dt;;
          rigid_rhs[pi * 2 + 0] -= force;
          rigid_rhs[pi * 2 + 1] -= ri.cross(force);
        }
      }
    }

    Vec new_vel = RigidStiffnessMatrix.colPivHouseholderQr().solve(MapVec(&rigid_rhs[0][0], part_num_ * 6));

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // update rigid motion, global position, and vertex offset
    OMP_FOR
    for (int p = 0; p < part_num_; ++p) {
      MapVec3 new_translational_vel(&new_vel[p * 6 + 0]);
      MapVec3 new_angular_vel(&new_vel[p * 6 + 3]);
      acceleration_[p] = (new_translational_vel - old_translation_vel[p]) / dt;
      angular_acceleration_[p] = (new_angular_vel - old_angular_vel[p]) / dt;

      translational_vel_[p] = new_translational_vel;
      angular_vel_[p] = new_angular_vel;
    }

    profiler.Start("update");
    Vec subspace_rhs = subspace_rhs0;
    //    if (0)
    for (int e = 0; e < interface_num_; ++e) {
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            MapVec(&subspace_rhs[basis_offset_[pi]], part_basis_size_[pi]) -=
              vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
              (element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] + angular_vel_[pj].cross(rj)) * weight * dt * dt);
          }
        }
      }
    }

    chol_solver_->SolveWithDecomposedMatrix(&subspace_rhs[0], &vel_q_[0]);
    // Compute error
    {
      std::vector<Vec3> rigid_rhs = rigid_rhs0;
      for (int t : interface_tet_) {
        Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = tet_[t * 4 + i];
          int pi = vert_part_id_[vi];
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          for (int j = 0; j < 4; ++j) {
            int vj = tet_[t * 4 + j];
            int pj = vert_part_id_[vj];
            Vec3 vel = vert_basis_[vj] * MapVec(&vel_q_[basis_offset_[pj]], part_basis_size_[pj]);
            vel = part_rotation_[pj] * (vel * dt * dt);
            Vec3 force = element_k.block<3, 3>(i * 3, j * 3) * vel;;
            rigid_rhs[pi * 2 + 0] -= force;
            rigid_rhs[pi * 2 + 1] -= ri.cross(force);
          }
        }
      }
      Vec error = Vec::Zero(total_basis_num_ + part_num_ * 6);
      MapVec e0 = MapVec(&error[0], part_num_ * 6);
      e0 = RigidStiffnessMatrix * new_vel - MapVec(&rigid_rhs[0][0], part_num_ * 6);

      Vec subspace_rhs = subspace_rhs0;
      for (int e = 0; e < interface_num_; ++e) {
        for (CubaturePoint & cubature : interface_cubature_[e]) {
          int t = cubature.first;
          double weight = cubature.second;
          int* verts = tet_ + t * 4;
          Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
          for (int i = 0; i < 4; ++i) {
            int vi = verts[i];
            int pi = vert_part_id_[vi];
            for (int j = 0; j < 4; ++j) {
              int vj = verts[j];
              int pj = vert_part_id_[vj];
              Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
              MapVec(&subspace_rhs[basis_offset_[pi]], part_basis_size_[pi]) -=
                vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
                (element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] + angular_vel_[pj].cross(rj)) * weight * dt * dt);
            }
          }
        }
      }
      MapVec e1 = MapVec(&error[part_num_ * 6], total_basis_num_);
      e1 = local_k * vel_q_ - subspace_rhs;
      P(iter, error.norm());//, vel_q_.norm(), new_vel.norm());
      //      P(e0.norm(), e1.norm());
      {
        Vec x = Vec::Zero(total_basis_num_ + part_num_ * 6);
        MapVec(&x[0], total_basis_num_) = vel_q_;
        MapVec(&x[total_basis_num_], part_num_ * 6) = new_vel;

        Vec rhs = x;
        MapVec(&rhs[0], total_basis_num_) = subspace_rhs0;
        MapVec(&rhs[total_basis_num_], part_num_ * 6) = MapVec(&rigid_rhs0[0][0], part_num_ * 6);
        Vec true_error = global_matrix * x - rhs;
        MapVec truth_e0(&true_error[0], total_basis_num_);
        MapVec truth_e1(&true_error[total_basis_num_], part_num_ * 6);
        P(iter, true_error.norm());//, x.norm());
        //        P(truth_e0.norm(), truth_e1.norm());
      }
      if (error.norm() < 1e-11) {
        break;
      }
    }
  }
  q_ += vel_q_ * dt;
  vel_q_ *= 0.98;
  //  P(vel_q_.dot(vel_q_));
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    center_of_mass_[p] += translational_vel_[p] * dt;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();
    translational_vel_[p] *= 0.98;
    angular_vel_[p] *= 0.98;

    quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
    part_rotation_transpose_[p] = part_rotation_[p].transpose();

    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];

    Mat3& rotation = part_rotation_[p];
    MapVec sub_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
      int v = vert_local_id2global_id_[p][local_v];
      Vec3 local_u = vert_basis_[v] * sub_q;
      MapVec3 map_x(X + v * 3);
      map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
    }
  }

  // compute global position and update vertex offset
  inv_fem_->UpdateOffset();
  UpdateObjMeshPosition();
  profiler.End("update");
  profiler.End("total");
  Build_VN();
}

void MultiDomainTet::SubspaceMultiBodySimulationWithCubature(double dt) {
  profiler.Start("total");
  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(all_cubature_tet_);
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  GetExternalForce(dt, &ext_force_);
  Vec new_rigid_vel = Vec::Zero(part_num_ * 6);
  profiler.Start("Ridid Motion");
  if (full_rigid_simulation_) {
    SimulateRigidMotionFull(dt, new_rigid_vel, &ext_force_);
  } else {
    SimulateRigidMotion(dt, new_rigid_vel, &ext_force_);
  }
  profiler.End("Ridid Motion");

  profiler.Start("rigid update");
  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapVec3 new_translational_vel(&new_rigid_vel[p * 6]);
    MapVec3 new_angular_vel(&new_rigid_vel[p * 6 + 3]);
    //        new_angular_vel .setZero();
    acceleration[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration[p] = (new_angular_vel - angular_vel_[p]) / dt;

    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
    //    translational_vel_[p] *= 0.96;
    //    angular_vel_[p] *= 0.96;
  }
  profiler.End("rigid update");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
#if 1
  // simulate internal dynamics
  // rhs = dt * R^T * F + M * vel
  Vec subspace_rhs = ComputeSubspaceRhs(dt, acceleration, angular_acceleration, internal_force_scaling_factor);
  ComputeSubspaceStiffnessMatrix(dt, internal_force_scaling_factor);
  // export full matrix
  {
    MatCol final_matrix = Mat::Zero(total_basis_num_ + part_num_ * 6, total_basis_num_ + part_num_ * 6);
    MatCol RigidStiffnessMatrix = ComputeFullRigidStiffnessMatrix(dt);
    MatCol subspace_k = MatCol::Zero(total_basis_num_, total_basis_num_);
    for (int i = 0; i < total_basis_num_; ++i) {
      Vec one = Vec::Zero(total_basis_num_);
      one[i] = 1.0f;
      chol_solver_->Multiply(&one[0], &subspace_k(0, i));
    }
    final_matrix.block(0, 0, total_basis_num_, total_basis_num_) = subspace_k;
    final_matrix.block(total_basis_num_, total_basis_num_, part_num_ * 6, part_num_ * 6) = RigidStiffnessMatrix;

    for (int e = 0; e < interface_num_; ++e) {
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            Mat3 skew_symmetrix_mat = GetSkewSymmetrixMatrix(rj);
            final_matrix.block(basis_offset_[pi], total_basis_num_ + pj * 6 + 0, part_basis_size_[pi], 3) +=
              vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
              (element_k.block<3, 3>(i * 3, j * 3) * weight * dt * dt);

            final_matrix.block(basis_offset_[pi], total_basis_num_ + pj * 6 + 3, part_basis_size_[pi], 3) -=
              vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
              (element_k.block<3, 3>(i * 3, j * 3) * (skew_symmetrix_mat * weight * dt * dt));
          }
        }
      }
    }

    for (int t : interface_tet_) {
      {
        Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = tet_[t * 4 + i];
          int pi = vert_part_id_[vi];
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          auto skew_symmetrix_mat = GetSkewSymmetrixMatrix(ri);
          for (int j = 0; j < 4; ++j) {
            int vj = tet_[t * 4 + j];
            int pj = vert_part_id_[vj];
            final_matrix.block(total_basis_num_ + pi * 6 + 0, basis_offset_[pj], 3, part_basis_size_[pj]) +=
              element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * dt * dt;

            final_matrix.block(total_basis_num_ + pi * 6 + 3, basis_offset_[pj], 3, part_basis_size_[pj]) +=
              skew_symmetrix_mat * element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * dt * dt;
          }
        }
      }
    }
    //    dj::WriteEigenMatrixToMatlab(final_matrix, "/tmp/full");
    //    L("written");
    //    exit(0);
  }
  profiler.Start("solver");
  // 1: use conjugate gradient, 0: use cholesky solver
#if 0
  //  int fixed_p = 16;
  //  MapVec(&subspace_rhs[0] + basis_offset_[fixed_p], part_basis_size_[fixed_p]).setZero();
  //  fixed_p = 5;
  //  MapVec(&subspace_rhs[0] + basis_offset_[fixed_p], part_basis_size_[fixed_p]).setZero();
  auto A = [&](double * x, double * result) {
    chol_solver_->Multiply(x, result);
    //    int fixed_p = 16;
    //    MapVec(result + basis_offset_[fixed_p], part_basis_size_[fixed_p]).setZero();
    //    fixed_p = 5;
    //    MapVec(result + basis_offset_[fixed_p], part_basis_size_[fixed_p]).setZero();
  };
  cg_solver_->Resize(total_basis_num_);
  auto info = cg_solver_->Solve(&subspace_rhs[0], &vel_q_[0], A, 3000, 1e-10);
  //  auto info = solver_->CG_Solve(&subspace_rhs[0], &vel_q_[0], 3000,k 1e-6);//, 3000, 1e-12);
  //  P(info.first, info.second);
#else
  int code = chol_solver_->Solve(&subspace_rhs[0], &vel_q_[0]);
  ASSERT(code == solver::BLOCK_MATRIX_GRAPH<double>::kSuccess);
  chol_solver_->ResetTopology();
#endif
  profiler.End("solver");
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#endif // LOCAL_DEFORMATION

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  profiler.Start("update");
  q_ += vel_q_ * dt;
  //  vel_q_ *= 0.99;
  //  P(vel_q_.dot(vel_q_));

  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    if (fixed_domain_ != p) {
      //      translational_vel_[p] += part_rotation_[p] * cm_offset[p] / dt;
    }
    center_of_mass_[p] += translational_vel_[p] * dt;
    quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
    quaternion_[p].Normalize();
    //    translational_vel_[p] *= 0.99;
    //    angular_vel_[p] *= 0.99;

    quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
    part_rotation_transpose_[p] = part_rotation_[p].transpose();

    current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];

    Mat3& rotation = part_rotation_[p];
    MapVec sub_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
      int v = vert_local_id2global_id_[p][local_v];
      Vec3 local_u = vert_basis_[v] * sub_q;
      MapVec3 map_x(X + v * 3);
      map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
    }
  }

  // compute global position and update vertex offset
  inv_fem_->UpdateOffset();
  UpdateObjMeshPosition();
  profiler.End("update");
  profiler.End("total");
  Build_VN();
}

MultiDomainTet::MatCol MultiDomainTet::ComputeGlobalStiffnessMatrix(double dt, MatCol & RigidStiffnessMatrix) {
  MatCol final_matrix = Mat::Zero(total_basis_num_ + part_num_ * 6, total_basis_num_ + part_num_ * 6);
  MatCol subspace_k = MatCol::Zero(total_basis_num_, total_basis_num_);
  for (int i = 0; i < total_basis_num_; ++i) {
    Vec one = Vec::Zero(total_basis_num_);
    one[i] = 1.0f;
    chol_solver_->Multiply(&one[0], &subspace_k(0, i));
  }
  final_matrix.block(0, 0, total_basis_num_, total_basis_num_) = subspace_k;
  final_matrix.block(total_basis_num_, total_basis_num_, part_num_ * 6, part_num_ * 6) = RigidStiffnessMatrix;

  for (int e = 0; e < interface_num_; ++e) {
    for (CubaturePoint & cubature : interface_cubature_[e]) {
      int t = cubature.first;
      double weight = cubature.second;
      int* verts = tet_ + t * 4;
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = vert_part_id_[vi];
        Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
        auto skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int pj = vert_part_id_[vj];
          Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
          Mat3 skew_symmetrix_mat_j = GetSkewSymmetrixMatrix(rj);
          final_matrix.block(basis_offset_[pi], total_basis_num_ + pj * 6 + 0, part_basis_size_[pi], 3) +=
            vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
            (element_k.block<3, 3>(i * 3, j * 3) * weight * dt * dt);

          final_matrix.block(basis_offset_[pi], total_basis_num_ + pj * 6 + 3, part_basis_size_[pi], 3) -=
            vert_basis_transpose_[vi] * part_rotation_transpose_[pi] *
            (element_k.block<3, 3>(i * 3, j * 3) * (skew_symmetrix_mat_j * weight * dt * dt));

          final_matrix.block(total_basis_num_ + pi * 6 + 0, basis_offset_[pj], 3, part_basis_size_[pj]) +=
            element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * weight * dt * dt;

          final_matrix.block(total_basis_num_ + pi * 6 + 3, basis_offset_[pj], 3, part_basis_size_[pj]) +=
            skew_symmetrix_mat_i * element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * weight * dt * dt;
        }
      }
    }
  }
  return std::move(final_matrix);
}

void MultiDomainTet::Simulate(double dt) {
  //  SubspaceMultiBodySimulation(dt); return;
  //    MultiBodyFullSimulationExplicit(dt); return;
  //    MultiBodyFullSimulation(dt); return;
  //  SimulateWithGlobalSubspace(dt); return;
  profiler.Start("asemble");
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  {
    //    L("add force at vertex", selected_vertex_);
    //      inv_fem_->rhs_[selected_vertex_][1] -= pose_conf.Get<double>("force");
    if (mouse_moved_) {
      Vec3 ui_force;
      int applied_vert = GetUIForce(&ui_force[0]);
      if (applied_vert >= 0) {
        inv_fem_->rhs_[applied_vert][0] += ui_force[0];
        inv_fem_->rhs_[applied_vert][1] += ui_force[1];
        inv_fem_->rhs_[applied_vert][2] += ui_force[2];
      }
    }
  }
  profiler.End("asemble");
  SparseMatrix* K = inv_fem_->tangent_stiffness_matrix_;
  for (int p = 0; p < part_num_; ++p) {
    tmp_v_matrix_[p] = Mat::Zero(tmp_v_matrix_[p].rows(), tmp_v_matrix_[p].cols());
  }
  for (int e = 0; e < (int)tmp_e_matrix_.size(); ++e) {
    tmp_e_matrix_[e] = Mat::Zero(tmp_e_matrix_[e].rows(), tmp_e_matrix_[e].cols());
  }

  profiler.Start("preprocess");
  for (int v0 = 0; v0 < vertex_num_; ++v0) {
    for (int j = 0; j < K->rowLength[v0 * 3]; j += 3) {
      int col = K->columnIndices[v0 * 3][j];
      int v1 = col / 3;
      // verify 3x3 block structure
      {
#if 0
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            ASSERT(j + c < K->rowLength[v0 * 3 + r] && K->columnIndices[v0 * 3 + r][j + c] == col + c);
          }
        }
#endif
      }
      int p0 = vert_part_id_[v0];
      int p1 = vert_part_id_[v1];
      if (p0 > p1) continue;
      int local_v0 = local_vert_idx_[v0];
      int local_v1 = local_vert_idx_[v1];
      double block[3][3] = {
        K->columnEntries[v0 * 3 + 0][j + 0], K->columnEntries[v0 * 3 + 0][j + 1], K->columnEntries[v0 * 3 + 0][j + 2],
        K->columnEntries[v0 * 3 + 1][j + 0], K->columnEntries[v0 * 3 + 1][j + 1], K->columnEntries[v0 * 3 + 1][j + 2],
        K->columnEntries[v0 * 3 + 2][j + 0], K->columnEntries[v0 * 3 + 2][j + 1], K->columnEntries[v0 * 3 + 2][j + 2],
      };
      // v_matrix
      if (p0 == p1) {
        for (int i = 0; i < basis_[p0].cols(); ++i) {
          for (int k = 0; k < 3; ++k) {
            tmp_v_matrix_[p0](local_v0 * 3 + k, i) += block[k][0] * basis_[p0](local_v1 * 3 + 0, i) +
                                                      block[k][1] * basis_[p0](local_v1 * 3 + 1, i) +
                                                      block[k][2] * basis_[p0](local_v1 * 3 + 2, i);
          }
        }
      } else {
        // e_matrix
        int e = topology_[p0][p1];
        for (int i = 0; i < basis_[p1].cols(); ++i) {
          for (int k = 0; k < 3; ++k) {
            tmp_e_matrix_[e](local_v0 * 3 + k, i) += block[k][0] * basis_[p1](local_v1 * 3 + 0, i) +
                                                     block[k][1] * basis_[p1](local_v1 * 3 + 1, i) +
                                                     block[k][2] * basis_[p1](local_v1 * 3 + 2, i);
          }
        }
      }
    }
  }


  const double dt_2 = dt * dt;
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    MapMat result(chol_solver_->Av[p].A, basis_[p].cols(), basis_[p].cols());
    result = basis_transpose_[p] * tmp_v_matrix_[p];
    result *= dt_2;
    result += Mat::Identity(basis_[p].cols(), basis_[p].cols());
  }


  OMP_FOR
  for (int e = 0; e < chol_solver_->e_number; ++e) {
    int p0 = chol_solver_->E[e * 2 + 0];
    int p1 = chol_solver_->E[e * 2 + 1];
    MapMat result(chol_solver_->Ae[e].A, basis_[p0].cols(), basis_[p1].cols());
    result = basis_transpose_[p0] * tmp_e_matrix_[e];
    result *= dt_2;
  }


  Vec subspace_force = Vec::Zero(total_basis_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    double* f = &(inv_fem_->rhs_[v][0]);
    int local_v = local_vert_idx_[v];
    for (int i = basis_offset_[p]; i < basis_offset_[p + 1]; ++i) {
      int basis_idx = i - basis_offset_[p];
      subspace_force[i] += basis_[p](local_v * 3 + 0, basis_idx) * f[0] +
                           basis_[p](local_v * 3 + 1, basis_idx) * f[1] +
                           basis_[p](local_v * 3 + 2, basis_idx) * f[2];
    }
  }
  subspace_force *= dt;
  subspace_force += vel_q_;
  profiler.End("preprocess");

  {
    // Verify
#if 0
    MatCol k_dot_u = MatCol::Zero(vertex_num_ * 3, total_basis_num_);
    for (int b = 0; b < total_basis_num_; ++b) {
      inv_fem_->tangent_stiffness_matrix_->MultiplyVector(global_basis_.col(b).data(), k_dot_u.col(b).data());
    }
#if 0
    //      int p = 0, local_v = 0;
    for (int p = 0; p < part_num_; ++p) {
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int global_idx = basis_offset_[p];
        //          P(global_idx);
        int global_v = vert_local_id2global_id_[p][local_v];
        //      P(basis_[p](local_v * 3 + 0, 0), global_basis_(global_v * 3 + 0, global_idx));
        //          int b = global_idx;
        for (int b = 0; b < part_basis_size_[p]; ++b) {
          ASSERT(dj::Abs(tmp_v_matrix_[p](local_v * 3 + 0, b) - k_dot_u(global_v * 3 + 0, global_idx + b)) < 1e-6, P(p, local_v, b));
          ASSERT(dj::Abs(tmp_v_matrix_[p](local_v * 3 + 1, b) - k_dot_u(global_v * 3 + 1, global_idx + b)) < 1e-6, P(p, local_v, b));
          ASSERT(dj::Abs(tmp_v_matrix_[p](local_v * 3 + 2, b) - k_dot_u(global_v * 3 + 2, global_idx + b)) < 1e-6, P(p, local_v, b));
        }
        //          P(tmp_v_matrix_[p](local_v * 3, 0), k_dot_u(global_v * 3, global_idx));
      }
    }
#endif
    MatCol reduced_k = global_basis_transpose_ * k_dot_u;
    MatCol my_k = MatCol::Zero(total_basis_num_, total_basis_num_);
    for (int i = 0; i < total_basis_num_; ++i) {
      std::vector<double> vec(total_basis_num_, 0);
      vec[i] = 1;
      solver_->Multiply(&vec[0], my_k.col(i).data());
    }
    reduced_k *= dt_2;
    reduced_k += MatCol::Identity(reduced_k.rows(), reduced_k.cols());

    // verify matrix
    if (0) {
      MatCol diff = my_k - reduced_k;
      P(diff.maxCoeff());
      P(diff.minCoeff());
      ASSERT(dj::Abs(diff.maxCoeff()) < 1e-5);
      ASSERT(dj::Abs(diff.minCoeff()) < 1e-5);
    }

    // verify rhs
    //    if (0) {
    MapVec rhs((double*)inv_fem_->rhs_, vertex_num_ * 3);
    Vec correct_reduced_force = global_basis_transpose_ * rhs;
    correct_reduced_force *= dt;
    correct_reduced_force += vel_q_;
    Vec force_diff = correct_reduced_force - subspace_force;
    double dot = force_diff.dot(force_diff);
    ASSERT(dot < 1e-6, P(dot));
    //    }
    Vec x = reduced_k.colPivHouseholderQr().solve(correct_reduced_force);
    P(x.dot(x));
#endif
  }

  profiler.Start("solver");
  //  solver_->CG_Solve(&subspace_force[0], &vel_q_[0]);//, 1000, 1e-12);
  chol_solver_->Solve(&subspace_force[0], &vel_q_[0]);
  chol_solver_->ResetTopology();
  profiler.End("solver");
  P(vel_q_.dot(vel_q_));
#if 0
  for (int p = 0; p < part_num_; ++p) {
    if (p != vert_compact_part_id_[selected_vertex_]) {
      for (int i = basis_offset_[p]; i < basis_offset_[p + 1]; ++i) {
        vel_q_[i] = 0;
      }
    }
  }
#endif
  q_ += vel_q_ * dt;
  vel_q_ *= 0.96;
  ComputeGlobalPositionFromSubspace();
  OMP_FOR
  for (int i = 0; i < vertex_num_; ++i) {
    inv_fem_->u_[i][0] = X[i * 3 + 0] - rest_pos_[i * 3 + 0];
    inv_fem_->u_[i][1] = X[i * 3 + 1] - rest_pos_[i * 3 + 1];
    inv_fem_->u_[i][2] = X[i * 3 + 2] - rest_pos_[i * 3 + 2];
  }
}

void MultiDomainTet::BuildTopology() {
  topology_ = std::vector<std::vector<int> >(part_num_, std::vector<int>(part_num_, 0));
  OMP_FOR
  for (int t = 0; t < tet_number; ++t) {
    int* v = tet_ + t * 4;
    int parts[4] = {
      vert_part_id_[v[0]],
      vert_part_id_[v[1]],
      vert_part_id_[v[2]],
      vert_part_id_[v[3]]
    };
    for (int i = 0; i < 4; ++i) {
      for (int j = i; j < 4; ++j) {
        int pi = parts[i], pj = parts[j];
        topology_[pi][pj] = 1;
        topology_[pj][pi] = 1;
      }
    }
  }
}

void MultiDomainTet::Render(int render_mode, double * pos) {
  UpdateRenderData();
  if (0) {
    glBegin(GL_LINES);
    for (int i = 0; i < int(surface_vertices_.size()); ++i) {
      int v = surface_vertices_[i];
      double* pos = X + v * 3;
      //      if (vert_tangent_[i][0] != vert_tangent_[i][0] ||
      //          vert_tangent_[i][1] != vert_tangent_[i][1] ||
      //          vert_tangent_[i][2] != vert_tangent_[i][2]) {
      //        continue;
      //      }
      //      vert_tangent_[i] *= 0.01;
      double pos1[3] = {
        pos[0] + vert_tangent_[i][0] * 0.01,
        pos[1] + vert_tangent_[i][1] * 0.01,
        pos[2] + vert_tangent_[i][2] * 0.01,
      };
      Vertex3v(pos);
      Vertex3v(&pos1[0]);
    }
    glEnd();
  }
  //  if (false && obj_) {
  if (obj_) {
    if (textured_surface_renderer_.size() != 0) {
      glEnable(GL_NORMALIZE);
      for (int g = 0; g < int(groups_.size()); ++g) {
        const ObjMesh::Group* group = obj_->getGroupHandle(g);
        auto materialHandle = obj_->getMaterialHandle(group->getMaterialIndex());
        Vec3d Ka = materialHandle->getKa();
        Vec3d Kd = materialHandle->getKd();
        float alpha = materialHandle->getAlpha();
        float ambient[4] = { (float)Ka[0], (float)Ka[1], (float)Ka[2], alpha };
        float diffuse[4] = { (float)Kd[0], (float)Kd[1], (float)Kd[2], alpha };
        float specular[4] = {1.0, 1.0, 1.0, 1.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
        textured_surface_renderer_[g].Render();
      }
      return;
    }
    if (1) {
      glEnable(GL_LIGHTING);
      glEnable(GL_NORMALIZE);
      for (int g = 0; g < int(groups_.size()); ++g) {
        //      for (int g = 0; g < 10; ++g) {
        const ObjMesh::Group* group = obj_->getGroupHandle(g);
        auto materialHandle = obj_->getMaterialHandle(group->getMaterialIndex());

        Vec3d Ka = materialHandle->getKa();
        Vec3d Kd = materialHandle->getKd();
        //        Vec3d Ks = materialHandle->getKs();

        //        float shininess = materialHandle->getShininess();
        float alpha = materialHandle->getAlpha();
        float ambient[4] = { (float)Ka[0], (float)Ka[1], (float)Ka[2], alpha };
        float diffuse[4] = { (float)Kd[0], (float)Kd[1], (float)Kd[2], alpha };
        //        float specular[4] = { (float)Ks[0], (float)Ks[1], (float)Ks[2], alpha };
        //        float diffuse[4] = {0.8, 0.0, 0.0, 1.0};
        float specular[4] = {1.0, 1.0, 1.0, 1.0};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        //                glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
        //        glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
        if (materialHandle->hasTextureFilename()) {
          ObjMeshRender::Texture* textureHandle = obj_render_->getTextureHandle(group->getMaterialIndex());
          glEnable(GL_TEXTURE_2D);
          //          glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_);
          ASSERT(materialHandle->hasBumpMapFilename(), P(textureHandle->getFullPath()));
          ObjMeshRender::Texture* bumpTextureHandle = obj_render_->getBumpTextureHandle(group->getMaterialIndex());
          glActiveTexture(GL_TEXTURE0 + bump_texture_id_);
          glBindTexture(GL_TEXTURE_2D, bumpTextureHandle->getTexture());
          glActiveTexture(GL_TEXTURE0 + texture_id_);
          glBindTexture(GL_TEXTURE_2D, textureHandle->getTexture());
          glEnable(GL_LIGHTING);
          glUseProgram(shader_);

          //          glBindBuffer(GL_ARRAY_BUFFER, tangent_vbo_);
          //          std::fill(vert_tangent_.begin(), vert_tangent_.end(), Vec3(1, 1, 0));
          //          glBufferData(GL_ARRAY_BUFFER, surface_vertices_.size() * sizeof(Vec3), &vert_tangent_[0][0], GL_DYNAMIC_DRAW);
          //          int uniloc = glGetAttribLocation(shader_, "tangent");
          //          ASSERT(uniloc >= 0, P(uniloc));
          //          glEnableVertexAttribArray(uniloc);
          //          glVertexAttribPointer(uniloc, 3, GL_DOUBLE, GL_FALSE, sizeof(Vec3), 0);
          //          glBindBuffer(GL_ARRAY_BUFFER, 0);
        } else {
          glDisable(GL_TEXTURE_2D);
        }
        //                  glDisable(GL_TEXTURE_2D);
        glBegin(GL_TRIANGLES);
        glColor3f(0, 0, 1);
        for (int t : groups_[g]) {
          int* verts = T + t * 3;
          double* tex = &tex_coord_[t * 6];
          for (int i = 0; i < 3; ++i) {
            glTexCoord2d(tex[i * 2 + 0], tex[i * 2 + 1]);
            Normal(VN + verts[i] * 3);
            //                        glColor3dv(&vert_tangent_[verts[i]][0]);
            //            Normal(TN + t * 3);
            Vertex3v(X + verts[i] * 3);
          }
        }
        glEnd();
        glUseProgram(0);
        glDisable(GL_TEXTURE_2D);
      }
      glDisable(GL_LIGHTING);
      return;
    } else {
      glEnable(GL_LIGHTING);
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
      glEnable(GL_NORMALIZE);
      obj_render_->render(OBJMESHRENDER_TRIANGLES, OBJMESHRENDER_MATERIAL | OBJMESHRENDER_TEXTURE | OBJMESHRENDER_SMOOTH);
      return;
    }
  }
  // interface tets
  if (0) {
    glDisable(GL_LIGHTING);
    glColor3fv(kOrage());
    glBegin(GL_LINES);
    for (int t : interface_tet_) {
      int* v = tet_ + t * 4;
      double* pos[4] = {
        X + v[0] * 3,
        X + v[1] * 3,
        X + v[2] * 3,
        X + v[3] * 3,
      };
      Vertex3v(pos[0]);
      Vertex3v(pos[1]);

      Vertex3v(pos[0]);
      Vertex3v(pos[2]);

      Vertex3v(pos[1]);
      Vertex3v(pos[2]);

      Vertex3v(pos[0]);
      Vertex3v(pos[3]);

      Vertex3v(pos[1]);
      Vertex3v(pos[3]);

      Vertex3v(pos[2]);
      Vertex3v(pos[3]);
    }
    glEnd();
  }

  // center of mass
  if (0) {
    glDisable(GL_LIGHTING);
    glColor3fv(kYellow());
    glPointSize(4.0);
    glBegin(GL_POINTS);
    for (int p = 0; p < part_num_; ++p) {
      Vertex3v(&center_of_mass_[p][0]);
      //    P(p, dj::Vec3d(&center_of_mass_[p][0]));
    }
    glEnd();
    for (int p = 0; p < part_num_; ++p) {
      glPushMatrix();
      glTranslated(center_of_mass_[p][0], center_of_mass_[p][1], center_of_mass_[p][2]);
      double angle;
      double axis[3];
      quaternion_[p].GetRotation(&angle, axis);
      glRotated(dj::Radian2Degree(angle), axis[0], axis[1], axis[2]);
      const float kScale = 0.05f;
      glScalef(kScale, kScale, kScale);
      DrawAxis();
      glPopMatrix();
      //    P(p, dj::Vec3d(&center_of_mass_[p][0]));
    }
    //  exit(1);
  }

  // render fixed vertex
  if (0) {
    glPointSize(6.0);
    glColor3fv(kBlue());
    glBegin(GL_POINTS);
    //    double cc[][3] = { { -0.0735189, 0.217201, 0.310038},
    //      { -0.0538111, 0.368427, 0.179364},
    //      {0.179715, 0.279858, 0.271521},
    //      {0.13692, 0.419744, 0.181553},
    //      {0.234812, 0.211822, 0.14103},
    //      {0.171691, 0.372852, 0.0850585},
    //      {0.234762, 0.148945, 0.00322043},
    //      {0.189236, 0.318458, -8.80534e-06},
    //      {0.175743, 0.140716, -0.368246},
    //      {0.101302, 0.285032, -0.230686},
    //      { -0.274315, 0.22093, 0.0980704},
    //      { -0.200936, 0.378716, 0.0658426},
    //      { -0.320562, 0.163922, 0.0136318},
    //      { -0.211874, 0.323375, -0.0145498},
    //      { -0.20524, 0.185171, -0.317437},
    //      { -0.119861, 0.295846, -0.20236},
    //      { -0.00271107, 0.493691, -0.0819799},
    //    };
    //    for (int p = 0; p < part_num_; ++p) {
    //      Vertex3v(&cc[p][0]);
    //    }
    if (0)
      for (int v = 0; v < vertex_num_; ++v) {
        if (vert_basis_[v].norm() < 1e-10) {
          Vertex3v(X + v * 3);
        }
      }

    for (int v : constrainted_vertex_) {
      Vertex3v(X + v * 3);
    }

    glEnd();
  }
  // render one vertex from each partition
  if (0) {
    std::set<int> rendered;
    QFont sansFont("Helvetica [Cronyx]", 18);
    sansFont.setFamily("sans serif");
    glColor3fv(kBlack());
    //  for (int v = 0; v < surface_; ++v)
    for (int v : surface_vertices_) {
      //      int b = closest_bone_[v];
      //      int b = vert_compact_part_id_[v];
      int b = vert_part_id_[v];
      //      int b = vert_bone_id_[v];
      if (rendered.count(b) == 0) {
        double pos[3] = {
          X[v * 3 + 0],
          X[v * 3 + 1],
          X[v * 3 + 2]
        };
        // handle retina display
#ifdef __APPLE__
        pos[0] *= 0.5; pos[1] *= 0.5; pos[2] *= 0.5;
#endif
        global::gl->renderText(
          pos[0],
          pos[1],
          pos[2],
          QString("%1").arg(b),
          sansFont
        );
        rendered.insert(b);
      }
    }

    rendered.clear();
    glBegin(GL_POINTS);
    //  for (int v = 0; v < vertex_num_; ++v) {
    for (int v : surface_vertices_) {
      int b = closest_bone_[v];
      if (rendered.count(b) == 0) {
        double* pos = X + v * 3;
        Vertex3v(pos);
        rendered.insert(b);
      }
    }
    glEnd();
  }

  const float* color_map[] = {
    kRed(),
    kGreen(),
    kBlue(),
    kYellow(),
    kOrage(),
    kChocolate(),
    kViolet(),
    kIndigo(),
  };
  const int color_num = sizeof(color_map) / sizeof(double*);
  if (render_mode & kParitionedMeshVertex) {
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    float dark[4] = {0, 0, 0, 1};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, dark);
    glPointSize(6.0);
    glBegin(GL_POINTS);
    for (int v = 0; v < vertex_num_; v++) {
      int bone = vert_part_id_[v];
      //      double color[] = {0, 0, 0};
      //      GetColor(0, skeleton_->bones_.size(), bone, color);
//      glColor3fv(color_map[bone % color_num]);

        if (dj::Abs(pos[v * 3 + 0]) > 0.95 && dj::Abs(pos[v * 3 + 2]) > 0.27) {
          glColor3fv(kRed());
        } else {
          glColor3fv(kBlack());
        }
      Vertex3v(pos + v * 3);
    }
    glEnd();
    glPopAttrib();
    glPopAttrib();
  }

  Super::Render(render_mode, pos);
  //  glPushMatrix();
  //  glTranslatef(0.5, 0, 0);
  //  Super::Render(render_mode, rest_pos_);
  //  glPopMatrix();
}

void MultiDomainTet::SetConstrainedVertex() {
  constrainted_vertex_.clear();
  for (int p = 0; p < part_num_; ++p) {
    double min_dist = 1e10;
    int min_idx = -1;
    for (int i = 0; i < vert_num_per_part_[p]; ++i) {
      int v = vert_local_id2global_id_[p][i];
      double dist = dj::Distance3(rest_pos_ + v * 3, &center_of_mass_[p][0]);
      if (dist < min_dist) {
        min_dist = dist;
        min_idx = v;
      }
    }
    ASSERT(min_idx >= 0);
    constrainted_vertex_.push_back(min_idx);
  }
}


template <class T>
void Clear(std::vector<T>& target) {
  for (int i = 0; i < int(target.size()); ++i) {
    target[i].clear();
  }
}

void MultiDomainTet::BuildParallelAccelerationStructure() {
  vv_cubature_.resize(part_num_ + chol_solver_->e_number0);
  vv_list_.resize(part_num_ + chol_solver_->e_number0);
  v_cubature_.resize(part_num_);
  v_list_.resize(part_num_);

  interface_vv_cubature_.resize(part_num_ + chol_solver_->e_number0);
  interface_vv_list_.resize(part_num_ + chol_solver_->e_number0);
  interface_v_cubature_.resize(part_num_);
  interface_v_list_.resize(part_num_);

  Clear(vv_cubature_);
  Clear(vv_list_);
  Clear(v_cubature_);
  Clear(v_list_);

  Clear(interface_vv_cubature_);
  Clear(interface_vv_list_);
  Clear(interface_v_cubature_);
  Clear(interface_v_list_);
  // accerlation structure for local deformation
  {
    std::map<std::pair<int, int>, int> vv_map;
    std::map<int, int> v_map;
    for (CubaturePoint & cubature : all_cubature_) {
      int t = cubature.first;
      double weight = cubature.second;
      for (int i = 0; i < 4; ++i) {
        int vi = tet_[t * 4 + i];
        int pi = vert_part_id_[vi];

        if (v_map.count(vi) > 0) {
          int idx = v_map[vi];
          ASSERT(idx < int(v_list_[pi].size()));
          ASSERT(idx < int(v_cubature_[pi].size()));
          ASSERT(v_list_[pi][idx] == vi);
          //          P(pi, v_cubature_[pi].size(), idx, v_list_[pi].size());
          v_cubature_[pi][idx].emplace_back(t, weight, i);
        } else {
          v_map[vi] = int(v_list_[pi].size());
          v_list_[pi].push_back(vi);
          v_cubature_[pi].push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i)));
        }

        for (int j = 0; j < 4; ++j) {
          int vj = tet_[t * 4 + j];
          int pj = vert_part_id_[vj];
          if (j == i) continue;
          if (pi > pj) continue;
          if (pi == pj && vi >= vj) continue;
          int vv_idx = (pi == pj) ? pi : topology_[pi][pj] + part_num_;
          auto v_pair = std::make_pair(vi, vj);
          if (vv_map.count(v_pair) > 0) {
            int idx = vv_map[v_pair];
            ASSERT(vv_list_[vv_idx][idx] == v_pair);
            vv_cubature_[vv_idx][idx].emplace_back(t, weight, i, j);
          } else {
            vv_map[v_pair] = int(vv_list_[vv_idx].size());
            vv_list_[vv_idx].emplace_back(vi, vj);
            vv_cubature_[vv_idx].push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i, j)));
          }
        }  // each col
      } // each row
    } // all cubature
  }
  // acceleration structure for local-rigid coupling
#if 1
  {
    std::map<std::pair<int, int>, int> vv_map;
    std::map<int, int> v_map;
    for (int e = 0; e < interface_num_; ++e) {
      for (CubaturePoint & cubature : interface_rigid_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        for (int i = 0; i < 4; ++i) {
          int vi = tet_[t * 4 + i];
          int pi = vert_part_id_[vi];

          if (v_map.count(vi) > 0) {
            int idx = v_map[vi];
            ASSERT(interface_v_list_[pi][idx] == vi);
            interface_v_cubature_[pi][idx].emplace_back(t, weight, i);
          } else {
            v_map[vi] = int(interface_v_list_[pi].size());
            interface_v_list_[pi].push_back(vi);
            interface_v_cubature_[pi].push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i)));
          }

          for (int j = 0; j < 4; ++j) {
            int vj = tet_[t * 4 + j];
            int pj = vert_part_id_[vj];
            if (j == i) continue;
            if (pi > pj) continue;
            if (pi == pj && vi >= vj) continue;
            int vv_idx = (pi == pj) ? pi : topology_[pi][pj] + part_num_;
            auto v_pair = std::make_pair(vi, vj);
            if (vv_map.count(v_pair) > 0) {
              int idx = vv_map[v_pair];
              ASSERT(interface_vv_list_[vv_idx][idx] == v_pair);
              interface_vv_cubature_[vv_idx][idx].emplace_back(t, weight, i, j);
            } else {
              vv_map[v_pair] = int(interface_vv_list_[vv_idx].size());
              interface_vv_list_[vv_idx].emplace_back(vi, vj);
              interface_vv_cubature_[vv_idx].push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i, j)));
            }
          }  // each col
        } // each row
      } // all cubature
    }
  }
#endif

  tmp_interface_rigid_k_.resize(chol_solver_->e_number0 * 8);

  domain_incident_interface_.resize(part_num_);
  Clear(domain_incident_interface_);

  for (int e = 0; e < interface_num_; ++e) {
    if (interface_domains_[e][0] > 2) continue;
    int p0 = interface_domains_[e][1];
    int p1 = interface_domains_[e][2];
    domain_incident_interface_[p0].emplace_back(e, 0);
    domain_incident_interface_[p1].emplace_back(e, 1);
  }
}

void MultiDomainTet::set_fixed_domain(int fixed_domain) {
  this->fixed_domain_ = fixed_domain;
}

void MultiDomainTet::UpdateTexturedSurfaceMeshVBO() {
  if (int(textured_surface_renderer_.size()) == 0) return;
  OMP_FOR
  for (int g = 0; g < int(groups_.size()); ++g) {
    int offset = group_triangle_offset_[g] * 9;
    for (int i = 0; i < int(groups_[g].size()); ++i) {
      int t = groups_[g][i];
      int* verts = T + t * 3;
      for (int i = 0; i < 3; ++i) {
        surface_vert_pos_[offset + i * 3 + 0] = X[verts[i] * 3 + 0];
        surface_vert_pos_[offset + i * 3 + 1] = X[verts[i] * 3 + 1];
        surface_vert_pos_[offset + i * 3 + 2] = X[verts[i] * 3 + 2];
        //        surface_vert_normal_[offset + i * 3 + 0] = TN[t * 3 + 0];
        //        surface_vert_normal_[offset + i * 3 + 1] = TN[t * 3 + 1];
        //        surface_vert_normal_[offset + i * 3 + 2] = TN[t * 3 + 2];
        surface_vert_normal_[offset + i * 3 + 0] = VN[verts[i] * 3 + 0];
        surface_vert_normal_[offset + i * 3 + 1] = VN[verts[i] * 3 + 1];
        surface_vert_normal_[offset + i * 3 + 2] = VN[verts[i] * 3 + 2];
      }
      offset += 9;
    }
  }
  for (int g = 0; g < int(groups_.size()); ++g) {
    int offset = group_triangle_offset_[g] * 9;
    textured_surface_renderer_[g].UpdatePosAndNormal(&surface_vert_pos_[offset], &surface_vert_normal_[offset], int(groups_[g].size()));
  }
}

void MultiDomainTet::EnalbeTexturedSurfaceRendering(const char *shader_source) {
  surface_vert_pos_.resize(triangle_num_ * 3 * 3);
  surface_vert_normal_.resize(triangle_num_ * 3 * 3);
  surface_vert_texture_coord_.resize(triangle_num_ * 3 * 2);
  int shader = SetupGLSL(shader_source);
  for (int i = 0; i < int(groups_.size()); ++i) {
    textured_surface_renderer_.emplace_back(shader);
  }

  auto tmp_tex_coord = tex_coord_;
  //  OMP_FOR
  for (int g = 0; g < int(groups_.size()); ++g) {
    int offset = group_triangle_offset_[g] * 6;
    double* group_tex_coord = &tex_coord_[offset];
    for (int i = 0; i < int(groups_[g].size()); ++i) {
      int t = groups_[g][i];
      for (int i = 0; i < 3; ++i) {
        tex_coord_[offset + i * 2 + 0] = tmp_tex_coord[t * 6 + i * 2 + 0];
        tex_coord_[offset + i * 2 + 1] = tmp_tex_coord[t * 6 + i * 2 + 1];
      }
      offset += 6;
    }

    const ObjMesh::Group* group = obj_->getGroupHandle(g);
    auto materialHandle = obj_->getMaterialHandle(group->getMaterialIndex());
    //    ASSERT(glGetError() == GL_NO_ERROR);
    if (materialHandle->hasTextureFilename()) {
      ObjMeshRender::Texture* textureHandle = obj_render_->getTextureHandle(group->getMaterialIndex());
      ObjMeshRender::Texture* bumpTextureHandle = obj_render_->getBumpTextureHandle(group->getMaterialIndex());
      std::vector<std::string> names({"texture", "bumpMap"});
      std::vector<unsigned int> texture_ids({textureHandle->getTexture(), bumpTextureHandle->getTexture()});
      textured_surface_renderer_[g].SetTextures(names, texture_ids);
      textured_surface_renderer_[g].UpdateTexCoord(group_tex_coord, int(groups_[g].size()));
    }
  }
  this->UpdateTexturedSurfaceMeshVBO();
}

void MultiDomainTet::SetFixedDomains(std::set<int> fixed_domains) {
  fixed_domains_ = std::vector<int>(fixed_domains.begin(), fixed_domains.end());
  std::sort(fixed_domains_.begin(), fixed_domains_.end());
  is_fixed_domain_.resize(part_num_);
  std::fill(is_fixed_domain_.begin(), is_fixed_domain_.end(), false);
  for (int fixed_p : fixed_domains) {
    ASSERT(fixed_p < part_num_ && fixed_p >= 0);
    is_fixed_domain_[fixed_p] = true;
  }

  non_constrained_dof_2_constrained_dof_.resize(total_basis_num_ + part_num_ * 6);
  constrained_dof_2_non_constrained_dof_.resize(total_basis_num_ + (part_num_ - fixed_domains_.size()) * 6);
  for (int p = 0, offset = 0; p < part_num_; ++p) {
    for (int i = 0; i < part_basis_size_[p]; ++i) {
      non_constrained_dof_2_constrained_dof_[basis_offset_[p] + p * 6 + i] = offset;
      constrained_dof_2_non_constrained_dof_[offset] = basis_offset_[p] + p * 6 + i;
      ++offset;
    }
    if (is_fixed_domain_[p]) {
      for (int i = 0; i < 6; ++i) {
        non_constrained_dof_2_constrained_dof_[basis_offset_[p + 1] + p * 6 + i] = INT_MIN;
      }
    } else {
      for (int i = 0; i < 6; ++i) {
        non_constrained_dof_2_constrained_dof_[basis_offset_[p + 1] + p * 6 + i] = offset;
        constrained_dof_2_non_constrained_dof_[offset] = basis_offset_[p + 1] + p * 6 + i;
        offset++;
      }
    }
  }

  {
    std::vector<int> edges(chol_solver_->E, chol_solver_->E + chol_solver_->e_number0 * 2);
    std::vector<int> block_size(part_num_);
    for (int p = 0; p < part_num_; ++p) {
      block_size[p] = is_fixed_domain_[p] ? part_basis_size_[p] : part_basis_size_[p] + 6;
    }
    constrained_chol_solver_ = new solver::BLOCK_MATRIX_GRAPH<double>(part_num_, edges, block_size);
  }
}

inline void MultiDomainTet::FullCholSolver2ConstrainedCholSolver(solver::BLOCK_MATRIX_GRAPH<double> *full_chol_solver,
                                                                 solver::BLOCK_MATRIX_GRAPH<double> *constrained_chol_solver) {
  OMP_FOR
  for (int p = 0; p < part_num_; ++p) {
    if (is_fixed_domain_[p]) {
      MapMat non_constrained_mat(full_chol_solver->Av[p].A, part_basis_size_[p] + 6, part_basis_size_[p] + 6);
      MapMat constrained_mat(constrained_chol_solver->Av[p].A, part_basis_size_[p], part_basis_size_[p]);
      constrained_mat = non_constrained_mat.block(0, 0, part_basis_size_[p], part_basis_size_[p]);
    } else {
      memcpy(constrained_chol_solver->Av[p].A, full_chol_solver->Av[p].A, sizeof(double) * full_chol_solver->Av[p].size);
    }
  }

  OMP_FOR
  for (int e = 0; e < full_chol_solver->e_number0; ++e) {
    int p0 = full_chol_solver->E[e * 2 + 0];
    int p1 = full_chol_solver->E[e * 2 + 1];
    if (!is_fixed_domain_[p0] && !is_fixed_domain_[p1]) {
      memcpy(constrained_chol_solver->Ae[e].A, full_chol_solver->Ae[e].A, sizeof(double) * full_chol_solver->Ae[e].size);
    } else if (is_fixed_domain_[p0] && !is_fixed_domain_[p1]) {
      MapMat non_constrained_mat(full_chol_solver->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
      MapMat constrained_mat(constrained_chol_solver->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1] + 6);
      constrained_mat = non_constrained_mat.block(0, 0, part_basis_size_[p0], part_basis_size_[p1] + 6);
    } else if (!is_fixed_domain_[p0] && is_fixed_domain_[p1]) {
      MapMat non_constrained_mat(full_chol_solver->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
      MapMat constrained_mat(constrained_chol_solver->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1]);
      constrained_mat = non_constrained_mat.block(0, 0, part_basis_size_[p0] + 6, part_basis_size_[p1]);
    } else {
      MapMat non_constrained_mat(full_chol_solver->Ae[e].A, part_basis_size_[p0] + 6, part_basis_size_[p1] + 6);
      MapMat constrained_mat(constrained_chol_solver->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1]);
      constrained_mat = non_constrained_mat.block(0, 0, part_basis_size_[p0], part_basis_size_[p1]);
    }
  }
}

void MultiDomainTet::ComputeTangentVector() {
  typedef Eigen::Vector2d Vec2;
  typedef Eigen::Map<Vec2> MapVec2;
  OMP_FOR
  for (int t = 0; t < triangle_num_; ++t) {
    int* vtx = T + t * 3;
    MapVec3 p0(X + vtx[0] * 3);
    MapVec3 p1(X + vtx[1] * 3);
    MapVec3 p2(X + vtx[2] * 3);
    MapVec2 p_tex0(&tex_coord_[t * 6 + 0]);
    MapVec2 p_tex1(&tex_coord_[t * 6 + 2]);
    MapVec2 p_tex2(&tex_coord_[t * 6 + 4]);
    Vec3 v1 = p1 - p0;
    Vec3 v2 = p2 - p0;
    Vec2 tex1 = p_tex1 - p_tex0;
    Vec2 tex2 = p_tex2 - p_tex0;
    double c =  (tex1[0] * tex2[1] - tex1[1] * tex2[0]) ;
    double inv = 1 / c ;
    if (inv != inv || !std::isfinite(inv)) {
      tri_tangent_[t][0] = 0;
      tri_tangent_[t][0] = 0;
      tri_tangent_[t][0] = 0;
    } else {
      tri_tangent_[t][0] = inv * (v1[0] * tex2[1] - v2[0] * tex1[1]) ;
      tri_tangent_[t][1] = inv * (v1[1] * tex2[1] - v2[1] * tex1[1]) ;
      tri_tangent_[t][2] = inv * (v1[2] * tex2[1] - v2[2] * tex1[1]) ;
      //      dj::Normalize3(&tri_tangent_[t][0]);
    }
    //    P(inv);
    //    PVEC(v1);
    //    PVEC(v2);
    //    PVEC(tex1);
    //    PVEC(tex2);
  }
  memset(&vert_tangent_[0][0], 0, sizeof(Vec3) * surface_vertices_.size());
  OMP_FOR
  for (int i = 0; i < int(surface_vertices_.size()); ++i) {
    int v = surface_vertices_[i];
    for (int t : incident_triangle_[v]) {
      vert_tangent_[i] += tri_tangent_[t];
    }
    if (vert_tangent_[i][0] !=  vert_tangent_[i][0] ||
        !std::isfinite(vert_tangent_[i][0]) ||
        !std::isfinite(vert_tangent_[i][1]) ||
        !std::isfinite(vert_tangent_[i][2]) ||
        vert_tangent_[i][1] !=  vert_tangent_[i][1] ||
        vert_tangent_[i][2] !=  vert_tangent_[i][2]) {
      vert_tangent_[i][0] = 0;
      vert_tangent_[i][1] = 0;
      vert_tangent_[i][2] = 1;
    }
    dj::Normalize3(&vert_tangent_[i][0]);
  }
}


void MultiDomainTet::EmbededObjMesh(ObjMesh * obj, ObjMeshRender * renderer, bool import_mapping, const char * folder) {
  obj_ = obj;
  obj_render_ = renderer;
  if (import_mapping) {
    ASSERT(folder);
    std::string vert_obj2tet_file = dj::Format("%z/vert_obj2tet.txt", folder);
    dj::Read1DVectorText<int>(vert_obj2tet_file.c_str(), vert_obj2tet_, int(surface_vertices_.size()));

    std::string vert_tet2obj_file = dj::Format("%z/vert_tet2obj.txt", folder);
    dj::Read1DVectorText<int>(vert_tet2obj_file.c_str(), vert_tet2obj_, int(surface_vertices_.size()));

    std::string tri_mapping_file = dj::Format("%z/tri_obj2tet.txt", folder);
    dj::Read1DVectorText<int>(tri_mapping_file.c_str(), tri_obj2tet_, triangle_num_);
  } else {
    BuildVertTriMapping(folder);
  }
  BuildRenderGroups();
  if (0) {
    dj::Vec3i verts(1614, 11554, 11571);
    for (int t = 0; t < triangle_num_; ++t) {
      dj::Vec3i v(T[t * 3 + 0], T[t * 3 + 1], T[t * 3 + 2]);
      std::sort(v.begin(), v.end());
      if (verts[0] == v[0] && verts[1] == v[1] && verts[2] == v[2]) {
        P(t);
        for (int i = 0; i < 3; ++i) {
          P(i, t, T[t * 3 + i], dj::Vec2d(&tex_coord_[0] + t * 6 + i * 2));
        }
        for (int t : incident_triangle_[1614]) {
          for (int i = 0; i < 3; ++i) {
            P(i, t, T[t * 3 + i], dj::Vec2d(&tex_coord_[0] + t * 6 + i * 2));
          }
        }
        exit(0);
      }
    }
    ASSERT(false);
  }
  {
    for (int t : {3799, 3268}) {
      for (int i = 0; i < 3; ++i) {
        if (T[t * 3 + i] == 1614) {
          tex_coord_[t * 6 + i * 2 + 0] = 0.189474;
          tex_coord_[t * 6 + i * 2 + 1] = 0.93;
        }
      }
    }
  }
  texture_id_ = 0;
  bump_texture_id_ = 1;
  shader_ = SetupGLSL(DATA_DIRECTORY "../shader/octopus");
  glUseProgram(shader_);
  {
    int uniloc = -1;
    //    uniloc = glGetUniformLocation(shader_, "test");
    //    ASSERT(uniloc >= 0, P(uniloc, glGetError(), shader_));
    //    glUniform1f(uniloc, 1.0f);

    uniloc = glGetUniformLocation(shader_, "bumpMap");
    //    ASSERT(uniloc >= 0, P(uniloc));
    glUniform1i(uniloc, bump_texture_id_);

    uniloc = glGetUniformLocation(shader_, "texture");
    //    ASSERT(uniloc >= 0, P(uniloc));
    glUniform1i(uniloc, texture_id_);
  }
  glGenBuffers(1, &tangent_vbo_);
  ASSERT(tangent_vbo_ > 0)
  vert_tangent_.resize(surface_vertices_.size());
  tri_tangent_.resize(triangle_num_);
  ComputeTangentVector();
  glUseProgram(0);
}


void MultiDomainTet::BuildRenderGroups() {
  groups_.resize(obj_->getNumGroups());
  int tri_count = 0;
  tex_coord_.resize(triangle_num_ * 2 * 3);
  for (int g = 0; g < int(obj_->getNumGroups()); ++g) {
    const ObjMesh::Group* group = obj_->getGroupHandle(g);
    for (int t = 0; t < int(group->getNumFaces()); ++t, ++tri_count) {
      const ObjMesh::Face* face = group->getFaceHandle(t);
      int tri_idx = tri_obj2tet_[tri_count];
      ASSERT(tri_idx >= 0 && tri_idx < triangle_num_);
      groups_[g].push_back(tri_idx);
      const ObjMesh::Vertex* verts[3] = {
        face->getVertexHandle(0),
        face->getVertexHandle(1),
        face->getVertexHandle(2),
      };
      int idx[3] = {
        vert_obj2tet_[verts[0]->getPositionIndex()],
        vert_obj2tet_[verts[1]->getPositionIndex()],
        vert_obj2tet_[verts[2]->getPositionIndex()],
      };
      Vec3d tex_uv[3] = {
        obj_->getTextureCoordinate(*verts[0]),
        obj_->getTextureCoordinate(*verts[1]),
        obj_->getTextureCoordinate(*verts[2]),
      };
      for (int i = 0; i < 3; ++i) {
        int v = idx[i];
        int j = 0;
        for (; j < 3; ++j) {
          if (T[tri_idx * 3 + j] == v) {
            break;
          }
        }
        ASSERT(j < 3);
        tex_coord_[tri_idx * 6 + j * 2 + 0] = tex_uv[i][0];
        tex_coord_[tri_idx * 6 + j * 2 + 1] = tex_uv[i][1];
      }
    }
  }

  group_triangle_offset_.resize(groups_.size() + 1);
  group_triangle_offset_[0] = 0;
  for (int i = 1; i < group_triangle_offset_.size(); ++i) {
    group_triangle_offset_[i] = group_triangle_offset_[i - 1] + int(groups_[i - 1].size());
  }

}

void MultiDomainTet::BuildVertTriMapping(const char * export_folder) {
  ASSERT(obj_->getNumVertices() == surface_vertices_.size());
  vert_tet2obj_ = std::vector<int>(obj_->getNumVertices(), -1);
  vert_obj2tet_ = std::vector<int>(obj_->getNumVertices(), -1);
  for (int v = 0; v < int(obj_->getNumVertices()); ++v) {
    Vec3d pos = obj_->getPosition(v);
    double min_dist = 1e10;
    int min_idx = 0;
    for (int i = 0; i < int(obj_->getNumVertices()); ++i) {
      double dist = dj::Distance3(&pos[0], X + i * 3);
      if (dist < min_dist) {
        min_idx = i;
        min_dist = dist;
      }
    }
    vert_tet2obj_[min_idx] = v;
    vert_obj2tet_[v] = min_idx;
    ASSERT(min_dist < 1e-6, P(v, min_dist));
  }
  std::unordered_set<int> set0(vert_tet2obj_.begin(), vert_tet2obj_.end());
  std::unordered_set<int> set1(vert_obj2tet_.begin(), vert_obj2tet_.end());
  ASSERT(set0.size() == obj_->getNumVertices());
  ASSERT(set1.size() == obj_->getNumVertices());
  for (int v = 0; v < int(obj_->getNumVertices()); ++v) {
    ASSERT(vert_tet2obj_[v] >= 0 && vert_tet2obj_[v] < int(obj_->getNumVertices()));
    ASSERT(vert_obj2tet_[v] >= 0 && vert_obj2tet_[v] < int(obj_->getNumVertices()));
  }
  ASSERT(triangle_num_ == int(obj_->getNumFaces()));
  tri_obj2tet_.resize(obj_->getNumFaces());
  std::map<dj::Vec3i, int> vert2tri;
  for (int t = 0; t < triangle_num_; ++t) {
    dj::Vec3i verts(T + t * 3);
    std::sort(verts.begin(), verts.end());
    vert2tri[verts] = t;
  }
  obj_->initTriangleLookup();
  for (int t = 0; t < int(obj_->getNumFaces()); ++t) {
  }

  int tri_count = 0;
  for (int g = 0; g < int(obj_->getNumGroups()); ++g) {
    const ObjMesh::Group* group = obj_->getGroupHandle(g);
    for (int t = 0; t < int(group->getNumFaces()); ++t, ++tri_count) {
      const ObjMesh::Face* face = group->getFaceHandle(t);
      dj::Vec3i tri_verts;
      tri_verts[0] = vert_obj2tet_[face->getVertexHandle(0)->getPositionIndex()];
      tri_verts[1] = vert_obj2tet_[face->getVertexHandle(1)->getPositionIndex()];
      tri_verts[2] = vert_obj2tet_[face->getVertexHandle(2)->getPositionIndex()];
      std::sort(tri_verts.begin(), tri_verts.end());
      auto iter = vert2tri.find(tri_verts);
      ASSERT(iter != vert2tri.end());
      tri_obj2tet_[tri_count] = iter->second;
    }
  }
  ASSERT(tri_count == int(obj_->getNumFaces()), P(tri_count, obj_->getNumFaces()));
  ASSERT(tri_count == triangle_num_, P(tri_count, triangle_num_));
  if (export_folder) {
    std::string vert_obj2tet_file = dj::Format(" % z / vert_obj2tet.txt", export_folder);
    dj::Write1DVectorText<int>(vert_obj2tet_file.c_str(), vert_obj2tet_, int(surface_vertices_.size()));

    std::string vert_tet2obj_file = dj::Format(" % z / vert_tet2obj.txt", export_folder);
    dj::Write1DVectorText<int>(vert_tet2obj_file.c_str(), vert_tet2obj_, int(surface_vertices_.size()));

    std::string tri_mapping_file = dj::Format(" % z / tri_obj2tet.txt", export_folder);
    dj::Write1DVectorText<int>(tri_mapping_file.c_str(), tri_obj2tet_, triangle_num_);
  }
}

void MultiDomainTet::UpdateObjMeshPosition() {
  if (!obj_) return;
  OMP_FOR
  for (int i = 0; i < int(surface_vertices_.size()); ++i) {
    int v = surface_vertices_[i];
    int obj_v = vert_tet2obj_[v];
    Vec3d pos(X[v * 3 + 0], X[v * 3 + 1], X[v * 3 + 2]);
    obj_->setPosition(obj_v, pos);
  }
}

#if  0
// subspace momentum compensation
std::vector<Vec3> cm_offset(part_num_, Vec3(0, 0, 0));
std::vector<Mat3> rotation(part_num_, Mat3::Zero());
if (0)
  for (int p = 0; p < part_num_; ++p) {
    if (p == fixed_domain_) continue;
    MapVec part_q(&q_[basis_offset_[p]], part_basis_size_[p]);
    Vec old_q = part_q;
    cm_offset[p] = momentum_matrix_[p] * old_q;
    //    P(cm_offset[p]);
    Mat3 lagragian_matrix = momentum_matrix_[p] * momentum_matrix_transpose_[p];
    Vec3 lambda = lagragian_matrix.ldlt().solve(cm_offset[p]);
    Vec offset = momentum_matrix_transpose_[p] * lambda;
    //    part_q = old_q - offset;
    if (0) {
      double norm = (momentum_matrix_[p].transpose() - momentum_matrix_transpose_[p]).norm();
      auto diff = momentum_matrix_[p] * part_q;
      auto diff0 = momentum_matrix_[p] * old_q;
      PMAT(lagragian_matrix);
      P(lambda);
      ASSERT(diff.norm() < 1e-8, P(diff.norm(), diff0.norm(), norm));
    }
    //    MapVec(&vel_q_[basis_offset_[p]], part_basis_size_[p]) -= offset / dt;
    //    cm_offset[p] /= mass_per_part_[p];
    //    P(cm_offset[p]);
  }
#endif
