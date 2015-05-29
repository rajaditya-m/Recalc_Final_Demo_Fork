#include <limits>
#include "mixed_sparse_matrix.h"
#include "mixed_multi_domain_tet.h"
#include "../basis_generation/pardio_solver/PardisoSolver.h"
#include "BLOCK_MATRIX_GRAPH.h"
#include "tet_mesh_simulator_bridge.h"
#include "global.h"
#include "profiler.h"
#include "matlab_io.h"

MixedSparseMatrix::MixedSparseMatrix(MixedMultiDomainTet *multi_tet) {
  this->multi_tet_ = multi_tet;
  //  profiler.Start("outline");
  if (multi_tet_->fixed_domains_.size() == 0) {
    SparseMatrixOutline outline(multi_tet_->total_basis_num_ + multi_tet_->subspace_domain_num_ * 6);
    for (int i = 0; i < multi_tet_->subspace_domain_num_; ++i) {
      int p = multi_tet_->subspace_domains_[i];
      int offset = multi_tet_->domain_offset_[p];
      int size = multi_tet_->part_basis_size_[p] + 6;
      outline.AddBlockMatrix(offset, offset, size, size);
    }

    for (int e = 0; e < multi_tet_->full_chol_solver_->e_number0; ++e) {
      int* E = multi_tet_->full_chol_solver_->E;
      int p0 = multi_tet_->subspace_domains_[E[e * 2 + 0]];
      int p1 = multi_tet_->subspace_domains_[E[e * 2 + 1]];
      int offset0 = multi_tet_->domain_offset_[p0];
      int offset1 = multi_tet_->domain_offset_[p1];
      int size0 = multi_tet_->part_basis_size_[p0] + 6;
      int size1 = multi_tet_->part_basis_size_[p1] + 6;
      outline.AddBlockMatrix(offset0, offset1, size0, size1);
      outline.AddBlockMatrix(offset1, offset0, size1, size0);
    }

    const int kFullOffset = multi_tet_->subspace_domain_basis_size_ + multi_tet_->subspace_domain_num_ * 6;
    for (int t : multi_tet_->full_tet_) {
      int* verts = multi_tet_->tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int local_vi = multi_tet_->full_vert_idx_[vi];
        ASSERT(local_vi != -1);
        int offset_i = kFullOffset + local_vi * 3;
        outline.AddBlockMatrix(offset_i, offset_i, 3, 3);
        for (int j = i + 1; j < 4; ++j) {
          int vj = verts[j];
          int local_vj = multi_tet_->full_vert_idx_[vj];
          int offset_j = kFullOffset + local_vj * 3;
          outline.AddBlockMatrix(offset_i, offset_j, 3, 3);
          outline.AddBlockMatrix(offset_j, offset_i, 3, 3);
        }
      }
    }

    for (int t : multi_tet_->full_subspace_interface_tet_) {
      int* verts = multi_tet_->tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = multi_tet_->vert_part_id_[vi];
        int offset_i, size_i;
        if (multi_tet_->is_subspace_domain_[pi]) {
          offset_i = multi_tet_->domain_offset_[pi];
          size_i = multi_tet_->part_basis_size_[pi] + 6;
        } else {
          int local_vi = multi_tet_->full_vert_idx_[vi];
          offset_i = kFullOffset + local_vi * 3;
          size_i = 3;
        }

        for (int j = i + 1; j < 4; ++j) {
          int vj = verts[j];
          int pj = multi_tet_->vert_part_id_[vj];
          if (multi_tet_->is_subspace_domain_[pi] && multi_tet_->is_subspace_domain_[pj]) {
            continue;
          } else {
            int offset_j, size_j;
            if (multi_tet_->is_subspace_domain_[pj]) {
              offset_j = multi_tet_->domain_offset_[pj];
              size_j = multi_tet_->part_basis_size_[pj] + 6;
            } else {
              int local_vj = multi_tet_->full_vert_idx_[vj];
              offset_j = kFullOffset + local_vj * 3;
              size_j = 3;
            }
            outline.AddBlockMatrix(offset_i, offset_j, size_i, size_j);
            outline.AddBlockMatrix(offset_j, offset_i, size_j, size_i);
          }
        }
      }
    }
    InitFromOutline(&outline);
  } else {
#if 1
    SparseMatrixOutline outline(multi_tet_->constrained_dof_ + int(multi_tet_->full_verts_.size() * 3));
    for (int i = 0; i < multi_tet_->subspace_domain_num_; ++i) {
      int p = multi_tet_->subspace_domains_[i];
      int offset = multi_tet_->constrained_domain_offset_[p];
      int size = multi_tet_->domain_size_[p];
      outline.AddBlockMatrix(offset, offset, size, size);
    }

    for (int e = 0; e < multi_tet_->full_chol_solver_->e_number0; ++e) {
      int* E = multi_tet_->full_chol_solver_->E;
      int p0 = multi_tet_->subspace_domains_[E[e * 2 + 0]];
      int p1 = multi_tet_->subspace_domains_[E[e * 2 + 1]];
      int offset0 = multi_tet_->constrained_domain_offset_[p0];
      int size0 = multi_tet_->domain_size_[p0];
      int offset1 = multi_tet_->constrained_domain_offset_[p1];
      int size1 = multi_tet_->domain_size_[p1];
      outline.AddBlockMatrix(offset0, offset1, size0, size1);
      outline.AddBlockMatrix(offset1, offset0, size1, size0);
    }

    const int kFullOffset = multi_tet_->constrained_dof_;
    for (int t : multi_tet_->full_tet_) {
      int* verts = multi_tet_->tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int local_vi = multi_tet_->full_vert_idx_[vi];
        ASSERT(local_vi != -1);
        int offset_i = kFullOffset + local_vi * 3;
        outline.AddBlockMatrix(offset_i, offset_i, 3, 3);
        for (int j = i + 1; j < 4; ++j) {
          int vj = verts[j];
          int local_vj = multi_tet_->full_vert_idx_[vj];
          int offset_j = kFullOffset + local_vj * 3;
          outline.AddBlockMatrix(offset_i, offset_j, 3, 3);
          outline.AddBlockMatrix(offset_j, offset_i, 3, 3);
        }
      }
    }

    for (int t : multi_tet_->full_subspace_interface_tet_) {
      int* verts = multi_tet_->tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = multi_tet_->vert_part_id_[vi];
        int offset_i, size_i;
        if (multi_tet_->is_subspace_domain_[pi]) {
          offset_i = multi_tet_->constrained_domain_offset_[pi];
          size_i = multi_tet_->domain_size_[pi];
        } else {
          int local_vi = multi_tet_->full_vert_idx_[vi];
          offset_i = kFullOffset + local_vi * 3;
          size_i = 3;
        }

        for (int j = i + 1; j < 4; ++j) {
          int vj = verts[j];
          int pj = multi_tet_->vert_part_id_[vj];
          if (multi_tet_->is_subspace_domain_[pi] && multi_tet_->is_subspace_domain_[pj]) {
            continue;
          } else {
            int offset_j, size_j;
            if (multi_tet_->is_subspace_domain_[pj]) {
              offset_j = multi_tet_->constrained_domain_offset_[pj];
              size_j = multi_tet_->domain_size_[pj];
            } else {
              int local_vj = multi_tet_->full_vert_idx_[vj];
              offset_j = kFullOffset + local_vj * 3;
              size_j = 3;
            }
            outline.AddBlockMatrix(offset_i, offset_j, size_i, size_j);
            outline.AddBlockMatrix(offset_j, offset_i, size_j, size_i);
          }
        }
      }
    }
    InitFromOutline(&outline);
#endif
  }

  //  profiler.Start("pardiso");
  pardiso_solver_ = new PardisoSolver(this, 8, 1, 1);
  //  profiler.End("pardiso");
  //  profiler.Start("idx");
  InitColumnIndex();
  //  profiler.End("idx");
  //  profiler.End("outline");
}

void MixedSparseMatrix::Solve(double *rhs, double *x) {
  pardiso_solver_->SolveLinearSystemDirectIterative(this, x, rhs);
}

void MixedSparseMatrix::UpdateMatrix(const double dt) {
  profiler.Start("copy");
  ResetToZero();
  const double dt_2 = dt * dt;
  for (int i = 0; i < multi_tet_->subspace_domain_num_; ++i) {
    int p = multi_tet_->subspace_domains_[i];
    int row_start = multi_tet_->constrained_domain_offset_[p];
    int size = multi_tet_->domain_size_[p];
    int col_start = domain_column_idx_(i, i);
//        ASSERT(columnIndices[row_start][col_start] == multi_tet_->domain_offset_[p],
//               P(row_start, col_start, multi_tet_->domain_offset_[p])
//               P(i, p, columnIndices[row_start][col_start]));
    AddSubMatrix(row_start, col_start, size, size,
                 multi_tet_->part_basis_size_[p] + 6,
                 &multi_tet_->full_chol_solver_->Av[i].A[0]);
//    P(i, columnIndices[102][3])P(p, row_start, col_start, size);
  }
//  if (0)
  for (int e = 0; e < multi_tet_->full_chol_solver_->e_number0; ++e) {
    int i0 = multi_tet_->full_chol_solver_->E[e * 2 + 0];
    int i1 = multi_tet_->full_chol_solver_->E[e * 2 + 1];
    int p0 = multi_tet_->subspace_domains_[i0];
    int p1 = multi_tet_->subspace_domains_[i1];
    int row_start0 = multi_tet_->constrained_domain_offset_[p0];
    int col_start0 = domain_column_idx_(i0, i1);
    int size0 = multi_tet_->domain_size_[p0];
    int row_start1 = multi_tet_->constrained_domain_offset_[p1];
    int size1 = multi_tet_->domain_size_[p1];
    int col_start1 = domain_column_idx_(i1, i0);
//    KK;
//    P(e);
//    P(row_start0, col_start0);
//    P(size0, size1);
//    P(e, i0, i1);
//    P(p0, p1);
//    P(row_start0, row_start1);
//    P(col_start0, col_start1);
//    P(size0, size1);
//     ASSERT(columnIndices[row_start0][col_start0] == multi_tet_->constrained_domain_offset_[p1],
//            P(row_start0, col_start0, p0, p1);
//         P(i0, i1);
//            P(columnIndices[row_start0][col_start0], multi_tet_->constrained_domain_offset_[p1]));
    AddSubMatrix(row_start0, col_start0, size0, size1,
                 multi_tet_->part_basis_size_[p1] + 6,
                 &multi_tet_->full_chol_solver_->Ae[e].A[0]);
    //    P(e, i0, i1);
    //    P(p0, p1);
//        P(row_start1, col_start1);
//        P(size1, size0);
//        ASSERT(columnIndices[row_start1][col_start1] == multi_tet_->constrained_domain_offset_[p0], P(columnIndices[row_start1][col_start1], multi_tet_->domain_offset_[p0]));
    AddSubMatrixTransposed(row_start1, col_start1, size1, size0,
                           multi_tet_->part_basis_size_[p1] + 6,
                           &multi_tet_->full_chol_solver_->Ae[e].A[0]);
  }
//  exit(0);
  profiler.End("copy");

  const int kFullOffset = multi_tet_->constrained_dof_;
  profiler.Start("assemble");
  profiler.Start("interface");
  for (int idx = 0; idx < int(multi_tet_->full_subspace_interface_tet_.size()); ++idx) {
    int t = multi_tet_->full_subspace_interface_tet_[idx];
    Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(multi_tet_->inv_fem_->element_k_ + t * 144);
    int* verts = multi_tet_->tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      int vi = verts[i];
      int pi = multi_tet_->vert_part_id_[vi];
      if (multi_tet_->is_subspace_domain_[pi]) {
        MultiDomainTet::Vec3 ri = MultiDomainTet::MapVec3(multi_tet_->X + vi * 3) - multi_tet_->center_of_mass_[pi];
        MultiDomainTet::Mat3 skew_symmetric_mat = GetSkewSymmetrixMatrix(ri);
        //        int idxi = multi_tet_->non_constrained_dof_2_constrained_dof_[multi_tet_->domain_offset_[pi]];
        int idxi = multi_tet_->constrained_domain_offset_[pi];
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int pj = multi_tet_->vert_part_id_[vj];
          if (multi_tet_->is_subspace_domain_[pj]) {
            continue;
          } else {
            MultiDomainTet::Mat3 sub_k = dt_2 * element_k.block<3, 3>(i * 3, j * 3);
            if (!multi_tet_->is_fixed_domain_[pi]) {
              Add3x3Matrix(multi_tet_->constrained_rigid_offset_[pi],
                           interface_tet_column_idx_[idx](i, j), sub_k);
            }
            MultiDomainTet::MatRow k = multi_tet_->vert_basis_transpose_[vi] * multi_tet_->part_rotation_transpose_[pi] * sub_k;
            AddSubMatrix(idxi, interface_tet_column_idx_[idx](i, j), multi_tet_->part_basis_size_[pi], 3, k.data());
            if (!multi_tet_->is_fixed_domain_[pi]) {
              sub_k = skew_symmetric_mat * sub_k;
              Add3x3Matrix(multi_tet_->constrained_rigid_offset_[pi] + 3,
                           interface_tet_column_idx_[idx](i, j), sub_k);
            }
          }
        }
      } else {
        int idxi = multi_tet_->constrained_dof_ + multi_tet_->full_vert_idx_[vi] * 3;
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int pj = multi_tet_->vert_part_id_[vj];
          MultiDomainTet::Mat3 sub_k = dt_2 * element_k.block<3, 3>(i * 3, j * 3);
          if (multi_tet_->is_subspace_domain_[pj])  {
            MultiDomainTet::Vec3 rj = MultiDomainTet::MapVec3(multi_tet_->X + vj * 3) - multi_tet_->center_of_mass_[pj];
            MultiDomainTet::Mat3 skew_symmetric_mat = GetSkewSymmetrixMatrix(-rj);
            //           ASSERT(columnIndices[idxi][interface_tet_column_idx_[idx](i, j)] == multi_tet_->full_vert_idx_[verts[j]] * 3 + kFullOffset);
            if (!multi_tet_->is_fixed_domain_[pj]) {
              Add3x3Matrix(idxi, interface_tet_column_idx_[idx](i, j) + multi_tet_->part_basis_size_[pj] + 0, sub_k);
            }

            MultiDomainTet::MatRow k = sub_k * multi_tet_->part_rotation_[pj] * multi_tet_->vert_basis_[vj];
            AddSubMatrix(idxi, interface_tet_column_idx_[idx](i, j), 3, multi_tet_->part_basis_size_[pj], k.data());
            if (!multi_tet_->is_fixed_domain_[pj]) {
              sub_k = sub_k * skew_symmetric_mat;
              Add3x3Matrix(idxi, interface_tet_column_idx_[idx](i, j) + multi_tet_->part_basis_size_[pj] + 3, sub_k);
            }
          } else {
            Add3x3Matrix(idxi, interface_tet_column_idx_[idx](i, j), sub_k);
          }
        }
      }
    }
  }
  profiler.End("interface");

  profiler.Start("inner tet");
  for (int idx = 0; idx < int(multi_tet_->full_tet_.size()); ++idx) {
    int t = multi_tet_->full_tet_[idx];
    int* verts = multi_tet_->tet_ + t * 4;
    Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(multi_tet_->inv_fem_->element_k_ + t * 144);
    for (int i = 0; i < 4; ++i) {
      int vi = multi_tet_->full_vert_idx_[verts[i]];
      int row = vi * 3 + kFullOffset;
      for (int j = 0; j < 4; ++j) {
        MultiDomainTet::Mat3 sub_k = dt_2 * element_k.block<3, 3>(i * 3, j * 3);
        Add3x3Matrix(row, full_tet_column_idx_[idx](i, j), sub_k);
      }
    }
  }
  profiler.End("inner tet");

  for (int i = 0; i < int(multi_tet_->full_verts_.size()); ++i) {
    int v = multi_tet_->full_verts_[i];
    int column_idx = diag_column_idx_[i];
    int row = kFullOffset + i * 3;
    columnEntries[row + 0][column_idx + 0] += multi_tet_->mass_[v];
    columnEntries[row + 1][column_idx + 1] += multi_tet_->mass_[v];
    columnEntries[row + 2][column_idx + 2] += multi_tet_->mass_[v];
  }
  profiler.End("assemble");
  if (0) {
    auto A = [&](double * x, double * rhs) {
      MultiplyVector(x, rhs);
    };
    P(multi_tet_->total_basis_num_ + multi_tet_->subspace_domain_num_ * 6);
    P(Getn());
    L("/tmp/test");
    dj::WriteImplicitMatrixToMatlab<double>("/tmp/test", A, Getn());
    L("exit");
    exit(0);
  }
}

void MixedSparseMatrix::InitColumnIndex() {
  domain_column_idx_ = Mati::Zero(multi_tet_->subspace_domain_num_, multi_tet_->subspace_domain_num_ + 1);
  for (int i = 0; i < domain_column_idx_.rows(); ++i) {
    int prev_size = 0;
    int pi = multi_tet_->subspace_domains_[i];
    for (int j = 0; j < domain_column_idx_.cols() - 1; ++j) {
      int pj = multi_tet_->subspace_domains_[j];
      if (multi_tet_->topology_[pi][pj] >= 0) {
        domain_column_idx_(i, j) = prev_size;
        //        P(i, j, domain_column_idx_(i, j));
        prev_size += multi_tet_->domain_size_[pj];
//        P(multi_tet_->domain_size_[pj], pj);
//        P(i, j, domain_column_idx_(i, j));
      } else {
        domain_column_idx_(i, j) = INT_MIN;
      }
    }
    domain_column_idx_(i, multi_tet_->subspace_domain_num_) = prev_size;
  }

  interface_tet_column_idx_.resize(multi_tet_->full_subspace_interface_tet_.size());
  full_tet_column_idx_.resize(multi_tet_->full_tet_.size());
  std::unordered_map<int, std::unordered_set<int>> topology; // for each vertex (or domain) list its neighbors
  std::map<std::pair<int, int>, int> column_indices;
  for (int t : multi_tet_->full_subspace_interface_tet_) {
    int* verts = multi_tet_->tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      int vi = verts[i];
      int pi = multi_tet_->vert_part_id_[vi];
      if (multi_tet_->is_subspace_domain_[pi]) {
        vi = pi - multi_tet_->part_num_;
      } else {
        vi = multi_tet_->full_vert_idx_[vi];
      }
      for (int j = i; j < 4; ++j) {
        int vj = verts[j];
        int pj = multi_tet_->vert_part_id_[vj];
        if (multi_tet_->is_subspace_domain_[pj] && multi_tet_->is_subspace_domain_[pi]) {
          continue;
        } else {
          if (multi_tet_->is_subspace_domain_[pj]) {
            vj = pj - multi_tet_->part_num_;
          } else {
            vj = multi_tet_->full_vert_idx_[vj];
          }
          topology[vi].insert(vj);
          topology[vj].insert(vi);
        }
      }
    }
  }
  for (int t : multi_tet_->full_tet_) {
    int* verts = multi_tet_->tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      int vi = multi_tet_->full_vert_idx_[verts[i]];
      topology[vi].insert(vi);
      for (int j = i + 1; j < 4; ++j) {
        int vj = multi_tet_->full_vert_idx_[verts[j]];
        topology[vi].insert(vj);
        topology[vj].insert(vi);
      }
    }
  }

  //  P(multi_tet_->subspace_domain_basis_size_ + multi_tet_->subspace_domain_num_ * 6);
  for (auto t : topology) {
    int v = t.first;
    auto& neighbor_set = t.second;
    std::vector<int> neighbors(neighbor_set.begin(), neighbor_set.end());
    std::sort(neighbors.begin(), neighbors.end());

    int column = 0;
    int i = 0;
    if (v < 0) {
      int p = v + multi_tet_->part_num_;
      for (; i < int(neighbors.size()) && neighbors[i] < 0; ++i) {
      }
      int offset = domain_column_idx_(multi_tet_->domain_index_[p], multi_tet_->subspace_domain_num_);
      column = offset;
    }

    for (; i < int(neighbors.size()); ++i) {
      int n  = neighbors[i];
      //      if (v == 319) {
      //        if (n >= 0) P(multi_tet_->full_verts_[n], multi_tet_->full_verts_[v]);
      //        P(v, column, n, multi_tet_->subspace_domain_basis_size_ + multi_tet_->subspace_domain_num_ * 6 + n * 3);
      //      }
      column_indices[std::make_pair(v, n)] = column;
      if (n < 0) {
        int p = n + multi_tet_->part_num_;
        column += multi_tet_->domain_size_[p];
      } else {
        column += 3;
      }
    }
  }
  //  for (auto pair : multi_tet_->adjacent_vertex_[8429]) {
  //    P(pair.first, pair.second);
  //  }
  //  exit(0);

  //  P(multi_tet_->full_subspace_interface_tet_.size());
  for (int idx = 0; idx < int(multi_tet_->full_subspace_interface_tet_.size()); ++idx) {
    int t = multi_tet_->full_subspace_interface_tet_[idx];
    int* verts = multi_tet_->tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      int vi = verts[i];
      int pi = multi_tet_->vert_part_id_[vi];
      if (multi_tet_->is_subspace_domain_[pi]) {
        vi = pi - multi_tet_->part_num_;
      } else {
        vi = multi_tet_->full_vert_idx_[vi];
      }
      for (int j = 0; j < 4; ++j) {
        int vj = verts[j];
        int pj = multi_tet_->vert_part_id_[vj];
        if (multi_tet_->is_subspace_domain_[pj]) {
          vj = pj - multi_tet_->part_num_;
        } else {
          vj = multi_tet_->full_vert_idx_[vj];
        }
        if (multi_tet_->is_subspace_domain_[pj] && multi_tet_->is_subspace_domain_[pi]) {
          interface_tet_column_idx_[idx](i, j) = domain_column_idx_(multi_tet_->domain_index_[pi], multi_tet_->domain_index_[pj]);
        } else {
          interface_tet_column_idx_[idx](i, j) = column_indices[std::make_pair(vi, vj)];
        }
      }
    }
  }

  for (int idx = 0; idx < int(multi_tet_->full_tet_.size()); ++idx) {
    int t = multi_tet_->full_tet_[idx];
    int* verts = multi_tet_->tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      int vi = multi_tet_->full_vert_idx_[verts[i]];
      for (int j = 0; j < 4; ++j) {
        int vj = multi_tet_->full_vert_idx_[verts[j]];
        full_tet_column_idx_[idx](i, j) = column_indices[std::make_pair(vi, vj)];
      }
    }
  }

  if (1) {
    diag_column_idx_.reserve(20000);
    diag_column_idx_.resize(multi_tet_->full_verts_.size());

    const int kFullOffset = multi_tet_->constrained_dof_;
    for (int i = 0; i < int(multi_tet_->full_verts_.size()); ++i) {
      int row = kFullOffset + i * 3;
      int idx = GetInverseIndex(row, row);
      ASSERT(idx != -1, P(i, row, idx, multi_tet_->full_verts_.size()));
      diag_column_idx_[i] = idx;
    }
  }
}

