#include <queue>
#include <functional>
#include <limits.h>
#include <utility>
#include "rainbow_color.h"
#include "global.h"
#include "opengl_helper.h"
#include "conjugate_gradient_solver.h"
#include "config_file.h"
#include "vector_lib.h"
#include "tet_mesh_simulator_bridge.h"
#include "triangle_mesh_render.h"
#include "isotropicHyperelasticFEM.h"
#include "textured_triangle_mesh_renderer.h"
#include "mixed_multi_domain_tet.h"
#include "biconjugate_gradient_stablized.h"
#include "objMeshRender.h"

#include "BLOCK_MATRIX_GRAPH.h"
#include "matlab_io.h"
#include "mixed_sparse_matrix.h"


const int MixedMultiDomainTet::kMaxVertexNum = 50000;
const int MixedMultiDomainTet::kMaxEdgeNum = 300000;
const int MixedMultiDomainTet::kMaxTriNum = 300000;
const int MixedMultiDomainTet::kMaxTetNum = 200000;

template <class T> FORCEINLINE
T DetMatrix(T *x) {
  return x[0] * (x[4] * x[8] - x[7] * x[5]) + x[3] * (x[7] * x[2] - x[1] * x[8]) + x[6] * (x[1] * x[5] - x[4] * x[2]);
}


inline void SwapIfGreater(int& v0, int& v1) {
  if (v0 > v1) dj::Swap(v0, v1);
}

MixedMultiDomainTet::MixedMultiDomainTet(const char *mesh_file, AffineTransformer<double> *transformer)
  : Super(mesh_file, 0, transformer) {
  original_vertex_num_ = vertex_num_;
  original_edge_num_ = edge_num_;
  original_tet_num_ = tet_number;

  // vertex data
  surface_vertices_.reserve(kMaxVertexNum);
  surface_vert_pos_.reserve(kMaxTriNum * 3 * 3);
  surface_vert_normal_.reserve(kMaxTriNum * 3 * 3);
  surface_vert_texture_coord_.reserve(kMaxTriNum * 3 * 3);
  is_cut_triangle_ = std::vector<char>(kMaxTriNum, false);
  vertex_id_ = buffer_.Malloc<int>(kMaxVertexNum * 3, -1);
  for (int v = 0; v < vertex_num_; ++v) {
    vertex_id_[v] = v;
  }
  is_surface_vert_ = buffer_.Malloc<char>(kMaxVertexNum * 2);
  global_vertex_id2surface_vertex_id_ = buffer_.Malloc<int>(kMaxVertexNum * 2);
  global_vertex_id2surface_vertex_id_ += kMaxVertexNum;
  is_surface_vert_ += kMaxVertexNum;
  double* new_mass = buffer_.Malloc<double>(kMaxVertexNum);
  double* new_rest_pos = buffer_.Malloc<double>(kMaxVertexNum * 3);
  double* new_X = buffer_.Malloc<double>(kMaxVertexNum * 3 * 2);
  double* new_VN = buffer_.Malloc<double>(kMaxVertexNum * 3 * 2);
  double* new_velocity = buffer_.Malloc<double>(kMaxVertexNum * 3 * 2);
  memcpy(new_rest_pos, rest_pos_, sizeof(double) * 3 * vertex_num_);
  buffer_.Free(rest_pos_);
  memcpy(new_X + kMaxVertexNum * 3, X, sizeof(double) * 3 * vertex_num_);
  buffer_.Free(X);
  buffer_.Free(VN);
  buffer_.Free(velocity_);
  memcpy(new_mass, mass_, sizeof(double) * vertex_num_);
  buffer_.Free(mass_);
  mass_ = new_mass;
  rest_pos_ = new_rest_pos;
  X = new_X + kMaxVertexNum * 3;
  VN = new_VN + kMaxVertexNum * 3;
  velocity_ = new_velocity + kMaxVertexNum * 3;
  embeded_vertex_.reserve(kMaxVertexNum);

  // edge data
  incident_tet_on_edge_.reserve(edge_num_ * 2);
  mid_point_weight_ = buffer_.Malloc<double>(kMaxEdgeNum);
  std::fill_n(mid_point_weight_, edge_num_, -1.0);
  edge_id_ = buffer_.Malloc<int>(kMaxEdgeNum);
  for (int e = 0; e < edge_num_; ++e) {
    edge_id_[e] = e;
  }
  int* new_edge = buffer_.Malloc<int>(kMaxEdgeNum * 2);
  memcpy(new_edge, edges_, sizeof(int) * 2 * edge_num_);
  original_edge_ = edges_;
  edges_ = new_edge;

  // triangle data
  surface_triangle_indices_.reserve(kMaxTriNum * 3);
  int* new_T = buffer_.Malloc<int>(kMaxTriNum * 3);
  double* new_TN = buffer_.Malloc<double>(kMaxTriNum * 3);
  memcpy(new_T, T, sizeof(int) * triangle_num_ * 3);
  buffer_.Free(T);
  buffer_.Free(TN);
  T = new_T;
  TN = new_TN;

  // tet data
  surface_tri_on_tet_.resize(kMaxTetNum);
  edge_on_tet_ = buffer_.Malloc<int>(kMaxTetNum * 6);
  tet_id_ = buffer_.Malloc<int>(kMaxTetNum);
  for (int t = 0; t < tet_number; ++t) {
    tet_id_[t] = t;
  }
  tet_type_ = buffer_.Malloc<int>(kMaxTetNum);
  std::fill_n(tet_type_, tet_number, kFullTet);
  int* new_tet = buffer_.Malloc<int>(kMaxTetNum * 4);
  memcpy(new_tet, tet_, sizeof(int) * 4 * tet_number);
  buffer_.Free(tet_);
  tet_ = new_tet;
  embeded_v_num_ = 0;

  //  for (int e = 0; e < edge_num_; ++e) {
  //    P(e, dj::Vec2i(edges_ + e * 2));
  //  }

  bi_cg_ = new BiconjugateGradientStablized<double>(kMaxVertexNum);
  delete cg_solver_;
  cg_solver_ = new ConjugateGradientSolver<double>(kMaxVertexNum * 3);
  cg_solver_->Resize(vertex_num_ * 3);
  BuildVertexAdjacencyList();
  BuildEdgeTetIndex();
  BuildTetSurfaceTriIndex();

  is_constrainted_ = std::vector<int>(kMaxVertexNum, 0);
  sparse_k_ = nullptr;
  position_changed_ = true;
  //  for (int v = 0; v < vertex_num_; ++v) {
  ////    if (rest_pos_[v * 3 + 2] > 0.65) {
  //    if (rest_pos_[v * 3 + 2] < 0.01) {
  ////    if (rest_pos_[v * 3 + 0] < 0.01) {
  //      is_constrainted_[v] = 1;
  //    }
  //  }

  //  is_constrainted_[0] = true;
  //  is_constrainted_[1] = true;
  //  is_constrainted_[3] = true;
  //  is_constrainted_[4] = true;
}

void MixedMultiDomainTet::Build_VN() {
  //  return;
  Build_TN();
  memset(VN, 0, sizeof(double)*vertex_num_ * 3);
  for (int i = 0; i < triangle_num_; i++) {
    if (is_cut_triangle_[i]) continue;
    int v0 = T[i * 3 + 0];
    int v1 = T[i * 3 + 1];
    int v2 = T[i * 3 + 2];

    VN[v0 * 3 + 0] += TN[i * 3 + 0];
    VN[v0 * 3 + 1] += TN[i * 3 + 1];
    VN[v0 * 3 + 2] += TN[i * 3 + 2];

    VN[v1 * 3 + 0] += TN[i * 3 + 0];
    VN[v1 * 3 + 1] += TN[i * 3 + 1];
    VN[v1 * 3 + 2] += TN[i * 3 + 2];

    VN[v2 * 3 + 0] += TN[i * 3 + 0];
    VN[v2 * 3 + 1] += TN[i * 3 + 1];
    VN[v2 * 3 + 2] += TN[i * 3 + 2];
  }
}

void MixedMultiDomainTet::BuildEdgeTetIndex() {
  for (int e = 0; e < edge_num_; ++e) {
    int v0 = edges_[e * 2 + 0];
    int v1 = edges_[e * 2 + 1];
    if (v0 > v1) std::swap(v0, v1);
    std::pair<int, int> vv(v0, v1);
    vertex_pair2edge_[vv] = e;
  }

  incident_tet_on_edge_.resize(edge_num_);
  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    int idx = 0;
    for (int i = 0; i < 4; ++i) {
      for (int j = i + 1; j < 4; ++j, ++idx) {
        int v0 = verts[i];
        int v1 = verts[j];
        if (v0 > v1) std::swap(v0, v1);
        std::pair<int, int> vv(v0, v1);
        int e = vertex_pair2edge_[vv];
        edge_on_tet_[t * 6 + idx] = e;
        incident_tet_on_edge_[e].insert(t);
      }
    }
  }
  //  P(incident_tet_on_edge_[0].size())
  //  P(incident_tet_on_edge_[1].size())
  //  P(incident_tet_on_edge_[2].size())
  //  P(incident_tet_on_edge_[3].size())
  //  P(incident_tet_on_edge_[4].size())
  //  P(incident_tet_on_edge_[5].size())
}

void MixedMultiDomainTet::BuildVertexAdjacencyList() {
  adjacent_vertex_.resize(kMaxVertexNum);
  for (int e = 0; e < edge_num_; ++e) {
    int v0 = edges_[e * 2 + 0];
    int v1 = edges_[e * 2 + 1];
    adjacent_vertex_[v0][v1] = e;
    adjacent_vertex_[v1][v0] = e;
    SwapIfGreater(v0, v1);
    if (v0 == 0 && v1 == 151) {
      ASSERT(false, L("found"));
    }
  }
}

inline void MixedMultiDomainTet::RenderTet(int tet_id) {
  //  if (tet_id != 39 && tet_id != -1) return;
  //  if (tet_id != 39) return;
  auto RenderTet = [](double * v0, double * v1, double * v2, double * v3) {
    //    return;
    //    if (v0[2] > 1 || v1[2] > 1 || v2[2] > 1 || v3[2] > 1) {
    //      P(dj::Vec3d(v0));
    //      P(dj::Vec3d(v1));
    //      P(dj::Vec3d(v2));
    //      P(dj::Vec3d(v3));
    //    }
    //    return;
    Vertex3v(v0); Vertex3v(v1);
    Vertex3v(v0); Vertex3v(v2);
    Vertex3v(v0); Vertex3v(v3);
    Vertex3v(v1); Vertex3v(v2);
    Vertex3v(v1); Vertex3v(v3);
    Vertex3v(v2); Vertex3v(v3);
  };
  auto RenderPrism = [](double * v0, double * v1, double * v2, double * v3, double * v4, double * v5) {
    //    return;
    Vertex3v(v0); Vertex3v(v1);
    Vertex3v(v0); Vertex3v(v2);
    Vertex3v(v1); Vertex3v(v2);
    Vertex3v(v3); Vertex3v(v4);
    Vertex3v(v3); Vertex3v(v5);
    Vertex3v(v4); Vertex3v(v5);
    Vertex3v(v0); Vertex3v(v3);
    Vertex3v(v1); Vertex3v(v4);
    Vertex3v(v2); Vertex3v(v5);
  };

  auto RenderFrustrum = [](double * v0, double * v1, double * v2, double * v3, double * e0, double * e1) {
    //    return;
    Vertex3v(v0); Vertex3v(v1);
    Vertex3v(v1); Vertex3v(v2);
    Vertex3v(v2); Vertex3v(v3);
    Vertex3v(v3); Vertex3v(v0);
    Vertex3v(e0); Vertex3v(e1);
    Vertex3v(e0); Vertex3v(v0);
    Vertex3v(e0); Vertex3v(v1);
    Vertex3v(e1); Vertex3v(v2);
    Vertex3v(e1); Vertex3v(v3);
  };

  //  const int kVerNumTable[16] = {
  //    -1, // 0000, 0
  //    +1, // 0001, 1
  //    +1, // 0010, 2
  //    +2, // 0011, 3
  //    +1, // 0100, 4
  //    +2, // 0101, 5
  //    +2, // 0110, 6
  //    +3, // 0111, 7
  //    +1, // 1000, 8
  //    +2, // 1001, 9
  //    +2, // 1010, 10
  //    +3, // 1011, 11
  //    +2, // 1100, 12
  //    +3, // 1101, 13
  //    +3, // 1110, 14
  //    +4, // 1111, 15
  //  };
  int idx0 = -1;
  int flag[2][2] = {
    -1, -1,
    -1, -1,
  };
  switch (tet_type_[tet_id]) {
    case 1:
      idx0 = 0;
    case 2:
      if (idx0 < 0) idx0 = 1;
    case 4:
      if (idx0 < 0) idx0 = 2;
    case 8: {
      if (idx0 < 0) idx0 = 3;
      //      KK;
      //      P(dj::Vec3d(X + tet_[tet_id * 4 + 0] * 3));
      //      P(dj::Vec3d(X + tet_[tet_id * 4 + 1] * 3));
      //      P(dj::Vec3d(X + tet_[tet_id * 4 + 2] * 3));
      //      P(dj::Vec3d(X + tet_[tet_id * 4 + 3] * 3));
      //      exit(0);
      int v = tet_[tet_id * 4 + idx0];
      Vec3 p[3];
      int count = 0;
      for (int i = 0; i < 4; ++i) {
        if (i == idx0) continue;
        int v0 = v;
        int v1 = tet_[tet_id * 4 + i];
        int old_v0 = vertex_id_[v0];
        int old_v1 = vertex_id_[v1];
        if (old_v0 > old_v1) std::swap(old_v0, old_v1);
        std::pair<int, int> pair(old_v0, old_v1);
        int e = vertex_pair2edge_[pair];
        double weight = mid_point_weight_[e];
        //        P(weight, tet_type_[tet_id], i, idx0);
        if (vertex_id_[v0] != original_edge_[e * 2 + 0]) {
          //          ASSERT(vertex_id[v1] == original_edge_[e * 2 + 0], P(v0, v1, dj::Vec2i(original_edge_ + e * 2)));
          weight = 1 - weight;
        }
        p[count] = weight * MapVec3(X + v0 * 3) + (1 - weight) * MapVec3(X + v1 * 3);
        if (0 && p[count][2] > 1) {
          PVEC(p[count]);
          P(weight, i);
          P(old_v0, old_v1, e);
          P(tet_id, tet_id_[tet_id]);
          P(dj::Vec3d(X + v0 * 3));
          P(dj::Vec3d(X + v1 * 3));
          //          exit(0);
        }
        ++count;
      }
      //      return;
      RenderTet(X + v * 3, &p[0][0], &p[1][0], &p[2][0]);
      break;

    }

    case 3: // 0011
      if (flag[0][0] < 0) {
        flag[0][0] = 0; flag[0][1] = 1;
        flag[1][0] = 2; flag[1][1] = 3;
      }
    case 5: // 0101
      if (flag[0][0] < 0) {
        flag[0][0] = 0; flag[0][1] = 2;
        flag[1][0] = 1; flag[1][1] = 3;
      }
    case 6: // 0110
      if (flag[0][0] < 0) {
        flag[0][0] = 1; flag[0][1] = 2;
        flag[1][0] = 0; flag[1][1] = 3;
      }
    case 9: // 1001
      if (flag[0][0] < 0) {
        flag[0][0] = 0; flag[0][1] = 3;
        flag[1][0] = 1; flag[1][1] = 2;
      }
    case 10: // 1010
      if (flag[0][0] < 0) {
        flag[0][0] = 1; flag[0][1] = 3;
        flag[1][0] = 0; flag[1][1] = 2;
      }
    case 12: { // 1100
      if (flag[0][0] < 0) {
        flag[0][0] = 2; flag[0][1] = 3;
        flag[1][0] = 0; flag[1][1] = 1;
      }
      double* q[2] = {
        X + tet_[tet_id * 4 + flag[0][0]] * 3,
        X + tet_[tet_id * 4 + flag[0][1]] * 3,
      };
      Vec3 p[4];
      int count = 0;
      for (int i = 0; i < 2; ++i) {
        int v0 = tet_[tet_id * 4 + flag[0][i]];
        for (int j = 0; j < 2; ++j, count++) {
          int v1;
          if (i == 0) {
            v1 = tet_[tet_id * 4 + flag[1][j]];
          } else {
            v1 = tet_[tet_id * 4 + flag[1][1 - j]];
          }
          int old_v0 = vertex_id_[v0];
          int old_v1 = vertex_id_[v1];
          if (old_v0 > old_v1) std::swap(old_v0, old_v1);
          std::pair<int, int> pair(old_v0, old_v1);
          int e = vertex_pair2edge_[pair];
          double weight = mid_point_weight_[e];
          //          ASSERT(weight >= 0 && weight <= 1, int* a = 0; *a = 0;);
          if (vertex_id_[v0] != original_edge_[e * 2 + 0]) {
            //                        ASSERT(vertex_id[v1] == original_edge_[e * 2 + 0], P(v0, v1, dj::Vec2i(original_edge_ + e * 2)));
            weight = 1 - weight;
          }
          //          P(count, old_v0, old_v1, v0);
          p[count] = weight * MapVec3(X + v0 * 3) + (1 - weight) * MapVec3(X + v1 * 3);
        }
      }
      //      P(dj::Vec3d(q[0]));
      //      P(dj::Vec3d(q[1]));
      //      PVEC(p[0]);
      //      PVEC(p[1]);
      //      PVEC(p[2]);
      //      PVEC(p[3]);
      RenderFrustrum(&p[0][0], &p[1][0], &p[2][0], &p[3][0], q[0], q[1]);
      //      exit(0);
      break;
    }

    case 7:
      idx0 = 3;
    case 11:
      if (idx0 < 0) idx0 = 2;
    case 13:
      if (idx0 < 0) idx0 = 1;
    case 14: {
      if (idx0 < 0) idx0 = 0;
      int v = tet_[tet_id * 4 + idx0];
      Vec3 p[3];
      double* q[3];
      int count = 0;
      for (int i = 0; i < 4; ++i) {
        if (i == idx0) continue;
        int v0 = v;
        int v1 = tet_[tet_id * 4 + i];
        int old_v0 = vertex_id_[v0];
        int old_v1 = vertex_id_[v1];
        if (old_v0 > old_v1) std::swap(old_v0, old_v1);
        std::pair<int, int> pair(old_v0, old_v1);
        int e = vertex_pair2edge_[pair];
        double weight = mid_point_weight_[e];
        if (vertex_id_[v0] != original_edge_[e * 2 + 0]) {
          //          ASSERT(vertex_id[v1] == original_edge_[e * 2 + 0], P(v0, v1, dj::Vec2i(original_edge_ + e * 2)));
          weight = 1 - weight;
        }
        p[count] = weight * MapVec3(X + v0 * 3) + (1 - weight) * MapVec3(X + v1 * 3);
        q[count] = X + v1 * 3;
        ++count;
      }
      //      return;
      RenderPrism(q[0], q[1], q[2],
                  &p[0][0], &p[1][0], &p[2][0]);
      break;
    }
    //      case 15: {
    //        break;
    //      }
    default: {
      //      return;
      for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
          Vertex3v(X + tet_[tet_id * 4 + i] * 3);
          Vertex3v(X + tet_[tet_id * 4 + j] * 3);
        }
      }

      break;
    }
  }
}

// apply fictitious force and gravity on subspace domains
void MixedMultiDomainTet::AddFictitiousForceAndGravity(const std::vector<MultiDomainTet::Vec3> &acceleration,
                                                       const std::vector<MultiDomainTet::Vec3> &angular_acceleration,
                                                       MultiDomainTet::Vec &subspace_rhs) {
  profiler.Start("fictitious_force");
  //  L("no fictitious force")
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    Mat3 rotation_transpose = part_rotation_[p].transpose();
    Vec3 acc = rotation_transpose * (acceleration[p] - gravity_);
    Vec3 angular_acc = rotation_transpose * angular_acceleration[p];
    MapVec subspace_force(&subspace_rhs[domain_offset_[p]], part_basis_size_[p]);
    MapVec map_vel_q(&part_vel_[p][0], part_basis_size_[p]);
    MapVec map_q(&part_q_[p][0], part_basis_size_[p]);
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
  profiler.End("fictitious_force");
}

void MixedMultiDomainTet::MixSimulaitonOneStepWithCubature(double dt) {
  profiler.Start("simulation");
  const int kFullOffset = subspace_domain_basis_size_ + subspace_domain_num_ * 6;
  const double dt_2 = dt * dt;
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  internal_force_scaling_factor = 1.0;

  profiler.Start("element stiffness");
  inv_fem_->inv_fem_force_model_->ComputePartialForceAndStiffnessMatrixWithSplitedMesh(&all_cubature_tet_[0],
                                                                                       int(all_cubature_tet_.size()),
                                                                                       tet_id_,
                                                                                       tet_, X);
  profiler.End("element stiffness");
  profiler.Start("ext. force");
  GetExternalForce(dt, &ext_force_);
  profiler.End("ext. force");
  Vec rhs_all = Vec::Zero(total_basis_num_ + subspace_domain_num_ * 6);
  // Compute rhs
  profiler.Start("rhs");
  {
    // external force
    profiler.Start("external f.");
    for (ExtForce & f : ext_force_) {
      int& v = f.first;
      Vec3& force = f.second;
      int p = vert_part_id_[v];
      if (is_subspace_domain_[p]) {
        MapVec subspace_force(&rhs_all[domain_offset_[p]], part_basis_size_[p]);
        subspace_force += vert_basis_transpose_[v] * (part_rotation_[p].transpose() * force);
        Vec3 r = MapVec3(X + v * 3) - center_of_mass_[p];
        double* rigid_part = &rhs_all[rigid_offset_[p]];
        // momentum
        MapVec3(rigid_part + 0) += force;
        // torque
        MapVec3(rigid_part + 3) += r.cross(force);
      } else {
        int local_v = full_vert_idx_[v];
        MapVec3 full_force(&rhs_all[kFullOffset + local_v * 3]);
        full_force += force;
      }
    }
    profiler.End("external f.");
    // internal force
    profiler.Start("internl f.");
    for (CubaturePoint & cubature : all_cubature_) {
      int& t = cubature.first;
      double& weight = cubature.second;
      double* force = inv_fem_->element_force_ + t * 12;
      for (int i4 = 0; i4 < 4; ++i4) {
        int v = tet_[t * 4 + i4];
        int p = vert_part_id_[v];
        if (is_subspace_domain_[p]) {
          Vec3 weighted_force = weight * MapVec3(force + i4 * 3);
          MapVec subspace_force(&rhs_all[domain_offset_[p]], part_basis_size_[p]);
          // negative sign before weight is because internal force from vega is in oppisite direction of actual force
          subspace_force -= vert_basis_transpose_[v] * (part_rotation_[p].transpose() * weighted_force);
          Vec3 r = MapVec3(X + v * 3) - center_of_mass_[p];
          double* rigid_part = &rhs_all[rigid_offset_[p]];
          // momentum
          MapVec3(rigid_part + 0) -= weighted_force;
          // torque
          MapVec3(rigid_part + 3) -= r.cross(weighted_force);
        } else {
          int local_v = full_vert_idx_[v];
          MapVec3 full_force(&rhs_all[kFullOffset + local_v * 3]);
          full_force -= MapVec3(force + i4 * 3);
        }
      }
    }
    profiler.End("internl f.");
    // subspace fictitious force
    profiler.Start("fictitious f.");
    AddFictitiousForceAndGravity(acceleration_, angular_acceleration_, rhs_all);
    profiler.End("fictitious f.");
    profiler.Start("rhs vel");
    rhs_all *= dt;
    // add gravity and velocity related term
    for (int p = 0; p < part_num_; ++p) {
      if (is_subspace_domain_[p]) {
        double* subspace_rhs = &rhs_all[domain_offset_[p]];
        double* rigid_rhs = &rhs_all[rigid_offset_[p]];
        MapVec(subspace_rhs, part_basis_size_[p]) += MapVec(&part_vel_[p][0], part_basis_size_[p]);
        MapVec3(rigid_rhs + 0) += mass_per_part_[p] * (translational_vel_[p] + gravity_ * dt);
        MapVec3(rigid_rhs + 3) += current_inertia_tensor_[p] * angular_vel_[p];
      } else {
        continue;
      }
    }
    for (int i = 0; i < int(full_verts_.size()); ++i) {
      int global_v = full_verts_[i];
      int offset = kFullOffset + i * 3;
      MapVec3(&rhs_all[offset]) += mass_[global_v] * (MapVec3(velocity_ + global_v * 3) + gravity_ * dt);
    }
    profiler.End("rhs vel");
  }
  profiler.End("rhs");

  //#define EXPLICIT
#ifdef EXPLICIT
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      double* subspace_rhs = &rhs_all[domain_offset_[p]];
      double* rigid_rhs = &rhs_all[rigid_offset_[p]];
      MapVec(&part_vel_[p][0], part_basis_size_[p]) = MapVec(subspace_rhs, part_basis_size_[p]);
      //      MapVec(&part_vel_[p][0], part_basis_size_[p]).setZero();
      translational_vel_[p] = MapVec3(rigid_rhs + 0) / mass_per_part_[p];
      angular_vel_[p] = current_inv_inertia_tensor_[p] * MapVec3(rigid_rhs + 3);
      //      angular_vel_[p].setZero();
    } else {
      continue;
    }
  }

  for (int i = 0; i < int(full_verts_.size()); ++i) {
    int global_v = full_verts_[i];
    int offset = kFullOffset + i * 3;
    MapVec3(velocity_ + global_v * 3) = MapVec3(&rhs_all[offset]) / mass_[global_v];
  }
#else
  MatCol rigid_matrix = MatCol::Zero(subspace_domain_num_ * 6, subspace_domain_num_ * 6);
  profiler.Start("rigid k");
  {
    for (int e = 0; e < int(interface_domains_.size()); ++e) {
      if (interface_domains_[e][0] > 2) continue;
      //      if (!is_subspace_domain_[p0] || !is_subspace_domain_[p1]) continue;
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int& t = cubature.first;
        double& weight = cubature.second;
        // the force is computed from vega, the direction of the force is opposite to the actual direction
        Eigen::Map<Eigen::Matrix<double, 12, 12>> element_k(inv_fem_->element_k_ + t * 144);
        int* verts = tet_ + t * 4;
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          if (!is_subspace_domain_[pi]) continue;
          int idxi = domain_index_[pi] * 6;
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          auto skew_symmetric_mat_i = GetSkewSymmetrixMatrix(ri);
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (!is_subspace_domain_[pj]) continue;
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            int idxj = domain_index_[pj] * 6;
            auto skew_symmetric_mat_j = GetSkewSymmetrixMatrix(rj);
            rigid_matrix.block<3, 3>(idxi + 0, idxj + 0) += element_k.block<3, 3>(i * 3, j * 3) * weight;
            rigid_matrix.block<3, 3>(idxi + 0, idxj + 3) -= element_k.block<3, 3>(i * 3, j * 3) * skew_symmetric_mat_j * weight;
            rigid_matrix.block<3, 3>(idxi + 3, idxj + 0) += skew_symmetric_mat_i * element_k.block<3, 3>(i * 3, j * 3) * weight;
            rigid_matrix.block<3, 3>(idxi + 3, idxj + 3) -= skew_symmetric_mat_i * element_k.block<3, 3>(i * 3, j * 3) * skew_symmetric_mat_j * weight;
          }
        }
      }
    }
    rigid_matrix *= (dt_2);
    // Diagonal terms
    for (int p = 0; p < subspace_domain_num_; ++p) {
      rigid_matrix(p * 6 + 0, p * 6 + 0) += mass_per_part_[subspace_domains_[p]];
      rigid_matrix(p * 6 + 1, p * 6 + 1) += mass_per_part_[subspace_domains_[p]];
      rigid_matrix(p * 6 + 2, p * 6 + 2) += mass_per_part_[subspace_domains_[p]];
      rigid_matrix.block<3, 3>(p * 6 + 3, p * 6 + 3) += current_inertia_tensor_[subspace_domains_[p]];
    }
  }
  profiler.End("rigid k");
  //  dj::WriteEigenMatrixToMatlab(rigid_matrix, "/tmp/rigid_k");
  //  exit(0);
  profiler.Start("local k");
  {
    chol_solver_->SetMatrixZero();
    profiler.Start("mul p");
    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      double* mat = chol_solver_->Av[i].A;
      Eigen::Map<MatRow> reduced_k(mat, part_basis_size_[p], part_basis_size_[p]);
      //    reduced_k.setZero();
      // all cubature vertex in this domain
      for (int j = 0; j < int(v_list_[p].size()); ++j) {
        int v = v_list_[p][j];
        int p = vert_part_id_[v];
        Mat3 sub_k = Mat3::Zero();
        // all cubature tet incident on this v
        for (int k = 0; k < int(v_cubature_[p][j].size()); ++k) {
          VVCubature& cubature = v_cubature_[p][j][k];
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
        reduced_k += vert_basis_transpose_[v]
                     * (part_rotation_transpose_[p] * sub_k * part_rotation_[p])
                     * vert_basis_[v];
      }
    }
    profiler.End("mul p");
    profiler.Start("mul e");
    OMP_FOR
    for (int i = 0; i < part_num_ + int(subspace_subspace_interface_.size()); ++i) {
      int e = i;
      int p0, p1;
      double* mat = NULL;
      if (e < part_num_) {
        if (!is_subspace_domain_[e]) continue; // ignor full domain
        p0 = p1 = e;
        mat = chol_solver_->Av[domain_index_[e]].A;
      } else {
        p0 = subspace_subspace_interface_[e - part_num_].first;
        p1 = subspace_subspace_interface_[e - part_num_].second;
        mat = chol_solver_->Ae[e - part_num_].A;
        e = part_num_ + topology_[p0][p1];
      }

      Eigen::Map<Mat> reduced_k(mat, part_basis_size_[p0], part_basis_size_[p1]);
      // list of vert-vert pair in this domain
      for (int i = 0; i < int(vv_list_[e].size()); ++i) {
        int vi = vv_list_[e][i].first;
        int vj = vv_list_[e][i].second;
        // list of cubature tets that contains this vert-vert pair
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

    OMP_FOR
    for (int j = 0; j < subspace_domain_num_; ++j) {
      int p = subspace_domains_[j];
      MapMat result(chol_solver_->Av[j].A, part_basis_size_[p], part_basis_size_[p]);
      result *= dt_2 * internal_force_scaling_factor;
      for (int i = 0; i < part_basis_size_[p]; ++i) {
        result(i, i) += 1;
      }
    }

    OMP_FOR
    for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
      int p0 = subspace_subspace_interface_[e].first;
      int p1 = subspace_subspace_interface_[e].second;
      MapMat result(chol_solver_->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1]);
      result *= (dt_2 * internal_force_scaling_factor);
    }
  }
  profiler.End("local k");
  if (0) {
    auto A = [&](double * a, double * b) {
      chol_solver_->Multiply(a, b);
    };
    dj::WriteImplicitMatrixToMatlab<double>("/tmp/truth", A, subspace_domain_basis_size_);
    exit(0);
  }
  // full chol
  profiler.Start("full chol");
  {
    full_chol_solver_->ResetTopology();
    full_chol_solver_->SetMatrixZero();

#if 0
    profiler.Start("int. rigid-local");
    OMP_FOR
    for (int e = 0; e < int(interface_domains_.size()); ++e) {
      //      int p0 = interface_domains_[e][1];
      //      int p1 = interface_domains_[e][2];
      if (interface_domains_[e][0] > 2) continue;
      //      P(domain_index_[p0]);
      //      P(domain_index_[p1]);
      //      int solver_edge_idx = full_chol_solver_->GetEdgeIdx(domain_index_[p0], domain_index_[p1]);
      //      KK;
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        const double factor = weight;// * dt_2;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          if (!is_subspace_domain_[pi]) continue;
          int idxi = domain_index_[pi];
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (!is_subspace_domain_[pj]) continue;
            if (pi >= pj) continue;
            //            printf("%06d\t%06d\t%06d\n", vi, vi, t);
            int idxj = domain_index_[pj];
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            Mat3 skew_symmetrix_mat_j = GetSkewSymmetrixMatrix(rj);
            MapMat full_mat = (pi == pj) ? MapMat(full_chol_solver_->Av[idxi].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6)
                              : MapMat(full_chol_solver_->Ae[full_chol_solver_->GetEdgeIdx(idxi, idxj)].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6);
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
    //    if (0)
    for (int e = 0; e < int(interface_domains_.size()); ++e) {
      //      int p0 = interface_domains_[e][1];
      //      int p1 = interface_domains_[e][2];
      if (interface_domains_[e][0] > 2) continue;
      //      if (!is_subspace_domain_[p0] || !is_subspace_domain_[p1]) continue;
      //      int solver_edge_idx = full_chol_solver_->GetEdgeIdx(domain_index_[p0], domain_index_[p1]);
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int t = cubature.first;
        double weight = cubature.second;
        const double factor = weight;// * dt_2;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          if (!is_subspace_domain_[pi]) continue;
          int idxi = domain_index_[pi];
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (!is_subspace_domain_[pj]) continue;
            if (pi != pj) continue;
            //              printf("%06d\t%06d\t%06d\n", vi, vj, t);
            int idxj = domain_index_[pj];
            MapMat full_mat = (pi == pj) ? MapMat(full_chol_solver_->Av[idxi].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6)
                              : MapMat(full_chol_solver_->Ae[full_chol_solver_->GetEdgeIdx(idxi, idxj)].A, part_basis_size_[pi] + 6, part_basis_size_[pj] + 6);

            Mat tmp_mat = (factor * element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj]) * vert_basis_[vj];
            full_mat.block(part_basis_size_[pi] + 0, 0, 3, part_basis_size_[pj]) += tmp_mat;
            full_mat.block(part_basis_size_[pi] + 3, 0, 3, part_basis_size_[pj]) += skew_symmetrix_mat_i * tmp_mat;
          }
        }
      }
    }
    profiler.End("domain rigid-local");
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#else
    //    if (0)
    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      MapMat full_mat =  MapMat(full_chol_solver_->Av[i].A, part_basis_size_[p] + 6, part_basis_size_[p] + 6);
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
          //            printf("%06d\t%06d\t%06d\n", v, v, t);
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

    OMP_FOR
    for (int idx = 0; idx < subspace_domain_num_ + full_chol_solver_->e_number0; ++idx) {
      int e = idx;
      int p0, p1;
      double* mat = nullptr;
      if (e < subspace_domain_num_) {
        p0 = p1 = subspace_domains_[e];
        mat = full_chol_solver_->Av[e].A;
        e = p0;
      } else {
        int d0 = full_chol_solver_->E[(e - subspace_domain_num_) * 2 + 0];
        int d1 = full_chol_solver_->E[(e - subspace_domain_num_) * 2 + 1];
        p0 = subspace_domains_[d0];
        p1 = subspace_domains_[d1];
        mat = full_chol_solver_->Ae[e - subspace_domain_num_].A;
        e = part_num_ + topology_[p0][p1];
        //            exit(0);
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
          //          if (p0 == p1) {
          //            printf("%06d\t%06d\t%06d\n", vi, vj, t);
          //            printf("%06d\t%06d\t%06d\n", vj, vi, t);
          //          }
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
          //          continue;
          const int p = p0;
          // dvi / dvj
          Mat tmp_mat = (sub_k * part_rotation_[p]) * vert_basis_[vj];
          full_mat.block(part_basis_size_[p] + 0, 0, 3, part_basis_size_[p]) += tmp_mat;
          full_mat.block(part_basis_size_[p] + 3, 0, 3, part_basis_size_[p]) += skew_symmetric_mat_i * tmp_mat;

          tmp_mat = (sub_k.transpose() * part_rotation_[p]) * vert_basis_[vi];
          full_mat.block(part_basis_size_[p] + 0, 0, 3, part_basis_size_[p]) += tmp_mat;
          full_mat.block(part_basis_size_[p] + 3, 0, 3, part_basis_size_[p]) += skew_symmetric_mat_j * tmp_mat;
        } else {
          //          continue;
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

    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      const int basis_size = part_basis_size_[p];
      MapMat mat(chol_solver_->Av[i].A, basis_size, basis_size);
      MapMat full_mat(full_chol_solver_->Av[i].A, basis_size + 6, basis_size + 6);
      full_mat.block(part_basis_size_[p], 0, 6, part_basis_size_[p]) *= dt_2;
      full_mat.block(0, part_basis_size_[p], part_basis_size_[p], 6) = full_mat.block(part_basis_size_[p], 0, 6, part_basis_size_[p]).transpose();
      full_mat.block(0, 0, basis_size, basis_size) = mat;
      full_mat.block<6, 6>(basis_size, basis_size) = rigid_matrix.block<6, 6>(i * 6, i * 6);
    }
    OMP_FOR
    for (int e = 0; e < full_chol_solver_->e_number0; ++e) {
      int i0 = full_chol_solver_->E[e * 2 + 0];
      int i1 = full_chol_solver_->E[e * 2 + 1];
      int p0 = subspace_domains_[i0];
      int p1 = subspace_domains_[i1];
      int basis_size0 = part_basis_size_[p0];
      int basis_size1 = part_basis_size_[p1];
      MapMat mat(chol_solver_->Ae[e].A, basis_size0, basis_size1);
      MapMat full_mat(full_chol_solver_->Ae[e].A, basis_size0 + 6, basis_size1 + 6);

      full_mat.block(0, part_basis_size_[p1], part_basis_size_[p0], 6) *= dt_2;
      full_mat.block(part_basis_size_[p0], 0, 6, part_basis_size_[p1]) *= dt_2;
      full_mat.block(0, 0, basis_size0, basis_size1) = mat;
      full_mat.block<6, 6>(basis_size0, basis_size1) = rigid_matrix.block<6, 6>(i0 * 6, i1 * 6);
    }

  }
  profiler.End("full chol");
  if (0) {
    auto A = [&](double * x, double * rhs) {
      full_chol_solver_->Multiply(x, rhs);
    };
    P(full_chol_solver_->size_);
    dj::WriteImplicitMatrixToMatlab<double>("/tmp/test", A, kFullOffset);
    //    dj::WriteImplicitMatrixToMatlab<double>("/tmp/truth", A, kFullOffset);
    exit(0);
  }

  auto K = [&](double * x, double * rhs) {
    memset(rhs, 0, sizeof(double) * (total_basis_num_ + subspace_domain_num_ * 6));
    //    memcpy(rhs, x, sizeof(double) * (total_basis_num_ + subspace_domain_num_ * 6));
    full_chol_solver_->Multiply(x, rhs);
    //        return;

    //    if (0)
    for (int t : full_subspace_interface_tet_) {
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      int* verts = tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int pi = vert_part_id_[vi];
        if (is_subspace_domain_[pi]) {
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          int idxi = domain_offset_[pi];
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (is_subspace_domain_[pj]) {
              continue;
            } else {
              //              continue;
              int idxj = kFullOffset + full_vert_idx_[vj] * 3;
              Vec3 implicit_force = element_k.block<3, 3>(i * 3, j * 3) * MapVec3(x + idxj) * dt_2;
              // momentum
              MapVec3(rhs + rigid_offset_[pi] + 0) += implicit_force;
              // torque
              MapVec3(rhs + rigid_offset_[pi] + 3) += ri.cross(implicit_force);
              //               local deformation
              MapVec(rhs + idxi, part_basis_size_[pi]) += vert_basis_transpose_[vi] * part_rotation_transpose_[pi] * implicit_force;
            }
          }
        } else {
          int idxi = kFullOffset + full_vert_idx_[vi] * 3;
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int pj = vert_part_id_[vj];
            if (is_subspace_domain_[pj])  {
              //              continue;
              int subspace_rhs = domain_offset_[pj];
              int rigid_idx = rigid_offset_[pj];
              Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
              MapVec3(rhs + idxi) += element_k.block<3, 3>(i * 3, j * 3) * (MapVec3(x + rigid_idx + 0) + MapVec3(x + rigid_idx + 3).cross(rj)) * dt_2;
              //              MapVec3(rhs + idxi) += element_k.block<3, 3>(i * 3, j * 3) * (MapVec3(x + rigid_idx + 3).cross(rj)) * dt_2;
              //              MapVec3(rhs + idxi) += element_k.block<3, 3>(i * 3, j * 3) * (MapVec3(x + rigid_idx + 0)) * dt_2;
              MapVec subspace_x(x + subspace_rhs, part_basis_size_[pj]);
              MapVec3(rhs + idxi) += element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] * subspace_x * dt_2;
            } else {
              int idxj = kFullOffset + full_vert_idx_[vj] * 3;
              MapVec3(rhs + idxi) += element_k.block<3, 3>(i * 3, j * 3) * MapVec3(x + idxj) * dt_2;
            }
          }
        }
      }
    }
    for (int t : full_tet_) {
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      int* verts = tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        int idxi = kFullOffset + full_vert_idx_[vi] * 3;
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j];
          int idxj = kFullOffset + full_vert_idx_[vj] * 3;
          MapVec3(rhs + idxi) += element_k.block<3, 3>(i * 3, j * 3) * MapVec3(x + idxj) * dt_2;
        }
      }
    }
    for (int i = 0; i < int(full_verts_.size()); ++i) {
      int v = full_verts_[i];
      int offset = kFullOffset + i * 3;
      MapVec3(rhs + offset) += mass_[v] * MapVec3(x + offset);
    }
  };
  if (0) {
    L("/tmp/truth");
    dj::WriteImplicitMatrixToMatlab<double>("/tmp/truth", K, total_basis_num_ + subspace_domain_num_ * 6);
    KK;
    //    dj::WriteVectorToMatlab(rhs_all.size(), &rhs_all[0], "/tmp/rhs");
    exit(0);
  }
  //  memset(&rhs_all[kFullOffset], 0, sizeof(double) * full_verts_.size() * 3);
  Vec new_vel = rhs_all;
  profiler.Start("solver");
#if 0
  cg_solver_->Resize(total_basis_num_ + subspace_domain_num_ * 6);
  auto info = cg_solver_->Solve(&rhs_all[0], &new_vel[0], K, 5000, 1e-8);
  if (0) {

    dj::WriteVectorToMatlab(rhs_all.size(), &new_vel[0], "/tmp/x");
    //    exit(0);
  }
  P(info.first, info.second);
#else
  sparse_k_->UpdateMatrix(dt);
  //  for (int i = 0; i < sparse_k_->rowLength[0]; ++i) {
  //    P(i, sparse_k_->columnIndices[0][i])
  //  }
  //  exit(0);
  profiler.Start("solve");
  if (fixed_domains_.size() == 0) {
    sparse_k_->Solve(&rhs_all[0], &new_vel[0]);
    //    auto A = [&](double * x, double * rhs) {
    //      sparse_k_->MultiplyVector(x, rhs);
    //    };
    //    L("/tmp/truth");
    //    dj::WriteImplicitMatrixToMatlab<double>("/tmp/truth", A, constrained_dof_ + int(full_verts_.size()) * 3);
    //    exit(0);
  } else {
    //    auto A = [&](double * x, double * rhs) {
    //      sparse_k_->MultiplyVector(x, rhs);
    //    };
    //    P(sparse_k_->columnIndices[125][0]);
    //    P(sparse_k_->rowLength[20]);
    //    for (int i = 0; i < sparse_k_->rowLength[20]; ++i) {
    //      P(i, sparse_k_->columnIndices[20][i], sparse_k_->columnEntries[20][i]);
    //    }
    //    for (int i = 0; i < sparse_k_->rowLength[0]; ++i) {
    //      P(i, sparse_k_->columnIndices[0][i], sparse_k_->columnEntries[0][i]);
    //    }
    //    Mat full = Mat::Zero(sparse_k_->Getn(), sparse_k_->Getn());
    //    for (int r = 0; r < sparse_k_->Getn(); ++r) {
    //      for (int c = 0; c < sparse_k_->rowLength[r]; ++c) {
    //        sparse_k_->columnEntries[r][c] = 1.0;
    //        full(r, sparse_k_->columnIndices[r][c]) = 1.0;
    //      }
    //    }
    //    int cnt = 0;
    //    for (int r = 0; r < full.rows(); ++r) {
    //      for (int c = 0; c < full.cols(); ++c) {
    //        if (full(r, c) != full(c, r)) {
    //          cnt++;
    //          P(r, c);
    //        }
    //      }
    //    }
    //    P(cnt);
    //    dj::WriteEigenMatrixToMatlab(full, "/tmp/full");
    //        sparse_k_->SaveToMatlabFormat("/tmp/save");
    //    exit(0);
    //    L("write /tmp/test");
    //    dj::WriteImplicitMatrixToMatlab<double>("/tmp/test", A, constrained_dof_ + int(full_verts_.size()) * 3);
    //    KK;
    //    exit(0);
    Vec constrained_rhs(constrained_dof_ + int(full_verts_.size()) * 3);
    MapVec non_constrained_rhs(&rhs_all[0], non_constrained_dof_);
    MapVec map_constrained_rhs(&constrained_rhs[0], constrained_dof_);
    Contract(constrained_dof_2_non_constrained_dof_, non_constrained_rhs, map_constrained_rhs);
    memcpy(&constrained_rhs[constrained_dof_], &rhs_all[non_constrained_dof_], sizeof(double) * int(full_verts_.size()) * 3);
    Vec constrained_x(constrained_dof_ + int(full_verts_.size()) * 3);
    sparse_k_->Solve(&constrained_rhs[0], &constrained_x[0]);
    MapVec map_constrained_x(&constrained_x[0], constrained_dof_);
    MapVec non_constrained_x(&new_vel[0], non_constrained_dof_);
    Expand(constrained_dof_2_non_constrained_dof_, map_constrained_x, non_constrained_x);
    memcpy(&new_vel[non_constrained_dof_], &constrained_x[constrained_dof_], sizeof(double) * int(full_verts_.size()) * 3);
  }
  profiler.End("solve");
#endif
  profiler.End("solver");
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      memcpy(&part_vel_[p][0], &new_vel[domain_offset_[p]], sizeof(double) * part_basis_size_[p]);
      translational_vel_[p] = MapVec3(&new_vel[rigid_offset_[p] + 0], 3);
      angular_vel_[p] = MapVec3(&new_vel[rigid_offset_[p] + 3], 3);
    } else {
      continue;
    }
  }
  for (int i = 0; i < int(full_verts_.size()); ++i) {
    int global_v = full_verts_[i];
    int offset = kFullOffset + i * 3;
    MapVec3(velocity_ + global_v * 3) = MapVec3(&new_vel[offset]);
  }
#endif

  profiler.Start("update");
  {
    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      //      memcpy(&part_vel_[p][0], &vel_q_[basis_offset_[p]], sizeof(double) * part_basis_size_[p]);
      //      memcpy(&part_vel_[p][0], &new_vel[basis_offset_[p] + p * 6], sizeof(double) * part_basis_size_[p]);
      MapVec map_q(&part_q_[p][0], part_basis_size_[p]);
      MapVec map_vel(&part_vel_[p][0], part_basis_size_[p]);
      map_q += map_vel * dt;
    }

    OMP_FOR
    for (int i = 0; i < int(full_verts_.size()); ++i) {
      int v = full_verts_[i];
      MapVec3 map_vel(velocity_ + v * 3);
      MapVec3 map_x(X + v * 3);
      map_x += map_vel * dt;
      //      map_vel *= 0.98;
    }

    for (int  i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      center_of_mass_[p] += translational_vel_[p] * dt;
      quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
      quaternion_[p].Normalize();
      //      translational_vel_[p] *= 0.98;
      //      angular_vel_[p] *= 0.98;

      quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
      part_rotation_transpose_[p] = part_rotation_[p].transpose();

      current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
      current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();

      Mat3& rotation = part_rotation_[p];
      MapVec sub_q(&part_q_[p][0], part_basis_size_[p]);
      OMP_FOR
      for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
        int v = vert_local_id2global_id_[p][local_v];
        Vec3 local_u = vert_basis_[v] * sub_q;
        MapVec3 map_x(X + v * 3);
        map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
      }
    }
    // compute global position and update vertex offset
    inv_fem_->UpdateOffset();
  }
  position_changed_ = true;
  profiler.End("update");
  profiler.End("simulation");
}

void MixedMultiDomainTet::SimulateCoupledFullRigid(double dt) {
  // solve rigid motion for reduced domains
  profiler.Start("total");
  const double dt_2 = dt * dt;
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  internal_force_scaling_factor = 1.0;
  inv_fem_->inv_fem_force_model_->ComputePartialForceAndStiffnessMatrixWithSplitedMesh(&all_cubature_tet_[0],
                                                                                       int(all_cubature_tet_.size()),
                                                                                       tet_id_,
                                                                                       tet_, X);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external collision force in world space to each vertex
  ext_force_.clear();
  profiler.Start("collision force");
  if (1) {
    const double kFloor = -0.20;
    const double kFloorStiffness = 25000;
    for (int v = 0; v < vertex_num_; ++v) {
      // Ground collision
      if (X[v * 3 + 1] < kFloor) {
        Vec3 force(0, 0, 0);
        force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
        //      force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
        ext_force_.emplace_back(v, force);
      }
    }
  }
  profiler.End("collision force");
  // user interaction force
  {
#if 1
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      ext_force_.emplace_back(v, force);
    }
#endif
  }

  Vec new_rigid_velocity(int(full_verts_.size()) * 3 + subspace_domain_num_ * 6);
  {
    double rigid_coupling_scaling = conf.Get<double>("rigid coupling scaling");
    // FIXME
    rigid_coupling_scaling = 1.0;
    //  double rigid_coupling_scaling = 1;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    std::vector<Vec3> net_force(subspace_domain_num_, Vec3(0, 0, 0));
    std::vector<Vec3> net_torque(subspace_domain_num_, Vec3(0, 0, 0));
    std::vector<Vec> subspace_interface_force(subspace_domain_num_);
    std::vector<Vec3> rigid_rhs(subspace_domain_num_ * 2 + full_verts_.size(), Vec3(0, 0, 0));
    // gravity on each part
    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      net_force[i] += mass_per_part_[p] * gravity_;
      subspace_interface_force[i] = Vec::Zero(part_basis_size_[p]);
    }
    // external force
    {
      for (ExtForce f : ext_force_) {
        const int& v = f.first;
        const Vec3& force = f.second;
        int p = vert_part_id_[v];
        if (is_subspace_domain_[p]) {
          int idx = domain_index_[p];
          net_force[idx][1] += force[1];
          MapVec3 map_x(X + v * 3);
          Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
          net_torque[idx] += r.cross(force);
        } else {
          int idx = full_vert_idx_[v];
          rigid_rhs[idx] += force;
        }
      }
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // elastic force from boundary tets
    profiler.Start("internal force");
    {
      // boundary force from subspace-subpsace interface
      for (int e = 0; e < int(interface_domains_.size()); ++e) {
        int p0 = interface_domains_[e][1];
        int p1 = interface_domains_[e][2];
        if (interface_domains_[e][0] > 2 || !is_subspace_domain_[p0] || !is_subspace_domain_[p1]) continue;
        // FIXME: ignor interface forces that come from interface of more than 2 domains
        Vec force0 = Vec::Zero(part_basis_size_[p0]);
        Vec force1 = Vec::Zero(part_basis_size_[p1]);
        for (CubaturePoint & cubature : interface_cubature_[e]) {
          int& t = cubature.first;
          double& weight = cubature.second;

          double* element_force = inv_fem_->element_force_ + t * 3 * 4;
          for (int local_v = 0; local_v < 4; ++local_v) {
            MapVec3 force(element_force + local_v * 3);
            int v = tet_[t * 4 + local_v];
            int p = vert_part_id_[v];
            // negative sign before weight is because internal force from vega is in oppisite direction of actual force
            if (p == p0) {
              force0 += vert_basis_transpose_[v] * (-weight * (part_rotation_transpose_[p] * force));
            } else {
              force1 += vert_basis_transpose_[v] * (-weight * (part_rotation_transpose_[p] * force));
            }
          }
        }
        ProjectInterfaceSubspaceForce(e, p0, p1, force0, force1);
        subspace_interface_force[domain_index_[p0]] += force0 * rigid_coupling_scaling;
        subspace_interface_force[domain_index_[p1]] += force1 * rigid_coupling_scaling;
      }

      OMP_FOR
      for (int i = 0; i < subspace_domain_num_; ++i) {
        int p = subspace_domains_[i];
        Mat3& rotation = part_rotation_[p];
        net_force[i] += rotation * momentum_matrix_[p] * subspace_interface_force[i];
        net_torque[i] += rotation * torque_matrix_[p] * subspace_interface_force[i];
      }

      // force from full-subspace interface
      if (1) {
        for (int i = 0; i < int(full_subspace_interface_tet_.size()); ++i) {
          int t = full_subspace_interface_tet_[i];
          int* vert = tet_ + t * 4;
          double* element_force = inv_fem_->element_force_ + t * 3 * 4;
          for (int i = 0; i < 4; ++i) {
            int v = vert[i];
            int p = vert_part_id_[v];
            if (is_subspace_domain_[p]) {
              Vec3 map_force(element_force + i * 3);
              map_force *= rigid_coupling_scaling;
              Vec3 r = MapVec3(X + v * 3) - center_of_mass_[p];
              // vega force is in negative direction
              net_force[domain_index_[p]] -= map_force;
              net_torque[domain_index_[p]] -= r.cross(map_force);
            } else {
              MapVec3 map_force(element_force + i * 3);
              int idx = full_vert_idx_[v];
              rigid_rhs[idx] -= map_force;
            }
          }
        }
      }
      // elastic force from tets within full domain
      if (1) {
        for (int t : full_tet_) {
          int* vert = tet_ + t * 4;
          double* element_force = inv_fem_->element_force_ + t * 3 * 4;
          for (int i = 0; i < 4; ++i) {
            int v = vert[i];
            MapVec3 map_force(element_force + i * 3);
            int idx = full_vert_idx_[v];
            rigid_rhs[idx] -= map_force;
          }
        }
      }
    }
    profiler.End("internal force");

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // multiply dt and add momentum terms
    OMP_FOR
    for (int idx = 0; idx < int(full_verts_.size()); ++idx) {
      int v = full_verts_[idx];
      rigid_rhs[idx] += gravity_ * mass_[v];
      rigid_rhs[idx] *= dt;
      rigid_rhs[idx] += MapVec3(velocity_ + v * 3) * mass_[v];
      //      rigid_rhs[idx].setZero();
    }
    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      int offset = int(full_verts_.size());
      rigid_rhs[offset + i * 2 + 0] = net_force[i] * dt + translational_vel_[p] * mass_per_part_[p];
      rigid_rhs[offset + i * 2 + 1] = net_torque[i] * dt + current_inertia_tensor_[p] * angular_vel_[p];
      //            rigid_rhs[offset + i * 2 + 0].setZero();
      //            rigid_rhs[offset + i * 2 + 1].setZero();
    }

    profiler.Start("rigid k");
    MatCol Rigid = MatCol::Zero(subspace_domain_num_ * 6, subspace_domain_num_ * 6);
    // Compute Stiffness matrix for solving rigid motion
    {
      // Off-diagonal terms from subspace-subspace interface
      {
        OMP_FOR
        for (int e = 0; e < int(interface_domains_.size()); ++e) {
          int p0 = interface_domains_[e][1];
          int p1 = interface_domains_[e][2];
          if (interface_domains_[e][0] > 2 || !is_subspace_domain_[p0] || !is_subspace_domain_[p1]) continue;

          Mat df_dt[2][2] = {
            {Mat::Zero(part_basis_size_[p0], 3), Mat::Zero(part_basis_size_[p0], 3)},
            {Mat::Zero(part_basis_size_[p1], 3), Mat::Zero(part_basis_size_[p1], 3)}
          }; // translation

          Mat df_dr[2][2] = {
            {Mat::Zero(part_basis_size_[p0], 3), Mat::Zero(part_basis_size_[p0], 3)},
            {Mat::Zero(part_basis_size_[p1], 3), Mat::Zero(part_basis_size_[p1], 3)}
          }; // rotation

          for (CubaturePoint & cubature : interface_cubature_[e]) {
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
                Vec3 r = -part_rotation_[pj] * vert_offset_from_mass_center_[vj];
                Mat3 skew_symmetrix_mat = GetSkewSymmetrixMatrix(r);
                int idx1 = (pj == p0) ? 0 : 1;
                df_dt[idx0][idx1] += vert_basis_transpose_[vi] *
                                     (weight * part_rotation_transpose_[pi] * element_k.block<3, 3>(i * 3, j * 3));
                df_dr[idx0][idx1] += vert_basis_transpose_[vi] *
                                     (weight * part_rotation_transpose_[pi] * (element_k.block<3, 3>(i * 3, j * 3) * skew_symmetrix_mat));

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

          tmp_interface_rigid_k_[e * 8 + 0] = part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][0];
          tmp_interface_rigid_k_[e * 8 + 1] = part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][0];
          tmp_interface_rigid_k_[e * 8 + 2] = part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][0];
          tmp_interface_rigid_k_[e * 8 + 3] = part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][0];
          tmp_interface_rigid_k_[e * 8 + 4] = part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][1];
          tmp_interface_rigid_k_[e * 8 + 5] = part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][1];
          tmp_interface_rigid_k_[e * 8 + 6] = part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][1];
          tmp_interface_rigid_k_[e * 8 + 7] = part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][1];


          int idx0 = domain_index_[p0] * 6;
          int idx1 = domain_index_[p1] * 6;
          Rigid.block<3, 3>(idx0 + 0, idx1 + 0) += part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][1];
          Rigid.block<3, 3>(idx0 + 0, idx1 + 3) += part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][1];
          Rigid.block<3, 3>(idx0 + 3, idx1 + 0) += part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][1];
          Rigid.block<3, 3>(idx0 + 3, idx1 + 3) += part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][1];
          Rigid.block<3, 3>(idx1 + 0, idx0 + 0) += part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][0];
          Rigid.block<3, 3>(idx1 + 0, idx0 + 3) += part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][0];
          Rigid.block<3, 3>(idx1 + 3, idx0 + 0) += part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][0];
          Rigid.block<3, 3>(idx1 + 3, idx0 + 3) += part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][0];
        }

        OMP_FOR
        for (int n = 0; n < subspace_domain_num_; ++n) {
          int p = subspace_domains_[n];
          int idx = domain_index_[p] * 6;
          for (std::pair<int, int>& edge_idx : domain_incident_interface_[p]) {
            int& e = edge_idx.first;
            int i = e * 8 + edge_idx.second * 4;
            Rigid.block<3, 3>(idx + 0, idx + 0) += tmp_interface_rigid_k_[i + 0];
            Rigid.block<3, 3>(idx + 0, idx + 3) += tmp_interface_rigid_k_[i + 1];
            Rigid.block<3, 3>(idx + 3, idx + 0) += tmp_interface_rigid_k_[i + 2];
            Rigid.block<3, 3>(idx + 3, idx + 3) += tmp_interface_rigid_k_[i + 3];
          }
        }
      }

      // Off-diagonal terms from subspace-full interface
      if (1) {
        for (int i = 0; i < int(full_subspace_interface_tet_.size()); ++i) {
          int t = full_subspace_interface_tet_[i];
          int* vert = tet_ + t * 4;
          Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
          for (int i4 = 0; i4 < 4; ++i4) {
            int vi = vert[i4];
            int pi = vert_part_id_[vi];
            if (!is_subspace_domain_[pi]) continue;
            int idx0 = domain_index_[pi] * 6;
            Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
            Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
            for (int j4 = 0; j4 < 4; ++j4) {
              int vj = vert[j4];
              int pj = vert_part_id_[vj];
              if (!is_subspace_domain_[pj]) continue;
              int idx1 = domain_index_[pj] * 6;
              Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
              Mat3 skew_symmetrix_mat_j = GetSkewSymmetrixMatrix(rj);
              Rigid.block<3, 3>(idx0 + 0, idx1 + 0) += element_k.block<3, 3>(i4 * 3, j4 * 3);
              Rigid.block<3, 3>(idx0 + 0, idx1 + 3) -= element_k.block<3, 3>(i4 * 3, j4 * 3) * skew_symmetrix_mat_j;
              Rigid.block<3, 3>(idx0 + 3, idx1 + 0) += skew_symmetrix_mat_i * element_k.block<3, 3>(i4 * 3, j4 * 3);
              Rigid.block<3, 3>(idx0 + 3, idx1 + 3) -= skew_symmetrix_mat_i * element_k.block<3, 3>(i4 * 3, j4 * 3) * skew_symmetrix_mat_j;
            }
          }
        }
      }

      Rigid *= (dt_2) * rigid_coupling_scaling;
      // Diagonal terms related to mass
      for (int i = 0; i < subspace_domain_num_; ++i) {
        int p = subspace_domains_[i];
        int idx = i * 6;
        Rigid(idx + 0, idx + 0) += mass_per_part_[p];
        Rigid(idx + 1, idx + 1) += mass_per_part_[p];
        Rigid(idx + 2, idx + 2) += mass_per_part_[p];
        Rigid.block<3, 3>(idx + 3, idx + 3) += current_inertia_tensor_[p];
      }
    }
    profiler.End("rigid k");

    const int offset = int(full_verts_.size()) * 3;
    // kx = K * x;
    auto K = [&](double * x, double * result) {
      memset(result, 0, sizeof(double) * bi_cg_->size());
      // explicity
      if (0) {
        for (int i = 0; i < int(full_verts_.size()); ++i) {
          int v = full_verts_[i];
          result[i * 3 + 0] = mass_[v] * x[i * 3 + 0];
          result[i * 3 + 1] = mass_[v] * x[i * 3 + 1];
          result[i * 3 + 2] = mass_[v] * x[i * 3 + 2];
        }
        for (int i = 0; i < subspace_domain_num_; ++i) {
          int p = subspace_domains_[i];
          MapVec3(result + offset + i * 6 + 0) = mass_per_part_[p] * MapVec3(x + offset + i * 6 + 0);
          MapVec3(result + offset + i * 6 + 3) = current_inertia_tensor_[p] * MapVec3(x + offset + i * 6 + 3);
          //                    MapVec3(result + offset + i * 6 + 0).setZero();
          //                    MapVec3(result + offset + i * 6 + 3).setZero();
        }
        return;
      }

      for (int t : full_tet_) {
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        int* verts = tet_ + t * 4;
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int idxi = full_vert_idx_[vi];
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j];
            int idxj = full_vert_idx_[vj];
            MapVec3(result + idxi * 3) += element_k.block<3, 3>(i * 3, j * 3) * MapVec3(x + idxj * 3);
          }
        }
      }

      for (int t : full_subspace_interface_tet_) {
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        int* verts = tet_ + t * 4;
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          int pi = vert_part_id_[vi];
          if (is_subspace_domain_[pi]) {
            Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
            int idxi = domain_index_[pi] * 6 + offset;
            for (int j = 0; j < 4; ++j) {
              int vj = verts[j];
              int pj = vert_part_id_[vj];
              if (is_subspace_domain_[pj]) continue;
              int idxj = full_vert_idx_[vj] * 3;
              Vec3 implicit_force = element_k.block<3, 3>(i * 3, j * 3) * MapVec3(x + idxj);
              MapVec3(result + idxi + 0) += implicit_force;
              MapVec3(result + idxi + 3) += ri.cross(implicit_force);
            }
          } else {
            int idxi = full_vert_idx_[vi] * 3;
            for (int j = 0; j < 4; ++j) {
              int vj = verts[j];
              int pj = vert_part_id_[vj];
              if (is_subspace_domain_[pj])  {
                int idxj = domain_index_[pj] * 6 + offset;
                Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
                MapVec3(result + idxi) += element_k.block<3, 3>(i * 3, j * 3) * (MapVec3(x + idxj + 0) + MapVec3(x + idxj + 3).cross(rj));
              } else {
                int idxj = full_vert_idx_[vj] * 3;
                MapVec3(result + idxi) += element_k.block<3, 3>(i * 3, j * 3) * MapVec3(x + idxj);
              }
            }
          }
        }
      }
      for (int i = 0; i < bi_cg_->size(); ++i) {
        result[i] *= dt_2;
        //        result[i] *= 0;
      }
      for (int i = 0; i < int(full_verts_.size()); ++i) {
        int v = full_verts_[i];
        result[i * 3 + 0] += mass_[v] * x[i * 3 + 0];
        result[i * 3 + 1] += mass_[v] * x[i * 3 + 1];
        result[i * 3 + 2] += mass_[v] * x[i * 3 + 2];
      }
      MapVec(result + offset, subspace_domain_num_ * 6) += Rigid * MapVec(x + offset, subspace_domain_num_ * 6);
      //            MapVec(result + offset, subspace_domain_num_ * 6).setZero();
    };

    profiler.Start("rigid solve");
    memcpy(&new_rigid_velocity[0], &rigid_rhs[0][0], sizeof(double) * bi_cg_->size());
    //    P(MapVec(&rigid_rhs[0][0], bi_cg_->size()).norm());
    //    for (int i = 0; i < bi_cg_->size() / 3; ++i) {
    //      P(dj::Vec3d(&rigid_rhs[i][0]));
    //    }
    //        P(bi_cg_->size(), full_verts_.size() * 3 + subspace_domain_num_ * 6);
    auto info = bi_cg_->Solve(&rigid_rhs[0][0], &new_rigid_velocity[0], K, 1000, 1e-10);
    //    new_rigid_velocity.setZero();
    P(info.first, info.second);
    profiler.End("rigid solve");
    //    exit(0);
  }

  profiler.Start("rigid update");
  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);

  // Compute translational acc and angular acc
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    int idx = i * 6 + int(full_verts_.size()) * 3;
    MapVec3 new_translational_vel(&new_rigid_velocity[idx]);
    MapVec3 new_angular_vel(&new_rigid_velocity[idx + 3]);
    //    new_angular_vel.setZero();
    acceleration[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration[p] = (new_angular_vel - angular_vel_[p]) / dt;

    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
  }
  profiler.End("rigid update");

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
#if 1
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // RHS
  // simulate internal dynamics
  Vec subspace_rhs = Vec::Zero(total_basis_num_);
  // collision and ui force
  for (std::pair<int, Vec3>& collision : ext_force_) {
    const int& v = collision.first;
    const Vec3& force = collision.second;
    int p = vert_part_id_[v];
    if (is_subspace_domain_[p]) {
      MapVec map_rhs(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
      map_rhs += vert_basis_transpose_[v] * (part_rotation_transpose_[p] * force);
    }
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  AddFictitiousForceAndGravity(acceleration, angular_acceleration, subspace_rhs);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // internal force
  profiler.Start("internal force");
  {
    for (CubaturePoint & cubature : subspace_cubature_) {
      int& t = cubature.first;
      double& weight = cubature.second;
      double* force = inv_fem_->element_force_ + t * 12;
      for (int i4 = 0; i4 < 4; ++i4) {
        int v = tet_[t * 4 + i4];
        int p = vert_part_id_[v];
        MapVec3 map_force(force + i4 * 3);
        if (is_subspace_domain_[p]) {
          MapVec subspace_force(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
          // negative sign before weight is because internal force from vega is in oppisite direction of actual force
          //          ASSERT(v >= 0 && v < vertex_num_);
          //          ASSERT(subspace_force.size() == vert_basis_transpose_[v].rows(), P(subspace_force.size(), vert_basis_transpose_[v].rows(), v, p));
          subspace_force -= vert_basis_transpose_[v] * (internal_force_scaling_factor * weight * (part_rotation_[p].transpose() * map_force));
          // implicit force
          if (0) {
            Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
            for (int j4 = 0; j4 < 4; ++j4) {
              int vj = tet_[t * 4 + j4];
              int pj = vert_part_id_[vj];
              if (is_subspace_domain_[pj]) continue;
              int idx = full_vert_idx_[vj] * 3;
              subspace_force -= vert_basis_transpose_[v] *
                                part_rotation_[p] *
                                (element_k.block<3, 3>(i4 * 3, j4 * 3) * (MapVec3(&new_rigid_velocity[idx]) * dt));
            }
          }
        }
      }
    }
    //    exit(0);
    subspace_rhs *= dt;
    // add velocity related term
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]) += MapVec(&part_vel_[p][0], part_basis_size_[p]);
    }
  }
  //  exit(0);
  profiler.End("internal force");

  // Compute reduced stiffness matrix
  // assemble stiffness matrix for reduced k
  profiler.Start("local k");
  {
    chol_solver_->SetMatrixZero();
    profiler.Start("mul p");
    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      double* mat = chol_solver_->Av[i].A;
      Eigen::Map<MatRow> reduced_k(mat, part_basis_size_[p], part_basis_size_[p]);
      //    reduced_k.setZero();
      // all cubature vertex in this domain
      for (int j = 0; j < int(v_list_[p].size()); ++j) {
        int v = v_list_[p][j];
        int p = vert_part_id_[v];
        Mat3 sub_k = Mat3::Zero();
        // all cubature tet incident on this v
        for (int k = 0; k < int(v_cubature_[p][j].size()); ++k) {
          VVCubature& cubature = v_cubature_[p][j][k];
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
        reduced_k += vert_basis_transpose_[v]
                     * (part_rotation_transpose_[p] * sub_k * part_rotation_[p])
                     * vert_basis_[v];
      }
    }
    profiler.End("mul p");
    profiler.Start("mul e");
    //    OMP_FOR
    for (int i = 0; i < part_num_ + int(subspace_subspace_interface_.size()); ++i) {
      int e = i;
      int p0, p1;
      double* mat = NULL;
      if (e < part_num_) {
        if (!is_subspace_domain_[e]) continue; // ignor full domain
        p0 = p1 = e;
        mat = chol_solver_->Av[domain_index_[e]].A;
      } else {
        p0 = subspace_subspace_interface_[e - part_num_].first;
        p1 = subspace_subspace_interface_[e - part_num_].second;
        // TODO delete
        {
          int test_p0 = subspace_domains_[chol_solver_->E[(e - part_num_) * 2 + 0]];
          int test_p1 = subspace_domains_[chol_solver_->E[(e - part_num_) * 2 + 1]];
          ASSERT(test_p0 == p0, P(p0, p1, e, part_num_) P(test_p0, p0, chol_solver_->E[(e - part_num_) * 2 + 0], chol_solver_->E[(e - part_num_) * 2 + 1]) );
          ASSERT(test_p1 == p1, P(test_p1, p1));
        }
        mat = chol_solver_->Ae[e - part_num_].A;
        e = part_num_ + topology_[p0][p1];
      }

      Eigen::Map<Mat> reduced_k(mat, part_basis_size_[p0], part_basis_size_[p1]);
      // list of vert-vert pair in this domain
      for (int i = 0; i < int(vv_list_[e].size()); ++i) {
        int vi = vv_list_[e][i].first;
        int vj = vv_list_[e][i].second;
        // list of cubature tets that contains this vert-vert pair
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

    OMP_FOR
    for (int j = 0; j < subspace_domain_num_; ++j) {
      int p = subspace_domains_[j];
      MapMat result(chol_solver_->Av[j].A, part_basis_size_[p], part_basis_size_[p]);
      result *= dt_2 * internal_force_scaling_factor;
      for (int i = 0; i < part_basis_size_[p]; ++i) {
        result(i, i) += 1;
      }
    }

    OMP_FOR
    for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
      int p0 = subspace_subspace_interface_[e].first;
      int p1 = subspace_subspace_interface_[e].second;
      MapMat result(chol_solver_->Ae[e].A, part_basis_size_[p0], part_basis_size_[p1]);
      result *= (dt_2 * internal_force_scaling_factor);
    }
  }
  profiler.End("local k");

  profiler.Start("local solver");
  memcpy(&vel_q_[0], &subspace_rhs[0], sizeof(double) * total_basis_num_);
  //  P(vel_q_.size(), total_basis_num_);
  //  P(subspace_rhs.norm());
  //  PVEC(subspace_rhs);
  int code = chol_solver_->Solve(&subspace_rhs[0], &vel_q_[0]);
  ASSERT(code == solver::BLOCK_MATRIX_GRAPH<double>::kSuccess)
  profiler.End("local solver");
#endif // LOCAL_DEFORMATION
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  profiler.Start("update");
  {
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      memcpy(&part_vel_[p][0], &vel_q_[basis_offset_[p]], sizeof(double) * part_basis_size_[p]);
      MapVec map_q(&part_q_[p][0], part_basis_size_[p]);
      MapVec map_vel(&part_vel_[p][0], part_basis_size_[p]);
      map_q += map_vel * dt;
    }

    for (int i = 0; i < int(full_verts_.size()); ++i) {
      int v = full_verts_[i];
      MapVec3 map_vel(velocity_ + v * 3);
      map_vel = MapVec3(&new_rigid_velocity[i * 3]);
      MapVec3 map_x(X + v * 3);
      map_x += map_vel * dt;
      map_vel *= 0.98;
    }

    OMP_FOR
    for (int  i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      center_of_mass_[p] += translational_vel_[p] * dt;
      quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
      quaternion_[p].Normalize();
      translational_vel_[p] *= 0.98;
      angular_vel_[p] *= 0.98;

      quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
      part_rotation_transpose_[p] = part_rotation_[p].transpose();

      current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
      current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();

      Mat3& rotation = part_rotation_[p];
      MapVec sub_q(&part_q_[p][0], part_basis_size_[p]);
      for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
        int v = vert_local_id2global_id_[p][local_v];
        Vec3 local_u = vert_basis_[v] * sub_q;
        MapVec3 map_x(X + v * 3);
        map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
      }
    }

    // compute global position and update vertex offset
    inv_fem_->UpdateOffset();
    UpdateEmbededVertexPosition();
  }
  profiler.End("update");
  profiler.End("total");
}

void MixedMultiDomainTet::SimulateMixed(double dt) {
  // solve rigid motion for reduced domains
  profiler.Start("total");
  const double dt_2 = dt * dt;
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  internal_force_scaling_factor = 1.0;
  inv_fem_->inv_fem_force_model_->ComputePartialForceAndStiffnessMatrixWithSplitedMesh(&all_cubature_tet_[0],
                                                                                       int(all_cubature_tet_.size()),
                                                                                       tet_id_,
                                                                                       tet_, X);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external collision force in world space to each vertex
  ext_force_.clear();
  profiler.Start("collision force");
  if (1) {
    const double kFloor = -0.50;
    const double kFloorStiffness = 25000;
    for (int v = 0; v < vertex_num_; ++v) {
      // Ground collision
      //    if (0)
      if (X[v * 3 + 1] < kFloor) {
        Vec3 force(0, 0, 0);
        force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
        //      force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
        ext_force_.emplace_back(v, force);
      }
    }
  }
  profiler.End("collision force");
  // user interaction force
  {
#if 1
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      ext_force_.emplace_back(v, force);
    }
#endif
  }
  Vec new_rigid_vel;
  SimulationRigidMotion(dt, new_rigid_vel, &ext_force_);
  profiler.Start("rigid update");
  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);
  std::vector<Vec3> old_vel(part_num_);

  // Compute translational acc and angular acc
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    int idx = i * 6;
    MapVec3 new_translational_vel(&new_rigid_vel[idx]);
    MapVec3 new_angular_vel(&new_rigid_vel[idx + 3]);
    //    new_angular_vel.setZero();
    acceleration[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration[p] = (new_angular_vel - angular_vel_[p]) / dt;

    // TODO delete
    old_vel[p] = translational_vel_[p];
    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
  }
  profiler.End("rigid update");

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local deformation
#if 1
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // RHS
  // simulate internal dynamics
  Vec subspace_rhs = Vec::Zero(total_basis_num_);
  // collision and ui force
  for (std::pair<int, Vec3>& collision : ext_force_) {
    const int& v = collision.first;
    const Vec3& force = collision.second;
    int p = vert_part_id_[v];
    if (is_subspace_domain_[p]) {
      MapVec map_rhs(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
      map_rhs += vert_basis_transpose_[v] * (part_rotation_transpose_[p] * force);
    } else {
      int local_v = local_vert_idx_[v];
      MapVec3 map_rhs(&subspace_rhs[basis_offset_[p] + local_v * 3]);
      map_rhs += force;
    }
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  AddFictitiousForceAndGravity(acceleration, angular_acceleration, subspace_rhs);

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // internal force
  Vec3 net_force_full(0, 0, 0);
  profiler.Start("internal force");
  {
    for (CubaturePoint & cubature : all_cubature_) {
      int& t = cubature.first;
      double& weight = cubature.second;
      double* force = inv_fem_->element_force_ + t * 12;
      for (int i4 = 0; i4 < 4; ++i4) {
        int v = tet_[t * 4 + i4];
        int p = vert_part_id_[v];
        MapVec3 map_force(force + i4 * 3);
        if (is_subspace_domain_[p]) {
          MapVec subspace_force(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
          // negative sign before weight is because internal force from vega is in oppisite direction of actual force
          subspace_force -= vert_basis_transpose_[v] * (internal_force_scaling_factor * weight * (part_rotation_[p].transpose() * map_force));
        } else {
          int local_v = local_vert_idx_[v];
          MapVec3 full_force(&subspace_rhs[basis_offset_[p] + local_v * 3]);
          full_force -= map_force;
          net_force_full -= map_force;
        }
      }
    }

    subspace_rhs *= dt;
    // add velocity related term
    for (int p = 0; p < part_num_; ++p) {
      if (is_subspace_domain_[p]) {
        MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]) += MapVec(&part_vel_[p][0], part_basis_size_[p]);
        //            MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]).setZero();
      } else {
        for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
          int global_v = vert_local_id2global_id_[p][local_v];
          int offset = basis_offset_[p] + local_v * 3;
          MapVec3(&subspace_rhs[offset]) += mass_[global_v] * (MapVec3(&part_vel_[p][local_v * 3]) + gravity_ * dt * conf.Get<double>("gravity scaling"));
          net_force_full += gravity_ * mass_[global_v];
          //        MapVec3(&subspace_rhs[offset]).setZero();
        }
      }
    }
  }
  net_force_full *= dt;
  KK
  P(dj::Vec3d(&net_force_full[0]));
  profiler.End("internal force");

  // Compute reduced stiffness matrix
  profiler.Start("mul p");
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    Mat& reduced_k = domain_k_[p];
    reduced_k.setZero();
    // all cubature vertex in this domain
    for (int j = 0; j < int(v_list_[p].size()); ++j) {
      int v = v_list_[p][j];
      int p = vert_part_id_[v];
      Mat3 sub_k = Mat3::Zero();
      // all cubature tet incident on this v
      for (int k = 0; k < int(v_cubature_[p][j].size()); ++k) {
        VVCubature& cubature = v_cubature_[p][j][k];
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
      reduced_k += vert_basis_transpose_[v]
                   * (part_rotation_transpose_[p] * sub_k * part_rotation_[p])
                   * vert_basis_[v];
    }
  }
  profiler.End("mul p");
  profiler.Start("mul e");
  OMP_FOR
  for (int i = 0; i < part_num_ + int(subspace_subspace_interface_.size()); ++i) {
    int e = i;
    int p0, p1;
    double* mat = NULL;
    if (e < part_num_) {
      if (!is_subspace_domain_[e]) continue; // ignor full domain
      p0 = p1 = e;
      mat = domain_k_[p0].data();
    } else {
      p0 = subspace_subspace_interface_[e - part_num_].first;
      p1 = subspace_subspace_interface_[e - part_num_].second;
      mat = interface_k_[e - part_num_].data();
      interface_k_[e - part_num_].setZero();
      e = part_num_ + topology_[p0][p1];
    }
    Eigen::Map<Mat> reduced_k(mat, part_basis_size_[p0], part_basis_size_[p1]);
    // list of vert-vert pair in this domain
    for (int i = 0; i < int(vv_list_[e].size()); ++i) {
      int vi = vv_list_[e][i].first;
      int vj = vv_list_[e][i].second;
      // list of cubature tets that contains this vert-vert pair
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

  // code for testing
  {
    //  for (int i = 0; i < subspace_domain_num_; ++i) {
    //    int p = subspace_domains_[i];
    //    Mat& result = domain_k_[p];
    //    result *= dt_2 * internal_force_scaling_factor;
    //    for (int i = 0; i < part_basis_size_[p]; ++i) {
    //      result(i, i) += 1;
    //    }
    //    P(p, result.norm());
    //  }
    //  exit(0);
    //  for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
    //    int p0 = subspace_subspace_interface_[e].first;
    //    int p1 = subspace_subspace_interface_[e].second;
    //    Mat& result = interface_k_[e];
    //    result *= (dt_2 * internal_force_scaling_factor);
    //        P(p0, p1, e, result.norm());
    //        Mat3 mm = result.block<3, 3>(0, 0);
    //        PMAT(mm);
    //  }
    //  exit(0);
    //    profiler.EndTimer("internal force");
  }

  // full-subspace interface k
  // pth domain
  for (int n = 0; n < subspace_domain_num_; ++n) {
    int p = subspace_domains_[n];
    // ith adjacent full domain
    for (int i = 0; i < int(full_vert_on_domain_[p].size()); ++i) {
      // jth vertex
      for (int j = 0; j < int(full_vert_on_domain_[p][i].second.size()); ++j) {
        Mat& vert_k = full_vert_subspace_matrix_block_[p][i][j];
        vert_k.setZero();
        // kth incident tet
        for (int k = 0; k < int(vert_incident_tet_domain_[p][i][j].size()); ++k) {
          const int t = vert_incident_tet_domain_[p][i][j][k].first;
          const int idx = vert_incident_tet_domain_[p][i][j][k].second;
          Eigen::Map<Eigen::Matrix<double, 12, 12> > tet_k(inv_fem_->element_k_ + t * 144);
          for (int l = 0; l < 4; ++l) {
            int v_subspace = tet_[t * 4  + l];
            if (l == idx || vert_part_id_[v_subspace] != p) continue;
            vert_k += vert_basis_transpose_[v_subspace] * (part_rotation_transpose_[p] * tet_k.block<3, 3>(l * 3, idx * 3));
          }
        }
      }
    }
  }

  auto StiffnessMatrix = [&](double * x, double * result) {
    memset(result, 0, sizeof(double) * total_basis_num_);
    // explicit
    if (0) {
      for (int p = 0; p < part_num_; ++p) {
        if (is_subspace_domain_[p]) {
          MapVec map_result(result + basis_offset_[p], part_basis_size_[p]);
          MapVec map_x(x + basis_offset_[p], part_basis_size_[p]);
          map_result = map_x;
        } else {
          for (int v = 0; v < vert_num_per_part_[p]; ++v) {
            int global_v = vert_local_id2global_id_[p][v];
            double mass = mass_[global_v];
            int offset = basis_offset_[p] + v * 3;
            result[offset + 0] = x[offset + 0] * mass;
            result[offset + 1] = x[offset + 1] * mass;
            result[offset + 2] = x[offset + 2] * mass;
          }
        }
      }
      return;
    }

    // subspace domain
    if (1) {
      for (int i = 0; i < subspace_domain_num_; ++i) {
        int p = subspace_domains_[i];
        MapVec map_result(&result[basis_offset_[p]], part_basis_size_[p]);
        MapVec map_x(&x[basis_offset_[p]], part_basis_size_[p]);
        map_result += domain_k_[p] * map_x;
        //      map_result.setZero();
      }
    }

    // subspace-subspace interface
    if (1) {
      for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
        int p0 = subspace_subspace_interface_[e].first;
        int p1 = subspace_subspace_interface_[e].second;
        {
          MapVec map_result(&result[basis_offset_[p0]], part_basis_size_[p0]);
          MapVec map_x(&x[basis_offset_[p1]], part_basis_size_[p1]);
          map_result += interface_k_[e] * map_x;
        }
        {
          MapVec map_result(&result[basis_offset_[p1]], part_basis_size_[p1]);
          Eigen::Map<Eigen::RowVectorXd> map_x(&x[basis_offset_[p0]], part_basis_size_[p0]);
          map_result += map_x * interface_k_[e];
        }
      }
    }

    // subspace-full inteface
    if (1) {
      for (int i = 0; i < subspace_domain_num_; ++i) {
        int p = subspace_domains_[i];
        // jth full domain
        for (int j = 0; j < int(full_vert_on_domain_[p].size()); ++j) {
          int full_domain = full_vert_on_domain_[p][j].first;
          // kth vertex
          for (int k = 0; k < int(full_vert_on_domain_[p][j].second.size()); ++k) {
            int v = full_vert_on_domain_[p][j].second[k];
            int local_v = local_vert_idx_[v];
            Mat& vert_k = full_vert_subspace_matrix_block_[p][j][k];
            {
              MapVec map_result(&result[basis_offset_[p]], part_basis_size_[p]);
              MapVec3 map_x(&x[basis_offset_[full_domain] + local_v * 3]);
              map_result += vert_k * map_x;
            }
            {
              MapVec3 map_result(&result[basis_offset_[full_domain] + local_v * 3]);
              Eigen::Map<Eigen::RowVectorXd> map_x(&x[basis_offset_[p]], part_basis_size_[p]);
              map_result += map_x * vert_k;
            }
          }
        }
      }
    }

    // tet complete inside full domain
    if (1) {
      for (int t : full_tet_) {
        // TODO delete
        //         if (t == 2) continue;
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > tet_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          const int v0 = verts[i];
          const int p0 = vert_part_id_[v0];
          const int local_v0 = local_vert_idx_[v0];
          for (int j = 0; j < 4; ++j) {
            const int v1 = verts[j];
            const int p1 = vert_part_id_[v1];
            const int local_v1 = local_vert_idx_[v1];
            MapVec3 map_result(result + basis_offset_[p0] + local_v0 * 3);
            MapVec3 map_x(x + basis_offset_[p1] + local_v1 * 3);
            map_result += tet_k.block<3, 3>(i * 3, j * 3) * map_x;
          }
        }
      }
    }
    // tet in full-subspace interface
    if (1) {
      for (int t : full_subspace_interface_tet_) {
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > tet_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          const int v0 = verts[i];
          const int p0 = vert_part_id_[v0];
          const int local_v0 = local_vert_idx_[v0];
          if (is_subspace_domain_[p0]) continue;
          for (int j = 0; j < 4; ++j) {
            const int v1 = verts[j];
            const int p1 = vert_part_id_[v1];
            if (is_subspace_domain_[p1]) continue;
            const int local_v1 = local_vert_idx_[v1];
            MapVec3 map_result(result + basis_offset_[p0] + local_v0 * 3);
            MapVec3 map_x(x + basis_offset_[p1] + local_v1 * 3);
            map_result += tet_k.block<3, 3>(i * 3, j * 3) * map_x;
          }
        }
      }
    }

    // multiply dt * dt and add mass term
    for (int p = 0; p < part_num_; ++p) {
      if (is_subspace_domain_[p]) {
        MapVec map_result(result + basis_offset_[p], part_basis_size_[p]);
        MapVec map_x(x + basis_offset_[p], part_basis_size_[p]);
        map_result *= dt_2;
        map_result += map_x;
        map_result.setZero();
      } else {
        for (int v = 0; v < vert_num_per_part_[p]; ++v) {
          int global_v = vert_local_id2global_id_[p][v];
          double mass = mass_[global_v];
          int offset = basis_offset_[p] + v * 3;
          result[offset + 0] = x[offset + 0] * mass + dt_2 * result[offset + 0];
          result[offset + 1] = x[offset + 1] * mass + dt_2 * result[offset + 1];
          result[offset + 2] = x[offset + 2] * mass + dt_2 * result[offset + 2];
          //          if (is_constrainted_[global_v]) {
          //            result[offset + 0] = 0;
          //            result[offset + 1] = 0;
          //            result[offset + 2] = 0;
          //          }

        }
      }
    }
  };

  // testing code
  if (0) {
    inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
    auto FullStiffnessMatrix = [&](double * x, double * result) {
      memset(result, 0, sizeof(double) * 3 * vertex_num_);
      //    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(x, result);
      //    if (0)
      for (int t = 0; t < tet_number; ++t) {
        int* verts = tet_ + t * 4;
        double (*k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          double* force = result + vi * 3;
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j] * 3;
            force[0] += k[i * 3 + 0][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 0][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 0][j * 3 + 2] * x[vj + 2];
            force[1] += k[i * 3 + 1][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 1][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 1][j * 3 + 2] * x[vj + 2];
            force[2] += k[i * 3 + 2][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 2][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 2][j * 3 + 2] * x[vj + 2];
          }
        }
      }
      for (int v = 0; v < vertex_num_; ++v) {
        result[v * 3 + 0] = mass_[v] * x[v * 3 + 0] + dt_2 * result[v * 3 + 0];
        result[v * 3 + 1] = mass_[v] * x[v * 3 + 1] + dt_2 * result[v * 3 + 1];
        result[v * 3 + 2] = mass_[v] * x[v * 3 + 2] + dt_2 * result[v * 3 + 2];
      }
    };

    int p = 2;
    Vec x = Vec::Zero(total_basis_num_);
    Vec result = Vec::Zero(total_basis_num_);
    double test[3][3];
    for (int e = 0; e < edge_num_; ++e) {
      int global_v = edges_[e * 2];
      int other_v = edges_[e * 2 + 1];
      //      global_v = 4196;
      //      other_v = 12648;
      int p0 = vert_part_id_[global_v];
      int p1 = vert_part_id_[other_v];
      if (p0 != p || p1 != p) continue;
      int v = local_vert_idx_[global_v];
      int local_other = local_vert_idx_[other_v];
      //        P(global_v, other_v);
      //        P(v, local_other);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          x[basis_offset_[p] + v * 3 + c] = 1;
          StiffnessMatrix(&x[0], &result[0]);
          test[r][c] = result[basis_offset_[p] + local_other * 3 + r];
          //          std::cout <<  result[basis_offset_[p] + local_other * 3 + r] << "\t";
          x[basis_offset_[p] + v * 3 + c] = 0;
        }
        //        std::cout << std::endl;
      }

      //        KK;
      double truth[3][3];
      {
        Vec x = Vec::Zero(vertex_num_ * 3);
        Vec rhs = Vec::Zero(vertex_num_ * 3);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            x[global_v * 3 + c] = 1;
            FullStiffnessMatrix(&x[0], &rhs[0]);
            truth[r][c] = rhs[other_v * 3 + r];
            //            std::cout <<  rhs[other_v * 3 + r] << "\t";
            ASSERT(dj::Abs(truth[r][c] - test[r][c]) < 1e-8, P(r, c, truth[r][c], test[r][c]) P(global_v, other_v));
            x[global_v * 3 + c] = 0;
          }
          //          std::cout << std::endl;
        }
      }
      //      break;
    }
    L("full domain k verified");
    exit(0);
  }
  profiler.Start("solver");
  cg_solver_->Resize(total_basis_num_);

  //  for (int v = 0; v < vertex_num_; ++v) {
  //    if (is_constrainted_[v]) {
  //      int p = vert_part_id_[v];
  //      int local_v = local_vert_idx_[v];
  //      int offset = basis_offset_[p] + local_v * 3;
  //      subspace_rhs[offset + 0] = 0;
  //      subspace_rhs[offset + 1] = 0;
  //      subspace_rhs[offset + 2] = 0;
  //    }
  //  }

  //  if (0)
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      MapVec map_result(&subspace_rhs[0] + basis_offset_[p], part_basis_size_[p]);
      map_result.setZero();
    } else {
      //      MapVec3 map_result(&subspace_rhs[0] + basis_offset_[p]);
      //      PVEC(map_result);
    }
  }
  memcpy(tmp_vertex_pos_, &subspace_rhs[0], sizeof(double) * subspace_rhs.size());
  //  memset(tmp_vertex_pos_, 0, sizeof(double) * subspace_rhs.size());
  auto info = cg_solver_->Solve(&subspace_rhs[0], tmp_vertex_pos_, StiffnessMatrix, 2000, 1e-20);
  P(info.first, info.second);
  profiler.End("solver");
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#else
#endif // LOCAL_DEFORMATION

  // Momentum conservation for full-subspace tets
  std::vector<Vec3> correct_force(full_subspace_interface_tet_.size());
  //  memset(&correct_force[0][0], 0, sizeof(Vec3) * full_subspace_interface_tet_.size());
  Vec3 net_f(0, 0, 0);
  Vec3 implicit_f(0, 0, 0);
  Vec3 implicit_f1(0, 0, 0);
  Vec3 net_f1(0, 0, 0);
  Mat K = Mat::Zero(3, 3);
  Vec3 net_imp_force(0, 0, 0);
  if (1) {
    for (int ii = 0; ii < int(full_subspace_interface_tet_.size()); ++ii) {
      int t = full_subspace_interface_tet_[ii];
      int* vert = tet_ + t * 4;
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      correct_force[ii].setZero();
      for (int i = 0; i < 4; ++i) {
        int vi = vert[i];
        int pi = vert_part_id_[vi];
        if (pi == 0) {
          net_f -= MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
        } else {
          net_f1 -= MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
        }
        if (is_subspace_domain_[pi]) {
          for (int j = 0; j < 4; ++j) {
            int vj = vert[j];
            int pj = vert_part_id_[vj];
            if (!is_subspace_domain_[pj]) continue;
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            correct_force[ii] -= element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] - rj.cross(angular_vel_[pj]));
            if (pi == 0) {
              implicit_f -= element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] - rj.cross(angular_vel_[pj])) * dt;
              K += element_k.block<3, 3>(i * 3, j * 3);
            }
          }
        } else {
          for (int j = 0; j < 4; ++j) {
            int vj = vert[j];
            int pj = vert_part_id_[vj];
            if (is_subspace_domain_[pj])  {
              correct_force[ii] -= element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] *
                                   MapVec(&tmp_vertex_pos_[basis_offset_[pj]], part_basis_size_[pj]);
              implicit_f1 -= element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] *
                             MapVec(&tmp_vertex_pos_[basis_offset_[pj]], part_basis_size_[pj]);
            } else {
              int local_v = local_vert_idx_[vj];
              correct_force[ii] -= element_k.block<3, 3>(i * 3, j * 3) * MapVec3(&tmp_vertex_pos_[basis_offset_[pj] + local_v * 3]);
              net_imp_force -= element_k.block<3, 3>(i * 3, j * 3) * MapVec3(&tmp_vertex_pos_[basis_offset_[pj] + local_v * 3]);
              implicit_f1 -= element_k.block<3, 3>(i * 3, j * 3) * MapVec3(&tmp_vertex_pos_[basis_offset_[pj] + local_v * 3]);
            }
          }
        }
      }
      //      correct_force[ii] *= 0.25 * dt * dt;
      correct_force[ii].setZero();
      //      correct_force[ii] *= dt * dt;
    }
    net_imp_force *= dt * dt;
    P(dj::Vec3d(&net_imp_force[0]));
#if 0
    PVEC(net_f1);
    PVEC(net_f);
    net_f *= dt;
    implicit_f1 *= dt * dt;
    implicit_f *= dt;
    PVEC(implicit_f1);
    PVEC(implicit_f);
    PVEC(net_f);
    Vec3 g = gravity_ * mass_per_part_[0] * dt + net_f;
    PVEC(g);
    //    Mat3 diff = K - Rigid.block<3, 3>(0, 0);
    //    PMAT(diff);
    //    implicit_f = K * translational_vel_[0] * dt * -dt;
    Vec3 f1 = gravity_ * mass_[4] * dt + implicit_f1 + net_f1 * dt;
    PVEC(f1);
    Vec3 m4 = mass_[4] * (MapVec3(velocity_ + basis_offset_[1] + local_vert_idx_[4] * 3) -
                          MapVec3(&part_vel_[1][local_vert_idx_[4] * 3]));
    PVEC(m4);
    ASSERT((f1 - m4).norm() < 1e-10);
    Vec3 test_rhs = g + old_vel[0] * mass_per_part_[0];
    Vec3 implicit_f01 = implicit_f + implicit_f1;
    PVEC(implicit_f01);
    PVEC(test_rhs);
    g += implicit_f;
    Vec3 net_int_force = g + f1;
    PVEC(net_int_force);
    PVEC(g);
    Vec3 dv = (translational_vel_[0] - old_vel[0]) * mass_per_part_[0];
    ASSERT((dv - g).norm() < 1e-8);
    Vec3 total_dv = dv + m4;
    PVEC(total_dv);
    PVEC(dv);
    PVEC(correct_force[0]);
    Vec3 correct_mom = correct_force[0] * 4;
    PVEC(correct_mom);
    KK;
    KK;
    //    exit(0);
#endif

    //    PVEC(translational_vel_[0]);
    //    PVEC(angular_vel_[0]);
    if (0) {
      Vec3 total_momentum(0, 0, 0);
      double total_mass = 0;
      for (int p = 0; p < part_num_; ++p) {
        if (is_subspace_domain_[p]) {
          total_mass += mass_per_part_[p];
          total_momentum += mass_per_part_[p] * translational_vel_[p];
        } else {
          for (int v = 0; v < vert_num_per_part_[p]; ++v) {
            int global = vert_local_id2global_id_[p][v];
            //            if (global != 4) continue;
            total_mass += mass_[global];
            //            P(global, dj::Vec3d(velocity_ + basis_offset_[p] + v * 3));
            total_momentum += MapVec3(tmp_vertex_pos_ + basis_offset_[p] + v * 3) * mass_[global];
          }
        }
      }
      static int i = 1;
      Vec3 g_mom = i * gravity_ * total_mass * dt;
      Vec rigid_mom = translational_vel_[0] * mass_per_part_[0];
      PVEC(rigid_mom);
      PVEC(total_momentum);
      PVEC(g_mom);

    }
    //    if (0)
    for (int ii = 0; ii < int(full_subspace_interface_tet_.size()); ++ii) {
      int t = full_subspace_interface_tet_[ii];
      int* vert = tet_ + t * 4;
      Vec3& force = correct_force[ii];
      //            P(t, force.norm());
      for (int i = 0; i < 4; ++i) {
        int vi = vert[i];
        int pi = vert_part_id_[vi];
        if (is_subspace_domain_[pi]) {
          translational_vel_[pi] -= (force) / mass_per_part_[pi];
          Vec3 r = MapVec3(X + vi * 3) - center_of_mass_[pi];
          angular_vel_[pi] -= current_inv_inertia_tensor_[pi] * r.cross(force);
        } else {
          int local_v = local_vert_idx_[vi];
          MapVec3(tmp_vertex_pos_ + basis_offset_[pi] + local_v * 3) -= force / mass_[vi];
        }
      }
    }

  }
  Vec3 total_momentum(0, 0, 0);
  double total_mass = 0;
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      total_mass += mass_per_part_[p];
      total_momentum += mass_per_part_[p] * translational_vel_[p];
    } else {
      for (int v = 0; v < vert_num_per_part_[p]; ++v) {
        int global = vert_local_id2global_id_[p][v];
        //        if (global != 4) continue;
        total_mass += mass_[global];
        //        P(global, dj::Vec3d(velocity_ + basis_offset_[p] + v * 3));
        total_momentum += MapVec3(tmp_vertex_pos_ + basis_offset_[p] + v * 3) * mass_[global];
      }
    }
  }

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  profiler.Start("update");
  {
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      memcpy(&part_vel_[p][0], tmp_vertex_pos_ + basis_offset_[p], sizeof(double) * part_basis_size_[p]);
      MapVec map_q(&part_q_[p][0], part_basis_size_[p]);
      MapVec map_vel(&part_vel_[p][0], part_basis_size_[p]);
      map_q += map_vel * dt;
    }

    static Vec3 pre_momentum(0, 0, 0);
    Vec3 net_momentum(0, 0, 0);
    for (int i = 0; i < int(full_simulation_domains_.size()); ++i) {
      int p = full_simulation_domains_[i];
      memcpy(&part_vel_[p][0], tmp_vertex_pos_ + basis_offset_[p], sizeof(double) * part_basis_size_[p]);
      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        int global_v = vert_local_id2global_id_[p][local_v];
        MapVec3 map_x(X + global_v * 3);
        MapVec3 map_vel(&part_vel_[p][local_v * 3]);
        net_momentum += mass_[global_v] * map_vel;
        map_x += map_vel * dt;
        //      map_vel *= 0.98;
      }
    }
    if (0) {
      Vec3 diff = net_force_full + net_imp_force - (net_momentum - pre_momentum);
      P(dj::Vec3d(&pre_momentum[0]));
      P(dj::Vec3d(&net_momentum[0]));
      P(dj::Vec3d(&net_force_full[0]));
      P(dj::Vec3d(&net_imp_force[0]));
      P(dj::Vec3d(&diff[0]));
      pre_momentum = net_momentum;
    }


    OMP_FOR
    for (int  i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      center_of_mass_[p] += translational_vel_[p] * dt;
      quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
      quaternion_[p].Normalize();
      //      translational_vel_[p] *= 0.98;
      //      angular_vel_[p] *= 0.98;

      quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
      part_rotation_transpose_[p] = part_rotation_[p].transpose();

      current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
      current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();

      Mat3& rotation = part_rotation_[p];
      MapVec sub_q(&part_q_[p][0], part_basis_size_[p]);
      //    Vec& sub_q = part_q_[p];//(&q_[basis_offset_[p]], part_basis_size_[p]);
      for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
        int v = vert_local_id2global_id_[p][local_v];
        Vec3 local_u = vert_basis_[v] * sub_q;
        MapVec3 map_x(X + v * 3);
        map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
      }
    }

    // compute global position and update vertex offset
    inv_fem_->UpdateOffset();
    UpdateEmbededVertexPosition();
  }
  profiler.End("update");
  profiler.End("total");
}

void MixedMultiDomainTet::SimulatePBD(double dt) {
  UpdateBasisOffSet(true);
  // solve rigid motion for reduced domains
  profiler.Start("total");
  const double dt_2 = dt * dt;
  double internal_force_scaling_factor = conf.Get<double>("internal force scaling");
  internal_force_scaling_factor = 1.0;
  inv_fem_->inv_fem_force_model_->ComputePartialForceAndStiffnessMatrixWithSplitedMesh(&all_cubature_tet_[0],
                                                                                       int(all_cubature_tet_.size()),
                                                                                       tet_id_,
                                                                                       tet_, X);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external collision force in world space to each vertex
  ext_force_.clear();
  profiler.Start("collision force");
  // Ground collision
  if (1) {
    const double kFloor = -0.20;
    const double kFloorStiffness = 25000;
    for (int v = 0; v < vertex_num_; ++v) {
      if (X[v * 3 + 1] < kFloor) {
        Vec3 force(0, 0, 0);
        force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
        //      force[1] = (kFloor - X[v * 3 + 1])  * kFloorStiffness;
        ext_force_.emplace_back(v, force);
      }
    }
  }
  profiler.End("collision force");
  // user interaction force
  {
#if 1
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      ext_force_.emplace_back(v, force);
    }
#endif
  }
  Vec new_rigid_vel;
  SimulationRigidMotion(dt, new_rigid_vel, &ext_force_);
  //  new_rigid_vel.setZero();
  profiler.Start("rigid update");
  std::vector<Vec3> acceleration(part_num_);
  std::vector<Vec3> angular_acceleration(part_num_);
  std::vector<Vec3> old_vel(part_num_);

  // Compute translational acc and angular acc
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    int idx = i * 6;
    MapVec3 new_translational_vel(&new_rigid_vel[idx]);
    MapVec3 new_angular_vel(&new_rigid_vel[idx + 3]);
    new_angular_vel.setZero();
    acceleration[p] = (new_translational_vel - translational_vel_[p]) / dt;
    angular_acceleration[p] = (new_angular_vel - angular_vel_[p]) / dt;

    // TODO delete
    old_vel[p] = translational_vel_[p];
    translational_vel_[p] = new_translational_vel;
    angular_vel_[p] = new_angular_vel;
  }
  profiler.End("rigid update");

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Local subspace deformation
#if 1
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // RHS
  // simulate internal dynamics
  Vec subspace_rhs = Vec::Zero(total_basis_num_);
  // collision and ui force
  for (std::pair<int, Vec3>& collision : ext_force_) {
    const int& v = collision.first;
    const Vec3& force = collision.second;
    int p = vert_part_id_[v];
    if (is_subspace_domain_[p]) {
      MapVec map_rhs(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
      map_rhs += vert_basis_transpose_[v] * (part_rotation_transpose_[p] * force);
    } else {
      //      int local_v = local_vert_idx_[v];
      //      MapVec3 map_rhs(&subspace_rhs[basis_offset_[p] + local_v * 3]);
      //      map_rhs += force;
    }
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  AddFictitiousForceAndGravity(acceleration, angular_acceleration, subspace_rhs);

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // internal force
  profiler.Start("internal force");
  {
    for (CubaturePoint & cubature : all_cubature_) {
      int& t = cubature.first;
      double& weight = cubature.second;
      double* force = inv_fem_->element_force_ + t * 12;
      for (int i4 = 0; i4 < 4; ++i4) {
        int v = tet_[t * 4 + i4];
        int p = vert_part_id_[v];
        MapVec3 map_force(force + i4 * 3);
        if (is_subspace_domain_[p]) {
          MapVec subspace_force(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]);
          // negative sign before weight is because internal force from vega is in oppisite direction of actual force
          subspace_force -= vert_basis_transpose_[v] * (internal_force_scaling_factor * weight * (part_rotation_[p].transpose() * map_force));
        } else {
          //          int local_v = local_vert_idx_[v];
          //          MapVec3 full_force(&subspace_rhs[basis_offset_[p] + local_v * 3]);
          //          full_force -= map_force;
        }
      }
    }

    subspace_rhs *= dt;
    // add velocity related term
    for (int p = 0; p < part_num_; ++p) {
      if (is_subspace_domain_[p]) {
        MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]) += MapVec(&part_vel_[p][0], part_basis_size_[p]);
        //            MapVec(&subspace_rhs[basis_offset_[p]], part_basis_size_[p]).setZero();
      } else {
        //        for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
        //          int global_v = vert_local_id2global_id_[p][local_v];
        //          int offset = basis_offset_[p] + local_v * 3;
        //          MapVec3(&subspace_rhs[offset]) += mass_[global_v] * (MapVec3(&part_vel_[p][local_v * 3]) + gravity_ * dt);
        //          //        MapVec3(&subspace_rhs[offset]).setZero();
        //        }
      }
    }
  }
  profiler.End("internal force");

  // Compute reduced stiffness matrix
  profiler.Start("mul p");
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    Mat& reduced_k = domain_k_[p];
    reduced_k.setZero();
    // all cubature vertex in this domain
    for (int j = 0; j < int(v_list_[p].size()); ++j) {
      int v = v_list_[p][j];
      int p = vert_part_id_[v];
      Mat3 sub_k = Mat3::Zero();
      // all cubature tet incident on this v
      for (int k = 0; k < int(v_cubature_[p][j].size()); ++k) {
        VVCubature& cubature = v_cubature_[p][j][k];
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
      reduced_k += vert_basis_transpose_[v]
                   * (part_rotation_transpose_[p] * sub_k * part_rotation_[p])
                   * vert_basis_[v];
    }
  }
  profiler.End("mul p");
  profiler.Start("mul e");
  OMP_FOR
  for (int i = 0; i < part_num_ + int(subspace_subspace_interface_.size()); ++i) {
    int e = i;
    int p0, p1;
    double* mat = NULL;
    if (e < part_num_) {
      if (!is_subspace_domain_[e]) continue; // ignor full domain
      p0 = p1 = e;
      mat = domain_k_[p0].data();
    } else {
      p0 = subspace_subspace_interface_[e - part_num_].first;
      p1 = subspace_subspace_interface_[e - part_num_].second;
      mat = interface_k_[e - part_num_].data();
      interface_k_[e - part_num_].setZero();
      e = part_num_ + topology_[p0][p1];
    }
    Eigen::Map<Mat> reduced_k(mat, part_basis_size_[p0], part_basis_size_[p1]);
    // list of vert-vert pair in this domain
    for (int i = 0; i < int(vv_list_[e].size()); ++i) {
      int vi = vv_list_[e][i].first;
      int vj = vv_list_[e][i].second;
      // list of cubature tets that contains this vert-vert pair
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

  // code for testing
  {
    //  for (int i = 0; i < subspace_domain_num_; ++i) {
    //    int p = subspace_domains_[i];
    //    Mat& result = domain_k_[p];
    //    result *= dt_2 * internal_force_scaling_factor;
    //    for (int i = 0; i < part_basis_size_[p]; ++i) {
    //      result(i, i) += 1;
    //    }
    //    P(p, result.norm());
    //  }
    //  exit(0);
    //  for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
    //    int p0 = subspace_subspace_interface_[e].first;
    //    int p1 = subspace_subspace_interface_[e].second;
    //    Mat& result = interface_k_[e];
    //    result *= (dt_2 * internal_force_scaling_factor);
    //        P(p0, p1, e, result.norm());
    //        Mat3 mm = result.block<3, 3>(0, 0);
    //        PMAT(mm);
    //  }
    //  exit(0);
    //    profiler.EndTimer("internal force");
  }

  // full-subspace interface k
  // pth domain
  //  for (int n = 0; n < subspace_domain_num_; ++n) {
  //    int p = subspace_domains_[n];
  //    // ith adjacent full domain
  //    for (int i = 0; i < int(full_vert_on_domain_[p].size()); ++i) {
  //      // jth vertex
  //      for (int j = 0; j < int(full_vert_on_domain_[p][i].second.size()); ++j) {
  //        Mat& vert_k = full_vert_subspace_matrix_block_[p][i][j];
  //        vert_k.setZero();
  //        // kth incident tet
  //        for (int k = 0; k < int(vert_incident_tet_domain_[p][i][j].size()); ++k) {
  //          const int t = vert_incident_tet_domain_[p][i][j][k].first;
  //          const int idx = vert_incident_tet_domain_[p][i][j][k].second;
  //          Eigen::Map<Eigen::Matrix<double, 12, 12> > tet_k(inv_fem_->element_k_ + t * 144);
  //          for (int l = 0; l < 4; ++l) {
  //            int v_subspace = tet_[t * 4  + l];
  //            if (l == idx || vert_part_id_[v_subspace] != p) continue;
  //            vert_k += vert_basis_transpose_[v_subspace] * (part_rotation_transpose_[p] * tet_k.block<3, 3>(l * 3, idx * 3));
  //          }
  //        }
  //      }
  //    }
  //  }

  auto StiffnessMatrix = [&](double * x, double * result) {
    memset(result, 0, sizeof(double) * total_basis_num_);
    // explicit
    if (0) {
      for (int p = 0; p < part_num_; ++p) {
        if (is_subspace_domain_[p]) {
          MapVec map_result(result + basis_offset_[p], part_basis_size_[p]);
          MapVec map_x(x + basis_offset_[p], part_basis_size_[p]);
          map_result = map_x;
        } else {
          for (int v = 0; v < vert_num_per_part_[p]; ++v) {
            int global_v = vert_local_id2global_id_[p][v];
            double mass = mass_[global_v];
            int offset = basis_offset_[p] + v * 3;
            result[offset + 0] = x[offset + 0] * mass;
            result[offset + 1] = x[offset + 1] * mass;
            result[offset + 2] = x[offset + 2] * mass;
          }
        }
      }
      return;
    }

    // subspace domain
    if (1) {
      for (int i = 0; i < subspace_domain_num_; ++i) {
        int p = subspace_domains_[i];
        MapVec map_result(&result[basis_offset_[p]], part_basis_size_[p]);
        MapVec map_x(&x[basis_offset_[p]], part_basis_size_[p]);
        map_result += domain_k_[p] * map_x;
        //      map_result.setZero();
      }
    }

    // subspace-subspace interface
    if (1) {
      for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
        int p0 = subspace_subspace_interface_[e].first;
        int p1 = subspace_subspace_interface_[e].second;
        {
          MapVec map_result(&result[basis_offset_[p0]], part_basis_size_[p0]);
          MapVec map_x(&x[basis_offset_[p1]], part_basis_size_[p1]);
          map_result += interface_k_[e] * map_x;
        }
        {
          MapVec map_result(&result[basis_offset_[p1]], part_basis_size_[p1]);
          Eigen::Map<Eigen::RowVectorXd> map_x(&x[basis_offset_[p0]], part_basis_size_[p0]);
          map_result += map_x * interface_k_[e];
        }
      }
    }

    // subspace-full inteface
    if (0) {
      for (int i = 0; i < subspace_domain_num_; ++i) {
        int p = subspace_domains_[i];
        // jth full domain
        for (int j = 0; j < int(full_vert_on_domain_[p].size()); ++j) {
          int full_domain = full_vert_on_domain_[p][j].first;
          // kth vertex
          for (int k = 0; k < int(full_vert_on_domain_[p][j].second.size()); ++k) {
            int v = full_vert_on_domain_[p][j].second[k];
            int local_v = local_vert_idx_[v];
            Mat& vert_k = full_vert_subspace_matrix_block_[p][j][k];
            {
              MapVec map_result(&result[basis_offset_[p]], part_basis_size_[p]);
              MapVec3 map_x(&x[basis_offset_[full_domain] + local_v * 3]);
              map_result += vert_k * map_x;
            }
            {
              MapVec3 map_result(&result[basis_offset_[full_domain] + local_v * 3]);
              Eigen::Map<Eigen::RowVectorXd> map_x(&x[basis_offset_[p]], part_basis_size_[p]);
              map_result += map_x * vert_k;
            }
          }
        }
      }
    }

    // tet complete inside full domain
    if (0) {
      for (int t : full_tet_) {
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > tet_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          const int v0 = verts[i];
          const int p0 = vert_part_id_[v0];
          const int local_v0 = local_vert_idx_[v0];
          for (int j = 0; j < 4; ++j) {
            const int v1 = verts[j];
            const int p1 = vert_part_id_[v1];
            const int local_v1 = local_vert_idx_[v1];
            MapVec3 map_result(result + basis_offset_[p0] + local_v0 * 3);
            MapVec3 map_x(x + basis_offset_[p1] + local_v1 * 3);
            map_result += tet_k.block<3, 3>(i * 3, j * 3) * map_x;
          }
        }
      }
    }
    // tet in full-subspace interface
    if (0) {
      for (int t : full_subspace_interface_tet_) {
        int* verts = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > tet_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          const int v0 = verts[i];
          const int p0 = vert_part_id_[v0];
          const int local_v0 = local_vert_idx_[v0];
          if (is_subspace_domain_[p0]) continue;
          for (int j = 0; j < 4; ++j) {
            const int v1 = verts[j];
            const int p1 = vert_part_id_[v1];
            if (is_subspace_domain_[p1]) continue;
            const int local_v1 = local_vert_idx_[v1];
            MapVec3 map_result(result + basis_offset_[p0] + local_v0 * 3);
            MapVec3 map_x(x + basis_offset_[p1] + local_v1 * 3);
            map_result += tet_k.block<3, 3>(i * 3, j * 3) * map_x;
          }
        }
      }
    }

    // multiply dt * dt and add mass term
    for (int p = 0; p < part_num_; ++p) {
      if (is_subspace_domain_[p]) {
        MapVec map_result(result + basis_offset_[p], part_basis_size_[p]);
        MapVec map_x(x + basis_offset_[p], part_basis_size_[p]);
        map_result *= dt_2;
        map_result += map_x;
        //        map_result.setZero();
      } else {
        //        for (int v = 0; v < vert_num_per_part_[p]; ++v) {
        //          int global_v = vert_local_id2global_id_[p][v];
        //          double mass = mass_[global_v];
        //          int offset = basis_offset_[p] + v * 3;
        //          result[offset + 0] = x[offset + 0] * mass + dt_2 * result[offset + 0];
        //          result[offset + 1] = x[offset + 1] * mass + dt_2 * result[offset + 1];
        //          result[offset + 2] = x[offset + 2] * mass + dt_2 * result[offset + 2];
        //          //          if (is_constrainted_[global_v]) {
        //          //            result[offset + 0] = 0;
        //          //            result[offset + 1] = 0;
        //          //            result[offset + 2] = 0;
        //          //          }
        //
        //        }
      }
    }
  };

  // testing code
  if (0) {
    inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
    auto FullStiffnessMatrix = [&](double * x, double * result) {
      memset(result, 0, sizeof(double) * 3 * vertex_num_);
      //    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(x, result);
      //    if (0)
      for (int t = 0; t < tet_number; ++t) {
        int* verts = tet_ + t * 4;
        double (*k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int vi = verts[i];
          double* force = result + vi * 3;
          for (int j = 0; j < 4; ++j) {
            int vj = verts[j] * 3;
            force[0] += k[i * 3 + 0][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 0][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 0][j * 3 + 2] * x[vj + 2];
            force[1] += k[i * 3 + 1][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 1][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 1][j * 3 + 2] * x[vj + 2];
            force[2] += k[i * 3 + 2][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 2][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 2][j * 3 + 2] * x[vj + 2];
          }
        }
      }
      for (int v = 0; v < vertex_num_; ++v) {
        result[v * 3 + 0] = mass_[v] * x[v * 3 + 0] + dt_2 * result[v * 3 + 0];
        result[v * 3 + 1] = mass_[v] * x[v * 3 + 1] + dt_2 * result[v * 3 + 1];
        result[v * 3 + 2] = mass_[v] * x[v * 3 + 2] + dt_2 * result[v * 3 + 2];
      }
    };

    int p = 2;
    Vec x = Vec::Zero(total_basis_num_);
    Vec result = Vec::Zero(total_basis_num_);
    double test[3][3];
    for (int e = 0; e < edge_num_; ++e) {
      int global_v = edges_[e * 2];
      int other_v = edges_[e * 2 + 1];
      //      global_v = 4196;
      //      other_v = 12648;
      int p0 = vert_part_id_[global_v];
      int p1 = vert_part_id_[other_v];
      if (p0 != p || p1 != p) continue;
      int v = local_vert_idx_[global_v];
      int local_other = local_vert_idx_[other_v];
      //        P(global_v, other_v);
      //        P(v, local_other);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          x[basis_offset_[p] + v * 3 + c] = 1;
          StiffnessMatrix(&x[0], &result[0]);
          test[r][c] = result[basis_offset_[p] + local_other * 3 + r];
          //          std::cout <<  result[basis_offset_[p] + local_other * 3 + r] << "\t";
          x[basis_offset_[p] + v * 3 + c] = 0;
        }
        //        std::cout << std::endl;
      }

      //        KK;
      double truth[3][3];
      {
        Vec x = Vec::Zero(vertex_num_ * 3);
        Vec rhs = Vec::Zero(vertex_num_ * 3);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            x[global_v * 3 + c] = 1;
            FullStiffnessMatrix(&x[0], &rhs[0]);
            truth[r][c] = rhs[other_v * 3 + r];
            //            std::cout <<  rhs[other_v * 3 + r] << "\t";
            ASSERT(dj::Abs(truth[r][c] - test[r][c]) < 1e-8, P(r, c, truth[r][c], test[r][c]) P(global_v, other_v));
            x[global_v * 3 + c] = 0;
          }
          //          std::cout << std::endl;
        }
      }
      //      break;
    }
    L("full domain k verified");
    exit(0);
  }
  profiler.Start("solver");
  cg_solver_->Resize(total_basis_num_);

  //  for (int v = 0; v < vertex_num_; ++v) {
  //    if (is_constrainted_[v]) {
  //      int p = vert_part_id_[v];
  //      int local_v = local_vert_idx_[v];
  //      int offset = basis_offset_[p] + local_v * 3;
  //      subspace_rhs[offset + 0] = 0;
  //      subspace_rhs[offset + 1] = 0;
  //      subspace_rhs[offset + 2] = 0;
  //    }
  //  }

  if (0)
    for (int p = 0; p < part_num_; ++p) {
      if (is_subspace_domain_[p]) {
        MapVec map_result(&subspace_rhs[0] + basis_offset_[p], part_basis_size_[p]);
        map_result.setZero();
      }
    }
  memset(tmp_vertex_pos_, 0, sizeof(double) * subspace_rhs.size());
  auto info = cg_solver_->Solve(&subspace_rhs[0], tmp_vertex_pos_, StiffnessMatrix, 2000, 1e-16);
  P(info.first, info.second);
  profiler.End("solver");
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#else
#endif // LOCAL_DEFORMATION

  // Momentum conservation for full-subspace tets
  std::vector<Vec3> correct_force(full_subspace_interface_tet_.size());
  //  memset(&correct_force[0][0], 0, sizeof(Vec3) * full_subspace_interface_tet_.size());
  Vec3 net_f(0, 0, 0);
  Vec3 implicit_f(0, 0, 0);
  Vec3 implicit_f1(0, 0, 0);
  Vec3 net_f1(0, 0, 0);
  Mat K = Mat::Zero(3, 3);
  if (0) {
    for (int ii = 0; ii < int(full_subspace_interface_tet_.size()); ++ii) {
      int t = full_subspace_interface_tet_[ii];
      int* vert = tet_ + t * 4;
      Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
      correct_force[ii].setZero();
      for (int i = 0; i < 4; ++i) {
        int vi = vert[i];
        int pi = vert_part_id_[vi];
        if (pi == 0) {
          net_f -= MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
        } else {
          net_f1 -= MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
        }
        if (is_subspace_domain_[pi]) {
          for (int j = 0; j < 4; ++j) {
            int vj = vert[j];
            int pj = vert_part_id_[vj];
            if (!is_subspace_domain_[pj]) continue;
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            correct_force[ii] -= element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] - rj.cross(angular_vel_[pj]));
            if (pi == 0) {
              implicit_f -= element_k.block<3, 3>(i * 3, j * 3) * (translational_vel_[pj] - rj.cross(angular_vel_[pj])) * dt;
              K += element_k.block<3, 3>(i * 3, j * 3);
            }
          }
        } else {
          for (int j = 0; j < 4; ++j) {
            int vj = vert[j];
            int pj = vert_part_id_[vj];
            if (is_subspace_domain_[pj])  {
              correct_force[ii] -= element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] *
                                   MapVec(&tmp_vertex_pos_[basis_offset_[pj]], part_basis_size_[pj]);
              implicit_f1 -= element_k.block<3, 3>(i * 3, j * 3) * part_rotation_[pj] * vert_basis_[vj] *
                             MapVec(&tmp_vertex_pos_[basis_offset_[pj]], part_basis_size_[pj]);
            } else {
              int local_v = local_vert_idx_[vj];
              correct_force[ii] -= element_k.block<3, 3>(i * 3, j * 3) * MapVec3(&tmp_vertex_pos_[basis_offset_[pj] + local_v * 3]);
              implicit_f1 -= element_k.block<3, 3>(i * 3, j * 3) * MapVec3(&tmp_vertex_pos_[basis_offset_[pj] + local_v * 3]);
            }
          }
        }
      }
      correct_force[ii] *= 0.25 * dt * dt * 1e-1;
      //      correct_force[ii] *= dt * dt;
    }
#if 0
    PVEC(net_f1);
    PVEC(net_f);
    net_f *= dt;
    implicit_f1 *= dt * dt;
    implicit_f *= dt;
    PVEC(implicit_f1);
    PVEC(implicit_f);
    PVEC(net_f);
    Vec3 g = gravity_ * mass_per_part_[0] * dt + net_f;
    PVEC(g);
    //    Mat3 diff = K - Rigid.block<3, 3>(0, 0);
    //    PMAT(diff);
    //    implicit_f = K * translational_vel_[0] * dt * -dt;
    Vec3 f1 = gravity_ * mass_[4] * dt + implicit_f1 + net_f1 * dt;
    PVEC(f1);
    Vec3 m4 = mass_[4] * (MapVec3(velocity_ + basis_offset_[1] + local_vert_idx_[4] * 3) -
                          MapVec3(&part_vel_[1][local_vert_idx_[4] * 3]));
    PVEC(m4);
    ASSERT((f1 - m4).norm() < 1e-10);
    Vec3 test_rhs = g + old_vel[0] * mass_per_part_[0];
    Vec3 implicit_f01 = implicit_f + implicit_f1;
    PVEC(implicit_f01);
    PVEC(test_rhs);
    g += implicit_f;
    Vec3 net_int_force = g + f1;
    PVEC(net_int_force);
    PVEC(g);
    Vec3 dv = (translational_vel_[0] - old_vel[0]) * mass_per_part_[0];
    ASSERT((dv - g).norm() < 1e-8);
    Vec3 total_dv = dv + m4;
    PVEC(total_dv);
    PVEC(dv);
    PVEC(correct_force[0]);
    Vec3 correct_mom = correct_force[0] * 4;
    PVEC(correct_mom);
    KK;
    KK;
    //    exit(0);
#endif

    PVEC(translational_vel_[0]);
    PVEC(angular_vel_[0]);
    if (0) {
      Vec3 total_momentum(0, 0, 0);
      double total_mass = 0;
      for (int p = 0; p < part_num_; ++p) {
        if (is_subspace_domain_[p]) {
          total_mass += mass_per_part_[p];
          total_momentum += mass_per_part_[p] * translational_vel_[p];
        } else {
          for (int v = 0; v < vert_num_per_part_[p]; ++v) {
            int global = vert_local_id2global_id_[p][v];
            //            if (global != 4) continue;
            total_mass += mass_[global];
            //            P(global, dj::Vec3d(velocity_ + basis_offset_[p] + v * 3));
            total_momentum += MapVec3(tmp_vertex_pos_ + basis_offset_[p] + v * 3) * mass_[global];
          }
        }
      }
      static int i = 1;
      Vec3 g_mom = i * gravity_ * total_mass * dt;
      Vec rigid_mom = translational_vel_[0] * mass_per_part_[0];
      PVEC(rigid_mom);
      PVEC(total_momentum);
      PVEC(g_mom);

    }
    //    if (0)
    for (int ii = 0; ii < int(full_subspace_interface_tet_.size()); ++ii) {
      int t = full_subspace_interface_tet_[ii];
      int* vert = tet_ + t * 4;
      Vec3& force = correct_force[ii];
      //            P(t, force.norm());
      for (int i = 0; i < 4; ++i) {
        int vi = vert[i];
        int pi = vert_part_id_[vi];
        if (is_subspace_domain_[pi]) {
          translational_vel_[pi] -= (force) / mass_per_part_[pi];
          Vec3 r = MapVec3(X + vi * 3) - center_of_mass_[pi];
          angular_vel_[pi] -= current_inv_inertia_tensor_[pi] * r.cross(force);
        } else {
          int local_v = local_vert_idx_[vi];
          MapVec3(tmp_vertex_pos_ + basis_offset_[pi] + local_v * 3) -= force / mass_[vi];
        }
      }
    }

  }
  Vec3 total_momentum(0, 0, 0);
  double total_mass = 0;
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      total_mass += mass_per_part_[p];
      total_momentum += mass_per_part_[p] * translational_vel_[p];
    } else {
      //      for (int v = 0; v < vert_num_per_part_[p]; ++v) {
      //        int global = vert_local_id2global_id_[p][v];
      //        //        if (global != 4) continue;
      //        total_mass += mass_[global];
      //        //        P(global, dj::Vec3d(velocity_ + basis_offset_[p] + v * 3));
      //        total_momentum += MapVec3(velocity_ + basis_offset_[p] + v * 3) * mass_[global];
      //      }
    }
  }
  profiler.Start("PBD");
  TetLimiting(dt, 1000);
  profiler.End("PBD");
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // update rigid motion, global position, and vertex offset
  //  OMP_FOR
  profiler.Start("update");
  {
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      memcpy(&part_vel_[p][0], tmp_vertex_pos_ + basis_offset_[p], sizeof(double) * part_basis_size_[p]);
      MapVec map_q(&part_q_[p][0], part_basis_size_[p]);
      MapVec map_vel(&part_vel_[p][0], part_basis_size_[p]);
      map_q += map_vel * dt;
    }

    //    for (int i = 0; i < int(full_simulation_domains_.size()); ++i) {
    //      int p = full_simulation_domains_[i];
    //      memcpy(&part_vel_[p][0], tmp_vertex_pos_ + basis_offset_[p], sizeof(double) * part_basis_size_[p]);
    //      for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
    //        int global_v = vert_local_id2global_id_[p][local_v];
    //        MapVec3 map_x(X + global_v * 3);
    //        MapVec3 map_vel(&part_vel_[p][local_v * 3]);
    //        map_x += map_vel * dt;
    //        //      map_vel *= 0.98;
    //      }
    //    }


    OMP_FOR
    for (int  i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      center_of_mass_[p] += translational_vel_[p] * dt;
      quaternion_[p] = quaternion_[p] + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[p][0], angular_vel_[p][1], angular_vel_[p][2]) * quaternion_[p];
      quaternion_[p].Normalize();
      //      translational_vel_[p] *= 0.98;
      //      angular_vel_[p] *= 0.98;

      quaternion_[p].Quaternion2Matrix(part_rotation_[p].data());
      part_rotation_transpose_[p] = part_rotation_[p].transpose();

      current_inertia_tensor_[p] = part_rotation_[p] * inertia_tensor_[p] * part_rotation_transpose_[p];
      current_inv_inertia_tensor_[p] = current_inertia_tensor_[p].inverse();

      Mat3& rotation = part_rotation_[p];
      MapVec sub_q(&part_q_[p][0], part_basis_size_[p]);
      //    Vec& sub_q = part_q_[p];//(&q_[basis_offset_[p]], part_basis_size_[p]);
      for (int local_v = 0; local_v < vert_num_per_part_[p]; local_v++) {
        int v = vert_local_id2global_id_[p][local_v];
        Vec3 local_u = vert_basis_[v] * sub_q;
        MapVec3 map_x(X + v * 3);
        map_x = rotation * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_[p];
      }
    }

    // compute global position and update vertex offset
    inv_fem_->UpdateOffset();
    UpdateEmbededVertexPosition();
  }
  profiler.End("update");
  profiler.End("total");
}




void MixedMultiDomainTet::SimulateFull(double dt) {
  Vec rhs = Vec::Zero(vertex_num_ * 3);
  inv_fem_->inv_fem_force_model_->ComputePartialForceAndStiffnessMatrixWithSplitedMesh(NULL, tet_number, tet_id_, tet_, X);
  //  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  // accumulate internal force
  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    double* force = inv_fem_->element_force_ + t * 12;
    for (int i = 0; i < 4; ++i) {
      rhs[verts[i] * 3 + 0] -= force[i * 3 + 0];
      rhs[verts[i] * 3 + 1] -= force[i * 3 + 1];
      rhs[verts[i] * 3 + 2] -= force[i * 3 + 2];
    }
  }

  //  P(dj::Vec3d(&rhs[5 * 3]));
  //  P(rhs.norm());
  double ui_force[3];
  int force_vert = GetUIForce(ui_force);
  if (force_vert >= 0) {
    //    P(force_vert, dj::Vec3d(ui_force));
    rhs[force_vert * 3 + 0] += ui_force[0];
    rhs[force_vert * 3 + 1] += ui_force[1];
    rhs[force_vert * 3 + 2] += ui_force[2];
  }

  //  rhs[2 * 3 + 0] += 10.5;
  for (int v = 0; v < vertex_num_; ++v) {
    //    P(v, dj::Vec3d(&rhs[v * 3]));
    rhs[v * 3 + 0] = rhs[v * 3 + 0] * dt + mass_[v] * velocity_[v * 3 + 0];
    rhs[v * 3 + 1] = rhs[v * 3 + 1] * dt + mass_[v] * velocity_[v * 3 + 1];
    rhs[v * 3 + 2] = rhs[v * 3 + 2] * dt + mass_[v] * velocity_[v * 3 + 2];

    if (is_constrainted_[v]) {
      rhs[v * 3 + 0] = 0;
      rhs[v * 3 + 1] = 0;
      rhs[v * 3 + 2] = 0;
    }
  }
  //  exit(0);

  //  std::vector<double> tmp_force(tet_number * 12);
  const double dt_2 = dt * dt;
  auto StiffnessMatrix = [&](double * x, double * result) {
    memset(result, 0, sizeof(double) * 3 * vertex_num_);
    //    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(x, result);
    //    if (0)
    for (int t = 0; t < tet_number; ++t) {
      int* verts = tet_ + t * 4;
      double (*k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
      for (int i = 0; i < 4; ++i) {
        int vi = verts[i];
        double* force = result + vi * 3;
        for (int j = 0; j < 4; ++j) {
          int vj = verts[j] * 3;
          force[0] += k[i * 3 + 0][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 0][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 0][j * 3 + 2] * x[vj + 2];
          force[1] += k[i * 3 + 1][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 1][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 1][j * 3 + 2] * x[vj + 2];
          force[2] += k[i * 3 + 2][j * 3 + 0] * x[vj + 0] +  k[i * 3 + 2][j * 3 + 1] * x[vj + 1] +  k[i * 3 + 2][j * 3 + 2] * x[vj + 2];
        }
      }
    }
    for (int v = 0; v < vertex_num_; ++v) {
      result[v * 3 + 0] = mass_[v] * x[v * 3 + 0] + dt_2 * result[v * 3 + 0];
      result[v * 3 + 1] = mass_[v] * x[v * 3 + 1] + dt_2 * result[v * 3 + 1];
      result[v * 3 + 2] = mass_[v] * x[v * 3 + 2] + dt_2 * result[v * 3 + 2];
      if (is_constrainted_[v]) {
        result[v * 3 + 0] = 0;
        result[v * 3 + 1] = 0;
        result[v * 3 + 2] = 0;
      }
    }
  };
  cg_solver_->Resize(vertex_num_ * 3);
  //  P(rhs.norm());
  auto info = cg_solver_->Solve(&rhs[0], velocity_, StiffnessMatrix, 2000, 1e-12);
  (void) info;
  //  P(info.first, info.second);
  //  MapVec vel(velocity_, vertex_num_ * 3);
  //  P(vel.norm());
  //  P(dj::Vec3d(velocity_ + 5 * 3), dj::Vec3d(&rhs[5 * 3]), dj::Vec3d(X + 5 * 3));
  //  OMP_FOR
  for (int v = 0; v < vertex_num_; ++v) {
    //    P(v, dj::Vec3d(velocity_ + v * 3));
    X[v * 3 + 0] += velocity_[v * 3 + 0] * dt;
    X[v * 3 + 1] += velocity_[v * 3 + 1] * dt;
    X[v * 3 + 2] += velocity_[v * 3 + 2] * dt;
    inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
    inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
    inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
    const double kDamping = 1.00;
    velocity_[v * 3 + 0] *= kDamping;
    velocity_[v * 3 + 1] *= kDamping;
    velocity_[v * 3 + 2] *= kDamping;
  }
  UpdateEmbededVertexPosition();
}

#if 1
void MixedMultiDomainTet::Cut(MultiDomainTet::Vec3 pa, MultiDomainTet::Vec3 pb, MultiDomainTet::Vec3 pc) {
  // TODO: update full_verts_, full_vert_idx_, bi_cg_, cg_ size
  profiler.Start("cut");
  std::unordered_set<int> affected_vertex_set;
  std::unordered_set<int> splitted_edge;
  std::unordered_set<int> splitted_tet;
  std::unordered_map<int, int> intersected_edge_per_tet_;
  Vec3 normal = (pb - pa).cross(pc - pa);
  // intersect edges with triangle
  for (int e = 0; e < original_edge_num_; ++e) {
    if (mid_point_weight_[e] >= 0) continue;
    int* vert = edges_ + e * 2;
    int p0 = vert_part_id_[vert[0]];
    int p1 = vert_part_id_[vert[1]];
    if (is_subspace_domain_[p0] || is_subspace_domain_[p1]) continue;
    double weight, u, v;
    bool intersect = dj::SegmentTriangleIntersection(X + vert[0] * 3,
                                                     X + vert[1] * 3,
                                                     &pa[0], &pb[0], &pc[0],
                                                     1e-5, u, v, weight);
    if (intersect) {
      mid_point_weight_[e] = 1 - weight;
      splitted_edge.insert(e);
      affected_vertex_set.insert(vert[0]);
      affected_vertex_set.insert(vert[1]);
      //      P(incident_tet_on_edge_[e].size());
      for (int tet : incident_tet_on_edge_[e]) {
        ++intersected_edge_per_tet_[tet];
        if (intersected_edge_per_tet_[tet] == 3) {
          splitted_tet.insert(tet);
        }
      }
      //      ASSERT(mid_point_weight_[e] > 0 && mid_point_weight_[e] < 1);
    }
  }
  //  P(splitted_tet.size());

  //  std::vector<int> affected_vertex(affected_vertex_set.begin(), affected_vertex_set.end());
  //  std::vector<std::unordered_map<int, int> > virtual_neighbor;
  std::vector<int> connected_neighbors;
  connected_neighbors.reserve(512);
  std::vector<int> tmp_neighor;
  tmp_neighor.reserve(512);
  std::unordered_set<int> disconnected_neighbors;
  std::queue<int> queue;
  std::vector<std::unordered_set<int> > connected_component;
  std::unordered_map<int, std::unordered_map<int, int> > virtual_nodes;
  int virtual_node_num = 0;
  // decide for each vertex if there is a scoop cut
  for (int v : affected_vertex_set) {
    connected_neighbors.clear();
    disconnected_neighbors.clear();
    connected_component.clear();
    // find all non-connected neighbors
    for (std::pair<int, int> neighbor : adjacent_vertex_[v]) {
      int& vert = neighbor.first;
      int& edge = neighbor.second;
      if (splitted_edge.count(edge) == 0) {
        connected_neighbors.push_back(vert);
      } else {
        disconnected_neighbors.insert(vert);
      }
    }
    if (disconnected_neighbors.size() == 0) continue;
    for (int i = 0; i < int(connected_neighbors.size()); ++i) {
      int n = connected_neighbors[i];
      tmp_neighor.clear();
      for (int disconnected_v : disconnected_neighbors) {
        if (adjacent_vertex_[n].count(disconnected_v) > 0) {
          int e = adjacent_vertex_[n][disconnected_v];
          if (splitted_edge.count(e) == 0) {
            tmp_neighor.push_back(disconnected_v);
            connected_neighbors.push_back(disconnected_v);
          }
        }
      }
      for (int reconnected_neighbor : tmp_neighor) {
        disconnected_neighbors.erase(reconnected_neighbor);
      }
    }

    //    P(disconnected_neighbors.size());
    // get disjoint scoops
    while (disconnected_neighbors.empty() == false) {
      auto iter = disconnected_neighbors.begin();
      int vert = *iter;
      ++iter;
      connected_component.emplace_back(disconnected_neighbors.begin(), iter);
      disconnected_neighbors.erase(vert);
      queue.push(vert);
      while (!queue.empty()) {
        int front = queue.front();
        queue.pop();
        for (std::pair<int, int> neighbor : adjacent_vertex_[front]) {
          if (disconnected_neighbors.count(neighbor.first) > 0 && splitted_edge.count(neighbor.second) == 0) {
            queue.push(neighbor.first);
            disconnected_neighbors.erase(neighbor.first);
            connected_component.back().insert(neighbor.first);
          }
        }
      }
    }

    // Create one virtual node for each scoop
    for (int i = 0; i < int(connected_component.size()); ++i) {
      int new_vertex_id = vertex_num_ + virtual_node_num;
      int p = vert_part_id_[v];
      int local_v = local_vert_idx_[v];
      int new_local_vert_id = vert_num_per_part_[p];
      virtual_node_num++;
      X[new_vertex_id * 3 + 0] = X[v * 3 + 0];
      X[new_vertex_id * 3 + 1] = X[v * 3 + 1];
      X[new_vertex_id * 3 + 2] = X[v * 3 + 2];

      part_vel_[p].push_back(part_vel_[p][local_v * 3 + 0]);
      part_vel_[p].push_back(part_vel_[p][local_v * 3 + 1]);
      part_vel_[p].push_back(part_vel_[p][local_v * 3 + 2]);

      velocity_[new_vertex_id * 3 + 0] = velocity_[v * 3 + 0];
      velocity_[new_vertex_id * 3 + 1] = velocity_[v * 3 + 1];
      velocity_[new_vertex_id * 3 + 2] = velocity_[v * 3 + 2];

      rest_pos_[new_vertex_id * 3 + 0] = rest_pos_[v * 3 + 0];
      rest_pos_[new_vertex_id * 3 + 1] = rest_pos_[v * 3 + 1];
      rest_pos_[new_vertex_id * 3 + 2] = rest_pos_[v * 3 + 2];
      mass_[new_vertex_id] = mass_[v];

      vertex_id_[new_vertex_id] = v;
      full_verts_.push_back(new_vertex_id);
      full_vert_idx_[new_vertex_id] = int(full_verts_.size()) - 1;
      vert_local_id2global_id_[p].push_back(new_vertex_id);
      local_vert_idx_.push_back(new_local_vert_id);
      vert_part_id_.push_back(p);
      vert_num_per_part_[p]++;

      for (int neighbor : connected_component[i]) {
        virtual_nodes[v][neighbor] = new_vertex_id;
      }
    }
  }

  if (0) {
    for (auto aa : virtual_nodes) {
      KK;
      P(aa.first);
      for (auto bb : aa.second) {
        P(bb.first, bb.second);
      }
    }
  }

  // table for surface triangle construction
  const int kFaceTable13[4][4][3] = {
    {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}},
    {{1, 2, 3}, {2, 0, 3}, {1, 3, 0}, {1, 0, 2}},
    {{2, 3, 1}, {2, 0, 3}, {0, 1, 3}, {2, 1, 0}},
    {{3, 1, 2}, {3, 2, 0}, {3, 0, 1}, {0, 2, 1}},
  };
  const int kEdgeTable[4][4] = {
    {INT_MIN, INT_MIN, INT_MIN, INT_MIN},
    {2, 3, 3, 2},
    {3, 1, 1, 3},
    {1, 2, 2, 1},
  };
  const int kFaceTable2[4][4][5] = {
    {INT_MIN, INT_MIN, INT_MIN, INT_MIN},
    {{1, 2, 3, 3, 2}, {0, 3, 2, 1, 0}, {3, 0, 1, 1, 2}, {2, 1, 0, 3, 0}},
    {{2, 3, 1, 3, 2}, {3, 2, 0, 3, 0}, {0, 1, 3, 1, 0}, {1, 0, 2, 1, 2}},
    {{3, 1, 2, 3, 2}, {2, 0, 3, 1, 2}, {1, 3, 0, 3, 0}, {0, 2, 1, 1, 0}},
  };

  // store which side of the cutting triangle the vertex lies in
  std::vector<int> vertex_side(vertex_num_, INT_MIN);
  for (int v : affected_vertex_set) {
    vertex_side[v] = (MapVec3(X + v * 3) - pa).dot(normal) >= 0;
  }

  std::unordered_map<std::pair<int, int>, int, IntPairHash> embeded_v_map;
  int new_tet_num = 0;
  //  P(intersected_edge_per_tet_[39]);exit(0);
  if (0) {
    int t = 39;
    P(dj::Vec4i(tet_ + t * 4));
    int* e = &edge_on_tet_[t * 6];
    for (int i = 0; i < 6; ++i) {
      P(mid_point_weight_[e[i]], e[i], dj::Vec2i(edges_ + e[i] * 2));
    }
    //    exit(0);
  }
  // split tetrahedra
  for (int tet : splitted_tet) {
    if (intersected_edge_per_tet_[tet] == 3) {
      int count = 0;
      int* verts = tet_ + tet * 4;
      for (int i = 0; i < 4; ++i) {
        if (vertex_side[verts[i]]) count++;
      }
      if (count != 1 && count != 3) continue;
    }
    int new_tet_id = tet_number + new_tet_num;
    new_tet_num++;

    int verts[4] = {INT_MIN, INT_MIN, INT_MIN, INT_MIN};
    int last_splitted_vertex_idx = INT_MIN;
    int new_nodes[4] = {INT_MIN, INT_MIN, INT_MIN, INT_MIN};
    // decide for each node if a virtual node is created and assign to
    // other nodes in the same tetrehedra
    for (int i = 0; i < 4; ++i) {
      int v = tet_[tet * 4 + i];
      if (virtual_nodes.count(v) > 0 &&
          (virtual_nodes[v].count(tet_[tet * 4 + 0]) > 0 ||
           virtual_nodes[v].count(tet_[tet * 4 + 1]) > 0 ||
           virtual_nodes[v].count(tet_[tet * 4 + 2]) > 0 ||
           virtual_nodes[v].count(tet_[tet * 4 + 3]) > 0)) {
        //        if (verbose) P(v, i);
        last_splitted_vertex_idx = i;
        for (int j = 0; j < 4; ++j) {
          if (virtual_nodes[v].count(tet_[tet * 4 + j]) > 0) {
            new_nodes[i] = virtual_nodes[v][tet_[tet * 4 + j]];
            break;
          }
        }
      } else {
        verts[i] = v;
      }
    }
    if (last_splitted_vertex_idx == INT_MIN) continue;
    int new_tet0[4] = {tet_[tet * 4 + 0], tet_[tet * 4 + 1], tet_[tet * 4 + 2], tet_[tet * 4 + 3]};
    int new_tet1[4] = {tet_[tet * 4 + 0], tet_[tet * 4 + 1], tet_[tet * 4 + 2], tet_[tet * 4 + 3]};
    int last_splitted_vert = tet_[tet * 4 + last_splitted_vertex_idx];
    new_tet1[last_splitted_vertex_idx] = new_nodes[last_splitted_vertex_idx];
    //    P(dj::Vec4i(new_tet0))
    //    P(dj::Vec4i(new_tet1))
    for (int i = last_splitted_vertex_idx - 1; i >= 0; --i) {
      if (verts[i] != INT_MIN) continue;
      int v = tet_[tet * 4 + i];
      if (virtual_nodes[last_splitted_vert].count(v) > 0) {
        new_tet0[i] = new_nodes[i];
        new_tet1[i] = v;
      } else {
        new_tet0[i] = v;
        new_tet1[i] = new_nodes[i];
      }
    }

    int tet_type = 0;
    int isolated_v_num = 0;
    for (int i = 0; i < 4; ++i) {
      int v = tet_[tet * 4 + i];
      int mark = (vertex_side[v] ^ vertex_side[last_splitted_vert]);
      isolated_v_num += mark;
      tet_type = tet_type | (mark << i);
    }

    int old_tet_verts[4] = {
      tet_[tet * 4 + 0],
      tet_[tet * 4 + 1],
      tet_[tet * 4 + 2],
      tet_[tet * 4 + 3],
    };
    tet_[tet * 4 + 0] = new_tet0[0];
    tet_[tet * 4 + 1] = new_tet0[1];
    tet_[tet * 4 + 2] = new_tet0[2];
    tet_[tet * 4 + 3] = new_tet0[3];
    tet_type_[tet] = (~tet_type) & 0xf;
    //    for (int i = 0; i < 4; ++i) {
    //      ASSERT(new_tet0[i] >= 0, P(tet, i));
    //      ASSERT(new_tet1[i] >= 0, P(tet, i));
    //    }
    tet_[new_tet_id * 4 + 0] = new_tet1[0];
    tet_[new_tet_id * 4 + 1] = new_tet1[1];
    tet_[new_tet_id * 4 + 2] = new_tet1[2];
    tet_[new_tet_id * 4 + 3] = new_tet1[3];
    tet_type_[new_tet_id] = tet_type;
    tet_id_[new_tet_id] = tet;
    full_tet_.push_back(new_tet_id);
    all_cubature_.emplace_back(new_tet_id, 1.0);
    all_cubature_tet_.push_back(new_tet_id);
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // update surface triangle
    if (isolated_v_num == 1 || isolated_v_num == 3) {
      int *tet0, *tet1;
      int type;
      if (isolated_v_num == 1) {
        tet0 = &new_tet0[0];
        tet1 = &new_tet1[0];
        type = tet_type;
      } else {
        tet0 = &new_tet1[0];
        tet1 = &new_tet0[0];
        type = (~tet_type) & 0xf;
      }
      int single_v = 0;
      for (; single_v < 4; ++single_v) {
        if (type & (1 << single_v)) {
          break;
        }
      }
      ASSERT(single_v < 4);
      int edge_vert[4][2] = {
        INT_MIN, INT_MIN, INT_MIN, INT_MIN,
        INT_MIN, INT_MIN, INT_MIN, INT_MIN,
      };
      //      P(dj::Vec2i(&edge_vert[0][0]));
      //      P(dj::Vec2i(&edge_vert[1][0]));
      //      P(dj::Vec2i(&edge_vert[2][0]));
      //      P(dj::Vec2i(&edge_vert[3][0]));
      for (int i = 1; i <= 3; ++i) {
        int idx = (single_v + i) % 4;
        std::pair<int, int> other_v_pair(tet0[single_v], tet0[idx]);
        std::pair<int, int> v_pair(tet1[single_v], tet1[idx]);
        SwapIfGreater(v_pair.first, v_pair.second);
        SwapIfGreater(other_v_pair.first, other_v_pair.second);
        auto iter = embeded_v_map.find(v_pair);
        if (iter == embeded_v_map.end()) {
          std::pair<int, int> orig_v_pair(vertex_id_[tet1[single_v]], vertex_id_[tet1[idx]]);
          if (orig_v_pair.first > orig_v_pair.second) {
            dj::Swap(orig_v_pair.first, orig_v_pair.second);
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end());
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet1[idx], tet1[single_v], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet0[idx], tet0[single_v], mid_point_weight_[e]);
          } else {
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end());
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet1[single_v], tet1[idx], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet0[single_v], tet0[idx], mid_point_weight_[e]);
          }
          //          KK;
          //          P(old_tet_verts[single_v], old_tet_verts[idx], -embeded_v_num_);
          //          P(tet0[single_v], tet0[idx]);
          embeded_vertex_.push_back(-embeded_v_num_ - 1);
          embeded_vertex_.push_back(-embeded_v_num_ - 2);
          embeded_v_map[v_pair] = -embeded_v_num_ - 1;
          embeded_v_map[other_v_pair] = -embeded_v_num_ - 2;
          edge_vert[idx][0] = -embeded_v_num_ - 1;
          edge_vert[idx][1] = -embeded_v_num_ - 2;
          embeded_v_num_ += 2;
        } else {
          edge_vert[idx][0] = iter->second;
          auto other_iter = embeded_v_map.find(other_v_pair);
          ASSERT(other_iter != embeded_v_map.end());
          edge_vert[idx][1] = other_iter->second;
        }
      }

      for (int i = 1; i <= 3; ++i) {
        int idx = (single_v + i) % 4;
        int tri_idx = surface_tri_on_tet_[tet][idx];
        if (tri_idx >= 0) {
          const int* tbl = &kFaceTable13[single_v][idx][0];
          SplitTriangle(tri_idx,
                        old_tet_verts[tbl[0]], old_tet_verts[tbl[1]], old_tet_verts[tbl[2]],
                        edge_vert[tbl[1]][0], edge_vert[tbl[1]][1],
                        edge_vert[tbl[2]][0], edge_vert[tbl[2]][1]);
        }
      }
      if (textured_surface_renderer_.size() != 0) {
        for (int tri = triangle_num_; tri <= triangle_num_ + 1; ++tri) {
          triangle_group_info_[tri].first = int(groups_.size()) - 1;
          triangle_group_info_[tri].second = int(groups_.back().size());
          groups_.back().push_back(tri);
        }
      }
      T[triangle_num_ * 3 + 0] = edge_vert[kFaceTable13[single_v][single_v][0]][0];
      T[triangle_num_ * 3 + 1] = edge_vert[kFaceTable13[single_v][single_v][1]][0];
      T[triangle_num_ * 3 + 2] = edge_vert[kFaceTable13[single_v][single_v][2]][0];
      is_cut_triangle_[triangle_num_] = true;
      triangle_num_++;

      T[triangle_num_ * 3 + 0] = edge_vert[kFaceTable13[single_v][single_v][0]][1];
      T[triangle_num_ * 3 + 2] = edge_vert[kFaceTable13[single_v][single_v][1]][1];
      T[triangle_num_ * 3 + 1] = edge_vert[kFaceTable13[single_v][single_v][2]][1];
      is_cut_triangle_[triangle_num_] = true;
      triangle_num_++;
    } else if (isolated_v_num == 2) {
      int *tet0, *tet1;
      int mark = tet_type & 0x1;
      if (mark == 0) {
        tet0 = &new_tet0[0];
        tet1 = &new_tet1[0];
      } else {
        tet0 = &new_tet1[0];
        tet1 = &new_tet0[0];
      }
      int other_v = 1;
      for (; other_v < 4 && (((tet_type >> other_v) & 0x1) != mark); other_v++) {
      }
      ASSERT(other_v < 4);
      int edge_vert[4][2];
      // find edge verts;
      for (int i = 0; i < 4; ++i) {
        int idx0 = (i < 2) ? 0 : other_v;
        int idx1 = kEdgeTable[other_v][i];

        std::pair<int, int> v_pair(tet0[idx0], tet0[idx1]);
        SwapIfGreater(v_pair.first, v_pair.second);
        std::pair<int, int> other_v_pair(tet1[idx0], tet1[idx1]);
        SwapIfGreater(other_v_pair.first, other_v_pair.second);
        auto iter = embeded_v_map.find(v_pair);
        if (iter == embeded_v_map.end()) {
          std::pair<int, int> orig_v_pair(vertex_id_[tet0[idx0]], vertex_id_[tet0[idx1]]);
          if (orig_v_pair.first > orig_v_pair.second) {
            dj::Swap(orig_v_pair.first, orig_v_pair.second);
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end());
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet0[idx1], tet0[idx0], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet1[idx1], tet1[idx0], mid_point_weight_[e]);
          } else {
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end(), int* a = NULL; *a = 0; P(orig_v_pair.first, orig_v_pair.second, idx0, idx1));
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet0[idx0], tet0[idx1], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet1[idx0], tet1[idx1], mid_point_weight_[e]);
          }

          embeded_vertex_.push_back(-embeded_v_num_ - 1);
          embeded_vertex_.push_back(-embeded_v_num_ - 2);
          embeded_v_map[v_pair] = -embeded_v_num_ - 1;
          embeded_v_map[other_v_pair] = -embeded_v_num_ - 2;
          edge_vert[i][0] = -embeded_v_num_ - 1;
          edge_vert[i][1] = -embeded_v_num_ - 2;
          embeded_v_num_ += 2;
        } else {
          edge_vert[i][0] = iter->second;
          auto other_iter = embeded_v_map.find(other_v_pair);
          ASSERT(other_iter != embeded_v_map.end());
          edge_vert[i][1] = other_iter->second;
        }
      }
      // split triangles
      for (int i = 0; i < 4; ++i) {
        if (surface_tri_on_tet_[tet][i] < 0) continue;
        int tri_verts[3] = {
          old_tet_verts[kFaceTable2[other_v][i][0]],
          old_tet_verts[kFaceTable2[other_v][i][1]],
          old_tet_verts[kFaceTable2[other_v][i][2]],
        };
        int *embeded_vert[2] = {
          &edge_vert[kFaceTable2[other_v][i][3]][0],
          &edge_vert[kFaceTable2[other_v][i][4]][0],
        };
        if (kFaceTable2[other_v][i][0] == 0 || kFaceTable2[other_v][i][0] == other_v) {
          SplitTriangle(surface_tri_on_tet_[tet][i],
                        tri_verts[0], tri_verts[1], tri_verts[2],
                        embeded_vert[0][0], embeded_vert[0][1],
                        embeded_vert[1][0], embeded_vert[1][1]);
        } else {
          SplitTriangle(surface_tri_on_tet_[tet][i],
                        tri_verts[0], tri_verts[1], tri_verts[2],
                        embeded_vert[0][1], embeded_vert[0][0],
                        embeded_vert[1][1], embeded_vert[1][0]);
        }
      }
      if (textured_surface_renderer_.size() != 0) {
        for (int tri = triangle_num_; tri < triangle_num_ + 4; ++tri) {
          triangle_group_info_[tri].first = int(groups_.size()) - 1;
          triangle_group_info_[tri].second = int(groups_.back().size());
          groups_.back().push_back(tri);
        }
      }

      T[triangle_num_ * 3 + 0] = edge_vert[0][0];
      T[triangle_num_ * 3 + 1] = edge_vert[1][0];
      T[triangle_num_ * 3 + 2] = edge_vert[2][0];
      is_cut_triangle_[triangle_num_] = true;
      triangle_num_++;
      T[triangle_num_ * 3 + 0] = edge_vert[0][1];
      T[triangle_num_ * 3 + 2] = edge_vert[1][1];
      T[triangle_num_ * 3 + 1] = edge_vert[2][1];
      is_cut_triangle_[triangle_num_] = true;
      triangle_num_++;

      T[triangle_num_ * 3 + 0] = edge_vert[2][0];
      T[triangle_num_ * 3 + 1] = edge_vert[3][0];
      T[triangle_num_ * 3 + 2] = edge_vert[0][0];
      is_cut_triangle_[triangle_num_] = true;
      triangle_num_++;
      T[triangle_num_ * 3 + 0] = edge_vert[2][1];
      T[triangle_num_ * 3 + 2] = edge_vert[3][1];
      T[triangle_num_ * 3 + 1] = edge_vert[0][1];
      is_cut_triangle_[triangle_num_] = true;
      triangle_num_++;
    } else {
      ASSERT(false);
    }
  }

  tet_number += new_tet_num;
  vertex_num_ += virtual_node_num;
  //   P(tet_number, vertex_num_, new_tet_num, virtual_node_num);
  if (0) {
    for (int t = 0; t < tet_number; ++t) {
      P(tet_type_[t]);
      P(dj::Vec4i(tet_ + t * 4));
    }
  }
  //  bi_cg_->Resize(int(full_verts_.size()) * 3 + int(subspace_domains_.size()) * 6);
  UpdateBasisOffSet();
  UpdateEmbededVertexPosition();
  UpdateSurfaceTriangleTopology();
  position_changed_ = true;
  UpdateTextureCoordinateVBO();
  UpdateTexturedSurfaceMeshVBO();
  profiler.Start("delete");
  delete sparse_k_;
  profiler.End("delete");
  profiler.Start("new solver");
  sparse_k_ = new MixedSparseMatrix(this);
  profiler.End("new solver");
  profiler.End("cut");
}
#else
void MixedMultiDomainTet::Cut(MultiDomainTet::Vec3 pa, MultiDomainTet::Vec3 pb, MultiDomainTet::Vec3 pc) {
  std::unordered_set<int> affected_vertex_set;
  std::unordered_set<int> splitted_edge;
  std::unordered_set<int> splitted_tet;
  std::unordered_map<int, int> intersected_edge_per_tet_;
  Vec3 normal = (pb - pa).cross(pc - pa);
  //  pvec(pa);
  //  pvec(pb);
  //  pvec(dir);
  //  pvec(dir0);
  //  pvec(normal);
  // intersect edges with scapel
  for (int e = 0; e < original_edge_num_; ++e) {
    if (mid_point_weight_[e] >= 0) continue;
    int* vert = edges_ + e * 2;
    double weight, u, v;
    bool intersect = dj::SegmentTriangleIntersection(X + vert[0] * 3,
                                                     X + vert[1] * 3,
                                                     &pa[0], &pb[0], &pc[0],
                                                     1e-5, u, v, weight);
    //    if (e == 9) P(weight, dj::Vec3d(X + v[0] * 2), dj::Vec3d(X + v[1] * 2));
    //    P(weight);
    //    P(vert[0], vert[1], e, weight);
    if (intersect) {
      if (vert[0] == 151 || vert[1] == 151) {
        P(weight, vert[0], vert[1], e)
      }
      //      P(e, weight);
      mid_point_weight_[e] = 1 - weight;
      splitted_edge.insert(e);
      affected_vertex_set.insert(vert[0]);
      affected_vertex_set.insert(vert[1]);
      //      P(incident_tet_on_edge_[e].size());
      for (int tet : incident_tet_on_edge_[e]) {
        ++intersected_edge_per_tet_[tet];
        if (intersected_edge_per_tet_[tet] == 3) {
          splitted_tet.insert(tet);
        }
      }
      ASSERT(mid_point_weight_[e] > 0 && mid_point_weight_[e] < 1);
    }
  }
  //  exit(0);
  //  P(splitted_tet.size());

  //  std::vector<int> affected_vertex(affected_vertex_set.begin(), affected_vertex_set.end());
  //  std::vector<std::unordered_map<int, int> > virtual_neighbor;
  std::vector<int> connected_neighbors;
  connected_neighbors.reserve(256);
  std::vector<int> tmp_neighor;
  tmp_neighor.reserve(512);
  std::unordered_set<int> disconnected_neighbors;
  std::queue<int> queue;
  std::vector<std::unordered_set<int> > connected_component;
  std::unordered_map<int, std::unordered_map<int, int> > virtual_nodes;
  int virtual_node_num = 0;
  // decide for each vertex if there is a scoop cut
  for (int v : affected_vertex_set) {
    disconnected_neighbors.clear();
    connected_component.clear();
    connected_neighbors.clear();
    // find all disconnected neighbors and connected neighbors
    for (std::pair<int, int> neighbor : adjacent_vertex_[v]) {
      int& vert = neighbor.first;
      int& edge = neighbor.second;
      if (splitted_edge.count(edge) == 0) {
        connected_neighbors.push_back(vert);
      } else {
        disconnected_neighbors.insert(vert);
      }
    }
    if (disconnected_neighbors.size() == 0) continue;
    for (int i = 0; i < int(connected_neighbors.size()); ++i) {
      int n = connected_neighbors[i];
      tmp_neighor.clear();
      for (int disconnected_v : disconnected_neighbors) {
        if (adjacent_vertex_[n].count(disconnected_v) > 0) {
          int e = adjacent_vertex_[n][disconnected_v];
          if (splitted_edge.count(e) == 0) {
            tmp_neighor.push_back(disconnected_v);
            connected_neighbors.push_back(disconnected_v);
          }
        }
      }
      for (int reconnected_neighbor : tmp_neighor) {
        disconnected_neighbors.erase(reconnected_neighbor);
      }
    }


    //    P(disconnected_neighbors.size());
    // get disjoint scoops
    while (disconnected_neighbors.empty() == false) {
      auto iter = disconnected_neighbors.begin();
      int vert = *iter;
      ++iter;
      connected_component.emplace_back(disconnected_neighbors.begin(), iter);
      disconnected_neighbors.erase(vert);
      queue.push(vert);
      while (!queue.empty()) {
        int front = queue.front();
        queue.pop();
        for (std::pair<int, int> neighbor : adjacent_vertex_[front]) {
          if (disconnected_neighbors.count(neighbor.first) > 0 && splitted_edge.count(neighbor.second) == 0) {
            queue.push(neighbor.first);
            disconnected_neighbors.erase(neighbor.first);
            connected_component.back().insert(neighbor.first);
          }
        }
      }
    }

    // Create one virtual node for each scoop
    for (int i = 0; i < int(connected_component.size()); ++i) {
      int new_vertex_id = vertex_num_ + virtual_node_num;
      //      int p = vert_part_id_[v];
      //      int local_v = vert_local_id2global_id_[p]
      //      int new_local_vert_id = vert_num_per_part_[p] + virtual_node_num;
      virtual_node_num++;
      X[new_vertex_id * 3 + 0] = X[v * 3 + 0];
      X[new_vertex_id * 3 + 1] = X[v * 3 + 1];
      X[new_vertex_id * 3 + 2] = X[v * 3 + 2];
      //      part_vel_[p].push_back

      velocity_[new_vertex_id * 3 + 0] = velocity_[v * 3 + 0];
      velocity_[new_vertex_id * 3 + 1] = velocity_[v * 3 + 1];
      velocity_[new_vertex_id * 3 + 2] = velocity_[v * 3 + 2];

      rest_pos_[new_vertex_id * 3 + 0] = rest_pos_[v * 3 + 0];
      rest_pos_[new_vertex_id * 3 + 1] = rest_pos_[v * 3 + 1];
      rest_pos_[new_vertex_id * 3 + 2] = rest_pos_[v * 3 + 2];
      mass_[new_vertex_id] = mass_[v];

      vertex_id_[new_vertex_id] = v;
      for (int neighbor : connected_component[i]) {
        virtual_nodes[v][neighbor] = new_vertex_id;
      }
    }
  }

  if (0)
    for (auto aa : virtual_nodes) {
      KK;
      P(aa.first);
      for (auto bb : aa.second) {
        P(bb.first, bb.second);
      }
    }

  std::vector<int> vertex_side(vertex_num_, INT_MIN);
  for (int v : affected_vertex_set) {
    vertex_side[v] = (MapVec3(X + v * 3) - pa).dot(normal) >= 0;
  }
  int new_tet_num = 0;
  //  P(intersected_edge_per_tet_[39]);exit(0);
  if (0) {
    int t = 39;
    P(dj::Vec4i(tet_ + t * 4));
    int* e = &edge_on_tet_[t * 6];
    for (int i = 0; i < 6; ++i) {
      P(mid_point_weight_[e[i]], e[i], dj::Vec2i(edges_ + e[i] * 2));
    }
    //    exit(0);
  }

  const int kFaceTable13[4][4][3] = {
    {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}},
    {{1, 2, 3}, {2, 0, 3}, {1, 3, 0}, {1, 0, 2}},
    {{2, 3, 1}, {2, 0, 3}, {0, 1, 3}, {2, 1, 0}},
    {{3, 1, 2}, {3, 2, 0}, {3, 0, 1}, {0, 2, 1}},
  };
  const int kEdgeTable[4][4] = {
    {INT_MIN, INT_MIN, INT_MIN, INT_MIN},
    {2, 3, 3, 2},
    {3, 1, 1, 3},
    {1, 2, 2, 1},
  };
  const int kFaceTable2[4][4][5] = {
    {INT_MIN, INT_MIN, INT_MIN, INT_MIN},
    {{1, 2, 3, 3, 2}, {0, 3, 2, 1, 0}, {3, 0, 1, 1, 2}, {2, 1, 0, 3, 0}},
    {{2, 3, 1, 3, 2}, {3, 2, 0, 3, 0}, {0, 1, 3, 1, 0}, {1, 0, 2, 1, 2}},
    {{3, 1, 2, 3, 2}, {2, 0, 3, 1, 2}, {1, 3, 0, 3, 0}, {0, 2, 1, 1, 0}},
  };

  std::unordered_map<std::pair<int, int>, int, IntPairHash> embeded_v_map;

  // split tetrahedra
  for (int tet : splitted_tet) {
    if (intersected_edge_per_tet_[tet] == 3) {
      int count = 0;
      int* verts = tet_ + tet * 4;
      for (int i = 0; i < 4; ++i) {
        if (vertex_side[verts[i]]) count++;
      }
      if (count != 1 && count != 3) continue;
    }
    bool verbose = false;
    if (tet == 8) verbose = true;
    int new_tet_id = tet_number + new_tet_num;
    new_tet_num++;

    int verts[4] = {INT_MIN, INT_MIN, INT_MIN, INT_MIN};
    int last_splitted_vertex_idx = INT_MIN;
    int new_nodes[4] = {INT_MIN, INT_MIN, INT_MIN, INT_MIN};
    // decide for each node if a virtual node is created and assign to
    // other nodes in the same tetrehedra
    for (int i = 0; i < 4; ++i) {
      int v = tet_[tet * 4 + i];
      if (virtual_nodes.count(v) > 0 &&
          (virtual_nodes[v].count(tet_[tet * 4 + 0]) > 0 ||
           virtual_nodes[v].count(tet_[tet * 4 + 1]) > 0 ||
           virtual_nodes[v].count(tet_[tet * 4 + 2]) > 0 ||
           virtual_nodes[v].count(tet_[tet * 4 + 3]) > 0)) {
        //        if (verbose) P(v, i);
        last_splitted_vertex_idx = i;
        for (int j = 0; j < 4; ++j) {
          if (virtual_nodes[v].count(tet_[tet * 4 + j]) > 0) {
            new_nodes[i] = virtual_nodes[v][tet_[tet * 4 + j]];
            break;
          }
        }
      } else {
        verts[i] = v;
      }
    }
    // TODO fix it
    if (last_splitted_vertex_idx == INT_MIN) continue;
    //    if (verbose) {
    //      P(dj::Vec4i(new_nodes));
    //      exit(0);
    //    }
    int new_tet0[4] = {tet_[tet * 4 + 0], tet_[tet * 4 + 1], tet_[tet * 4 + 2], tet_[tet * 4 + 3]};
    int new_tet1[4] = {tet_[tet * 4 + 0], tet_[tet * 4 + 1], tet_[tet * 4 + 2], tet_[tet * 4 + 3]};
    int last_splitted_vert = tet_[tet * 4 + last_splitted_vertex_idx];
    new_tet1[last_splitted_vertex_idx] = new_nodes[last_splitted_vertex_idx];
    //    P(dj::Vec4i(new_tet0))
    //    P(dj::Vec4i(new_tet1))
    for (int i = last_splitted_vertex_idx - 1; i >= 0; --i) {
      if (verts[i] != INT_MIN) continue;
      int v = tet_[tet * 4 + i];
      if (virtual_nodes[last_splitted_vert].count(v) > 0) {
        new_tet0[i] = new_nodes[i];
        new_tet1[i] = v;
      } else {
        new_tet0[i] = v;
        new_tet1[i] = new_nodes[i];
      }
    }
    int tet_type = 0;
    //    ASSERT(vertex_side[last_splitted_vert] != INT_MIN);
    int isolated_v_num = 0;
    for (int i = 0; i < 4; ++i) {
      int v = tet_[tet * 4 + i];
      int mark = (vertex_side[v] ^ vertex_side[last_splitted_vert]);
      isolated_v_num += mark;
      tet_type = tet_type | (mark << i);
    }
    //    KK;
    int old_tet_verts[4] = {
      tet_[tet * 4 + 0],
      tet_[tet * 4 + 1],
      tet_[tet * 4 + 2],
      tet_[tet * 4 + 3],
    };

    tet_[tet * 4 + 0] = new_tet0[0];
    tet_[tet * 4 + 1] = new_tet0[1];
    tet_[tet * 4 + 2] = new_tet0[2];
    tet_[tet * 4 + 3] = new_tet0[3];
    //    P(tet_type, tet);
    tet_type_[tet] = (~tet_type) & 0xf;

    tet_[new_tet_id * 4 + 0] = new_tet1[0];
    tet_[new_tet_id * 4 + 1] = new_tet1[1];
    tet_[new_tet_id * 4 + 2] = new_tet1[2];
    tet_[new_tet_id * 4 + 3] = new_tet1[3];
    tet_type_[new_tet_id] = tet_type;
    tet_id_[new_tet_id] = tet;
    //    if (verbose) P(tet_type0, tet_type1);
    //    P(dj::Vec4i(new_tet0))
    //    P(dj::Vec4i(new_tet1))
    //    P(dj::Vec4i(tet_))
    //    P(dj::Vec4i(tet_ + 4))

    // update surface triangle
    if (isolated_v_num == 1 || isolated_v_num == 3) {
      int *tet0, *tet1;
      int type;
      if (isolated_v_num == 1) {
        tet0 = &new_tet0[0];
        tet1 = &new_tet1[0];
        type = tet_type;
      } else {
        tet0 = &new_tet1[0];
        tet1 = &new_tet0[0];
        type = (~tet_type) & 0xf;
      }
      int single_v = 0;
      for (; single_v < 4; ++single_v) {
        if (type & (1 << single_v)) {
          break;
        }
      }
      ASSERT(single_v < 4);
      int edge_vert[4][2] = {
        INT_MIN, INT_MIN, INT_MIN, INT_MIN,
        INT_MIN, INT_MIN, INT_MIN, INT_MIN,
      };
      //      P(dj::Vec2i(&edge_vert[0][0]));
      //      P(dj::Vec2i(&edge_vert[1][0]));
      //      P(dj::Vec2i(&edge_vert[2][0]));
      //      P(dj::Vec2i(&edge_vert[3][0]));
      for (int i = 1; i <= 3; ++i) {
        int idx = (single_v + i) % 4;
        std::pair<int, int> other_v_pair(tet0[single_v], tet0[idx]);
        std::pair<int, int> v_pair(tet1[single_v], tet1[idx]);
        SwapIfGreater(v_pair.first, v_pair.second);
        SwapIfGreater(other_v_pair.first, other_v_pair.second);
        auto iter = embeded_v_map.find(v_pair);
        if (iter == embeded_v_map.end()) {
          std::pair<int, int> orig_v_pair(vertex_id_[tet1[single_v]], vertex_id_[tet1[idx]]);
          if (orig_v_pair.first > orig_v_pair.second) {
            dj::Swap(orig_v_pair.first, orig_v_pair.second);
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end());
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet1[idx], tet1[single_v], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet0[idx], tet0[single_v], mid_point_weight_[e]);
          } else {
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end());
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet1[single_v], tet1[idx], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet0[single_v], tet0[idx], mid_point_weight_[e]);
          }
          //          KK;
          //          P(old_tet_verts[single_v], old_tet_verts[idx], -embeded_v_num_);
          //          P(tet0[single_v], tet0[idx]);
          embeded_vertex_.push_back(-embeded_v_num_ - 1);
          embeded_vertex_.push_back(-embeded_v_num_ - 2);
          //          if (v_pair == std::make_pair(170, 186)) {
          //            P(other_v_pair.first, other_v_pair.second); ASSERT(false);
          //          }
          embeded_v_map[v_pair] = -embeded_v_num_ - 1;
          embeded_v_map[other_v_pair] = -embeded_v_num_ - 2;
          edge_vert[idx][0] = -embeded_v_num_ - 1;
          edge_vert[idx][1] = -embeded_v_num_ - 2;
          //          P(idx);
          //          KK; KK;
          embeded_v_num_ += 2;
        } else {
          edge_vert[idx][0] = iter->second;
          auto other_iter = embeded_v_map.find(other_v_pair);
          ASSERT(other_iter != embeded_v_map.end(),
                 P(vertex_id_[other_v_pair.first]);
                 P(vertex_id_[other_v_pair.second]);
                 P(vertex_id_[v_pair.first]);
                 P(vertex_id_[v_pair.second]););
          edge_vert[idx][1] = other_iter->second;
        }
      }

      for (int i = 1; i <= 3; ++i) {
        int idx = (single_v + i) % 4;
        int tri_idx = surface_tri_on_tet_[tet][idx];
        if (tri_idx >= 0) {
          const int* tbl = &kFaceTable13[single_v][idx][0];
          SplitTriangle(tri_idx,
                        old_tet_verts[tbl[0]], old_tet_verts[tbl[1]], old_tet_verts[tbl[2]],
                        edge_vert[tbl[1]][0], edge_vert[tbl[1]][1],
                        edge_vert[tbl[2]][0], edge_vert[tbl[2]][1]);
        }
      }
      //      P(dj::Vec3i(&kFaceTable13[single_v][single_v][0]));
      T[triangle_num_ * 3 + 0] = edge_vert[kFaceTable13[single_v][single_v][0]][0];
      T[triangle_num_ * 3 + 1] = edge_vert[kFaceTable13[single_v][single_v][1]][0];
      T[triangle_num_ * 3 + 2] = edge_vert[kFaceTable13[single_v][single_v][2]][0];
      //      P(edge_vert[2][0], edge_vert[2][1]);
      //      P(dj::Vec3i(T + triangle_num_ * 3));
      triangle_num_++;
      T[triangle_num_ * 3 + 0] = edge_vert[kFaceTable13[single_v][single_v][0]][1];
      T[triangle_num_ * 3 + 2] = edge_vert[kFaceTable13[single_v][single_v][1]][1];
      T[triangle_num_ * 3 + 1] = edge_vert[kFaceTable13[single_v][single_v][2]][1];
      //      P(dj::Vec3i(T + triangle_num_ * 3));
      triangle_num_++;
    } else if (isolated_v_num == 2) {
#if 1
      int *tet0, *tet1;
      int mark = tet_type & 0x1;
      if (mark == 0) {
        tet0 = &new_tet0[0];
        tet1 = &new_tet1[0];
      } else {
        tet0 = &new_tet1[0];
        tet1 = &new_tet0[0];
      }
      int other_v = 1;
      for (; other_v < 4 && (((tet_type >> other_v) & 0x1) != mark); other_v++) {
      }
      ASSERT(other_v < 4);
      int edge_vert[4][2];
      // find edge verts;
      for (int i = 0; i < 4; ++i) {
        int idx0 = (i < 2) ? 0 : other_v;
        int idx1 = kEdgeTable[other_v][i];

        std::pair<int, int> v_pair(tet0[idx0], tet0[idx1]);
        SwapIfGreater(v_pair.first, v_pair.second);
        std::pair<int, int> other_v_pair(tet1[idx0], tet1[idx1]);
        SwapIfGreater(other_v_pair.first, other_v_pair.second);
        auto iter = embeded_v_map.find(v_pair);
        if (iter == embeded_v_map.end()) {
          std::pair<int, int> orig_v_pair(vertex_id_[tet0[idx0]], vertex_id_[tet0[idx1]]);
          if (orig_v_pair.first > orig_v_pair.second) {
            dj::Swap(orig_v_pair.first, orig_v_pair.second);
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            ASSERT(e_iter != vertex_pair2edge_.end());
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet0[idx1], tet0[idx0], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet1[idx1], tet1[idx0], mid_point_weight_[e]);
          } else {
            auto e_iter = vertex_pair2edge_.find(orig_v_pair);
            //            ASSERT(e_iter != vertex_pair2edge_.end(), int* a = NULL; *a = 0; P(orig_v_pair.first, orig_v_pair.second, idx0, idx1));
            int e = e_iter->second;
            embeded_edge_.emplace_back(tet0[idx0], tet0[idx1], mid_point_weight_[e]);
            embeded_edge_.emplace_back(tet1[idx0], tet1[idx1], mid_point_weight_[e]);
          }

          //          P(other_v, tet_type);
          //          P(old_tet_verts[idx0], old_tet_verts[idx1], -embeded_v_num_ - 1);
          //          P(tet1[idx0], tet1[idx1]);
          embeded_vertex_.push_back(-embeded_v_num_ - 1);
          embeded_vertex_.push_back(-embeded_v_num_ - 2);
          embeded_v_map[v_pair] = -embeded_v_num_ - 1;
          embeded_v_map[other_v_pair] = -embeded_v_num_ - 2;
          edge_vert[i][0] = -embeded_v_num_ - 1;
          edge_vert[i][1] = -embeded_v_num_ - 2;
          embeded_v_num_ += 2;
        } else {
          edge_vert[i][0] = iter->second;
          auto other_iter = embeded_v_map.find(other_v_pair);
          ASSERT(other_iter != embeded_v_map.end());
          edge_vert[i][1] = other_iter->second;
        }
      }
      // split triangles
      for (int i = 0; i < 4; ++i) {
        if (surface_tri_on_tet_[tet][i] < 0) continue;
        int tri_verts[3] = {
          old_tet_verts[kFaceTable2[other_v][i][0]],
          old_tet_verts[kFaceTable2[other_v][i][1]],
          old_tet_verts[kFaceTable2[other_v][i][2]],
        };
        int *embeded_vert[2] = {
          &edge_vert[kFaceTable2[other_v][i][3]][0],
          &edge_vert[kFaceTable2[other_v][i][4]][0],
        };
        //        P(dj::Vec3i(tri_verts));
        //        P(dj::Vec2i(embeded_vert[0]));
        //        P(dj::Vec2i(embeded_vert[1]));
        if (kFaceTable2[other_v][i][0] == 0 || kFaceTable2[other_v][i][0] == other_v) {
          SplitTriangle(surface_tri_on_tet_[tet][i],
                        tri_verts[0], tri_verts[1], tri_verts[2],
                        embeded_vert[0][0], embeded_vert[0][1],
                        embeded_vert[1][0], embeded_vert[1][1]);
        } else {
          SplitTriangle(surface_tri_on_tet_[tet][i],
                        tri_verts[0], tri_verts[1], tri_verts[2],
                        embeded_vert[0][1], embeded_vert[0][0],
                        embeded_vert[1][1], embeded_vert[1][0]);
        }
      }

      T[triangle_num_ * 3 + 0] = edge_vert[0][0];
      T[triangle_num_ * 3 + 1] = edge_vert[1][0];
      T[triangle_num_ * 3 + 2] = edge_vert[2][0];
      triangle_num_++;
      T[triangle_num_ * 3 + 0] = edge_vert[0][1];
      T[triangle_num_ * 3 + 2] = edge_vert[1][1];
      T[triangle_num_ * 3 + 1] = edge_vert[2][1];
      triangle_num_++;

      T[triangle_num_ * 3 + 0] = edge_vert[2][0];
      T[triangle_num_ * 3 + 1] = edge_vert[3][0];
      T[triangle_num_ * 3 + 2] = edge_vert[0][0];
      triangle_num_++;
      T[triangle_num_ * 3 + 0] = edge_vert[2][1];
      T[triangle_num_ * 3 + 2] = edge_vert[3][1];
      T[triangle_num_ * 3 + 1] = edge_vert[0][1];
      triangle_num_++;
#endif
    } else {
      ASSERT(false, int* a = 0; *a = 0;);
    }
  }

  tet_number += new_tet_num;
  vertex_num_ += virtual_node_num;
  P(tet_number, vertex_num_);
  if (0)
    for (int t = 0; t < tet_number; ++t) {
      P(tet_type_[t], tet_id_[t], t, dj::Vec4i(tet_ + t * 4));
    }
  //    exit(0);
  UpdateEmbededVertexPosition();
}
#endif

/**
 * @brief MixedMultiDomainTet::SetFullSimulationDomains
 * @param domains the set of domain ids that will be simulated in full space
 */
void MixedMultiDomainTet::SetFullSimulationDomains(std::unordered_set<int> full_domains) {
  full_simulation_domains_ = std::vector<int>(full_domains.begin(), full_domains.end());
  is_reduced_vertex_.resize(kMaxVertexNum);
  std::fill(is_reduced_vertex_.begin(), is_reduced_vertex_.end(), false);
  for (int v : all_cubature_vert_) {
    is_reduced_vertex_[v] = true;
  }

  // get list of subspace domains
  subspace_domains_.clear();
  is_subspace_domain_.resize(part_num_);
  domain_index_ = std::vector<int>(part_num_, 0);
  int count = 0;
  for (int p = 0; p < part_num_; ++p) {
    // full domain
    if (full_domains.count(p) == 0) {
      subspace_domains_.push_back(p);
      is_subspace_domain_[p] = true;
      domain_index_[p] = count;
      count++;
    } else {
      // subspace domain
      domain_index_[p] = INT_MIN;
    }
  }
  subspace_domain_num_ = int(subspace_domains_.size());
  // verify
  {
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      ASSERT(i == domain_index_[p]);
    }
  }

  // get the list of subspace-subpsace and full-subspace interface
  for (int p0 = 0; p0 < part_num_; ++p0) {
    for (int p1 = p0 + 1; p1 < part_num_; ++p1) {
      if (topology_[p0][p1] < 0) continue;
      if (full_domains.count(p0) == 0 && full_domains.count(p1) == 0) {
        subspace_subspace_interface_.emplace_back(p0, p1);
      } else  if (full_domains.count(p0) == 0 || full_domains.count(p1) == 0) { // one subspace domain, one full domain
        if (full_domains.count(p0) > 0) {
          full_subspace_interface_.emplace_back(p0, p1);
        } else {
          full_subspace_interface_.emplace_back(p1, p0);
        }
      }
    }
  }

  full_vert_on_domain_.resize(part_num_);
  full_vert_subspace_matrix_block_.resize(part_num_);
  std::unordered_map<std::pair<int, int>, std::unordered_set<int>, IntPairHash> full_subspace_interface_vert;
  for (int e = 0; e < edge_num_; ++e) {
    int v0 = edges_[e * 2 + 0];
    int v1 = edges_[e * 2 + 1];
    int p0 = vert_part_id_[v0];
    int p1 = vert_part_id_[v1];
    std::pair<int, int> interface_pair = std::make_pair(p0, p1);
    if (full_domains.count(p0) == 0 && full_domains.count(p1) > 0) {
    } else if (full_domains.count(p0) > 0 && full_domains.count(p1) == 0) {
      std::swap(interface_pair.first, interface_pair.second);
    } else {
      continue;
    }
    int full_domain = interface_pair.second;
    if (full_subspace_interface_vert.count(interface_pair) == 0) {
      std::unordered_set<int> empty_set;
      full_subspace_interface_vert[interface_pair] = empty_set;
    }
    int full_vert = (p0 == full_domain) ? v0 : v1;
    full_subspace_interface_vert[interface_pair].insert(full_vert);
  }

  vert_incident_tet_domain_.resize(part_num_);
  auto iter = full_subspace_interface_vert.begin();
  for (; iter != full_subspace_interface_vert.end(); ++iter) {
    const std::pair<int, int>& interface = iter->first;
    const std::unordered_set<int>& verts_set = iter->second;
    std::vector<int> verts(verts_set.begin(), verts_set.end());
    full_vert_on_domain_[interface.first].emplace_back(interface.second, verts);

    std::vector<std::vector<std::pair<int, int> > > tet_incident_on_verts;
    tet_incident_on_verts.resize(verts.size());
    for (int i = 0; i < int(verts.size()); ++i) {
      int v = verts[i];
      for (int t : incident_tet_[v]) {
        bool is_adjacent2subspace_domain = false;
        for (int j = 0; j < 4; ++j) {
          int v1 = tet_[t * 4 + j];
          int p = vert_part_id_[v1];
          if (p == interface.first) {
            is_adjacent2subspace_domain = true;
            break;
          }
        }
        if (is_adjacent2subspace_domain) {
          int j = 0;
          for (; j < 4; ++j) {
            if (tet_[t * 4 + j] == v) {
              tet_incident_on_verts[i].emplace_back(t, j);
              break;
            }
          }
          ASSERT(j < 4);
        }
      }
    }
    vert_incident_tet_domain_[interface.first].push_back(tet_incident_on_verts);
  }
  //p = 1	i = 0	j = 0	k = 2
  //    P(vert_incident_tet_domain_.size());
  //  P(vert_incident_tet_domain_[1].size());
  //  P(vert_incident_tet_domain_[1][0].size());
  //  P(vert_incident_tet_domain_[1][0][0].size());
  //  P(vert_incident_tet_domain_[1][0][0][2].first);
  //  P(vert_incident_tet_domain_[1][0][0][2].second);

  // For each subspace domain
  for (int p = 0; p < part_num_; ++p) {
    // has at least one adjacent full domain
    if (full_vert_on_domain_[p].size() == 0) continue;
    // for each adjacent full domain
    for (int i = 0; i < int(full_vert_on_domain_[p].size()); ++i) {
      std::vector<int>& verts = full_vert_on_domain_[p][i].second;
      full_vert_subspace_matrix_block_[p].push_back(std::vector<Mat>(verts.size(), Mat::Zero(part_basis_size_[p], 3)));
    }
  }

  std::unordered_map<int, double> cubature_map;
  for (CubaturePoint & cubature : all_cubature_) {
    cubature_map[cubature.first] = cubature.second;
  }

  std::unordered_set<int> all_cubature_tet_set(all_cubature_tet_.begin(), all_cubature_tet_.end());
  std::unordered_set<int> all_cubature_vert_set(all_cubature_vert_.begin(), all_cubature_vert_.end());
  for (int p  : full_simulation_domains_) {
    for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
      int global_v = vert_local_id2global_id_[p][local_v];
      is_reduced_vertex_[global_v] = false;
      all_cubature_vert_set.insert(global_v);
    }
  }

  std::vector<std::unordered_set<int>> interface_cubature_tet_set(interface_domains_.size());
  for (int e = 0; e < int(interface_domains_.size()); ++e) {
    for (int i = 0; i < int(interface_cubature_[e].size()); ++i) {
      interface_cubature_tet_set[e].insert(interface_cubature_[e][i].first);
    }
  }

  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    int p[4] = {
      vert_part_id_[verts[0]],
      vert_part_id_[verts[1]],
      vert_part_id_[verts[2]],
      vert_part_id_[verts[3]],
    };
    // the tet contains at least one node from full domain
    if (full_domains.count(p[0]) > 0 ||  full_domains.count(p[1]) > 0 ||  full_domains.count(p[2]) > 0 ||  full_domains.count(p[3]) > 0) {
      all_cubature_tet_set.insert(t);
      cubature_map[t] = 1.0;
      // at least one node from the tet is from subspace domain
      if (full_domains.count(p[0]) == 0 || full_domains.count(p[1]) == 0 || full_domains.count(p[2]) == 0 || full_domains.count(p[3]) == 0) {
        full_subspace_interface_tet_.push_back(t);
        {
          std::set<int> domain_set(p + 0, p + 4);
          std::vector<int> interface_id(domain_set.begin(), domain_set.end());
          std::sort(interface_id.begin(), interface_id.end());
          Vec5i interface_domain;
          interface_domain[0] = int(interface_id.size());
          while (int(interface_id.size()) < 4) {
            interface_id.push_back(-1);
          }
          for (int i = 0; i < int(interface_id.size()); ++i) {
            interface_domain[i + 1] = interface_id[i];
          }
          bool found = false;
          for (int e = 0; e < int(interface_domains_.size()); ++e) {
            if (interface_domains_[e] == interface_domain) {
              interface_cubature_tet_set[e].insert(t);
              found = true;
              break;
            }
          }
          ASSERT(found);
        }
      }
    }
  }

  //  if (0)
  for (int e = 0; e < int(interface_domains_.size()); ++e) {
    bool is_full_interface = false;
    for (int i = 0; i < interface_domains_[e][0]; ++i) {
      if (!is_subspace_domain_[interface_domains_[e][i + 1]]) {
        is_full_interface = true;
        break;
      }
    }
    if (is_full_interface) {
      interface_cubature_[e].clear();
      interface_cubature_[e].reserve(interface_cubature_tet_set[e].size());
      for (int t : interface_cubature_tet_set[e]) {
        interface_cubature_[e].emplace_back(t, 1.0);
      }
      interface_rigid_cubature_[e] = interface_cubature_[e];
    }
  }


  all_cubature_tet_ = std::vector<int>(all_cubature_tet_set.begin(), all_cubature_tet_set.end());
  all_cubature_vert_ = std::vector<int>(all_cubature_vert_set.begin(), all_cubature_vert_set.end());
  all_cubature_.clear();
  all_cubature_.reserve(cubature_map.size());
  for (auto cubature : cubature_map) {
    all_cubature_.emplace_back(cubature.first, cubature.second);
  }
  BuildParallelAccelerationStructure();
  UpdateBasisOffSet();

  domain_k_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      domain_k_[p] = Mat::Zero(part_basis_size_[p], part_basis_size_[p]);
    }
  }

  interface_k_.resize(subspace_subspace_interface_.size());
  for (int i = 0; i < int(subspace_subspace_interface_.size()); ++i) {
    std::pair<int, int>& interface = subspace_subspace_interface_[i];
    int p0 = interface.first;
    int p1 = interface.second;
    ASSERT(p0 < p1);
    interface_k_[i] = Mat::Zero(part_basis_size_[p0], part_basis_size_[p1]);
  }

  full_tet_.reserve(kMaxTetNum);
  for (int t = 0; t < tet_number; ++t) {
    int* verts = tet_ + t * 4;
    bool is_full_tet = true;
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      int p = vert_part_id_[v];
      if (is_subspace_domain_[p]) {
        is_full_tet = false;
        break;
      }
    }
    if (is_full_tet) {
      full_tet_.push_back(t);
    }
  }

  part_vel_.resize(part_num_);
  part_q_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      part_vel_[p] = std::vector<double>(part_basis_size_[p], 0);
      part_q_[p] = std::vector<double>(part_basis_size_[p], 0);
    } else {
      part_vel_[p].reserve(kMaxVertexNum / 10);
      part_vel_[p].resize(vert_num_per_part_[p] * 3);
      std::fill(part_vel_[p].begin(), part_vel_[p].end(), 0);
    }
  }

  // remove full domain from subspace domains adacency list
  if (1) {
    for (int p = 0; p < part_num_; ++p) {
      if (!is_subspace_domain_[p]) {
        domain_incident_interface_[p].clear();
      } else {
        int i = 0;
        while (i < int(domain_incident_interface_[p].size())) {
          int e = domain_incident_interface_[p][i].first;
          int idx = 2 - domain_incident_interface_[p][i].second;
          int other_p = interface_domains_[e][idx];
          if (!is_subspace_domain_[other_p]) {
            domain_incident_interface_[p].erase(domain_incident_interface_[p].begin() + i);
          } else {
            ++i;
          }
        }
      }
    }
  }

  local_vert_idx_.reserve(kMaxVertexNum);
  vert_part_id_.reserve(kMaxVertexNum);
  all_cubature_.reserve(1.2 * tet_number);
  full_tet_.reserve(kMaxTetNum);
  for (int p : full_simulation_domains_) {
    vert_local_id2global_id_[p].reserve(kMaxVertexNum);
    part_vel_[p].reserve(kMaxVertexNum);
  }
  P(full_tet_.size(), full_subspace_interface_tet_.size(), tet_number, full_verts_.size());

  full_verts_.clear();
  full_verts_.reserve(vertex_num_);
  full_vert_idx_.resize(kMaxVertexNum);
  std::fill_n(full_vert_idx_.begin(), vertex_num_, -1);
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    if (!is_subspace_domain_[p]) {
      for (auto neighbor : adjacent_vertex_[v]) {
        int v_neighbor = neighbor.first;
        int p_neighbor = vert_part_id_[v_neighbor];
        if (is_subspace_domain_[p_neighbor]) {
          full_verts_.push_back(v);
          full_vert_idx_[v] = int(full_verts_.size() - 1);
          break;
        }
      }
    }
  }
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    if (!is_subspace_domain_[p]) {
      bool is_on_full_subspace_interface = false;
      for (auto neighbor : adjacent_vertex_[v]) {
        int v_neighbor = neighbor.first;
        int p_neighbor = vert_part_id_[v_neighbor];
        if (is_subspace_domain_[p_neighbor]) {
          is_on_full_subspace_interface = true;
          break;
        }
      }
      if (!is_on_full_subspace_interface) {
        full_verts_.push_back(v);
        full_vert_idx_[v] = int(full_verts_.size() - 1);
      }
    }
  }

  bi_cg_->Resize(int(full_verts_.size()) * 3 + int(subspace_domains_.size()) * 6);


  // setup cholesky solver
  {
    std::vector<int> edges(subspace_subspace_interface_.size() * 2);
    for (int e = 0; e < int(subspace_subspace_interface_.size()); ++e) {
      edges[e * 2 + 0] = domain_index_[subspace_subspace_interface_[e].first];
      edges[e * 2 + 1] = domain_index_[subspace_subspace_interface_[e].second];
    }
    std::vector<int> block_size(subspace_domain_num_);
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      ASSERT(domain_index_[p] == i);
      block_size[i] = part_basis_size_[p];
    }
    delete chol_solver_;
    chol_solver_ = new solver::BLOCK_MATRIX_GRAPH<double>(subspace_domain_num_, edges, block_size);
    delete full_chol_solver_;
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      block_size[i] = part_basis_size_[p] + 6;
    }
    full_chol_solver_ = new solver::BLOCK_MATRIX_GRAPH<double>(subspace_domain_num_, edges, block_size);
  }
  // cubatures used in local deformation
  {
    subspace_cubature_.clear();
    for (CubaturePoint & cubature : all_cubature_) {
      int t = cubature.first;
      int* verts = tet_ + t * 4;
      for (int i = 0; i < 4; ++i) {
        int p = vert_part_id_[verts[i]];
        if (is_subspace_domain_[p]) {
          subspace_cubature_.push_back(cubature);
          break;
        }
      }
    }
  }

  is_fixed_domain_ = std::vector<int>(part_num_, false);
  domain_offset_.resize(part_num_);
  rigid_offset_.resize(part_num_);
  domain_size_.resize(part_num_);
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      domain_offset_[p] = basis_offset_[p] + domain_index_[p] * 6;
      rigid_offset_[p] = domain_offset_[p] + part_basis_size_[p];
      domain_size_[p] = part_basis_size_[p] + 6;
    } else {
      domain_offset_[p] = subspace_domain_basis_size_ + subspace_domain_num_ * 6;
      domain_size_[p] = INT_MIN;
      rigid_offset_[p] = INT_MIN;
    }
  }
  constrained_domain_offset_ = domain_offset_;
  constrained_rigid_offset_ = rigid_offset_;

  non_constrained_dof_ = subspace_domain_basis_size_ + subspace_domain_num_ * 6;
  constrained_dof_ = subspace_domain_basis_size_ + subspace_domain_num_ * 6;
  sparse_k_ = new MixedSparseMatrix(this);
}

void MixedMultiDomainTet::SimulationRigidMotion(const double dt,
                                                MultiDomainTet::Vec & new_rigid_velocity,
                                                const std::vector<MultiDomainTet::ExtForce> *ext_force) {
  double rigid_coupling_scaling = conf.Get<double>("rigid coupling scaling");
  // FIXME
  rigid_coupling_scaling = 1.0;
  //  double rigid_coupling_scaling = 1;
  //  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(all_cubature_tet_);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  std::vector<Vec3> net_force(subspace_domain_num_, Vec3(0, 0, 0));
  std::vector<Vec3> net_torque(subspace_domain_num_, Vec3(0, 0, 0));
  std::vector<Vec> subspace_interface_force(subspace_domain_num_);
  // gravity on each part
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    net_force[i] += mass_per_part_[p] * gravity_;
    subspace_interface_force[i] = Vec::Zero(part_basis_size_[p]);
  }
  // external force
  {
    if (ext_force == nullptr) { ext_force = &ext_force_; }
    for (ExtForce f : *ext_force) {
      const int& v = f.first;
      const Vec3& force = f.second;
      int p = vert_part_id_[v];
      if (!is_subspace_domain_[p]) continue;
      int idx = domain_index_[p];

      net_force[idx][1] += force[1];

      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_[vert_part_id_[v]];
      net_torque[idx] += r.cross(force);
    }
  }


  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // elastic force from boundary tets
  profiler.Start("boundary force");
  {
    // boundary force from subspace-subpsace interface
    for (int e = 0; e < int(interface_domains_.size()); ++e) {
      int p0 = interface_domains_[e][1];
      int p1 = interface_domains_[e][2];
      if (interface_domains_[e][0] > 2 || !is_subspace_domain_[p0] || !is_subspace_domain_[p1]) continue;
      // FIXME: ignor interface forces that come from interface of more than 2 domains
      Vec force0 = Vec::Zero(part_basis_size_[p0]);
      Vec force1 = Vec::Zero(part_basis_size_[p1]);
      for (CubaturePoint & cubature : interface_cubature_[e]) {
        int& t = cubature.first;
        double& weight = cubature.second;

        double* element_force = inv_fem_->element_force_ + t * 3 * 4;
        for (int local_v = 0; local_v < 4; ++local_v) {
          MapVec3 force(element_force + local_v * 3);
          int v = tet_[t * 4 + local_v];
          int p = vert_part_id_[v];
          // negative sign before weight is because internal force from vega is in oppisite direction of actual force
          if (p == p0) {
            force0 += vert_basis_transpose_[v] * (-weight * (part_rotation_transpose_[p] * force));
          } else {
            force1 += vert_basis_transpose_[v] * (-weight * (part_rotation_transpose_[p] * force));
          }
        }
      }
      ProjectInterfaceSubspaceForce(e, p0, p1, force0, force1);
      subspace_interface_force[domain_index_[p0]] += force0 * rigid_coupling_scaling;
      subspace_interface_force[domain_index_[p1]] += force1 * rigid_coupling_scaling;
    }

    OMP_FOR
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      Mat3& rotation = part_rotation_[p];
      net_force[i] += rotation * momentum_matrix_[p] * subspace_interface_force[i];
      net_torque[i] += rotation * torque_matrix_[p] * subspace_interface_force[i];
    }

    // force from full-subspace interface
    if (1) {
      //    for (int t : full_subspace_interface_tet_) {
      for (int i = 0; i < int(full_subspace_interface_tet_.size()); ++i) {
        int t = full_subspace_interface_tet_[i];
        int* vert = tet_ + t * 4;

        double* element_force = inv_fem_->element_force_ + t * 3 * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i = 0; i < 4; ++i) {
          int v = vert[i];
          int p = vert_part_id_[v];
          if (is_subspace_domain_[p]) {
            Vec3 map_force(element_force + i * 3);
            map_force *= rigid_coupling_scaling;
            Vec3 r = MapVec3(X + v * 3) - center_of_mass_[p];
            // vega force is in negative direction
            net_force[domain_index_[p]] -= map_force;
            net_torque[domain_index_[p]] -= r.cross(map_force);

            // implicit force
            if (0) {
              for (int j = 0; j < 4; ++j) {
                if (j == i) continue;
                int vj = vert[j];
                int pj = vert_part_id_[vj];
                if (!is_subspace_domain_[pj]) {
                  int local_v = local_vert_idx_[vj];
                  MapVec3 vel(&part_vel_[pj][local_v * 3]);
                  Vec3 d_force = element_k.block<3, 3>(i * 3, j * 3) * vel;
                  d_force *= dt;
                  net_force[domain_index_[p]] -= d_force;
                  net_torque[domain_index_[p]] -= r.cross(d_force);
                }
              }
            }
          }
        }
      }
    }
  }
  profiler.End("boundary force");


  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // compute rigid transformation
  std::vector<Vec3> rhs(subspace_domain_num_ * 2);
  OMP_FOR
  for (int i = 0; i < subspace_domain_num_; ++i) {
    int p = subspace_domains_[i];
    rhs[i * 2 + 0] = net_force[i] * dt + translational_vel_[p] * mass_per_part_[p];
    rhs[i * 2 + 1] = net_torque[i] * dt + current_inertia_tensor_[p] * angular_vel_[p];
  }

  profiler.Start("rigid k");
  MatCol Rigid = MatCol::Zero(subspace_domain_num_ * 6, subspace_domain_num_ * 6);
  // Compute Stiffness matrix for solving rigid motion
  {
    // Off-diagonal terms from subspace-subspace interface
    {
      OMP_FOR
      for (int e = 0; e < int(interface_domains_.size()); ++e) {
        int p0 = interface_domains_[e][1];
        int p1 = interface_domains_[e][2];
        if (interface_domains_[e][0] > 2 || !is_subspace_domain_[p0] || !is_subspace_domain_[p1]) continue;

        Mat df_dt[2][2] = {
          {Mat::Zero(part_basis_size_[p0], 3), Mat::Zero(part_basis_size_[p0], 3)},
          {Mat::Zero(part_basis_size_[p1], 3), Mat::Zero(part_basis_size_[p1], 3)}
        }; // translation

        Mat df_dr[2][2] = {
          {Mat::Zero(part_basis_size_[p0], 3), Mat::Zero(part_basis_size_[p0], 3)},
          {Mat::Zero(part_basis_size_[p1], 3), Mat::Zero(part_basis_size_[p1], 3)}
        }; // rotation

        for (CubaturePoint & cubature : interface_cubature_[e]) {
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
              df_dt[idx0][idx1] += vert_basis_transpose_[vi] *
                                   (weight * part_rotation_transpose_[pi] * element_k.block<3, 3>(i * 3, j * 3));
              df_dr[idx0][idx1] += vert_basis_transpose_[vi] *
                                   (weight * part_rotation_transpose_[pi] * (element_k.block<3, 3>(i * 3, j * 3) * skew_symmetrix_mat));

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
              //          ASSERT(f0.norm() < 1e10, P(f0.norm()));
              //          ASSERT(f1.norm() < 1e10, P(f1.norm()));
            }
          }
        }

        tmp_interface_rigid_k_[e * 8 + 0] = part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][0];
        tmp_interface_rigid_k_[e * 8 + 1] = part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][0];
        tmp_interface_rigid_k_[e * 8 + 2] = part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][0];
        tmp_interface_rigid_k_[e * 8 + 3] = part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][0];
        tmp_interface_rigid_k_[e * 8 + 4] = part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][1];
        tmp_interface_rigid_k_[e * 8 + 5] = part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][1];
        tmp_interface_rigid_k_[e * 8 + 6] = part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][1];
        tmp_interface_rigid_k_[e * 8 + 7] = part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][1];


        int idx0 = domain_index_[p0] * 6;
        int idx1 = domain_index_[p1] * 6;
        Rigid.block<3, 3>(idx0 + 0, idx1 + 0) += part_rotation_[p0] * momentum_matrix_[p0] * df_dt[0][1];
        Rigid.block<3, 3>(idx0 + 0, idx1 + 3) += part_rotation_[p0] * momentum_matrix_[p0] * df_dr[0][1];
        Rigid.block<3, 3>(idx0 + 3, idx1 + 0) += part_rotation_[p0] * torque_matrix_[p0]   * df_dt[0][1];
        Rigid.block<3, 3>(idx0 + 3, idx1 + 3) += part_rotation_[p0] * torque_matrix_[p0]   * df_dr[0][1];
        Rigid.block<3, 3>(idx1 + 0, idx0 + 0) += part_rotation_[p1] * momentum_matrix_[p1] * df_dt[1][0];
        Rigid.block<3, 3>(idx1 + 0, idx0 + 3) += part_rotation_[p1] * momentum_matrix_[p1] * df_dr[1][0];
        Rigid.block<3, 3>(idx1 + 3, idx0 + 0) += part_rotation_[p1] * torque_matrix_[p1]   * df_dt[1][0];
        Rigid.block<3, 3>(idx1 + 3, idx0 + 3) += part_rotation_[p1] * torque_matrix_[p1]   * df_dr[1][0];
      }

      OMP_FOR
      for (int n = 0; n < subspace_domain_num_; ++n) {
        int p = subspace_domains_[n];
        int idx = domain_index_[p] * 6;
        for (std::pair<int, int>& edge_idx : domain_incident_interface_[p]) {
          int& e = edge_idx.first;
          int i = e * 8 + edge_idx.second * 4;
          //      if (!is_subspace_domain_[interface_domains_[e][2 - edge_idx.second]]) continue;
          Rigid.block<3, 3>(idx + 0, idx + 0) += tmp_interface_rigid_k_[i + 0];
          Rigid.block<3, 3>(idx + 0, idx + 3) += tmp_interface_rigid_k_[i + 1];
          Rigid.block<3, 3>(idx + 3, idx + 0) += tmp_interface_rigid_k_[i + 2];
          Rigid.block<3, 3>(idx + 3, idx + 3) += tmp_interface_rigid_k_[i + 3];
        }
      }
    }

    // Off-diagonal terms from subspace-full interface
    if (1) {
      for (int i = 0; i < int(full_subspace_interface_tet_.size()); ++i) {
        int t = full_subspace_interface_tet_[i];
        int* vert = tet_ + t * 4;
        Eigen::Map<Eigen::Matrix<double, 12, 12> > element_k(inv_fem_->element_k_ + t * 144);
        for (int i4 = 0; i4 < 4; ++i4) {
          int vi = vert[i4];
          int pi = vert_part_id_[vi];
          if (!is_subspace_domain_[pi]) continue;
          int idx0 = domain_index_[pi] * 6;
          Vec3 ri = MapVec3(X + vi * 3) - center_of_mass_[pi];
          Mat3 skew_symmetrix_mat_i = GetSkewSymmetrixMatrix(ri);
          for (int j4 = 0; j4 < 4; ++j4) {
            int vj = vert[j4];
            int pj = vert_part_id_[vj];
            if (!is_subspace_domain_[pj]) continue;
            int idx1 = domain_index_[pj] * 6;
            Vec3 rj = MapVec3(X + vj * 3) - center_of_mass_[pj];
            Mat3 skew_symmetrix_mat_j = GetSkewSymmetrixMatrix(rj);
            Rigid.block<3, 3>(idx0 + 0, idx1 + 0) += element_k.block<3, 3>(i4 * 3, j4 * 3);
            Rigid.block<3, 3>(idx0 + 0, idx1 + 3) -= element_k.block<3, 3>(i4 * 3, j4 * 3) * skew_symmetrix_mat_j;
            Rigid.block<3, 3>(idx0 + 3, idx1 + 0) += skew_symmetrix_mat_i * element_k.block<3, 3>(i4 * 3, j4 * 3);
            Rigid.block<3, 3>(idx0 + 3, idx1 + 3) -= skew_symmetrix_mat_i * element_k.block<3, 3>(i4 * 3, j4 * 3) * skew_symmetrix_mat_j;
          }
        }
      }
    }

    //    PMATCOL(Rigid);
    //    PVEC(rhs[0]);
    Rigid *= (dt * dt) * rigid_coupling_scaling;

    // Diagonal terms
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      int idx = i * 6;
      Rigid(idx + 0, idx + 0) += mass_per_part_[p];
      Rigid(idx + 1, idx + 1) += mass_per_part_[p];
      Rigid(idx + 2, idx + 2) += mass_per_part_[p];
      Rigid.block<3, 3>(idx + 3, idx + 3) += current_inertia_tensor_[p];
    }
  }
  profiler.End("rigid k");

  profiler.Start("rigid solve");
  MapVec map_rhs(&rhs[0][0], subspace_domain_num_ * 6);
  new_rigid_velocity = Rigid.colPivHouseholderQr().solve(map_rhs);
  profiler.End("rigid solve");
}

void MixedMultiDomainTet::RenderSurface() {
  UpdateRenderData();
#if 1
  if (textured_surface_renderer_.size() != 0) {
    glEnable(GL_NORMALIZE);
    for (int g = 0; g < int(groups_.size()) - 1; ++g) {
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
    //    if (0)
    if (groups_.back().size() != 0) {
      //      textured_surface_renderer_.back().Render();
      float specular[4] = {0.0, 0.0, 0.0, 1.0};
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 0);
      float diffuse[] = {138 / 255.0f, 7 / 255.0f, 7 / 255.0f, 1.0};
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
      glBegin(GL_TRIANGLES);
      for (int t : groups_.back()) {
        int* verts = T + t * 3;
        //        Normal(VN + verts[0] * 3);
        Normal(TN + t * 3);
        Vertex3v(X + verts[0] * 3);

        //        Normal(VN + verts[1] * 3);
        //        Normal(TN + t * 3);
        Vertex3v(X + verts[1] * 3);

        //        Normal(VN + verts[2] * 3);
        //        Normal(TN + t * 3);
        Vertex3v(X + verts[2] * 3);
      }
      glEnd();
    }
  } else if (surface_renderer_ != nullptr) {
    surface_renderer_->Render();
  } else {
    glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (int t = 0; t < triangle_num_; t++) {
      int* verts = T + t * 3;
      Normal(VN + verts[0] * 3);
      Vertex3v(X + verts[0] * 3);


      Normal(VN + verts[1] * 3);
      Vertex3v(X + verts[1] * 3);

      Normal(VN + verts[2] * 3);
      Vertex3v(X + verts[2] * 3);
    }
    glEnd();
  }
#else
  glEnable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  for (int t = 0; t < triangle_num_; t++) {
    int* verts = &surface_triangle_indices_[0] + t * 3;
    //    ASSERT(verts[0] >= 0 && verts[0] < surface_vertex_num_, P(verts[0], surface_vertex_num_));
    //    ASSERT(verts[1] >= 0 && verts[1] < surface_vertex_num_, P(verts[1], surface_vertex_num_));
    //    ASSERT(verts[2] >= 0 && verts[2] < surface_vertex_num_, P(verts[2], surface_vertex_num_));
    Normal(&surface_vert_normal_[verts[0] * 3]);
    Vertex3v(&surface_vert_pos_[verts[0] * 3]);

    Normal(&surface_vert_normal_[verts[1] * 3]);
    Vertex3v(&surface_vert_pos_[verts[1] * 3]);

    Normal(&surface_vert_normal_[verts[2] * 3]);
    Vertex3v(&surface_vert_pos_[verts[2] * 3]);
  }
  glEnd();
  if (0) {
    glDisable(GL_LIGHTING);
    glPointSize(6.0);
    glColor3f(0, 0, 0);
    glBegin(GL_POINTS);
    //  for (int v : surface_vertices_) {
    //    Vertex3v(X + v * 3);
    //  }
    for (int i = 0; i < surface_vertex_num_; ++i) {
      //    Vertex3v(&surface_vert_pos_[i * 3]);
      Vertex3v(X + surface_vertices_[i] * 3);
    }
    glEnd();
  }
#endif
  //  exit(0);
}


int MixedMultiDomainTet::Intersect(double * ray_start, double * ray_end, double * clicked_world_pos, double * selected_pos) {
  (void) clicked_world_pos;
  double start[3] = {ray_start[0], ray_start[1], ray_start[2]};
  double end[3] = {ray_end[0], ray_end[1], ray_end[2]};
  double dir[3] = {
    ray_end[0] - ray_start[0],
    ray_end[1] - ray_start[1],
    ray_end[2] - ray_start[2],
  };
  const double kThreshold = 0.001;
  for (int tri_idx = 0; tri_idx < triangle_num_; ++tri_idx) {
    if (dj::Dot3(dir, TN + tri_idx * 3) >= 0) continue;
    double u, v, t;
    bool intersect = dj::SegmentTriangleIntersection<double>(start, end,
                                                             X + T[tri_idx * 3 + 0] * 3,
                                                             X + T[tri_idx * 3 + 1] * 3,
                                                             X + T[tri_idx * 3 + 2] * 3,
                                                             kThreshold, u, v, t);
    if (intersect) {
      selected_pos[0] = ray_start[0] * (1 - t) + ray_end[0] * t;
      selected_pos[1] = ray_start[1] * (1 - t) + ray_end[1] * t;
      selected_pos[2] = ray_start[2] * (1 - t) + ray_end[2] * t;
      return tri_idx;
    }
  }
  return -1;
}

void MixedMultiDomainTet::EnalbeTexturedSurfaceRendering(const char* shader_file) {
  surface_vert_pos_.resize(kMaxTriNum * 3 * 3);
  surface_vert_normal_.resize(kMaxTriNum * 3 * 3);
  surface_vert_texture_coord_.resize(kMaxTriNum * 3 * 2);
  tex_coord_.resize(kMaxTriNum * 3 * 2);
  int shader = SetupGLSL(shader_file);
  for (int i = 0; i < int(groups_.size()); ++i) {
    textured_surface_renderer_.emplace_back(shader);
  }

  triangle_group_info_.resize(kMaxTriNum);
  group_triangle_offset_.resize(groups_.size() + 2);
  group_triangle_offset_.back() = group_triangle_offset_[group_triangle_offset_.size() - 2];
  groups_.resize(groups_.size() + 1); // one more group for cut face
  groups_.back().reserve(10000);
  triangle_texture_coordinate_ = tex_coord_;
  for (int g = 0; g < int(groups_.size()) - 1; ++g) {
    int offset = group_triangle_offset_[g] * 6;
    double* group_tex_coord = &tex_coord_[offset];
    groups_[g].reserve(groups_[g].size() * 4);
    for (int i = 0; i < int(groups_[g].size()); ++i) {
      int t = groups_[g][i];
      triangle_group_info_[t].first = g;
      triangle_group_info_[t].second = i;
      for (int i = 0; i < 3; ++i) {
        tex_coord_[offset + i * 2 + 0] = triangle_texture_coordinate_[t * 6 + i * 2 + 0];
        tex_coord_[offset + i * 2 + 1] = triangle_texture_coordinate_[t * 6 + i * 2 + 1];
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

void MixedMultiDomainTet::UpdateRenderData() {
  if (position_changed_) {
    UpdateEmbededVertexPosition();
    Build_TN();
    if (textured_surface_renderer_.size() != 0) {
      for (int t = 0; t < triangle_num_; t++) {
        if (is_cut_triangle_[t]) continue;
        int* verts = T + t * 3;
        VN[verts[0] * 3 + 0] += TN[t * 3 + 0];
        VN[verts[0] * 3 + 1] += TN[t * 3 + 1];
        VN[verts[0] * 3 + 2] += TN[t * 3 + 2];

        VN[verts[1] * 3 + 0] += TN[t * 3 + 0];
        VN[verts[1] * 3 + 1] += TN[t * 3 + 1];
        VN[verts[1] * 3 + 2] += TN[t * 3 + 2];

        VN[verts[2] * 3 + 0] += TN[t * 3 + 0];
        VN[verts[2] * 3 + 1] += TN[t * 3 + 1];
        VN[verts[2] * 3 + 2] += TN[t * 3 + 2];
      }
      UpdateTexturedSurfaceMeshVBO();
    }  else if (surface_renderer_ != nullptr) {
      for (int t = 0; t < triangle_num_; t++) {
        int* verts = T + t * 3;
        VN[verts[0] * 3 + 0] += TN[t * 3 + 0];
        VN[verts[0] * 3 + 1] += TN[t * 3 + 1];
        VN[verts[0] * 3 + 2] += TN[t * 3 + 2];

        VN[verts[1] * 3 + 0] += TN[t * 3 + 0];
        VN[verts[1] * 3 + 1] += TN[t * 3 + 1];
        VN[verts[1] * 3 + 2] += TN[t * 3 + 2];

        VN[verts[2] * 3 + 0] += TN[t * 3 + 0];
        VN[verts[2] * 3 + 1] += TN[t * 3 + 1];
        VN[verts[2] * 3 + 2] += TN[t * 3 + 2];
      }
      UpdateSurfaceMeshVBO();
    }
    position_changed_ = false;
  }
}

void MixedMultiDomainTet::UpdateTextureCoordinateVBO() {
  if (textured_surface_renderer_.size() == 0) return;
  group_triangle_offset_[0] = 0;
  for (int i = 1; i < int(group_triangle_offset_.size()); ++i) {
    group_triangle_offset_[i] = group_triangle_offset_[i - 1] + groups_[i - 1].size();
  }
  for (int g = 0; g < int(groups_.size()) - 1; ++g) {
    const ObjMesh::Group* group = obj_->getGroupHandle(g);
    auto materialHandle = obj_->getMaterialHandle(group->getMaterialIndex());
    //    ASSERT(glGetError() == GL_NO_ERROR);
    if (!materialHandle->hasTextureFilename()) { continue; }
    int offset = group_triangle_offset_[g] * 6;
    double* group_tex_coord = &tex_coord_[offset];
    for (int i = 0; i < int(groups_[g].size()); ++i) {
      int t = groups_[g][i];
      for (int i = 0; i < 3; ++i) {
        tex_coord_[offset + i * 2 + 0] = triangle_texture_coordinate_[t * 6 + i * 2 + 0];
        tex_coord_[offset + i * 2 + 1] = triangle_texture_coordinate_[t * 6 + i * 2 + 1];
      }
      offset += 6;
    }
    textured_surface_renderer_[g].UpdateTexCoord(group_tex_coord, int(groups_[g].size()));
  }
}

void MixedMultiDomainTet::TetLimiting(double dt, int iteration) {
  // store old pos
  memcpy(tmp_X, X, sizeof(double) * vertex_num_ * 3);

  for (ExtForce & f : ext_force_) {
    int v = f.first;
    int p = vert_part_id_[v];
    if (is_subspace_domain_[p]) continue;
    MapVec3(velocity_ + v * 3) += f.second * (dt / mass_[v]);
  }

  // apply gravity and update initial position
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    if (is_subspace_domain_[p]) continue;
    MapVec3(velocity_ + v * 3) += (gravity_ * dt);
    MapVec3(X + v * 3) += MapVec3(velocity_ + v * 3) * dt;
  }

  for (int l = 0; l < iteration; ++l) {
    memset(tmp_vertex_pos_, 0, sizeof(vertex_num_) * 3);
    for (unsigned int i = 0; i < full_tet_.size() + full_subspace_interface_tet_.size(); ++i) {
      int t = (i < full_tet_.size()) ? full_tet_[i] : full_subspace_interface_tet_[i - full_tet_.size()];
      int t9 = t * 9;
      const int p0 = tet_[t * 4 + 0] * 3;
      const int p1 = tet_[t * 4 + 1] * 3;
      const int p2 = tet_[t * 4 + 2] * 3;
      const int p3 = tet_[t * 4 + 3] * 3;
      typedef Eigen::Matrix3d Mat3Col;
      double Ds[3][3] = {
        X[p1 + 0] - X[p0 + 0],
        X[p2 + 0] - X[p0 + 0],
        X[p3 + 0] - X[p0 + 0],

        X[p1 + 1] - X[p0 + 1],
        X[p2 + 1] - X[p0 + 1],
        X[p3 + 1] - X[p0 + 1],

        X[p1 + 2] - X[p0 + 2],
        X[p2 + 2] - X[p0 + 2],
        X[p3 + 2] - X[p0 + 2]
      };
      Mat3 F;
      dj::MulMatrix3x3<double>(Ds, inv_Dm + t9, F.data());
      Eigen::JacobiSVD<Mat3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::Vector3d S = svd.singularValues();
      S(0) = 1;
      S(1) = 1;
      S(2) = 1;
      if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) S(2) = -S(2);
      F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
      double new_ds[9];
      dj::MulMatrix3x3<double>(F.data(), Dm + t9, new_ds);

      double shift[3];
      shift[0] = (X[p0 + 0] + X[p1 + 0] + X[p2 + 0] + X[p3 + 0]
                  - new_ds[0] - new_ds[1] - new_ds[2]) * 0.25;
      shift[1] = (X[p0 + 1] + X[p1 + 1] + X[p2 + 1] + X[p3 + 1]
                  - new_ds[3] - new_ds[4] - new_ds[5]) * 0.25;
      shift[2] = (X[p0 + 2] + X[p1 + 2] + X[p2 + 2] + X[p3 + 2]
                  - new_ds[6] - new_ds[7] - new_ds[8]) * 0.25;

      tmp_vertex_pos_[p0 + 0] += shift[0] - X[p0 + 0];
      tmp_vertex_pos_[p0 + 1] += shift[1] - X[p0 + 1];
      tmp_vertex_pos_[p0 + 2] += shift[2] - X[p0 + 2];

      tmp_vertex_pos_[p1 + 0] += (shift[0] + new_ds[0]) - X[p1 + 0];
      tmp_vertex_pos_[p1 + 1] += (shift[1] + new_ds[3]) - X[p1 + 1];
      tmp_vertex_pos_[p1 + 2] += (shift[2] + new_ds[6]) - X[p1 + 2];

      tmp_vertex_pos_[p2 + 0] += (shift[0] + new_ds[1]) - X[p2 + 0];
      tmp_vertex_pos_[p2 + 1] += (shift[1] + new_ds[4]) - X[p2 + 1];
      tmp_vertex_pos_[p2 + 2] += (shift[2] + new_ds[7]) - X[p2 + 2];

      tmp_vertex_pos_[p3 + 0] += (shift[0] + new_ds[2]) - X[p3 + 0];
      tmp_vertex_pos_[p3 + 1] += (shift[1] + new_ds[5]) - X[p3 + 1];
      tmp_vertex_pos_[p3 + 2] += (shift[2] + new_ds[8]) - X[p3 + 2];
    } // For each tet

    for (int v = 0; v < vertex_num_; ++v) {
      int p = vert_part_id_[v];
      if (is_subspace_domain_[p]) continue;
      const double kWeight = 0.01;
      MapVec3 tmp_pos(tmp_vertex_pos_ + v * 3);
      tmp_pos *= kWeight;
      MapVec3(X + v * 3) += tmp_pos;
    }
  }
  for (int v = 0; v < vertex_num_; ++v) {
    int p = vert_part_id_[v];
    if (is_subspace_domain_[p]) continue;
    velocity_[v * 3 + 0] = (X[v * 3 + 0] - tmp_X[v * 3 + 0]) / dt;
    velocity_[v * 3 + 1] = (X[v * 3 + 1] - tmp_X[v * 3 + 1]) / dt;
    velocity_[v * 3 + 2] = (X[v * 3 + 2] - tmp_X[v * 3 + 2]) / dt;
  }
} // For each iteration

inline void MixedMultiDomainTet::UpdateEmbededVertexPosition() {
  OMP_FOR
  for (int i = 0; i < embeded_v_num_; ++i) {
    int v = embeded_vertex_[i];
    int& v0 = embeded_edge_[i].v0;
    int& v1 = embeded_edge_[i].v1;
    double& weight = embeded_edge_[i].weight;
    X[v * 3 + 0] = weight * (X[v0 * 3 + 0]) + (1 - weight) * (X[v1 * 3 + 0]);
    X[v * 3 + 1] = weight * (X[v0 * 3 + 1]) + (1 - weight) * (X[v1 * 3 + 1]);
    X[v * 3 + 2] = weight * (X[v0 * 3 + 2]) + (1 - weight) * (X[v1 * 3 + 2]);
    //    P(dj::Vec3d(X + v * 3), weight, v0, v1);
  }
  //  exit(0);
}

void MixedMultiDomainTet::UpdateBasisOffSet(bool using_pbd) {
  basis_offset_[0] = 0;
  total_basis_num_ = 0;
  subspace_domain_basis_size_ = 0;
  // WARNING: old code uses this version
#if 0
  for (int p = 0; p < part_num_; ++p) {
    // full domain
    if (is_subspace_domain_[p] == false) {
      part_basis_size_[p] = (using_pbd) ? 0 : vert_num_per_part_[p] * 3;
    } else {
      subspace_domain_basis_size_ += part_basis_size_[p];
    }
    total_basis_num_ += part_basis_size_[p];
    basis_offset_[p + 1] = basis_offset_[p] + part_basis_size_[p];
  }
#else
  (void) using_pbd;
  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p] == false) {
      part_basis_size_[p] = vert_num_per_part_[p] * 3;
      basis_offset_[p + 1] = basis_offset_[p];
    } else {
      subspace_domain_basis_size_ += part_basis_size_[p];
      basis_offset_[p + 1] = basis_offset_[p] + part_basis_size_[p];
    }
    total_basis_num_ += part_basis_size_[p];
  }
#endif
}

void MixedMultiDomainTet::UpdateSurfaceTriangleTopology() {
  if (surface_renderer_ == nullptr) return;
  memset(is_surface_vert_ - embeded_v_num_, 0, sizeof(char) * (embeded_v_num_ + vertex_num_));
  OMP_FOR
  for (int t = 0; t < triangle_num_; ++t) {
    ASSERT(T[t * 3 + 0] >= -embeded_v_num_ && T[t * 3 + 0] < vertex_num_);
    ASSERT(T[t * 3 + 1] >= -embeded_v_num_ && T[t * 3 + 1] < vertex_num_);
    ASSERT(T[t * 3 + 2] >= -embeded_v_num_ && T[t * 3 + 2] < vertex_num_);
    is_surface_vert_[T[t * 3 + 0]] = 1;
    is_surface_vert_[T[t * 3 + 1]] = 1;
    is_surface_vert_[T[t * 3 + 2]] = 1;
  }
  surface_vertex_num_ = 0;
  for (int i = -embeded_v_num_; i < vertex_num_; ++i) {
    if (is_surface_vert_[i]) {
      global_vertex_id2surface_vertex_id_[i] = surface_vertex_num_;
      surface_vertices_[surface_vertex_num_] = i;
      surface_vertex_num_++;
    }
  }
  surface_triangle_indices_.resize(triangle_num_ * 3);
  OMP_FOR
  for (int t = 0; t < triangle_num_; ++t) {
    surface_triangle_indices_[t * 3 + 0] = global_vertex_id2surface_vertex_id_[T[t * 3 + 0]];
    surface_triangle_indices_[t * 3 + 1] = global_vertex_id2surface_vertex_id_[T[t * 3 + 1]];
    surface_triangle_indices_[t * 3 + 2] = global_vertex_id2surface_vertex_id_[T[t * 3 + 2]];
  }
  surface_renderer_->UpdateTriangles(triangle_num_, &surface_triangle_indices_[0]);
  surface_renderer_->set_vertex_num(surface_vertex_num_);
  //  for (int i = 0; i < surface_vertex_num_; ++i) {
  //    surface_vert_pos_[i * 3 + 0] = X[surface_vertices_[i] * 3 + 0];
  //    surface_vert_pos_[i * 3 + 1] = X[surface_vertices_[i] * 3 + 1];
  //    surface_vert_pos_[i * 3 + 2] = X[surface_vertices_[i] * 3 + 2];
  //  }
}

void MixedMultiDomainTet::BuildTetSurfaceTriIndex() {
  surface_tri_on_tet_.resize(kMaxTetNum);
  std::map<dj::Vec3i, int> vert2tri;
  for (int t = 0; t < triangle_num_; ++t) {
    dj::Vec3i verts(T[t * 3 + 0], T[t * 3 + 1], T[t * 3 + 2]);
    std::sort(&verts[0], &verts[0] + 3);
    vert2tri[verts] = t;
  }

  for (int t = 0; t < tet_number; ++t) {
    for (int i = 0; i < 4; ++i) {
      dj::Vec3i verts;
      // get the other 3 vertices
      for (int j = 0, count = 0; j < 4; ++j) {
        if (j == i) continue;
        verts[count] = tet_[t * 4 + j];
        ++count;
      }
      std::sort(&verts[0], &verts[0] + 3);
      if (vert2tri.count(verts) > 0) {
        surface_tri_on_tet_[t][i] = vert2tri[verts];
      } else {
        surface_tri_on_tet_[t][i] = INT_MIN;
      }
    }
  }
}

inline void MixedMultiDomainTet::SplitTriangle(int tri_idx, int v0, int v1, int v2, int v10, int v11, int v20, int v21) {
  if (textured_surface_renderer_.size() != 0) {
    double old_tex_coord[3][2] = {DBL_MAX, DBL_MAX, DBL_MAX,
                                  DBL_MAX, DBL_MAX, DBL_MAX
                                 };
    const int verts[3] = {v0, v1, v2};
    int* old_tri_verts = T + tri_idx * 3;
    for (int ii = 0; ii < 3; ++ii) {
      int jj = 0;
      for (; jj < 3; ++jj) {
        if (old_tri_verts[jj] == verts[ii]) {
          old_tex_coord[ii][0] = triangle_texture_coordinate_[tri_idx * 6 + jj * 2 + 0];
          old_tex_coord[ii][1] = triangle_texture_coordinate_[tri_idx * 6 + jj * 2 + 1];
          break;
        }
      }
      ASSERT(jj < 3);
    }
    double edge_weight[2] = {
      embeded_edge_[-v10 - 1].weight,
      embeded_edge_[-v20 - 1].weight,
    };
    ASSERT((vertex_id_[embeded_edge_[-v10 - 1].v0] == v0 && vertex_id_[embeded_edge_[-v10 - 1].v1] == v1) ||
           (vertex_id_[embeded_edge_[-v10 - 1].v0] == v1 && vertex_id_[embeded_edge_[-v10 - 1].v1] == v0));
    ASSERT((vertex_id_[embeded_edge_[-v20 - 1].v0] == v0 && vertex_id_[embeded_edge_[-v20 - 1].v1] == v2) ||
           (vertex_id_[embeded_edge_[-v20 - 1].v0] == v2 && vertex_id_[embeded_edge_[-v20 - 1].v1] == v0));
    if (v0 != vertex_id_[embeded_edge_[-v10 - 1].v0]) edge_weight[0] = 1 - edge_weight[0];
    if (v0 != vertex_id_[embeded_edge_[-v20 - 1].v0]) edge_weight[1] = 1 - edge_weight[1];
    double new_tex_coord[2][2] = {
      {
        old_tex_coord[0][0] * edge_weight[0] + (1 - edge_weight[0]) * old_tex_coord[1][0],
        old_tex_coord[0][1] * edge_weight[0] + (1 - edge_weight[0]) * old_tex_coord[1][1]
      },
      {
        old_tex_coord[0][0] * edge_weight[1] + (1 - edge_weight[1]) * old_tex_coord[2][0],
        old_tex_coord[0][1] * edge_weight[1] + (1 - edge_weight[1]) * old_tex_coord[2][1]
      }
    };
    int tri = tri_idx;
    triangle_texture_coordinate_[tri * 6 + 0] = old_tex_coord[0][0];
    triangle_texture_coordinate_[tri * 6 + 1] = old_tex_coord[0][1];
    triangle_texture_coordinate_[tri * 6 + 2] = new_tex_coord[0][0];
    triangle_texture_coordinate_[tri * 6 + 3] = new_tex_coord[0][1];
    triangle_texture_coordinate_[tri * 6 + 4] = new_tex_coord[1][0];
    triangle_texture_coordinate_[tri * 6 + 5] = new_tex_coord[1][1];

    int group = triangle_group_info_[tri_idx].first;
    tri = triangle_num_;
    triangle_texture_coordinate_[tri * 6 + 0] = old_tex_coord[1][0];
    triangle_texture_coordinate_[tri * 6 + 1] = old_tex_coord[1][1];
    triangle_texture_coordinate_[tri * 6 + 2] = new_tex_coord[1][0];
    triangle_texture_coordinate_[tri * 6 + 3] = new_tex_coord[1][1];
    triangle_texture_coordinate_[tri * 6 + 4] = new_tex_coord[0][0];
    triangle_texture_coordinate_[tri * 6 + 5] = new_tex_coord[0][1];
    triangle_group_info_[tri].first = group;
    triangle_group_info_[tri].second = int(groups_[group].size());
    groups_[group].push_back(tri);

    tri++;
    triangle_texture_coordinate_[tri * 6 + 0] = old_tex_coord[1][0];
    triangle_texture_coordinate_[tri * 6 + 1] = old_tex_coord[1][1];
    triangle_texture_coordinate_[tri * 6 + 2] = old_tex_coord[2][0];
    triangle_texture_coordinate_[tri * 6 + 3] = old_tex_coord[2][1];
    triangle_texture_coordinate_[tri * 6 + 4] = new_tex_coord[1][0];
    triangle_texture_coordinate_[tri * 6 + 5] = new_tex_coord[1][1];
    triangle_group_info_[tri].first = group;
    triangle_group_info_[tri].second = int(groups_[group].size());
    groups_[group].push_back(tri);
  }

  //  KK;
  //  P(tri_idx, v0, v1, v2);
  //  P(v10, v11);
  //  P(v20, v21);
  T[tri_idx * 3 + 0] = v0;
  T[tri_idx * 3 + 1] = v10;
  T[tri_idx * 3 + 2] = v20;
  int idx = triangle_num_ * 3;
  T[idx + 0] = v1;
  T[idx + 1] = v21;
  T[idx + 2] = v11;

  T[idx + 3] = v1;
  T[idx + 4] = v2;
  T[idx + 5] = v21;
  triangle_num_ += 2;
  //  P(dj::Vec3i(tri_idx * 3 + T), tri_idx);
  //  P(dj::Vec3i(idx + T), idx / 3);
  //  P(dj::Vec3i(idx + 3 + T), idx / 3 + 1);
}

void MixedMultiDomainTet::SetFixedDomains(std::set<int> fixed_domains) {
  fixed_domains_ = std::vector<int>(fixed_domains.begin(), fixed_domains.end());
  std::sort(fixed_domains_.begin(), fixed_domains_.end());
  const unsigned int kNonConstrainedDof = subspace_domain_basis_size_ + subspace_domain_num_ * 6;
  const unsigned int kConstraiendDof = kNonConstrainedDof - int(fixed_domains_.size()) * 6;
  non_constrained_dof_ = kNonConstrainedDof;
  constrained_dof_ = kConstraiendDof;
  is_fixed_domain_.resize(part_num_);
  std::fill(is_fixed_domain_.begin(), is_fixed_domain_.end(), false);
  for (int fixed_p : fixed_domains) {
    ASSERT(fixed_p < part_num_ && fixed_p >= 0);
    ASSERT(is_subspace_domain_[fixed_p], L("try to set a full domained to be fixed"));
    is_fixed_domain_[fixed_p] = true;
  }

  non_constrained_dof_2_constrained_dof_.resize(kNonConstrainedDof);
  constrained_dof_2_non_constrained_dof_.resize(kConstraiendDof);
  for (int idx = 0, offset = 0; idx < subspace_domain_num_; ++idx) {
    int p = subspace_domains_[idx];
    for (int i = 0; i < part_basis_size_[p]; ++i) {
      non_constrained_dof_2_constrained_dof_[domain_offset_[p] + i] = offset;
      constrained_dof_2_non_constrained_dof_[offset] = domain_offset_[p] + i;
      ++offset;
    }
    if (is_fixed_domain_[p]) {
      for (int i = 0; i < 6; ++i) {
        non_constrained_dof_2_constrained_dof_[rigid_offset_[p] + i] = INT_MIN;
      }
    } else {
      for (int i = 0; i < 6; ++i) {
        non_constrained_dof_2_constrained_dof_[rigid_offset_[p] + i] = offset;
        constrained_dof_2_non_constrained_dof_[offset] = rigid_offset_[p] + i;
        offset++;
      }
    }
  }

  for (int p = 0; p < part_num_; ++p) {
    if (is_subspace_domain_[p]) {
      if (is_fixed_domain_[p]) {
        domain_size_[p] = part_basis_size_[p];
        constrained_rigid_offset_[p] = INT_MIN;
      } else {
        constrained_rigid_offset_[p] = non_constrained_dof_2_constrained_dof_[rigid_offset_[p]];
        domain_size_[p] = part_basis_size_[p] + 6;
      }
      constrained_domain_offset_[p] = non_constrained_dof_2_constrained_dof_[domain_offset_[p]];
    }
  }
  if (0) {
    std::vector<int> edges(chol_solver_->E, full_chol_solver_->E + full_chol_solver_->e_number0 * 2);
    std::vector<int> block_size(subspace_domain_num_);
    for (int i = 0; i < subspace_domain_num_; ++i) {
      int p = subspace_domains_[i];
      block_size[i] = domain_size_[p];
    }
    constrained_chol_solver_ = new solver::BLOCK_MATRIX_GRAPH<double>(subspace_domain_num_, edges, block_size);
  }
  delete sparse_k_;
  sparse_k_ = new MixedSparseMatrix(this);
  //      exit(0);
}

MixedMultiDomainTet::~MixedMultiDomainTet() {
  delete sparse_k_;
}

#include "open_gl_qt.h"
void MixedMultiDomainTet::Render() {
  //  MultiDomainTet::Render(Tet::kDefaultRendering, X); return;
  //  static int i = 0;
  //  if (i == 0) {
  //    P(tet_type_[12], dj::Vec4i(tet_ + 12 * 4));
  //    ++i;
  //  }
  if (0) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (int t = 0; t < tet_number; ++t) {
      //        for (int t : {156}) {
      //    int* verts = tet_ + t * 4;
      //    for (int i = 0; i < 4; ++i) {
      //      for (int j = i + 1; j < 4; ++j) {
      //        Vertex3v(X + verts[i] * 3);
      //        Vertex3v(X + verts[j] * 3);
      //      }
      //    }
      //    t = 1;
      //            if (tet_id_[t] != 503) continue;
      RenderTet(t);
    }
    glEnd();
  }
  // render edge
  if (0) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (int e = 0; e < edge_num_; ++e) {
      int* verts = edges_ + e * 2;
      if (verts[0] != 151 && verts[1] != 151 && verts[0] != 166 && verts[1] != 166) continue;
      Vertex3v(X + verts[0] * 3);
      Vertex3v(X + verts[1] * 3);
    }
    glEnd();
  }

  if (0) {
    RenderSurface();
  }
  glPointSize(6.0);
  glBegin(GL_POINTS);
  for (int v = 0; v < vertex_num_; ++v) {
    if (is_constrainted_[v]) {
      glColor3fv(kGreen());
    } else {
      glColor3fv(kRed());
    }
    Vertex3v(X + v * 3);
  }
  glEnd();
  // cut plane
  if (0) {
    Eigen::Vector3d a(-1.08, 1.55, 0.15);
    Eigen::Vector3d b(1.08, 1.55, 0.15);
    Eigen::Vector3d dir(0.0, -5, 0);
    Eigen::Vector3d c = b + dir;
    Eigen::Vector3d d = a + dir;
    glColor3fv(kBlue());
    glBegin(GL_QUADS);
    Vertex3v(&a[0]);
    Vertex3v(&b[0]);
    Vertex3v(&c[0]);
    Vertex3v(&d[0]);
    glEnd();

  }

  // vertex
  if (0) {
    QFont sansFont("Helvetica [Cronyx]", 18);
    sansFont.setFamily("sans serif");
    //    int t = 503;
    //    int tet_v[] = {
    //      tet_[t * 4 + 0],
    //      tet_[t * 4 + 1],
    //      tet_[t * 4 + 2],
    //      tet_[t * 4 + 3],
    //      155,
    //    };
    //    for (int v = 0; v < vertex_num_; ++v) {
    //          for (int v : tet_v) {
    for (int v : {26, 27}) {
      glColor3fv(kBlack());
      global::gl->renderText(
        X[v * 3 + 0], X[v * 3 + 1], X[v * 3 + 2],
        QString("%1").arg(v),
        sansFont
      );
    }

  }

  // tet
  if (0) {
    QFont sansFont("Helvetica [Cronyx]", 18);
    sansFont.setFamily("sans serif");
    for (int t = 0; t < tet_number; ++t) {
      //      if (tet_id_[t] != 4) continue;
      //    for (int t = original_tet_num_; t < tet_number; ++t) {
      int* verts = tet_ + t * 4;
      Vec3 center(0, 0, 0);
      for (int i = 0; i < 4; ++i) {
        center += MapVec3(&X[verts[i] * 3]);
      }
      center /= 4;
      glColor3fv(kOrage());
      global::gl->renderText(center[0], center[1], center[2],
                             QString("%1").arg(t), sansFont);
    }
  }

  // edge
  if (0) {
    QFont sansFont("Helvetica [Cronyx]", 18);
    sansFont.setFamily("sans serif");
    for (int e = 0; e < edge_num_; ++e) {
      int* verts = edges_ + e * 2;
      //      if (!(verts[0] == 151 || verts[1] == 151)) continue;
      Vec3 center = MapVec3(&X[verts[0] * 3]) + MapVec3(&X[verts[1] * 3]);
      center /= 2;
      glColor3fv(kOrage());
      global::gl->renderText(center[0], center[1], center[2],
                             QString("%1").arg(e), sansFont);
    }
  }

  //  Tet::Render();
}

void MixedMultiDomainTet::UpdateSurfaceMeshVBO() {
  if (surface_renderer_ == nullptr) return;
  OMP_FOR
  for (int i = 0; i < surface_vertex_num_; ++i) {
    int global_v = surface_vertices_[i];
    surface_vert_pos_[i * 3 + 0] = X[global_v * 3 + 0];
    surface_vert_pos_[i * 3 + 1] = X[global_v * 3 + 1];
    surface_vert_pos_[i * 3 + 2] = X[global_v * 3 + 2];
    surface_vert_normal_[i * 3 + 0] = VN[global_v * 3 + 0];
    surface_vert_normal_[i * 3 + 1] = VN[global_v * 3 + 1];
    surface_vert_normal_[i * 3 + 2] = VN[global_v * 3 + 2];
  }
  surface_renderer_->UpdateVertexVBO(&surface_vert_pos_[0], &surface_vert_normal_[0]);
}

void MixedMultiDomainTet::EnalbeSurfaceRenderingWithVBO(const char *shader_file) {

  surface_renderer_ = new TriangleMeshRender<double>(int(surface_vertices_.size()),
                                                     triangle_num_,
                                                     &surface_triangle_indices_[0],
                                                     shader_file);
  UpdateSurfaceTriangleTopology();
  UpdateSurfaceMeshVBO();
  position_changed_ = true;
}



#if 0
double ParallelogramSegmentIntersection(MultiDomainTet::Vec3 & pa, MultiDomainTet::Vec3 & dir0, MultiDomainTet::Vec3 & dir1,
                                        MultiDomainTet::Vec3 & normal,
                                        MultiDomainTet::MapVec3 & a, MultiDomainTet::MapVec3 & b) {
  typedef MultiDomainTet::Vec3 Vec3;
  Vec3 a2b = b - a;
  Vec3 pa2a = a - pa;
  double dots[2] = {
    normal.dot(pa2a),
    normal.dot(a2b),
  };
  //  PVEC(pa);
  //  PVEC(dir0);
  //  PVEC(dir1);
  //  PVEC(normal);
  //  PVEC(a);
  //  PVEC(b);
  //  P(dots[0]);
  //  P(dots[1]);
  if (dots[0] * dots[1] > 0 || dj::Abs(dots[1]) < 1e-8) return -1.0;
  double ratio = dots[0] / dots[1];
  //  P(ratio);
  if (ratio < -0.98 || ratio > -0.02) return -1.0; // out side of segment
  Vec3 intersection = a - ratio * a2b;
  //  PVEC(intersection)
  Vec3 pa2intersection = intersection - pa;
  dots[0] = pa2intersection.dot(dir0) / dj::Square(dir0.norm());
  dots[1] = pa2intersection.dot(dir1) / dj::Square(dir1.norm());
  const double kMinThreshold = -1e-5;
  const double kMaxThreshold = 1 + 1e-5;
  P(dots[0], dots[1], dir0.norm(), dir1.norm());
  if (dots[0] > kMinThreshold && dots[0] < kMaxThreshold && dots[1] > kMinThreshold && dots[1] < kMaxThreshold) {
    //    P(-ratio);
    return -ratio;
  } else {
    return -1;
  }
}
#endif

