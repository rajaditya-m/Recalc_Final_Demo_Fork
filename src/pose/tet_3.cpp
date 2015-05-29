#include "tet.h"
#include <unordered_set>
#include <set>
#include <numeric>
#include <iostream>
#include <fstream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <thread>
#include <set>
#include "string_formatter.h"
#include "VECTOR_TOOLS.h"
#include "opengl_helper.h"
#include "binary_file_io.h"
#include "fixed_vector.h"
#include "open_gl_qt.h"
#include "vector_io.h"
#include "config_file.h"
#include "timer.h"
#include "affine_transformer.h"
#include "skeleton_node.h"
#include "rainbow_color.h"
#include "buffer_manager.h"
#include "text_file_io.h"
#include "triangle_mesh_render.h"
#include "vector_lib.h"
#include "fixed_matrix_utility.h"
#include "textured_triangle_mesh_renderer.h"
#include "conjugate_gradient_solver.h"
#include "tet_gen_mesh_io.h"
#include "svd.h"
#include "global.h"
#include "skeleton.h"
#include "svd.h"
#include "tet.h"
#include "print_macro.h"
#include "opengl_helper.h"
#include "tet_mesh_simulator_bridge.h"
#include "pose_sampler.h"
extern  void SVDTest(double (*A)[3]);
extern  void ArmaSVD(double (*A)[3]);
using dj::Mat3d;
using dj::Wmat3d;
using dj::Vec3d;
typedef double(*Matrix3x3)[3];
const double kMinStiffness = 0.10f;
const double kMaxStiffness = 9.0f;
namespace MultiThread {
const static int kThreadNum = 4;
};

namespace EdgeLimitingParameter {
const static double kStiffness = .1;
const static double kMaxStretch = 1.0 + EdgeLimitingParameter::kStiffness;
const static double kMinStretch = 1.0 - EdgeLimitingParameter::kStiffness;
};

int GetRealIndex(int idx) {
  std::map<int, int> map;
  const int kCommon = 0;
  map[0] = kCommon;
  map[2] = kCommon;
  map[7] = kCommon;
  map[11] = kCommon;
  map[12] = kCommon;
  map[13] = kCommon;
  map[14] = kCommon;
  map[15] = kCommon;
  map[16] = kCommon;
  map[17] = kCommon;
  map[18] = kCommon;
  map[19] = kCommon;
  map[20] = kCommon;
  map[21] = kCommon;
  map[23] = kCommon;
  map[26] = kCommon;

  map[11] = 8;
  map[6] = 3;
  //  map[1] = 11;
  //  map[5] = 11;
  //  map[9] = 11;
  //  map[10] = 11;
  //  map[11] = 11;
  //  map[13] = 11;
  //  map[16] = 11;
  if (map.count(idx) > 0) {
    return map[idx];
  } else {
    return idx;
  }
}

std::set<int> inverted;

int Tet::Select(double* ray_start, double* ray_end,
                double* clicked_world_pos,
                double*selected_pos) {
  (void) clicked_world_pos;
  double start[3] = {ray_start[0], ray_start[1], ray_start[2]};
  double end[3] = {ray_end[0], ray_end[1], ray_end[2]};
  const double kThreshold = 0.001;
  int selected_vertex = -1;
  std::vector<std::pair<double, int> > intersect_triangles;
  for (int tri_idx = 0; tri_idx < triangle_num_; ++tri_idx) {
    double u, v, t;
    bool intersect = dj::SegmentTriangleIntersection<double>(start, end,
                                                             X + T[tri_idx * 3 + 0] * 3,
                                                             X + T[tri_idx * 3 + 1] * 3,
                                                             X + T[tri_idx * 3 + 2] * 3,
                                                             kThreshold, u, v, t);
    //    Real length;
    //    Real distance = dj::PointLineDistance<real>(start, end, &vert_[v][0], &length);
    if (intersect) {
      intersect_triangles.push_back(std::pair<double, int>(t, tri_idx));
    }
  }

  if (intersect_triangles.size() > 0) {
    std::sort(intersect_triangles.begin(), intersect_triangles.end());
    int tri_idx =  intersect_triangles.front().second;
    double min_dist[3] = {
      dj::PointLineDistance<double>(start, end, &X[T[tri_idx * 3  + 0] * 3]),
      dj::PointLineDistance<double>(start, end, &X[T[tri_idx * 3  + 1] * 3]),
      dj::PointLineDistance<double>(start, end, &X[T[tri_idx * 3  + 2] * 3]),
    };
    int min_idx = 0;
    if (min_dist[1] < min_dist[0]) {
      min_idx = 1;
    }
    if (min_dist[2] < min_dist[min_idx]) {
      min_idx = 2;
    }
    selected_vertex = T[tri_idx * 3 + min_idx];
    selected_pos[0] = X[selected_vertex * 3 + 0];
    selected_pos[1] = X[selected_vertex * 3 + 1];
    selected_pos[2] = X[selected_vertex * 3 + 2];
    selected_vertex_ = selected_vertex;
    return selected_vertex;
  } else {
    return -1;
  }
}

bool Tet::GetVertexPosition(double *pos) {
  if (IsSelected()) {
    pos[0] = X[ui_selected_vertex_ * 3 + 0];
    pos[1] = X[ui_selected_vertex_ * 3 + 1];
    pos[2] = X[ui_selected_vertex_ * 3 + 2];
    return true;
  } else {
    return false;
  }
}

int Tet::GetUIForce(double *force) {
  if (ui_selected_vertex_ >= 0 && mouse_moved_) {
    force[0] = last_mouse_pos_[0] - X[ui_selected_vertex_ * 3 + 0];
    force[1] = last_mouse_pos_[1] - X[ui_selected_vertex_ * 3 + 1];
    force[2] = last_mouse_pos_[2] - X[ui_selected_vertex_ * 3 + 2];
    force[0] *= ui_force_scaling_factor_;
    force[1] *= ui_force_scaling_factor_;
    force[2] *= ui_force_scaling_factor_;
    return ui_selected_vertex_;
  } else {
    return -1;
  }
}

void Tet::UpdateRenderData() {
  if (position_changed_) {
    Build_VN();
    if (surface_renderer_ != nullptr) {
      UpdateSurfaceMeshVBO();
    } else if (textured_surface_renderer_.size() != 0) {
      UpdateTexturedSurfaceMeshVBO();
    }
    position_changed_ = false;
  }
}


void Tet::UpdateVegaVertexOffset() {
  for (int v = 0; v < vertex_num_; ++v) {
    inv_fem_->u_[v][0] = X[v * 3 + 0] - rest_pos_[v * 3 + 0];
    inv_fem_->u_[v][1] = X[v * 3 + 1] - rest_pos_[v * 3 + 1];
    inv_fem_->u_[v][2] = X[v * 3 + 2] - rest_pos_[v * 3 + 2];
  }
}

void Tet::GetJointPositions(double *joint_pos) {
  for (int i = 0; i < int(joints_.size()) / 2; ++i) {
    Vec3d end_points[2] = {
      Vec3d(skeleton_->GetNode(joints_[i * 2 + 0].first)->world_pos()),
      Vec3d(skeleton_->GetNode(joints_[i * 2 + 1].first)->world_pos())
    };
    Vec3d pos = end_points[0] * joints_[i * 2].second + end_points[1] * joints_[i * 2 + 1].second;
    joint_pos[i * 3 + 0] = pos[0];
    joint_pos[i * 3 + 1] = pos[1];
    joint_pos[i * 3 + 2] = pos[2];
  }
}

Tet::Tet(const char *filename, double _limit_threshold, AffineTransformer<double>* affine_transformer, bool initialize_fem_module, bool load_interface)
  : vertex_num_(0) {
  P(filename);
  higher_sampler_ = NULL;
  lower_sampler_ = NULL;
  skeleton_ = NULL;
  tet_limit_threshold = _limit_threshold;
  edge_num_ = 0;
  std::vector<double> vert_pos;
  std::vector<int> tets;
  TetrahedralMeshIO* mesh_io = TetGenMeshIO::Instance();
  mesh_io->Read(filename, vert_pos, tets);
  mesh_io->OrientTetrahedron(vert_pos, tets);
  //    Read_Original_File(filename);
  vertex_num_ = int(vert_pos.size()) / 3;
  tet_number = int(tets.size()) / 4;
  if (affine_transformer) {
    affine_transformer->Transform(&vert_pos[0], vertex_num_);
  }
  AllocateVertexData(vertex_num_);
  AllocateTetData(tet_number);
  memcpy(X, &vert_pos[0], sizeof(double) * vertex_num_ * 3);
  memcpy(rest_pos_, X, sizeof(double) * vertex_num_ * 3);
  memcpy(tet_, &tets[0], sizeof(int) * tet_number * 4);
  // create one tet
  {
#if 0
    tet_number = 1;
    vertex_num_ = 4;
    int v;
    v = 0;
    X[v * 3 + 0] = 0;
    X[v * 3 + 1] = 0;
    X[v * 3 + 2] = 0;

    v = 1;
    X[v * 3 + 0] = 1;
    X[v * 3 + 1] = 0;
    X[v * 3 + 2] = 0;

    v = 2;
    X[v * 3 + 0] = 0;
    X[v * 3 + 1] = 1;
    X[v * 3 + 2] = 0;

    v = 3;
    X[v * 3 + 0] = 0;
    X[v * 3 + 1] = 0;
    X[v * 3 + 2] = 1;

    tet_[0] = 0;
    tet_[1] = 1;
    tet_[2] = 2;
    tet_[3] = 3;
#endif
  }
  memcpy(rest_pos_, X, sizeof(double) * vertex_num_ * 3);

  Build_Boundary_Triangles();
  BuildIncidentTet();
  ComputeRestDeformationGradientAndVolume();
  Build_Edge_From_Tets();

  selected_vertex_ = conf.Get<int>("selected_vertex");
  int tmp_pos_size = dj::Max(edge_num_ * 2 * 3, tet_number * 4 * 3);
  tmp_vertex_pos_ = buffer_.Malloc<double>(tmp_pos_size);


  per_vertex_weight_ = 1.0f / surface_vertices_.size();
  const double kStiffness = 2000000000;
  edge_stiffness_.resize(edge_num_);
  std::fill(edge_stiffness_.begin(), edge_stiffness_.end(), kStiffness);
  cg_solver_ = new ConjugateGradientSolver<double>(vertex_num_ * 3);
  // TODO change density
  //  density_ = 1000;
  //  density_ = 1100;
  density_ = conf.Get<double>("density");
  P(density_);
  render_mode_ = 0;
  InitializeMass();

  if(load_interface) {
      std::string ifn = conf.Get<std::string>("interface_file");
      std::string interface_file_path = dj::Format("%z/%z",GetDataFolder(),ifn.c_str());
      std::ifstream fileIn(interface_file_path.c_str());
      std::string line;
      std::istringstream iss;
      int v1,v2;
      while(std::getline(fileIn,line)) {
          iss.clear();
          iss.str("");
          iss.str(line);
          iss >> v1 >> v2;
          interface_vertices_.push_back(std::pair<int,int>(v1,v2));
      }
  }

  //  std::fill_n(mass_, vertex_num_, 1e-2f);
  if (initialize_fem_module) {
    inv_fem_ = new TetMeshSimulatorBridge(this,interface_vertices_);
  } else {
    inv_fem_ = NULL;
  }
  is_constrainted_ = std::vector<int>(vertex_num_, 0);
  P(triangle_num_, tet_number, vertex_num_, edge_num_);
  internal_force_scaling_factor_ = 1.0;
  ui_force_scaling_factor_ = conf.Get<double>("ui force scaling");
  surface_renderer_ = nullptr;
  Build_VN();
  position_changed_ = true;
}

Tet::~Tet() {
  delete surface_renderer_;
}

void Tet::MoveSkeletonJoint(int joint_idx, int angle_idx, double angle_change_in_radian) {
  GetJointPositions(&joint_vel_[0]);
  double* angle = skeleton_->GetNode(joint_idx)->rotation_angle();
  angle[angle_idx] += angle_change_in_radian;
  skeleton_->UpdatePosition();
  GetJointPositions(&joint_pos_[0]);
  for (int i = 0; i < (int) joint_vel_.size(); ++i) {
    joint_vel_[i] = (joint_pos_[i] - joint_vel_[i]) / global::time_step;
  }
}

void Tet::SaveConstrainedVertex(const char *file_name) {
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  out << vertex_num_ << std::endl;
  out << constrainted_vertex_.size() << std::endl;
  for (int i = 0; i < (int) constrainted_vertex_.size(); ++i) {
    out << constrainted_vertex_[i] << std::endl;
  }
  out.close();
}

void Tet::SaveBoneVertex(const char *file_name) {
  std::ofstream out(file_name);
  ASSERT(out.is_open());
  std::vector<int> bone_vertex;
  for (int v = 0; v < vertex_num_; ++v) {
    if (attached_joint_of_vertex_[v] >= 0) {
      bone_vertex.push_back(v);
    }
  }
  out << vertex_num_ << " " << bone_vertex.size() << std::endl;
  for (int i = 0; i < int(bone_vertex.size()); ++i) {
    out << bone_vertex[i] << std::endl;
  }
  out.close();
}



void Tet::SaveMeshPartition(const char *file_name) {
  std::ofstream out(file_name);
  out << vertex_num_ << std::endl;
  for (int v = 0; v < vertex_num_; ++v) {
    out << closest_bone_[v] << std::endl;
  }
  out.close();
}

void Tet::LoadMeshPartition(const char *file_name) {
  std::map<int, int> map;
  // map for 400k armadillo partition
  if (0) {
    map[2] = 3;
    map[9] = 2;
    map[19] = 1;
    map[16] = 5;
    map[4] = 7;
    map[11] = 6;
    map[36] = 8;
    map[22] = 17;
    map[25] = 14;
    map[34] = 13;
    map[37] = 9;
    map[28] = 11;
    map[30] = 16;
  }
  // map for 150k armadillo partition
  if (0) {
    map[0] = 3;
    map[3] = 2;
    map[5] = 1;
    map[4] = 5;
    map[2] = 6;
    map[1] = 7;
    map[9] = 13;
    map[7] = 14;
    map[6] = 17;
    map[8] = 16;
    map[12] = 9;
    map[11] = 8;
    map[10] = 11;
  }
  TextFileReader in(file_name);
  int v_num = 0;
  in >> v_num;
  ASSERT(v_num == vertex_num_, P(v_num, vertex_num_));
  for (int v = 0; v < vertex_num_; ++v) {
    in >> closest_bone_[v];
    //        ASSERT(map.count(closest_bone_[v]) > 0, P(closest_bone_[v]));
    //        closest_bone_[v] = map[closest_bone_[v]];
  }
  ReinitializeOffsetToSkeleton();
}

void Tet::EnalbeSurfaceRenderingWithVBO(const char *shader_source) {
  surface_renderer_ = new TriangleMeshRender<double>(int(surface_vertices_.size()),
                                                     triangle_num_,
                                                     &surface_triangle_indices_[0],
                                                     shader_source);
}

void Tet::EnalbeTexturedSurfaceRendering(const char *shader_source) {
  textured_surface_renderer_.emplace_back(shader_source);
  surface_vert_pos_.resize(triangle_num_ * 3 * 3);
  surface_vert_normal_.resize(triangle_num_ * 3 * 3);
  surface_vert_texture_coord_.resize(triangle_num_ * 3 * 2);
}

void Tet::Reset() {
  memcpy(X, rest_pos_, sizeof(double) * 3 * vertex_num_);
  memset(velocity_, 0, sizeof(double) * 3 * vertex_num_);
}


void Tet::UpdateSurfaceMeshVBO() {
  if (surface_renderer_ == nullptr) return;
  OMP_FOR
  for (int i = 0; i < int(surface_vertices_.size()); ++i) {
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

void Tet::UpdateTexturedSurfaceMeshVBO() {
  if (int(textured_surface_renderer_.size()) == 0) return;
  OMP_FOR
  for (int t = 0; t < triangle_num_; ++t) {
    int* verts = T + t * 3;
    for (int i = 0; i < 3; ++i) {
      surface_vert_pos_[t * 9 + i * 3 + 0] = X[verts[i] * 3 + 0];
      surface_vert_pos_[t * 9 + i * 3 + 1] = X[verts[i] * 3 + 1];
      surface_vert_pos_[t * 9 + i * 3 + 2] = X[verts[i] * 3 + 2];

      surface_vert_normal_[t * 9 + i * 3 + 0] = VN[verts[i] * 3 + 0];
      surface_vert_normal_[t * 9 + i * 3 + 1] = VN[verts[i] * 3 + 1];
      surface_vert_normal_[t * 9 + i * 3 + 2] = VN[verts[i] * 3 + 2];
    }
  }
  textured_surface_renderer_[0].UpdatePosAndNormal(&surface_vert_pos_[0], &surface_vert_normal_[0], triangle_num_);
}

void Tet::MassSpringSimulation(double dt) {
  using namespace dj;
  typedef double Real;
  const Real kGravity[3] = {0, -9.8, 0};
  (void) kGravity;
  std::vector<Real> force_gradient(edge_num_ * 9, 0);
  // Compute force gradient for each edge
  //   OMP_FOR
  for (int e = 0; e < edge_num_; ++e) {
    int e9 = e * 9;
    int* v = &edges_[e * 2];
    Real* pos[] = {
      &X[v[0] * 3],
      &X[v[1] * 3],
    };
    Real direction[3];
    dj::SubVec3(pos[1], pos[0], direction);
    Real length = dj::Normalize3(direction);
    if (length < rest_length_[e]) {
      Real (*matrix)[3] = (Real (*)[3]) &force_gradient[e9];
      matrix[0][0] = direction[0] * direction[0];
      matrix[0][1] = direction[0] * direction[1];
      matrix[0][2] = direction[0] * direction[2];

      matrix[1][0] = direction[1] * direction[0];
      matrix[1][1] = direction[1] * direction[1];
      matrix[1][2] = direction[1] * direction[2];

      matrix[2][0] = direction[2] * direction[0];
      matrix[2][1] = direction[2] * direction[1];
      matrix[2][2] = direction[2] * direction[2];
      {
        Real * tmp = &force_gradient[e9];
        for (int i = 0; i < 9; ++i) {
          tmp[i] *= edge_stiffness_[e];
        }
      }
    } else {
      Real factor = edge_stiffness_[e] * rest_length_[e] / length;
      Real (*matrix)[3] = (Real (*)[3]) &force_gradient[e9];
      matrix[0][0] = direction[0] * direction[0] - 1;
      matrix[0][1] = direction[0] * direction[1];
      matrix[0][2] = direction[0] * direction[2];

      matrix[1][0] = direction[1] * direction[0];
      matrix[1][1] = direction[1] * direction[1] - 1;
      matrix[1][2] = direction[1] * direction[2];

      matrix[2][0] = direction[2] * direction[0];
      matrix[2][1] = direction[2] * direction[1];
      matrix[2][2] = direction[2] * direction[2] - 1;
      {
        Real * tmp = &force_gradient[e9];
        for (int i = 0; i < 9; ++i) {
          tmp[i] *= factor;
        }
      }
      matrix[0][0] += edge_stiffness_[e];
      matrix[1][1] += edge_stiffness_[e];
      matrix[2][2] += edge_stiffness_[e];
    }
  }

  const Real t2 = dt * dt;
  Real damping_ = 0.02;
  const Real kDampingFactor = damping_ * dt;
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    memcpy(result, x, sizeof(Real) * vertex_num_ * 3);
    OMP_FOR
    for (int v = 0; v < vertex_num_; ++v) {
      if (attached_joint_of_vertex_[v] >= 0) {
        result[v * 3 + 0] = 0;
        result[v * 3 + 1] = 0;
        result[v * 3 + 2] = 0;
        continue;
      }
      int v3 = v * 3;
      const double factor = kDampingFactor / mass_[v];
      result[v3 + 0] += factor * x[v3 + 0];
      result[v3 + 1] += factor * x[v3 + 1];
      result[v3 + 2] += factor * x[v3 + 2];
      Real tmp_result[3] = {0, 0, 0};
      for (int num = 0; num < (int) incident_edge_[v].size(); ++num) {
        int e = incident_edge_[v][num];
        Real* fg = &force_gradient[e * 9];
        int* verts = &edges_[e * 2];
        Real* pos[] = {
          x + verts[0] * 3,
          x + verts[1] * 3,
        };
        Real diff[3];
        dj::SubVec3(pos[1], pos[0], diff);
        Real sign = (verts[0] == v) ? 1 : -1;
        Real tmp[3];
        dj::MulMatrix3x3Vec<Real>((Real (*)[3]) fg, diff, tmp);
        tmp_result[0] += tmp[0] * sign;
        tmp_result[1] += tmp[1] * sign;
        tmp_result[2] += tmp[2] * sign;
      }
      result[v3 + 0] += -t2 * tmp_result[0] / mass_[v];
      result[v3 + 1] += -t2 * tmp_result[1] / mass_[v];
      result[v3 + 2] += -t2 * tmp_result[2] / mass_[v];
    }
  };

  memset(tmp_X, 0, sizeof(Real) * 3 * vertex_num_);
  for (int joint = 0; joint < (int)attached_tet_of_joints_.size(); ++joint) {
    int v = attached_tet_of_joints_[joint];
    tmp_X[v * 3 + 0] = joint_vel_[joint * 3 + 0];
    tmp_X[v * 3 + 1] = joint_vel_[joint * 3 + 1];
    tmp_X[v * 3 + 2] = joint_vel_[joint * 3 + 2];
  }
  StiffnessMatrix(tmp_X, tmp_vertex_pos_);
  // Assemble right hand side
  //  memcpy(tmp_vertex_pos_, velocity_, sizeof(Real) * 3 * vertex_num_);
  //  OMP_FOR
  for (int v_idx = 0; v_idx < vertex_num_; ++v_idx) {
    if (attached_joint_of_vertex_[v_idx] >= 0) {
      tmp_vertex_pos_[v_idx * 3 + 0] = 0;
      tmp_vertex_pos_[v_idx * 3 + 1] = 0;
      tmp_vertex_pos_[v_idx * 3 + 2] = 0;
      continue;
    }
    Real* force = &tmp_vertex_pos_[v_idx * 3 + 0];
    force[0] += velocity_[v_idx * 3 + 0];
    force[1] += velocity_[v_idx * 3 + 1];
    force[2] += velocity_[v_idx * 3 + 2];
    //    force[0] += dt * kGravity[0];
    //    force[1] += dt * kGravity[1];
    //    force[2] += dt * kGravity[2];
    //    vel_[v_idx] *= kDamping_;
    //    force[0] += vel_[v_idx][0];
    //    force[1] += vel_[v_idx][1];
    //    force[2] += vel_[v_idx][2];

    for (int num = 0; num < (int) incident_edge_[v_idx].size(); ++num) {
      int e = incident_edge_[v_idx][num];
      int* verts = &edges_[e * 2 + 0];
      Real* pos[] = {
        &X[verts[0] * 3 + 0],
        &X[verts[1] * 3 + 0],
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      Real sign = (verts[1] == v_idx) ? 1 : -1;
      Real length = dj::Normalize3(diff);

      force[0] += -edge_stiffness_[e] * sign * diff[0] * (length - rest_length_[e]) * dt / mass_[v_idx];
      force[1] += -edge_stiffness_[e] * sign * diff[1] * (length - rest_length_[e]) * dt / mass_[v_idx];
      force[2] += -edge_stiffness_[e] * sign * diff[2] * (length - rest_length_[e]) * dt / mass_[v_idx];

      //      if (attached_joint_of_vertex_[verts[0]] >= 0) {
      //        Real diff[3] = { -pos[0][0], -pos[0][1], -pos[0][2]};
      //        //        dj::SubVec3(pos[1], pos[0], diff);
      //        Real sign = (verts[0] == v_idx) ? 1 : -1;
      //        ASSERT(sign == -1);
      //        Real tmp[3];
      //        Real* fg = &force_gradient[e * 9];
      //        dj::MulMatrix3x3Vec<Real>((real (*)[3]) fg, diff, tmp);
      //        force[0] -= tmp[0] * sign;
      //        force[1] -= tmp[1] * sign;
      //        force[2] -= tmp[2] * sign;
      //      }
      //      if (attached_joint_of_vertex_[verts[1]] >= 0) {
      //        Real diff[3] = {pos[1][0], pos[1][1], pos[1][2]};
      //        //        dj::SubVec3(pos[1], pos[0], diff);
      //        Real sign = (verts[0] == v_idx) ? 1 : -1;
      //        ASSERT(sign == 1);
      //        Real tmp[3];
      //        Real* fg = &force_gradient[e * 9];
      //        dj::MulMatrix3x3Vec<Real>((real (*)[3]) fg, diff, tmp);
      //        force[0] -= tmp[0] * sign;
      //        force[1] -= tmp[1] * sign;
      //        force[2] -= tmp[2] * sign;
      //      }
    }
  }
  cg_solver_->Solve(tmp_vertex_pos_, velocity_, StiffnessMatrix, 500, 1e-10);
  //  vel_[0] = Vec3(0, 0, 0);
  for (int i = 0; i < (int) attached_tet_of_joints_.size(); ++i) {
    int v = attached_tet_of_joints_[i];
    velocity_[v * 3 + 0] = joint_vel_[i * 3 + 0];
    velocity_[v * 3 + 1] = joint_vel_[i * 3 + 1];
    velocity_[v * 3 + 2] = joint_vel_[i * 3 + 2];
  }

  OMP_FOR
  for (int v_idx = 0; v_idx < vertex_num_; ++v_idx) {
    //    if (attached_joint_of_vertex_[v_idx] >= 0) {
    //      int joint = attached_joint_of_vertex_[v_idx];
    //      X[v_idx * 3 + 0] = joint_pos_[joint * 3 + 0];
    //      X[v_idx * 3 + 1] = joint_pos_[joint * 3 + 1];
    //      X[v_idx * 3 + 2] = joint_pos_[joint * 3 + 2];
    //      continue;;
    //    }
    X[v_idx * 3 + 0] += velocity_[v_idx * 3 + 0] * dt;
    X[v_idx * 3 + 1] += velocity_[v_idx * 3 + 1] * dt;
    X[v_idx * 3 + 2] += velocity_[v_idx * 3 + 2] * dt;
    //    vert_[v_idx] += vel_[v_idx] * dt;
    if (X[v_idx * 3 + 1] < 0) {
      velocity_[v_idx * 3 + 1] -= X[v_idx * 3 + 1] / dt;
      X[v_idx * 3 + 1] = 0;
    }
  }
}

void Tet::AllocateVertexData(int v_num) {
  X = buffer_.Malloc<double>(v_num * 3);
  velocity_ = buffer_.Malloc<double>(v_num * 3);
  rest_pos_ = buffer_.Malloc<double>(v_num * 3);
  mass_ = buffer_.Malloc<double>(v_num);
  tmp_X = buffer_.Malloc<double>(v_num * 3);
  strain_limiting_offset_ = buffer_.Malloc<double>(v_num * 3);
  old_X = buffer_.Malloc<double>(v_num * 3);
  weight = buffer_.Malloc<int>(v_num);
  closest_bone_ = buffer_.Malloc<int>(v_num, false);
  initial_offset_ = buffer_.Malloc<double>(v_num * 3);
  VN = buffer_.Malloc<double>(v_num * 3);
}

void Tet::AllocateTriangleData(int tri_num) {
  T = buffer_.Malloc<int>(tri_num * 3);
  triangle_vector_ = buffer_.Malloc<double>(tri_num * 6);
  triangle_bounding_box_ = buffer_.Malloc<double>(tri_num * 6);
  TN = buffer_.Malloc<double>(tri_num * 3);
  triangle_area_ = buffer_.Malloc<double>(tri_num);
}

void Tet::AllocateEdgeData(int e_num) {
  edges_ = buffer_.Malloc<int>(e_num * 2);
  rest_length_ = buffer_.Malloc<double>(e_num);
}

void Tet::AllocateTetData(int tet_num) {
  volume_ = buffer_.Malloc<double>(tet_num);
  tet_ = buffer_.Malloc<int>(tet_num * 4);
  tet_weights_ = buffer_.Malloc<double>(tet_num);

  Dm = buffer_.Malloc<double>(tet_num * 9);
  inv_Dm = buffer_.Malloc<double>(tet_num * 9);
  Area_Coeffs = buffer_.Malloc<double>(tet_num * 9);
}


void Tet::ComputeTriVecNormalAndBoundingBox(double threshold) {
  OMP_FOR
  for (int t = 0; t < triangle_num_; ++t) {
    int t3 = t * 3;
    int t6 = t * 6;
    double* tri_vec = triangle_vector_ + t6;
    double* box = triangle_bounding_box_ + t6;
    double* pos[] = {
      T[t3 + 0] * 3 + X,
      T[t3 + 1] * 3 + X,
      T[t3 + 2] * 3 + X
    };

    // Compute normal and triangle vector
    tri_vec[0] = pos[1][0] - pos[0][0];
    tri_vec[1] = pos[1][1] - pos[0][1];
    tri_vec[2] = pos[1][2] - pos[0][2];

    tri_vec[3] = pos[2][0] - pos[0][0];
    tri_vec[4] = pos[2][1] - pos[0][1];
    tri_vec[5] = pos[2][2] - pos[0][2];

    dj::Cross3(tri_vec + 0, tri_vec + 3, TN + t3);
    triangle_area_[t] = dj::Normalize3(TN + t3);

    // Building bounding box
    // min & max x
    box[0] = pos[0][0];
    box[1] = pos[0][0];
    // min & max y
    box[2] = pos[0][1];
    box[3] = pos[0][1];
    // min & max z
    box[4] = pos[0][2];
    box[5] = pos[0][2];

    for (int i = 0, idx = 0; i < 3; ++i, idx += 2) {
      if (pos[1][i] < box[idx]) box[idx] = pos[1][i];
      else if (pos[1][i] > box[idx + 1]) box[idx + 1] = pos[1][i];

      if (pos[2][i] < box[idx]) box[idx] = pos[2][i];
      else if (pos[2][i] > box[idx + 1]) box[idx + 1] = pos[2][i];

      box[idx] -= threshold;
      box[idx + 1] += threshold;
    }
  }
}

void Tet::Simulate(double time_step) {
  static int cnt = 0;
  if (0) {
    if (cnt < 1300) {
      int moving_node = 18;
      const double rotate_vel = dj::Degree2Radian(70.0f);
      MoveSkeletonJoint(moving_node, 1, -rotate_vel * time_step);
      cnt++;
    } else {
      L("done") ;
      std::fill(joint_vel_.begin(), joint_vel_.end(), 0);
    }
  }
  //  pose_sampler->Move();
  //  GetConstrainVertexVelocity();
  //  exit(0);
  //  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(time_step);
  inv_fem_->Simulate(time_step);
  //  MassSpringSimulation(time_step);
  //    EnforceSkeletonConstraint(time_step);
  return;
}


void Tet::ExportTetMesh(const char * file_name) {
  std::ofstream out(file_name, std::ios::binary);
  if (!out.is_open()) {
    std::cerr << "Tet::ExportTetMesh() => failed to open output file " << file_name << std::endl;
    return;
  }
  // Vertex data
  out.write((char*)&vertex_num_, sizeof(int));
  out.write((char*)X, sizeof(double)* vertex_num_ * 3);
  out.write((char*)&mass_, sizeof(double)* vertex_num_);
  dj::Write1DVectorBinary<int>(out, surface_vertices_);
  dj::Write2DVector<int>(out, adjacent_vertex_);
  dj::Write2DVector<int>(out, incident_triangle_);
  dj::Write2DVector<int>(out, incident_edge_);
  dj::Write2DVector<int>(out, incident_tet_);
  dj::Write2DVector<int>(out, rank_in_tet_);
  // Tet data
  out.write((char*)&tet_number, sizeof(int));
  out.write((char*)tet_, sizeof(int)* tet_number * 4);
  out.write((char*)volume_, sizeof(double)* tet_number);
  out.write((char*)Dm, sizeof(double)* tet_number * 9);
  out.write((char*)inv_Dm, sizeof(double)* tet_number * 9);
  // Triangle data
  out.write((char*)&triangle_num_, sizeof(int));
  out.write((char*)T, sizeof(int)* triangle_num_ * 3);
  out.write((char*)triangle_area_, sizeof(double)* triangle_num_);
  // Edge data
  out.write((char*)&edge_num_, sizeof(int));
  out.write((char*)edges_, sizeof(int)* edge_num_ * 2);
  out.write((char*)rest_length_, sizeof(double)* edge_num_);

  out.close();
  //    double min[3], max[3];
  //  dj::GetRange(X, vertex_num_,
  //               min[0],
  //               min[1],
  //               min[2],
  //               max[0],
  //               max[1],
  //               max[2]
  //               );
  //  P(Vec3f(min));
  //  P(Vec3f(max));

  //  out.write();
}

void Tet::ImportTetMesh(const char * file_name) {
  BinaryFileReader in(file_name);
  // Vertex data
  in.Read(&vertex_num_, 1);
  AllocateVertexData(vertex_num_);
  in.Read(X, vertex_num_ * 3);
  in.Read(mass_, vertex_num_);
  in.Read1DVector(surface_vertices_);
  in.Read2DVector(adjacent_vertex_);
  in.Read2DVector(incident_triangle_);
  in.Read2DVector(incident_edge_);
  in.Read2DVector(incident_tet_);
  in.Read2DVector(rank_in_tet_);
  // Tet data
  in.Read(&tet_number, 1);
  AllocateTetData(tet_number);
  in.Read(tet_, tet_number * 4);
  in.Read(volume_, tet_number);
  in.Read(Dm, tet_number * 9);
  in.Read(inv_Dm, tet_number * 9);
  // Triangle data
  in.Read(&triangle_num_, 1);
  AllocateTriangleData(triangle_num_);
  in.Read(T, triangle_num_ * 3);
  in.Read(triangle_area_, triangle_num_);
  // Edge data
  in.Read(&edge_num_, 1);
  AllocateEdgeData(edge_num_);
  in.Read(edges_, edge_num_ * 2);
  in.Read(rest_length_, edge_num_);
}

void Tet::InitializeMass() {
  //    std::fill_n(mass_, vertex_num_, 1);
  //    return;
  memset(mass_, 0, sizeof(double)* vertex_num_);
  //  double min_volume = 1e10;
  double total_volume = 0;
  for (int t = 0; t < tet_number; ++t) {
    mass_[tet_[t * 4 + 0]] += volume_[t] * density_ / 4;
    mass_[tet_[t * 4 + 1]] += volume_[t] * density_ / 4;
    mass_[tet_[t * 4 + 2]] += volume_[t] * density_ / 4;
    mass_[tet_[t * 4 + 3]] += volume_[t] * density_ / 4;
    //    min_volume = std::min(min_volume, volume_[t]);
    total_volume += volume_[t];
  }
  //  L("using uniform mass");
  //  std::fill_n(mass_, vertex_num_, 1.0 / vertex_num_);
  P(total_volume, total_volume * density_);
  if (0) {
    std::ofstream out("D:/mass.txt");
    out << vertex_num_ << std::endl;
    for (int v = 0; v < vertex_num_; ++v) {
      out << mass_[v] << "\n";
    }
    out.close();
    L("mass matrix exported");
    exit(0);
  }
  //  P(min_volume, total_volume, total_volume * density_);
  //  double total_mass = 0;
  //  double min_mass = 1e10;
  //  double max_mass = -100;
  //  for (int v = 0; v < vertex_num_; ++v) {
  //    total_mass += mass_[v];
  //    min_mass = std::min(min_mass, mass_[v]);
  //    max_mass = std::max(max_mass, mass_[v]);
  //  }
  //  P(total_mass);
  //  P(min_mass, max_mass);
  //  P(total_mass / vertex_num_);
}

std::vector<int> Tet::GetIncidentTetStat() {
  std::vector<int> stat(1000, 0);
  int max = 0;
  for (int v = 0; v < vertex_num_; ++v) {
    int num = (int) incident_tet_[v].size();
    max = dj::Max(max, num);
    ASSERT(num < 1000);
    stat[num]++;
  }
  stat.resize(max + 1);
  int sum = std::accumulate(stat.begin(), stat.end(), 0);
  ASSERT(sum == vertex_num_);
  return stat;
}

void Tet::FindBodyAttachedTetOfJoints(Skeleton *skeleton) {
  (void)skeleton;
  using namespace body;
  joints_.clear();
  for (int i = 0; i < skeleton_->Size(); ++i) {
    SkeletonNode* node = skeleton_->GetNode(i);
    SkeletonNode* parent = node->parent();
    if (node->parent() == NULL) continue;
    int parent_idx = parent - skeleton_->GetNode(0);
    double length = dj::Distance3(node->world_pos(), parent->world_pos());
    int segment_num = std::floor(length / kMinJointLength - 1e-8) + 1;
    segment_num = dj::Max(segment_num, 2);
    for (int j = 0; j < segment_num; j++) {
      double alpha = j * 1.0f / (segment_num - 1);
      double beta = 1 - alpha;
      joints_.push_back(std::make_pair(i, alpha));
      joints_.push_back(std::make_pair(parent_idx, beta));
    }
  }

  // remove duplicate joints
  std::set<Vec3d> unique_joints;
  for (int i = 0; i < int(joints_.size()) / 2;) {
    Vec3d end_points[2] = {
      Vec3d(skeleton_->GetNode(joints_[i * 2 + 0].first)->world_pos()),
      Vec3d(skeleton_->GetNode(joints_[i * 2 + 1].first)->world_pos())
    };
    Vec3d pos = end_points[0] * joints_[i * 2].second + end_points[1] * joints_[i * 2 + 1].second;
    if (unique_joints.count(pos) > 0) {
      joints_.erase(joints_.begin() + i * 2, joints_.begin() + (i + 1) * 2);
    } else {
      ++i;
    }
    unique_joints.insert(pos);
  }

  attached_joint_of_vertex_.resize(vertex_num_);
  std::fill(attached_joint_of_vertex_.begin(), attached_joint_of_vertex_.end(), -1);
  joint_vel_.resize((joints_.size() / 2) * 3);
  std::fill(joint_vel_.begin(), joint_vel_.end(), 0);

  joint_pos_.resize((joints_.size() / 2) * 3);
  std::fill(joint_pos_.begin(), joint_pos_.end(), 0);

  attached_tet_of_joints_.resize(joints_.size() / 2);
  barycentric_of_joint_.resize(attached_tet_of_joints_.size() * 4);
  is_constrainted_.resize(vertex_num_);
  std::fill(is_constrainted_.begin(), is_constrainted_.end(), 0);
  {
#if 0
    int i = 0;
    std::ofstream out(DATA_DIRECTORY  "joints.txt");
    out << unique_joints.size() << " 3 0 0" << std::endl;
    for (auto & p : unique_joints) {
      out << i << " "
          << p[0] << " "
          << p[1] << " "
          << p[2] << " " << std::endl;
      ++i;
    }
    out.close();
    exit(0);
#endif
  }
  std::set<int> constrained_vert_set_;
  // Find the tets that contain the joints
  OMP_FOR
  for (int i = 0; i < int(joints_.size()) / 2; ++i) {
    Vec3d end_points[2] = {
      Vec3d(skeleton_->GetNode(joints_[i * 2 + 0].first)->world_pos()),
      Vec3d(skeleton_->GetNode(joints_[i * 2 + 1].first)->world_pos())
    };
    Vec3d pos = end_points[0] * joints_[i * 2].second + end_points[1] * joints_[i * 2 + 1].second;
    attached_tet_of_joints_[i] = -1;
#if 0
    //    bool found = false;
    double min_distnace = 1e10;
    for (int v = 0; v < vertex_num_; ++v) {
      double dist = dj::Distance3(X + v * 3, &pos[0]);
      min_distnace = dj::Min(min_distnace, dist);
      if (dist < 1e-5) {
        attached_tet_of_joints_[i] = v;
        attached_joint_of_vertex_[v] = i;
        break;
      }
    }
    //    ASSERT(attached_vert_of_joints_[i] >= 0, P(min_distnace));
    //    P(attached_vert_of_joints_[i]);
#else
    double double_pos[4] = {pos[0], pos[1], pos[2]};
    for (int t = 0; t < tet_number; ++t) {
      double* v[] = {
        X + tet_[t * 4 + 0] * 3,
        X + tet_[t * 4 + 1] * 3,
        X + tet_[t * 4 + 2] * 3,
        X + tet_[t * 4 + 3] * 3,
      };
      double double_v[4][3] = {
        {v[0][0], v[0][1], v[0][2]},
        {v[1][0], v[1][1], v[1][2]},
        {v[2][0], v[2][1], v[2][2]},
        {v[3][0], v[3][1], v[3][2]}
      };
      double barycentric[4];
      //      if (dj::IsInsideTet(v[0], v[1], v[2], v[3], pos, 0, &barycentric_of_joint_[i * 4])) {
      if (dj::IsInsideTet<double>(double_v[0], double_v[1], double_v[2], double_v[3],
                                  double_pos, 0, barycentric)) {// &)) {
        barycentric_of_joint_[i * 4 + 0] = barycentric[0];
        barycentric_of_joint_[i * 4 + 1] = barycentric[1];
        barycentric_of_joint_[i * 4 + 2] = barycentric[2];
        barycentric_of_joint_[i * 4 + 3] = barycentric[3];
        attached_tet_of_joints_[i] = t;
        is_constrainted_[tet_[t * 4 + 0]] = 1;
        is_constrainted_[tet_[t * 4 + 1]] = 1;
        is_constrainted_[tet_[t * 4 + 2]] = 1;
        is_constrainted_[tet_[t * 4 + 3]] = 1;
        constrained_vert_set_.insert(tet_[t * 4 + 0]);
        constrained_vert_set_.insert(tet_[t * 4 + 1]);
        constrained_vert_set_.insert(tet_[t * 4 + 2]);
        constrained_vert_set_.insert(tet_[t * 4 + 3]);
        break;
      }
    }
    ASSERT(attached_tet_of_joints_[i] != -1);
#endif
  }
  constrainted_vertex_.clear();
  constrainted_vertex_.insert(constrainted_vertex_.end(), constrained_vert_set_.begin(), constrained_vert_set_.end());
  constrainted_vert_vel_.resize(constrainted_vertex_.size() * 3);
  std::fill(constrainted_vert_vel_.begin(), constrainted_vert_vel_.end(), 0);
  P(constrainted_vertex_.size());
}



void Tet::AttachSkeleton(Skeleton * skeleton, bool find_attached_tet) {
  skeleton_ = skeleton;
  // Compute all the joints
  if (find_attached_tet) {
    FindBodyAttachedTetOfJoints(skeleton);
  }

  skeleton_->UpdatePosition();
  skeleton_->InitializeBoneRotation();
  for (int v = 0; v < vertex_num_; ++v) {
    closest_bone_[v] = skeleton_->AttachedToSkeleton(X + v * 3, initial_offset_ + v * 3);
  }
}

void Tet::AttachSkeleton(Skeleton * skeleton, Skeleton * rest_skeleton, const char *file_name) {
  BinaryFileReader in(file_name);
  int tmp_vertex_num;
  in.Read(&tmp_vertex_num, 1);
  //  P(tmp_vertex_num, vertex_num_);
  ASSERT(tmp_vertex_num == vertex_num_);
  in.Read(tmp_X, vertex_num_ * 3);
  dj::Swap(X, tmp_X);
  AttachSkeleton(rest_skeleton, false);
  dj::Swap(X, tmp_X);
  AttachSkeleton(skeleton, true);
}

void Tet::ReinitializeOffsetToSkeleton() {
  OMP_FOR
  for (int v = 0; v < vertex_num_; ++v) {
    double* pos = X + v * 3;
    double* offset = initial_offset_ + v * 3;
    skeleton_->AttachedToBone(pos, closest_bone_[v], offset);
  }
}

void Tet::Build_Edge_From_Tets() {
  std::vector<int> RE(tet_number * 18);
  for (int t = 0; t < tet_number; t++) {
    int v0 = tet_[t * 4 + 0];
    int v1 = tet_[t * 4 + 1];
    int v2 = tet_[t * 4 + 2];
    int v3 = tet_[t * 4 + 3];

    RE[t * 18 + 0] = v0;
    RE[t * 18 + 1] = v1;
    RE[t * 18 + 2] = t;
    RE[t * 18 + 3] = v1;
    RE[t * 18 + 4] = v2;
    RE[t * 18 + 5] = t;
    RE[t * 18 + 6] = v2;
    RE[t * 18 + 7] = v0;
    RE[t * 18 + 8] = t;

    RE[t * 18 + 0 + 9] = v0;
    RE[t * 18 + 1 + 9] = v3;
    RE[t * 18 + 2 + 9] = t;
    RE[t * 18 + 3 + 9] = v1;
    RE[t * 18 + 4 + 9] = v3;
    RE[t * 18 + 5 + 9] = t;
    RE[t * 18 + 6 + 9] = v2;
    RE[t * 18 + 7 + 9] = v3;
    RE[t * 18 + 8 + 9] = t;
  }

  for (int i = 0; i < tet_number * 6; i++)
    if (RE[i * 3 + 0] > RE[i * 3 + 1])	Swap(RE[i * 3 + 0], RE[i * 3 + 1]);

  //Quick_Sort
  Quick_Sort_RE(&RE[0], 0, tet_number * 6 - 1);

  std::vector<int> tmp_edge;
  for (int i = 0; i < tet_number * 6; i++) {
    if (i != 0 && RE[i * 3] == RE[(i - 1) * 3] && RE[i * 3 + 1] == RE[(i - 1) * 3 + 1]) {
    } else {
      //Add the edge to E
      tmp_edge.push_back(RE[i * 3 + 0]);
      tmp_edge.push_back(RE[i * 3 + 1]);
      if (tmp_edge[edge_num_ * 2 + 0] > tmp_edge[edge_num_ * 2 + 1]) {
        dj::Swap(tmp_edge[edge_num_ * 2 + 0], tmp_edge[edge_num_ * 2 + 1]);
      }
      edge_num_++;
    }
  }
  AllocateEdgeData(edge_num_);
  memcpy(edges_, &tmp_edge[0], sizeof(int) * edge_num_ * 2);

  for (int e = 0; e < edge_num_; e++) {
    rest_length_[e] = dj::Distance3(&X[edges_[e * 2 + 0] * 3], &X[edges_[e * 2 + 1] * 3]);
  }
  adjacent_vertex_.resize(vertex_num_);
  incident_edge_.resize(vertex_num_);
  for (int e = 0; e < edge_num_; e++) {
    int& v1 = edges_[e * 2 + 0];
    int& v2 = edges_[e * 2 + 1];
    adjacent_vertex_[v1].push_back(v2);
    incident_edge_[v1].push_back(e);
    adjacent_vertex_[v2].push_back(v1);
    incident_edge_[v2].push_back(e);
  }
  //    for (int i = 0; i < vertex_num_; ++i) {
  //      std::sort(adjacent_vertex_[i].begin(), adjacent_vertex_[i].end());
  //    }
}


void Tet::Quick_Sort_RE(int a[], int l, int r) {
  int j;
  if (l < r) {
    j = Quick_Sort_Partition_RE(a, l, r);
    Quick_Sort_RE(a, l, j - 1);
    Quick_Sort_RE(a, j + 1, r);
  }
}

int Tet::Quick_Sort_Partition_RE(int a[], int l, int r) {
  int pivot[3], i, j, t[3];
  pivot[0] = a[l * 3 + 0];
  pivot[1] = a[l * 3 + 1];
  pivot[2] = a[l * 3 + 2];
  i = l;
  j = r + 1;
  while (1) {
    do ++i;
    while ((a[i * 3] < pivot[0] || (a[i * 3] == pivot[0] && a[i * 3 + 1] <= pivot[1])) && i <= r);
    do --j;
    while (a[j * 3] > pivot[0] || (a[j * 3] == pivot[0] && a[j * 3 + 1] > pivot[1]));
    if (i >= j) break;
    //Swap i and j
    t[0] = a[i * 3 + 0];
    t[1] = a[i * 3 + 1];
    t[2] = a[i * 3 + 2];
    a[i * 3 + 0] = a[j * 3 + 0];
    a[i * 3 + 1] = a[j * 3 + 1];
    a[i * 3 + 2] = a[j * 3 + 2];
    a[j * 3 + 0] = t[0];
    a[j * 3 + 1] = t[1];
    a[j * 3 + 2] = t[2];
  }
  //Swap l and j
  t[0] = a[l * 3 + 0];
  t[1] = a[l * 3 + 1];
  t[2] = a[l * 3 + 2];
  a[l * 3 + 0] = a[j * 3 + 0];
  a[l * 3 + 1] = a[j * 3 + 1];
  a[l * 3 + 2] = a[j * 3 + 2];
  a[j * 3 + 0] = t[0];
  a[j * 3 + 1] = t[1];
  a[j * 3 + 2] = t[2];
  return j;
}

void Tet::Write_File(const char * file_name) {
  std::fstream output;
  output.open(file_name, std::ios::out | std::ios::binary);
  if (!output.is_open())	printf("Error, file not open.\n");
  Write_Binary(output, vertex_num_);
  Write_Binaries(output, X, vertex_num_ * 3);
  output.close();
  printf("Write file %s successfully.\n", file_name);
}

void Tet::Read_File(const char * file_name) {
  std::fstream input;
  input.open(file_name, std::ios::in | std::ios::binary);
  if (!input.is_open())	printf("Error, file not open.\n");
  Read_Binary(input, vertex_num_);
  Read_Binaries(input, X, vertex_num_ * 3);
  input.close();
  printf("Write file %s successfully.\n", file_name);
}

void Tet::ComputeRestDeformationGradientAndVolume() {
  memcpy(old_X, X, sizeof(double)*vertex_num_ * 3);
  for (int t = 0; t < tet_number; t++) {
    int p0 = tet_[t * 4 + 0] * 3;
    int p1 = tet_[t * 4 + 1] * 3;
    int p2 = tet_[t * 4 + 2] * 3;
    int p3 = tet_[t * 4 + 3] * 3;

    Dm[t * 9 + 0] = X[p1 + 0] - X[p0 + 0];
    Dm[t * 9 + 3] = X[p1 + 1] - X[p0 + 1];
    Dm[t * 9 + 6] = X[p1 + 2] - X[p0 + 2];
    Dm[t * 9 + 1] = X[p2 + 0] - X[p0 + 0];
    Dm[t * 9 + 4] = X[p2 + 1] - X[p0 + 1];
    Dm[t * 9 + 7] = X[p2 + 2] - X[p0 + 2];
    Dm[t * 9 + 2] = X[p3 + 0] - X[p0 + 0];
    Dm[t * 9 + 5] = X[p3 + 1] - X[p0 + 1];
    Dm[t * 9 + 8] = X[p3 + 2] - X[p0 + 2];

    volume_[t] = fabsf(Matrix_Inverse(&Dm[t * 9], &inv_Dm[t * 9])) / 6;

    //      if (Vol[t] < limit_threshold)	printf("too small %d: %ef\n", t, Vol[t]);
    //if(t==5120)	printf("TTT: %.16f\n", ttt);

    //printf("V: %f\n", Vol[t]);

    double e01[3], e02[3], e03[3], e12[3], e13[3];
    for (int i = 0; i < 3; i++) {
      e01[i] = X[p1 + i] - X[p0 + i];
      e02[i] = X[p2 + i] - X[p0 + i];
      e03[i] = X[p3 + i] - X[p0 + i];
      e12[i] = X[p2 + i] - X[p1 + i];
      e13[i] = X[p3 + i] - X[p1 + i];
    }
    double AN012[3], AN023[3], AN031[3], AN132[3];
    Cross(e01, e02, AN012);
    Cross(e02, e03, AN023);
    Cross(e03, e01, AN031);
    Cross(e13, e12, AN132);

    Area_Coeffs[t * 9 + 0] = (AN012[0] + AN132[0] + AN031[0]) / 6;
    Area_Coeffs[t * 9 + 3] = (AN012[1] + AN132[1] + AN031[1]) / 6;
    Area_Coeffs[t * 9 + 6] = (AN012[2] + AN132[2] + AN031[2]) / 6;
    Area_Coeffs[t * 9 + 1] = (AN012[0] + AN132[0] + AN023[0]) / 6;
    Area_Coeffs[t * 9 + 4] = (AN012[1] + AN132[1] + AN023[1]) / 6;
    Area_Coeffs[t * 9 + 7] = (AN012[2] + AN132[2] + AN023[2]) / 6;
    Area_Coeffs[t * 9 + 2] = (AN023[0] + AN132[0] + AN031[0]) / 6;
    Area_Coeffs[t * 9 + 5] = (AN023[1] + AN132[1] + AN031[1]) / 6;
    Area_Coeffs[t * 9 + 8] = (AN023[2] + AN132[2] + AN031[2]) / 6;
  }
}


void Tet::Read_Original_File(const char * prefix_name) {
  ASSERT(false, L("can't be used"));
  std::vector<double> verts;
  std::vector<int> tets;
  TetGenMeshIO::Instance()->Read(prefix_name, verts, tets);
  vertex_num_ = int(verts.size()) / 3;
  X = buffer_.Malloc<double>(vertex_num_);
  memcpy(X, &verts[0], sizeof(double) * verts.size());
  tet_number = int(tets.size()) / 4;
  tet_ = buffer_.Malloc<int>(tet_number);
  memcpy(tet_, &tets[0], sizeof(int) * tets.size());
#if 0
  tet_number = 2;

  tet_[0] = 0;
  tet_[1] = 1;
  tet_[2] = 2;
  tet_[3] = 3;

  tet_[4] = 0;
  tet_[5] = 2;
  tet_[6] = 1;
  tet_[7] = 4;

  vertex_num_ = 5;
  int v;
  double l = 0.3f;
  v = 0;
  X[v * 3 + 0] = l;
  X[v * 3 + 1] = 0;
  X[v * 3 + 2] = 0;

  v = 1;
  X[v * 3 + 0] = -l;
  X[v * 3 + 1] = 0;
  X[v * 3 + 2] = 0;

  v = 2;
  X[v * 3 + 0] = 0;
  X[v * 3 + 1] = 0;
  X[v * 3 + 2] = l;

  v = 3;
  X[v * 3 + 0] = 0;
  X[v * 3 + 1] = l;
  X[v * 3 + 2] = 0;

  v = 4;
  X[v * 3 + 0] = 0;
  X[v * 3 + 1] = -l;
  X[v * 3 + 2] = 0;
#endif
}

void Tet::Write_Original_File(const char *name) {
  char filename[1024];

  sprintf(filename, "%s.node", name);
  FILE *fp = fopen(filename, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }
  fprintf(fp, "%d %d %d %d\n", vertex_num_, 3, 0, 0);
  for (int i = 0; i < vertex_num_; i++)
    fprintf(fp, "%d %f %f %f\n", i + 1, X[i * 3 + 0], X[i * 3 + 1], X[i * 3 + 2]);
  fclose(fp);

  sprintf(filename, "%s.ele", name);
  fp = fopen(filename, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }
  fprintf(fp, "%d %d %d\n", tet_number, 4, 0);

  for (int i = 0; i < tet_number; i++)
    fprintf(fp, "%d %d %d %d %d\n", i + 1, tet_[i * 4 + 0] + 1, tet_[i * 4 + 1] + 1, tet_[i * 4 + 2] + 1, tet_[i * 4 + 3] + 1);
  fclose(fp);
}

void Tet::Load_SMF_File(char * name) {
  (void) name;
  //  FILE *fp = fopen(name, "r+");
  //  if (fp == NULL)	{
  //    printf("ERROR: file %s not open.\n", name);
  //    return;
  //  }

  //  vertex_num_ = 0;
  //  triangle_num_ = 0;

  //  fscanf(fp, "begin\n");
  //  while (1) {
  //    if (feof(fp))	break;
  //    char tag[1024];
  //    fscanf(fp, "%s", tag);

  //    if (tag[0] == 'v') {
  //      fscanf(fp, "%f %f %f\n", &X[vertex_num_ * 3 + 0], &X[vertex_num_ * 3 + 1], &X[vertex_num_ * 3 + 2]);
  //      vertex_num_++;
  //    }
  //    if (tag[0] == 'f') {
  //      fscanf(fp, "%d %d %d\n", &T[triangle_num_ * 3 + 0], &T[triangle_num_ * 3 + 1], &T[triangle_num_ * 3 + 2]);
  //      T[triangle_num_ * 3 + 0]--;
  //      T[triangle_num_ * 3 + 1]--;
  //      T[triangle_num_ * 3 + 2]--;
  //      triangle_num_++;
  //    }
  //  }
  //  fclose(fp);
}

void Tet::Output_Mesh(char * file_name) {
  std::fstream output;
  output.open(file_name, std::ios::out | std::ios::binary);
  if (!output.is_open())	printf("Error, file not open.\n");
  Write_Binary(output, vertex_num_);
  for (int i = 0; i < vertex_num_ * 3; i++)
    Write_Binary<double>(output, (double)(X[i]));
  Write_Binary(output, triangle_num_);
  Write_Binaries(output, T, triangle_num_ * 3);
  output.close();
  printf("Write file %s successfully.\n", file_name);
}

void Tet::Input_Mesh(char * file_name) {
  std::fstream input;
  input.open(file_name, std::ios::in | std::ios::binary);
  if (!input.is_open()) printf("Error, Triangle Mesh file %s not open.\n", file_name);
  Read_Binary(input, vertex_num_);

  double *temp_X = new double[vertex_num_ * 3];
  Read_Binaries(input, temp_X, vertex_num_ * 3);
  for (int i = 0; i < vertex_num_; i++) {
    X[i * 3 + 0] = temp_X[i * 3 + 0];
    X[i * 3 + 1] = temp_X[i * 3 + 1];
    X[i * 3 + 2] = temp_X[i * 3 + 2];
  }
  delete[]temp_X;
  Read_Binary(input, triangle_num_);
  Read_Binaries(input, T, triangle_num_ * 3);
  input.close();
  printf("here??? %d, %d\n", vertex_num_, triangle_num_);

  double min_area = 99999;
  int min_area_id = 0;

  for (int i = 0; i < triangle_num_; i++) {

    double normal[3];
    double area = Normal(&X[T[i * 3 + 0] * 3], &X[T[i * 3 + 1] * 3], &X[T[i * 3 + 2] * 3], normal);

    //printf("area %d: %f\n", i, area);
    if (area < min_area) {
      min_area = area;
      min_area_id = i;
    }
  }

  printf("min_area %d: %.14f\n", min_area_id, min_area);

}

void Tet::Load_PLY_File(char * name) {
  FILE *fp = fopen(name, "r+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", name);
    return;
  }

  vertex_num_ = 437645;
  triangle_num_ = 871414;


  for (int i = 0; i < vertex_num_; i++) {
    float x, y, z;
    fscanf(fp, "%f %f %f\n", &x, &y, &z);
    X[i * 3 + 0] = x;
    X[i * 3 + 1] = y;
    X[i * 3 + 2] = z;
  }

  int temp;

  double min_area = 99999;
  int min_area_id = 0;

  for (int i = 0; i < triangle_num_; i++) {
    //printf("start\n");
    fscanf(fp, "%d %d %d %d\n", &temp, &T[i * 3 + 0], &T[i * 3 + 1], &T[i * 3 + 2]);
    //printf("%d: %d, %d, %d\n", i, T[i*3], T[i*3+1], T[i*3+2]);

    double normal[3];
    double area = Normal(&X[T[i * 3 + 0] * 3], &X[T[i * 3 + 1] * 3], &X[T[i * 3 + 2] * 3], normal);

    //printf("area %d: %f\n", i, area);
    if (area < min_area) {
      min_area = area;
      min_area_id = i;
    }
  }
  printf("min_area %d: %f (%d, %d)\n", min_area_id, min_area, vertex_num_, triangle_num_);

  fclose(fp);
}

void Tet::Output_POLY_File(const char * name) {
#if 0
  TextFileWriter out(name);
  out << vertex_num_ << " 3 0 0" << std::endl;
  for (int i = 0; i < vertex_num_; i++) {
    out << i + 1 << " " << X[i * 3 + 0] << " " << X[i * 3 + 1] << " " << X[i * 3 + 2] << std::endl;
  }
  out << std::endl;

  out << triangle_num_;
  for (int i = 0; i < triangle_num_; i++) {
    out << "1 0 1" << std::endl;
    out << "3 " << T[i * 3 + 0] << " " << T[i * 3 + 1] << " " << T[i * 3 + 2] << std::endl;
  }
  out << std::endl;
  out << "0" << std::endl;
  out << std::endl;
  out << "0" << std::endl;
#else
  FILE *fp = fopen(name, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", name);
    return;
  }

  fprintf(fp, "%d %d %d %d\n", vertex_num_, 3, 0, 0);
  for (int i = 0; i < vertex_num_; i++)
    fprintf(fp, "%d %f %f %f\n", i + 1, X[i * 3 + 0], X[i * 3 + 1], X[i * 3 + 2]);
  fprintf(fp, "\n");

  fprintf(fp, "%d %d\n", triangle_num_, 1);
  for (int i = 0; i < triangle_num_; i++) {
    fprintf(fp, "%d %d %d\n", 1, 0, 1);

    fprintf(fp, "%d %d %d %d\n", 3, T[i * 3 + 0] + 1, T[i * 3 + 1] + 1, T[i * 3 + 2] + 1);
    //printf("%d: %d, %d, %d\n", i, T[i*3], T[i*3+1], T[i*3+2]);

  }

  fprintf(fp, "\n");
  fprintf(fp, "0\n");
  fprintf(fp, "\n");
  fprintf(fp, "0\n");
  fclose(fp);
#endif
}

void Tet::Save_OBJ_File(const char * name) {
#if 1
  FILE *fp = fopen(name, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", name);
    return;
  }

  for (int i = 0; i < vertex_num_; i++)
    fprintf(fp, "v %f %f %f\n", X[i * 3 + 0], X[i * 3 + 1], X[i * 3 + 2]);
  for (int i = 0; i < triangle_num_; i++)
    fprintf(fp, "f %d %d %d\n", T[i * 3 + 0] + 1, T[i * 3 + 1] + 1, T[i * 3 + 2] + 1);
  fclose(fp);
#else
  std::ofstream out(name);
  ASSERT(out.is_open(), P(name));
  std::vector<int> global_vert2surface_vert(vertex_num_, -1);
  int i = 0;
  for (int v : surface_vertices_) {
    global_vert2surface_vert[v] = i;
    out << "v " << X[v * 3 + 0] << " " << X[v * 3 + 1] << " " << X[v * 3 + 2] << "\n";
    ++i;
  }
  for (int i = 0; i < triangle_num_; i++) {
    int verts[3] = {
      global_vert2surface_vert[T[i * 3 + 0]] + 1,
      global_vert2surface_vert[T[i * 3 + 1]] + 1,
      global_vert2surface_vert[T[i * 3 + 2]] + 1,
    };
    out << "f " << verts[0] << " " << verts[1] << " " << verts[2] << "\n";
  }
  out.close();
#endif
}

void Tet::Save_PLY_File(char * name) {
  FILE *fp = fopen(name, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", name);
    return;
  }

  fprintf(fp, "ply\n");
  fprintf(fp, "format ascii 1.0\n");
  fprintf(fp, "comment generated by www\n");
  fprintf(fp, "element vertex %d\n", vertex_num_);
  fprintf(fp, "property double x\n");
  fprintf(fp, "property double y\n");
  fprintf(fp, "property double z\n");
  fprintf(fp, "element face %d\n", triangle_num_);
  fprintf(fp, "property list uchar int vertex_indices\n");
  fprintf(fp, "end_header\n");

  for (int i = 0; i < vertex_num_; i++)
    fprintf(fp, "%f %f %f\n", X[i * 3 + 0], X[i * 3 + 1], X[i * 3 + 2]);
  for (int i = 0; i < triangle_num_; i++) {
    fprintf(fp, "%d %d %d %d\n", 3, T[i * 3 + 0], T[i * 3 + 1], T[i * 3 + 2]);
  }
  fclose(fp);
}

void Tet::Output_STL_File(char * name) {
  char filename[1024];
  FILE *fp = fopen(name, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }

  fprintf(fp, "solid\n");
  for (int t = 0; t < triangle_num_; t++) {
    int v0 = T[t * 3 + 0];
    int v1 = T[t * 3 + 1];
    int v2 = T[t * 3 + 2];
    fprintf(fp, "facet normal %f %f %f\n", TN[t * 3 + 0], TN[t * 3 + 1], TN[t * 3 + 2]);
    fprintf(fp, "outer loop\n");
    fprintf(fp, "vertex %f %f %f\n", X[v0 * 3 + 0], X[v0 * 3 + 1], X[v0 * 3 + 2]);
    fprintf(fp, "vertex %f %f %f\n", X[v1 * 3 + 0], X[v1 * 3 + 1], X[v1 * 3 + 2]);
    fprintf(fp, "vertex %f %f %f\n", X[v2 * 3 + 0], X[v2 * 3 + 1], X[v2 * 3 + 2]);
    fprintf(fp, "endloop\n");
  }
  fprintf(fp, "endsolid\n");
  fclose(fp);
}

void Tet::SavePosition(const char * file_name) {
  std::ofstream out(file_name, std::ios::out | std::ios::binary);
  if (!out.is_open()) {
    std::cerr << "failed to open file " << file_name << std::endl;
    exit(0);
  }
  out.write((char*)&vertex_num_, sizeof(int));
  out.write((char*)X, sizeof(double) * 3 * vertex_num_);
  out.close();
}

void Tet::LoadPosition(const char * file_name) {
  std::ifstream in(file_name, std::ios::in | std::ios::binary);
  if (!in.is_open()) {
    std::cerr << "failed to open file " << file_name << std::endl;
    exit(0);
  }
  int new_vertex_num;
  in.read((char*)&new_vertex_num, sizeof(int));
  ASSERT(new_vertex_num == vertex_num_);
  in.read((char*)X, sizeof(double) * 3 * vertex_num_);
  in.close();
}

void Tet::LoadPosition(std::vector<double> &pos) {
  ASSERT(int(pos.size()) == vertex_num_ * 3);
  memcpy(X, &pos[0], sizeof(double) * 3 * vertex_num_);
}




void Tet::Render(int render_mode, double * pos) {
  if (render_mode == kDefaultRendering) render_mode = kRenderMode[render_mode_];
  if (pos == NULL) pos = X;
  if (0) {
    glDisable(GL_LIGHTING);
    if (selected_vertex_ >= 0 && selected_vertex_ < vertex_num_) {
      glPointSize(8);
      glColor3fv(kRed());
      glBegin(GL_POINTS);
      Vertex3v(X + selected_vertex_ * 3);
      glEnd();
    }
    glPointSize(6);
    glBegin(GL_POINTS);
    glColor3fv(kGreen());
    for (int v : constrainted_vertex_) {
      Vertex3v(X + v * 3);
    }
    glEnd();
  }

  // render bone point
  if (skeleton_) {
    glColor3fv(kBlue());
    glBegin(GL_POINTS);
    for (int i = 0; i < int(joints_.size()) / 2; ++i) {
      Vec3d end_points[2] = {
        Vec3d(skeleton_->GetNode(joints_[i * 2 + 0].first)->world_pos()),
        Vec3d(skeleton_->GetNode(joints_[i * 2 + 1].first)->world_pos())
      };
      Vec3d pos = end_points[0] * joints_[i * 2].second + end_points[1] * joints_[i * 2 + 1].second;
      Vertex3v(&pos[0]);
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
  if (render_mode & kNoRendering) return;
  if ((render_mode & kSmoothLighting) || (render_mode & kFlatLighting)) {
    if (textured_surface_renderer_.size() != 0) {
      textured_surface_renderer_[0].Render();
    } else if (surface_renderer_ != nullptr) {
      surface_renderer_->Render();
    } else {
      //    Build_VN();
      glPushAttrib(GL_ENABLE_BIT);
      glPushAttrib(GL_LIGHTING_BIT);
      glColor3fv(kBlue());
      float diffuse_color[4] = { 1, 1, 1, 0.6f};
      glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
      float dark[4] = {0, 0, 0, 1};
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, dark);
      if (render_mode & kFlatLighting) {
        glShadeModel(GL_FLAT);
      } else {
        glShadeModel(GL_SMOOTH);
      }

      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.0, 1.0);
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < triangle_num_; i++) {
        double *p0 = &pos[T[i * 3 + 0] * 3];
        double *p1 = &pos[T[i * 3 + 1] * 3];
        double *p2 = &pos[T[i * 3 + 2] * 3];
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
        glNormal3f(VN[T[i * 3 + 0] * 3 + 0], VN[T[i * 3 + 0] * 3 + 1], VN[T[i * 3 + 0] * 3 + 2]);
        Vertex3v(p0);
        glNormal3f(VN[T[i * 3 + 1] * 3 + 0], VN[T[i * 3 + 1] * 3 + 1], VN[T[i * 3 + 1] * 3 + 2]);
        Vertex3v(p1);
        glNormal3f(VN[T[i * 3 + 2] * 3 + 0], VN[T[i * 3 + 2] * 3 + 1], VN[T[i * 3 + 2] * 3 + 2]);
        Vertex3v(p2);
      }
      glEnd();

      glPopAttrib();
      glPopAttrib();
      glDisable(GL_CULL_FACE);
    }
  }

  //  std::set<int> tt({131038,131045,131046,131055,131056,131061});
  if (render_mode & kTets) {
    //    glBegin(GL_TRIANGLES);
    glBegin(GL_LINES);
    for (int t = 0; t < tet_number; ++t) {
      //    for (int t : tt) {
      int* verts = tet_ + t * 4;
      int b = closest_bone_[verts[0]];
      //      if (tt.count(t) == 0) {
      if (closest_bone_[verts[1]] != b || closest_bone_[verts[2]] != b || closest_bone_[verts[3]] != b) {
        continue;
      }
      //      } else {

      //      }
      //      glColor3fv(kRed());
      glColor3fv(color_map[b * 71 % color_num]);
      Vertex3v(pos + verts[0] * 3);
      Vertex3v(pos + verts[1] * 3);
      Vertex3v(pos + verts[2] * 3);

      Vertex3v(pos + verts[0] * 3);
      Vertex3v(pos + verts[1] * 3);
      Vertex3v(pos + verts[3] * 3);

      Vertex3v(pos + verts[1] * 3);
      Vertex3v(pos + verts[2] * 3);
      Vertex3v(pos + verts[3] * 3);

      Vertex3v(pos + verts[0] * 3);
      Vertex3v(pos + verts[2] * 3);
      Vertex3v(pos + verts[3] * 3);
    }
    glEnd();
  }

  if (render_mode & kWireFrame) {
    glColor3fv(kGreen());
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LINE_WIDTH);
    glDisable(GL_LIGHTING);
    glLineWidth(0.2f);
    //    glLineStipple(3, 0xAAAA);
    //    glEnable(GL_LINE_STIPPLE);
    //    for (int i = 0; i < triangle_num_; i++) {
    //      double *p0 = &X[T[i * 3 + 0] * 3];
    //      double *p1 = &X[T[i * 3 + 1] * 3];
    //      double *p2 = &X[T[i * 3 + 2] * 3];
    //      glBegin(GL_LINES);
    //      glVertex3fv(p0);
    //      glVertex3fv(p1);
    //      glVertex3fv(p2);
    //      glEnd();
    //    }
    glBegin(GL_LINES);
    for (int i = 0; i < edge_num_; ++i) {
      double* v0 = pos + 3 * edges_[i * 2 + 0];
      double* v1 = pos + 3 * edges_[i * 2 + 1];
      Vertex3v(v0);
      Vertex3v(v1);
    }
    glEnd();

    //    glColor3fv(kRed());
    //    glPointSize(5.0);
    //    glBegin(GL_POINTS);
    //    for (int i = 0; i < vertex_num_; ++i) {
    //      glColor3fv(kBlue());
    //      for (int j = 0; j < incident_tet_[i].size(); ++j) {
    //        if (incident_tet_[i][j] == 4537) {
    //          glColor3fv(kRed());
    //        }
    //      }
    ////      if (skeleton_tag_[i] >= 0)
    //        glVertex3fv(X + 3 * i);
    //    }
    //    glEnd();
    glPopAttrib();
    glPopAttrib();
  }
  if (render_mode & kSurfaceEdge) {
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LINE_WIDTH);
    glLineStipple(3, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);
    glDisable(GL_LIGHTING);
    glColor3fv(kBlack());
    glBegin(GL_LINES);
    for (int i = 0; i < triangle_num_; ++i) {
      int* v = T + i * 3;
      double* v_pos[3] = {
        v[0] * 3 + pos,
        v[1] * 3 + pos,
        v[2] * 3 + pos
      };
      Vertex3v(v_pos[0]);
      Vertex3v(v_pos[1]);
      Vertex3v(v_pos[1]);
      Vertex3v(v_pos[2]);
      Vertex3v(v_pos[2]);
      Vertex3v(v_pos[0]);
    }
    glEnd();
    glPopAttrib();
    glPopAttrib();
  }

  if (render_mode & kTransparent) {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    glDrawBuffer(GL_NONE);
    Render(kSmoothLighting, pos);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDrawBuffer(GL_BACK);
    glEnable(GL_LIGHTING);
    glColor3f(0.0, 0.0, 0.0);
    Render(kSmoothLighting, pos);
    glDisable(GL_BLEND);
  }

}


void Tet::Build_Boundary_Triangles() {
  std::vector<int> dummy(tet_number * 4 * 4);
  int* temp_T = &dummy[0];
  for (int i = 0; i < tet_number; i++) {
    temp_T[i * 16 + 0] = tet_[i * 4 + 0];
    temp_T[i * 16 + 1] = tet_[i * 4 + 1];
    temp_T[i * 16 + 2] = tet_[i * 4 + 2];
    temp_T[i * 16 + 3] = i;

    temp_T[i * 16 + 4] = tet_[i * 4 + 0];
    temp_T[i * 16 + 5] = tet_[i * 4 + 2];
    temp_T[i * 16 + 6] = tet_[i * 4 + 3];
    temp_T[i * 16 + 7] = i;

    temp_T[i * 16 + 8] = tet_[i * 4 + 0];
    temp_T[i * 16 + 9] = tet_[i * 4 + 3];
    temp_T[i * 16 + 10] = tet_[i * 4 + 1];
    temp_T[i * 16 + 11] = i;

    temp_T[i * 16 + 12] = tet_[i * 4 + 1];
    temp_T[i * 16 + 13] = tet_[i * 4 + 3];
    temp_T[i * 16 + 14] = tet_[i * 4 + 2];
    temp_T[i * 16 + 15] = i;
  }

  for (int i = 0; i < tet_number * 4; i++) {
    std::sort(temp_T + i * 4, temp_T + i * 4 + 3);
    //    if (temp_T[i * 4 + 1] < temp_T[i * 4 + 0]) {
    //      Swap(temp_T[i * 4 + 0], temp_T[i * 4 + 1]);
    //      //      temp_T[i * 4 + 3] *= -1;
    //    }
    //    if (temp_T[i * 4 + 2] < temp_T[i * 4 + 0]) {
    //      Swap(temp_T[i * 4 + 0], temp_T[i * 4 + 2]);
    //      //      temp_T[i * 4 + 3] = (temp_T[i * 4 + 3] + 1) % 2;
    //      //      temp_T[i * 4 + 3] *= -1;
    //    }
    //    if (temp_T[i * 4 + 2] < temp_T[i * 4 + 1]) {
    //      Swap(temp_T[i * 4 + 1], temp_T[i * 4 + 2]);
    //      //      temp_T[i * 4 + 3] = (temp_T[i * 4 + 3] + 1) % 2;
    //      //      temp_T[i * 4 + 3] *= -1;
    //    }
  }

  QuickSort(&temp_T[0], 0, tet_number * 4 - 1);

  incident_triangle_.resize(vertex_num_);
  triangle_num_ = 0;
  std::vector<int> tmp_triangle;
  std::unordered_set<int> surface_vertex_set;
  struct IntPairHash {
    std::size_t operator()(const std::pair<int, int>& p) const {
      static_assert(sizeof(int) * 2 == sizeof(long long), "expect long type to be twice the size of int.");
      long long val = p.first;
      val = (val << 32) | p.second;
      return std::hash<long long>()(val);
    }
  };
  // each edge (v0, v1) is mapped to its two incident triangles (assume manifold mesh)
  std::unordered_set<int> surface_tet_set;
  for (int i = 0; i < tet_number * 4; i++) {
    //    P(dj::Vec4i(temp_T + i * 4));
    if (i != tet_number * 4 - 1
        && temp_T[i * 4 + 0] == temp_T[i * 4 + 4] && temp_T[i * 4 + 1] == temp_T[i * 4 + 5]
        && temp_T[i * 4 + 2] == temp_T[i * 4 + 6]) {
      i++;
      continue;
    }

    int* verts = temp_T + i * 4;

    tmp_triangle.push_back(verts[0]);
    tmp_triangle.push_back(verts[1]);
    tmp_triangle.push_back(verts[2]);
    //    P(dj::Vec3i(T + triangle_num_ * 3))

    int forth_v = -1;
    int tet = temp_T[i * 4 + 3];
    for (int ii = 0; ii < 4; ++ii) {
      if (tet_[tet * 4 + ii] != temp_T[i * 4 + 0] &&
          tet_[tet * 4 + ii] != temp_T[i * 4 + 1] &&
          tet_[tet * 4 + ii] != temp_T[i * 4 + 2]) {
        forth_v = tet_[tet * 4 + ii];
        break;
      }
    }
    ASSERT(forth_v != -1);
    double vec[3][3];
    dj::SubVec3(X + forth_v * 3, X + verts[0] * 3, vec[2]);
    dj::SubVec3(X + verts[1] * 3, X + verts[0] * 3, vec[0]);
    dj::SubVec3(X + verts[2] * 3, X + verts[1] * 3, vec[1]);
    double cross[3];
    dj::Cross3(vec[0], vec[1], cross);
    if (dj::Dot3(cross, vec[2]) > 0) {
      dj::Swap(tmp_triangle[triangle_num_ * 3 + 2], tmp_triangle[triangle_num_ * 3 + 1]);
    }


    incident_triangle_[verts[0]].push_back(triangle_num_);
    incident_triangle_[verts[1]].push_back(triangle_num_);
    incident_triangle_[verts[2]].push_back(triangle_num_);
    surface_vertex_set.insert(tmp_triangle[triangle_num_ * 3 + 0]);
    surface_vertex_set.insert(tmp_triangle[triangle_num_ * 3 + 1]);
    surface_vertex_set.insert(tmp_triangle[triangle_num_ * 3 + 2]);

    int tet_idx = temp_T[i * 4 + 3];
    surface_tet_set.insert(tet_idx);
    triangle_num_++;
  }
  AllocateTriangleData(triangle_num_);
  memcpy(T, &tmp_triangle[0], sizeof(int) * triangle_num_ * 3);
  surface_vertices_.clear();
  surface_vertices_.insert(surface_vertices_.end(), surface_vertex_set.begin(), surface_vertex_set.end());
  surface_tets_.clear();
  surface_tets_.insert(surface_tets_.end(), surface_tet_set.begin(), surface_tet_set.end());

  is_boundary_tri_on_tet_.resize(surface_tets_.size() * 4);
  //  std::fill(is_boundary_tri_on_tet_.begin(), is_boundary_tri_on_tet_.end(), true);

  for (int i = 0; i < int(surface_tets_.size()); ++i) {
    int tet = surface_tets_[i];
    int is_boundary_vertex[4] = {
      surface_vertex_set.find(tet_[tet * 4 + 0]) != surface_vertex_set.end(),
      surface_vertex_set.find(tet_[tet * 4 + 1]) != surface_vertex_set.end(),
      surface_vertex_set.find(tet_[tet * 4 + 2]) != surface_vertex_set.end(),
      surface_vertex_set.find(tet_[tet * 4 + 3]) != surface_vertex_set.end()
    };
    is_boundary_tri_on_tet_[i * 4 + 0] = is_boundary_vertex[1] && is_boundary_vertex[2] && is_boundary_vertex[3];
    is_boundary_tri_on_tet_[i * 4 + 1] = is_boundary_vertex[2] && is_boundary_vertex[3] && is_boundary_vertex[0];
    is_boundary_tri_on_tet_[i * 4 + 2] = is_boundary_vertex[3] && is_boundary_vertex[0] && is_boundary_vertex[1];
    is_boundary_tri_on_tet_[i * 4 + 3] = is_boundary_vertex[0] && is_boundary_vertex[1] && is_boundary_vertex[2];
  }

  std::vector<int> global_vert2surface_vert(vertex_num_, -1);
  surface_vert_pos_.resize(surface_vertices_.size()  * 3);
  for (int i = 0; i < int(surface_vertices_.size()); ++i) {
    int v = surface_vertices_[i];
    global_vert2surface_vert[v] = i;
  }
  surface_triangle_indices_.resize(triangle_num_ * 3);
  surface_vert_normal_.resize(triangle_num_ * 3);
  for (int i = 0; i < triangle_num_; i++) {
    int verts[3] = {
      global_vert2surface_vert[T[i * 3 + 0]],
      global_vert2surface_vert[T[i * 3 + 1]],
      global_vert2surface_vert[T[i * 3 + 2]],
    };
    surface_triangle_indices_[i * 3 + 0] = verts[0];
    surface_triangle_indices_[i * 3 + 1] = verts[1];
    surface_triangle_indices_[i * 3 + 2] = verts[2];
  }
}

void Tet::QuickSort(int a[], int l, int r) {
  int j;
  if (l < r) {
    j = QuickSort_Partition(a, l, r);
    QuickSort(a, l, j - 1);
    QuickSort(a, j + 1, r);
  }
}

int Tet::QuickSort_Partition(int a[], int l, int r) {
  int pivot[4], i, j, t[4];
  pivot[0] = a[l * 4 + 0];
  pivot[1] = a[l * 4 + 1];
  pivot[2] = a[l * 4 + 2];
  pivot[3] = a[l * 4 + 3];
  i = l;
  j = r + 1;

  while (1) {
    do ++i;
    while (
      (a[i * 4 + 0] < pivot[0] || (a[i * 4 + 0] == pivot[0] && a[i * 4 + 1] < pivot[1]) || (a[i * 4 + 0] == pivot[0] && a[i * 4 + 1] == pivot[1] && a[i * 4 + 2] <= pivot[2]))
      && i <= r);
    do --j;
    while (a[j * 4 + 0] > pivot[0] || (a[j * 4 + 0] == pivot[0] && a[j * 4 + 1] > pivot[1]) || (a[j * 4 + 0] == pivot[0] && a[j * 4 + 1] == pivot[1] && a[j * 4 + 2] > pivot[2]));
    if (i >= j) break;
    //Swap i and j
    t[0] = a[i * 4 + 0];
    t[1] = a[i * 4 + 1];
    t[2] = a[i * 4 + 2];
    t[3] = a[i * 4 + 3];
    a[i * 4 + 0] = a[j * 4 + 0];
    a[i * 4 + 1] = a[j * 4 + 1];
    a[i * 4 + 2] = a[j * 4 + 2];
    a[i * 4 + 3] = a[j * 4 + 3];
    a[j * 4 + 0] = t[0];
    a[j * 4 + 1] = t[1];
    a[j * 4 + 2] = t[2];
    a[j * 4 + 3] = t[3];
  }
  //Swap l and j
  t[0] = a[l * 4 + 0];
  t[1] = a[l * 4 + 1];
  t[2] = a[l * 4 + 2];
  t[3] = a[l * 4 + 3];
  a[l * 4 + 0] = a[j * 4 + 0];
  a[l * 4 + 1] = a[j * 4 + 1];
  a[l * 4 + 2] = a[j * 4 + 2];
  a[l * 4 + 3] = a[j * 4 + 3];
  a[j * 4 + 0] = t[0];
  a[j * 4 + 1] = t[1];
  a[j * 4 + 2] = t[2];
  a[j * 4 + 3] = t[3];
  return j;
}

void Tet::Build_TN() {
  //  memset(TN, 0, sizeof(double)* triangle_num_ * 3);
  OMP_FOR
  for (int i = 0; i < triangle_num_; i++) {
    double*p0 = &X[T[i * 3 + 0] * 3];
    double*p1 = &X[T[i * 3 + 1] * 3];
    double*p2 = &X[T[i * 3 + 2] * 3];
    Normal(p0, p1, p2, &TN[i * 3]);
  }
}

void Tet::Build_VN() {
  //  return;
  Build_TN();
#if 1
  memset(VN, 0, sizeof(double)*vertex_num_ * 3);
  OMP_FOR
  for (int v = 0; v < vertex_num_; ++v) {
    double* normal = VN + v * 3;
    for (auto tri : incident_triangle_[v]) {
      normal[0] += TN[tri * 3 + 0];
      normal[1] += TN[tri * 3 + 1];
      normal[2] += TN[tri * 3 + 2];
    }
    //    dj::Normalize3(normal);
  }
  //  return;
#else
  for (int i = 0; i < triangle_num_; i++) {
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
#endif

  //  OMP_FOR
  //  for (int i = 0; i < vertex_num_; i++) {
  //    double length2, inv_length;
  //    length2 = VN[i * 3 + 0] * VN[i * 3 + 0] + VN[i * 3 + 1] * VN[i * 3 + 1] + VN[i * 3 + 2] * VN[i * 3 + 2];
  //    if (length2 < 1e-16f)	continue;
  //    inv_length = 1.0f / sqrtf(length2);

  //    VN[i * 3 + 0] *= inv_length;
  //    VN[i * 3 + 1] *= inv_length;
  //    VN[i * 3 + 2] *= inv_length;
  //  }
}

double Tet::ComputeAvgEdgeLength() {
  double avg = 0;
  for (int e = 0; e < edge_num_; ++e) {
    int* v = edges_ + e * 2;
    avg += dj::Distance3(X + v[0] * 3, X + v[1] * 3);
  }
  return avg / edge_num_;
}


void Tet::BuildIncidentTet() {
  incident_tet_.resize(vertex_num_);
  rank_in_tet_.resize(vertex_num_);
  for (int i = 0; i < tet_number; ++i) {
    int v0 = tet_[i * 4 + 0];
    int v1 = tet_[i * 4 + 1];
    int v2 = tet_[i * 4 + 2];
    int v3 = tet_[i * 4 + 3];
    rank_in_tet_[v0].push_back(0);
    rank_in_tet_[v1].push_back(1);
    rank_in_tet_[v2].push_back(2);
    rank_in_tet_[v3].push_back(3);
    incident_tet_[v0].push_back(i);
    incident_tet_[v1].push_back(i);
    incident_tet_[v2].push_back(i);
    incident_tet_[v3].push_back(i);
  }
  for (int i = 0; i < tet_number; ++i) {
    int* v = tet_ +  4 * i;
    int max_incident_tet = (int) incident_tet_[v[0]].size();
    max_incident_tet = dj::Max<int>(max_incident_tet, (int) incident_tet_[v[1]].size());
    max_incident_tet = dj::Max<int>(max_incident_tet, (int) incident_tet_[v[2]].size());
    max_incident_tet = dj::Max<int>(max_incident_tet, (int) incident_tet_[v[3]].size());
    tet_weights_[i] = 1.0f / max_incident_tet;
  }
}



void Tet::ApplySkeletonTransformationToVertex() {
  double* skeleton_transformation = &skeleton_->skeleton_transformation_[0];
  for (int v = 0; v < vertex_num_; ++v) {
    double* pos = X + v * 3;
    int attached_bone = closest_bone_[v];
    dj::MulMatrix3x3Vec<double>((double(*)[3]) (skeleton_transformation + attached_bone * 12), initial_offset_ + v * 3, pos);
    pos[0] += skeleton_transformation[attached_bone * 12 + 9];
    pos[1] += skeleton_transformation[attached_bone * 12 + 10];
    pos[2] += skeleton_transformation[attached_bone * 12 + 11];
  }
}

void Tet::GetConstrainVertexVelocity() {
  double* skeleton_transformation = &skeleton_->skeleton_transformation_[0];
  for (int i = 0; i < (int) constrainted_vertex_.size(); ++i) {
    int v = constrainted_vertex_[i];
    double pos[3] = {
      X[v * 3 + 0],
      X[v * 3 + 1],
      X[v * 3 + 2]
    };
    double old_pos[3] = {pos[0], pos[1], pos[2]};
    int attached_bone = closest_bone_[v];
    dj::MulMatrix3x3Vec<double>((double(*)[3]) (skeleton_transformation + attached_bone * 12), initial_offset_ + v * 3, pos);
    pos[0] += skeleton_transformation[attached_bone * 12 + 9];
    pos[1] += skeleton_transformation[attached_bone * 12 + 10];
    pos[2] += skeleton_transformation[attached_bone * 12 + 11];
    constrainted_vert_vel_[i * 3 + 0] = (pos[0] - old_pos[0]) / global::time_step;
    constrainted_vert_vel_[i * 3 + 1] = (pos[1] - old_pos[1]) / global::time_step;
    constrainted_vert_vel_[i * 3 + 2] = (pos[2] - old_pos[2]) / global::time_step;
  }
}



const int Tet::kRenderMode[] = {
  Tet::kSmoothLighting,
  Tet::kSurfaceEdge | Tet::kSmoothLighting,
  Tet::kParitionedMeshVertex,
  Tet::kParitionedMeshVertex | Tet::kSurfaceEdge | Tet::kSmoothLighting,
  //  Tet::kWireFrame | Tet::kSmoothLighting,
  Tet::kFlatLighting,
  Tet::kFlatLighting | Tet::kSurfaceEdge,
  Tet::kSurfaceEdge,
  Tet::kWireFrame,
  Tet::kTransparent,
  Tet::kTets,
  Tet::kTets | Tet::kWireFrame,
  Tet::kNoRendering,
};
const int Tet::kRenderModeNum = sizeof(Tet::kRenderMode) / sizeof(int);
