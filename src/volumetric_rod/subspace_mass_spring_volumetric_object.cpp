#include "subspace_mass_spring_volumetric_object.h"
#include <vector>
#include <set>
#include "tet_collider.h"
#include "vector_io.h"
#include "config_file.h"
#include "affine_transformer.h"
#include "random.h"
#include "conjugate_gradient_solver.h"
#include "opengl_helper.h"
#include "rainbow_color.h"
#include "text_file_io.h"
#include "print_macro.h"

int v = 3019;
std::set<int> force_vert({
  2999,
  //                             3399,
  //                           1799,
  //  3399,
  //  4599,
  //  1399,
  //  4599,
  //  3219,
  //  2959,
  //  2979,
  //  2999,
  //  3019,
  //  3039,
  //  3059,
  //  3079,
  //  3099,
});
//   13019

int ring = 1;
SubspaceMassSpringVolumetricObject::Vec3 force(0, -400000, 0);
//SubspaceMassSpringVolumetricObject::Vec3 force(-600000, 0, 0);

bool IsColumnMajor(SubspaceMassSpringVolumetricObject::MatrixX& mat)
{
  typedef SubspaceMassSpringVolumetricObject::MatrixX::Scalar Float;
  Float* data = mat.data();
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (data[j * mat.rows() + i] != mat(i, j)) {
        return false;
      }
    }
  }
  return true;
}

std::vector<std::vector<SubspaceMassSpringVolumetricObject::Vec3> > frames;


SubspaceMassSpringVolumetricObject::SubspaceMassSpringVolumetricObject(TetrahedralMeshIO *mesh_io,
                                                                       const char *file_name,
                                                                       AffineTransformer<Real>* transformer)
  : Super(mesh_io, file_name, transformer)
{
  Construct();
}

//SubspaceMassSpringVolumetricObject::SubspaceMassSpringVolumetricObject(TetrahedralMeshIO *mesh_io,
//                                                                       std::vector<const char *> file_names,
//                                                                       std::vector<AffineTransformer<Real> *> transformers)
//  : Super(mesh_io, file_names, transformers) {
//  Construct();
//}


SubspaceMassSpringVolumetricObject::SubspaceMassSpringVolumetricObject(TetrahedralMesh *mesh)
  : Super(mesh)
{
  Construct();
}

void SubspaceMassSpringVolumetricObject::LoadSamplePose(const char *file_name, std::vector<SubspaceMassSpringVolumetricObject::VectorX>& poses)
{
  TextFileReader in(file_name);
  int pos_num;
  in >> pos_num;
  int q_num;
  in >> q_num;
  ASSERT(q_num == global_subspace_coord_num_);
  poses.resize(pos_num);
  for (int i = 0; i < pos_num; ++i) {
    poses[i] = VectorX::Zero(global_subspace_coord_num_);
    for (int j = 0; j < global_subspace_coord_num_; ++j) {
      in >> poses[i][j];
    }
  }
}

void SubspaceMassSpringVolumetricObject::ViewSamlePose(const char *pos_file)
{
  static int i = -1;
  static std::vector<VectorX> poses;
  if (i == -1) {
    LoadSamplePose(pos_file, poses);
    i = 0;
  }
  //  P(i, poses.size());
  q_ = poses[i];
  ReconstructObjectFromSubspace();
  i = (i + 1) % poses.size();
}

void SubspaceMassSpringVolumetricObject::GenerateSamplePose(const char* pose_file, int pos_num)
{
  EXECUTE_TIMES(1);
  const Real kPi = 3.1415926;
  const Real kMaxForce = 700000;
  const Real kMinForce = 600000;
  std::set<int> applied_verts;
  std::map<int, Vec3> applied_force;
  Random rand;
  for (int i = 0; i < pos_num;) {
    int v = rand.GetRandom(kVertexNumPerSegment - 1, 0);
    if (vert_[v][1] < 1) continue;
    //    if (applied_verts.count(v) > 0) continue;
    Real theta = rand.GetRandom<Real>(0, 2 * kPi);
    Real phi = rand.GetRandom<Real>(-kPi / 2, kPi / 2);
    Vec3 direction(0, 0, 0);
    direction[0] = sin(theta) * cos(phi);
    direction[1] = sin(theta) * sin(phi);
    direction[2] = cos(theta);
    Real force_mag = rand.GetRandom<Real>(kMaxForce, kMinForce);
    P(force_mag);
    direction *= force_mag;
    applied_force[v] = direction;
    applied_verts.insert(v);
    ++i;
  }
  L("random force genrated");
  std::ofstream out(pose_file);
  out << pos_num << " " << global_subspace_coord_num_ << std::endl;
  int pose_idx = 0;
  for (int v : applied_verts) {
    P(pose_idx);
    force_vert.clear();
    force_vert.insert(v);
    force = applied_force[v];
    q_ = VectorX::Zero(global_subspace_coord_num_);
    vel_q_ = VectorX::Zero(global_subspace_coord_num_);
    for (int i = 0; i < 100; ++i) {
      SubspaceSimulation(global::time_step);
    }
    for (int i = 0; i < global_subspace_coord_num_; ++i) {
      out << q_[i] << " ";
    }
    ++pose_idx;
  }
  out.close();
}

SubspaceMassSpringVolumetricObject::~SubspaceMassSpringVolumetricObject()
{
  //  return;
  //  std::ofstream out(DATA_DIRECTORY"frames.txt");
  //  out << frames.size() << " ";
  //  out << v_num_ << " ";
  //  for (int i = 0; i < frames.size(); ++i) {
  //    for (int j = 0; j < v_num_; ++j) {
  //      out << frames[i][j][0] << " "
  //          << frames[i][j][1] << " "
  //          << frames[i][j][2] << " ";
  //    }
  //  }
  //  out.close();
  delete collider_;
  delete conf_;
}

bool SubspaceMassSpringVolumetricObject::VerifySubspaceDiagonality(MatrixX &subspace)
{
  //  MatrixX mass_matrix = MatrixX::Zero(subspace.rows(), subspace.rows());
  auto tmp_subspace = subspace;
  for (int i = 0; i < subspace.rows(); ++i) {
    for (int j = 0; j < subspace.cols(); ++j) {
      tmp_subspace(i, j) *= mass_[i / 3];
    }
    //    mass_matrix(i * 3 + 0, i * 3 + 0) = mass_[i];
    //    mass_matrix(i * 3 + 1, i * 3 + 1) = mass_[i];
    //    mass_matrix(i * 3 + 2, i * 3 + 2) = mass_[i];
  }
  MatrixX zero = subspace.transpose() * tmp_subspace - MatrixX::Identity(subspace.cols(), subspace.cols());
  MatrixX::Scalar sum = 0;
  for (int i = 0; i < zero.rows(); i++) {
    for (int j = 0; j < zero.rows(); j++) {
      sum += zero(i, j) * zero(i, j);
    }
  }
  P(sum);
  return sum < 1e-6;
}

void SubspaceMassSpringVolumetricObject::Construct()
{
  conf_ = new ConfigFile(DATA_DIRECTORY "spring_param.conf");
  gravity_ = Vec3(global::gravity);
  stiffness_ = conf_->Get<Real>("stiffness");
  damping_ = conf_->Get<Real>("damping");
  P(stiffness_, damping_);
  prev_vert_.resize(v_num_);
  vel_ = std::vector<Vec3>(v_num_, Vec3(0, 0, 0));
  local_u_ = vel_;
  tmp_vert_buffer_.resize(v_num_);
  cg_solver_ = new ConjugateGradientSolver<Real>(v_num_ * 3);
  //  mass_.resize(v_num_);
  if (v_num_ <= 0) {
    std::cerr << CURRENT_LINE << " mass not initialized." << std::endl;
    exit(0);
  }
  rest_edge_length_.resize(e_num_);
  ComputeRestEdgeLength();
  //  ComputeTetrahedraVolume();
  //  ComputeLumpedMass();
  prev_vert_ = vert_;
  collider_ = new TetCollider(&vert_[0][0], v_num_, &tet_[0][0], tet_num_, &edge_[0][0], e_num_);
  if (1) {
    const Real kMass = 1;
    std::fill(mass_.begin(), mass_.end(), kMass);
    total_mass_ = kMass * v_num_;
  } else {
    TextFileReader in(DATA_DIRECTORY "vertex_mass.txt");
    for (int v = 0; v < kVertexNumPerSegment; ++v) {
      in >> mass_[v];
    }
    //    memcpy(&mass_[kVertexNumPerSegment], &mass_[0], sizeof(Real) * kVertexNumPerSegment);
  }
  segment_num_ = v_num_ / kVertexNumPerSegment;
}

void SubspaceMassSpringVolumetricObject::LoadMassMatrix(const char *mass_file)
{
  TextFileReader in(mass_file);
  int row, col;
  in >> row >> col;
  ASSERT(row == kVertexNumPerSegment);
}


void SubspaceMassSpringVolumetricObject::LoadSubspace(const char* file_name)
{
  int row, col;
#if 1
  TextFileReader file(file_name);
  file >> row >> col;
  subspace_coordinate_per_segment_ = col;
  subspace_.resize(row, col);
  for (int r = 0; r < row; ++r) {
    for (int c = 0; c < col; ++c) {
      file >> subspace_(r, c);
    }
  }
  int NUM = conf_->Get<int>("basis_num");
  subspace_coordinate_per_segment_ = NUM;
  MatrixX tmp(kVertexNumPerSegment * 3, NUM);
  for (int i = 0; i < NUM; ++i) {
    //    Vec3 center;
    //    AffineTransformer<Real>::Centerize(subspace_.col(i).data(), v_num_, &center[0]);
    tmp.col(i) = subspace_.col(i);
    //    AffineTransformer<Real>::GetCenter(subspace_.col(i).data(), v_num_, &center[0]);
    //    P(center);
  }
  subspace_ = tmp;
  //  Real axis[] = {1, 1, 1};
  //  for (int i = 0; i < NUM; ++i) {
  //    AffineTransformer<Real>::Rotate(subspace_.col(i).data(), kVertexNumPerSegment, axis, 90);
  //  }
  //  ASSERT(VerifySubspaceDiagonality(subspace_));
  //  exit(0);
  col = NUM;
#else
  full_subspace_ = MatrixX::Zero(v_num_ * 3, 20 * 2);
  {
    TextFileReader file(DATA_DIRECTORY "subspace/basis_1st.txt");
    file >> row >> col;
    subspace_coordinate_per_segment_ = col;
    subspace_.resize(row, col);
    for (int r = 0; r < row; ++r) {
      for (int c = 0; c < col; ++c) {
        file >> subspace_(r, c);
        full_subspace_(r, c) = subspace_(r, c);
      }
    }
    int NUM = conf_->Get<Real>("basis_num");
    subspace_coordinate_per_segment_ = NUM;
    MatrixX tmp(kVertexNumPerSegment * 3, NUM);
    //    for (int i = 0; i < NUM; ++i) {
    //      tmp.col(i) = subspace_.col(i);
    //    }
    for (int r = 0; r < row; ++r) {
      for (int c = 0; c < col; ++c) {
        file >> subspace_(r, c);
        full_subspace_(r, c) = subspace_(r, c);
      }
    }
    subspace_ = tmp;
    col = NUM;
  }

  {
    //  TextFileReader file(file_name);
    TextFileReader file(DATA_DIRECTORY "subspace/basis_2nd.txt");
    file >> row >> col;
    subspace_coordinate_per_segment_ = col;
    subspace_.resize(row, col);
    for (int r = 0; r < row; ++r) {
      for (int c = 0; c < col; ++c) {
        file >> subspace_(r, c);
        full_subspace_(kVertexNumPerSegment * 3 + r, 20 + c) = subspace_(r, c);
      }
    }
    int NUM = conf_->Get<Real>("basis_num");
    subspace_coordinate_per_segment_ = NUM;
    MatrixX tmp(kVertexNumPerSegment * 3, NUM);
    for (int i = 0; i < NUM; ++i) {
      tmp.col(i) = subspace_.col(i);
    }
    subspace_ = tmp;
    col = NUM;
  }
#endif
  ASSERT(IsColumnMajor(subspace_));
  ASSERT(VerifySubspaceDiagonality(subspace_));
  //  for (int i = 0; i < subspace_coordinate_per_segment_; ++i) {
  //    Vec3 center(0, 0, 0);
  //    AffineTransformer<Real>::GetCenter(subspace_.col(i).data(), kVertexNumPerSegment, &center[0]);
  //    AffineTransformer<Real>::Centerize(subspace_.col(i).data(), kVertexNumPerSegment, &center[0]);
  //    P(i, center);
  //  }
  subspace_transpose_ = subspace_.transpose();
  ConstructFullSubspace();
  ASSERT(VerifySubspaceDiagonality(full_subspace_));
  global_subspace_coord_num_ = subspace_coordinate_per_segment_ * segment_num_;
  InitializeSubspaceData();
  BuildMomentumMatrix();
  P(row, col);
  P(subspace_coordinate_per_segment_);
  P(segment_num_);
  //  P(rest_pos_[0]);
  //  ASSERT(segment_num_ == 1);
}

void SubspaceMassSpringVolumetricObject::BuildMomentumMatrix()
{
  //  reduced_linear_momentum_matrix_ = MatrixX::Zero(3, global_subspace_coord_num_);
  //  reduced_angular_momentum_matrix_[0] = MatrixX::Zero(3, global_subspace_coord_num_);
  //  reduced_angular_momentum_matrix_[1] = MatrixX::Zero(3, global_subspace_coord_num_);
  //  reduced_angular_momentum_matrix_[2] = MatrixX::Zero(3, global_subspace_coord_num_);
  reduced_linear_momentum_matrix_ = MatrixX::Zero(3, subspace_coordinate_per_segment_);
  for (int v = 0; v < kVertexNumPerSegment; ++v) {
    for (int xyz = 0; xyz < 3; ++xyz) {
      reduced_linear_momentum_matrix_.row(xyz) += mass_[v] * subspace_.row(v * 3 + xyz);
      //      reduced_angular_momentum_matrix_[0].row(xyz) += mass_[v] * rest_pos_[v][0] * full_subspace_.row(v * 3 + xyz);
      //      reduced_angular_momentum_matrix_[1].row(xyz) += mass_[v] * rest_pos_[v][1] * full_subspace_.row(v * 3 + xyz);
      //      reduced_angular_momentum_matrix_[2].row(xyz) += mass_[v] * rest_pos_[v][2] * full_subspace_.row(v * 3 + xyz);
    }
  }
}

void SubspaceMassSpringVolumetricObject::SubspaceSimulation(Real dt)
{
  //  return;
  if (0) {
    static int c = 0;
    q_ = VectorX::Zero(q_.size());
    q_(c) = 1;
    //    MapVectorX v(&vert_[0][0], v_num_ * 3);
    c = (c + 1) % full_subspace_.cols();
    ReconstructObjectFromSubspace();
    return;
  }
  using namespace dj;
  std::vector<Real> force_gradient(e_num_ * 9, 0);
  ComputeForceGradient(force_gradient);
  ComputeForce(tmp_vert_buffer_, dt);
  const Real t2 = dt * dt;
  const Real kDampingFactor = damping_ * dt;
  (void) kDampingFactor;
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    //    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    //    OMP_FOR
    //    for (int v = 0; v < v_num_; ++v) {
    //      int v3 = v * 3;
    //      const float factor = kDampingFactor / mass_[v];
    //      result[v3 + 0] += factor * x[v3 + 0];
    //      result[v3 + 1] += factor * x[v3 + 1];
    //      result[v3 + 2] += factor * x[v3 + 2];
    //    }

    for (int e = 0; e < e_num_; ++e) {
      Real* fg = &force_gradient[e * 9];
      int* verts = &edge_[e][0];
      Real* pos[] = {
        x + verts[0] * 3,
        x + verts[1] * 3,
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      Real tmp[3];
      dj::MulMatrix3x3Vec<Real>((Real (*)[3]) fg, diff, tmp);
#if 1
      result[verts[0] * 3 + 0] += -t2 * tmp[0];
      result[verts[0] * 3 + 1] += -t2 * tmp[1];
      result[verts[0] * 3 + 2] += -t2 * tmp[2];

      result[verts[1] * 3 + 0] -= -t2 * tmp[0];
      result[verts[1] * 3 + 1] -= -t2 * tmp[1];
      result[verts[1] * 3 + 2] -= -t2 * tmp[2];
#else
      result[verts[0] * 3 + 0] += -t2 * tmp[0] / mass_[verts[0]];
      result[verts[0] * 3 + 1] += -t2 * tmp[1] / mass_[verts[0]];
      result[verts[0] * 3 + 2] += -t2 * tmp[2] / mass_[verts[0]];

      result[verts[1] * 3 + 0] -= -t2 * tmp[0] / mass_[verts[1]];
      result[verts[1] * 3 + 1] -= -t2 * tmp[1] / mass_[verts[1]];
      result[verts[1] * 3 + 2] -= -t2 * tmp[2] / mass_[verts[1]];
#endif
    }
  };

  // Compute reduced stiffness matrix
  MatrixX k_dot_u(v_num_ * 3, global_subspace_coord_num_);
  for (int col = 0; col < full_subspace_.cols(); ++col) {
    StiffnessMatrix(full_subspace_.col(col).data(), k_dot_u.col(col).data());
  }
  //  P(k_dot_u); exit(0);
  MatrixX reduced_stiffness(global_subspace_coord_num_, global_subspace_coord_num_);
  reduced_stiffness = full_subspace_transpose_ * k_dot_u + MatrixX::Identity(global_subspace_coord_num_, global_subspace_coord_num_);
  //  Super::Simulate(dt);

  // Compute reduced rhs
  VectorX rhs(global_subspace_coord_num_);
  MapVectorX unreduced_force(&tmp_vert_buffer_[0][0], v_num_ * 3);
  //  VectorX k_dot_translation(v_num_ * 3);
  //  StiffnessMatrix(&rest_pos_[0][0], k_dot_translation.data());
  //  unreduced_force -= k_dot_translation;
  rhs = full_subspace_transpose_ * unreduced_force + vel_q_;
  //  VectorX det_q = vel_q_;
  //  VectorX  prev_vel_q = vel_q_;
  vel_q_ = reduced_stiffness.colPivHouseholderQr().solve(rhs);
  //  WriteVectorToMatlab(global_subspace_coord_num_, vel_q_.data(), std::cout);
  // Remove translational mode
  //  VectorX det_q = vel_q_ - prev_vel_q;
  //  if (0)
  //    for (int seg = 0; seg < segment_num_; ++seg) {
  //      MapVectorX sub_det_q(det_q.data() + seg * subspace_coordinate_per_segment_, subspace_coordinate_per_segment_);
  //      VectorX rhs = reduced_linear_momentum_matrix_ * sub_det_q;// - net_ext_force_[seg];
  //      MatrixX M = reduced_linear_momentum_matrix_ * reduced_linear_momentum_matrix_.transpose();
  //      VectorX lambda =  M.colPivHouseholderQr().solve(rhs);
  //      sub_det_q = sub_det_q - reduced_linear_momentum_matrix_.transpose() * lambda;
  //    }
  //  vel_q_ = det_q + prev_vel_q;
  //  vel_q_ *= 0.95;
  q_ += dt * vel_q_;
  //  L("Q");
  //  std::cout << q_ << std::endl;

  // Apply rigid motion
  if (0) {
    Vector3 net_force(0, 0, 0);
    for (int seg = 0; seg < segment_num_; ++seg) {
      net_force += net_ext_force_[seg];
    }
    translational_vel_[0] += net_force[0] * dt / total_mass_;
    translational_vel_[1] += net_force[1] * dt / total_mass_;
    translational_vel_[2] += net_force[2] * dt / total_mass_;
    for (int seg = 0; seg < segment_num_; ++seg) {
      segment_center_[seg][0] += translational_vel_[0] * dt;
      segment_center_[seg][1] += translational_vel_[1] * dt;
      segment_center_[seg][2] += translational_vel_[2] * dt;
    }
  }
  //  for (int i = 0; i < subspace_coordinate_per_segment_; ++i) {
  //    q_[i + subspace_coordinate_per_segment_] = q_[i];
  //  }
  ReconstructObjectFromSubspace();
  //  UpdatePosition();
}

void SubspaceMassSpringVolumetricObject::ImplicitStepOneSegment(Real dt)
{
  using namespace dj;
  std::vector<Real> force_gradient(e_num_ * 9, 0);
  ComputeForceGradient(force_gradient);
  ComputeForce(tmp_vert_buffer_, dt);
  const Real t2 = dt * dt;
  const Real kDampingFactor = damping_ * dt;

  auto StiffnessMatrix = [&](Real * x, Real * result) {
    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      const Real factor = kDampingFactor / mass_[v];
      result[v3 + 0] += factor * x[v3 + 0];
      result[v3 + 1] += factor * x[v3 + 1];
      result[v3 + 2] += factor * x[v3 + 2];
      Real tmp_result[3] = {0, 0, 0};
      for (int num = 0; num < (int) incident_edge_[v].size(); ++num) {
        int e = incident_edge_[v][num];
        Real* fg = &force_gradient[e * 9];
        int* verts = &edge_[e][0];
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

  MatrixX k_dot_u(v_num_ * 3, subspace_coordinate_per_segment_);
  Real* data = k_dot_u.data();
  Real* fullsubspace_data = full_subspace_.data();
  for (int i = 0; i < subspace_coordinate_per_segment_; ++i) {
    StiffnessMatrix(fullsubspace_data + i * (v_num_ * 3), data + i * (v_num_ * 3));
  }
  MatrixX reduced_stiffness = full_subspace_transpose_ * k_dot_u;
  // Assemble right hand side
  //  memcpy(&tmp_vert_buf_[0][0], &vel_[0][0], sizeof(Vec3) * v_num_);
  VectorX net_ext_force(3);
  net_ext_force[0] = 0;
  net_ext_force[1] = 0;
  net_ext_force[2] = 0;
  VectorX rhs(subspace_coordinate_per_segment_);
  Eigen::Map<VectorX> unreduced_force(&tmp_vert_buffer_[0][0], v_num_ * 3);
  rhs = full_subspace_transpose_ * unreduced_force + vel_q_;
  VectorX prev_vel_q_ = vel_q_;
  vel_q_ = reduced_stiffness.colPivHouseholderQr().solve(rhs);
  MatrixX U = MatrixX::Zero(3, vel_q_.size());
  for (int i = 0; i < full_subspace_.rows(); i += 3) {
    U.row(0) += full_subspace_.row(i + 0);
    U.row(1) += full_subspace_.row(i + 1);
    U.row(2) += full_subspace_.row(i + 2);
  }
  VectorX det_q = vel_q_ - prev_vel_q_;
  {
    VectorX rhs = U * det_q - net_ext_force;
    MatrixX M = U * U.transpose();
    VectorX lambda =  M.colPivHouseholderQr().solve(rhs);
    det_q = det_q - U.transpose() * lambda ;
  }
  vel_q_ = det_q + prev_vel_q_;
  //  vel_q_ *= (1 - damping_);
  q_ += vel_q_ * dt;
  Eigen::Map<VectorX> vert(&vert_[0][0], v_num_ * 3);
  vert = full_subspace_ * q_;
}

void SubspaceMassSpringVolumetricObject::ComputeForceGradient(std::vector<Real> &force_gradient)
{
  // Compute force gradient for each edge
  // OMP_FOR
  for (int e = 0; e < e_num_; ++e) {
    int e9 = e * 9;
    int* v = &edge_[e][0];
    Real* pos[] = {
      &vert_[v[0]][0],
      &vert_[v[1]][0],
    };
    Real direction[3];
    dj::SubVec3(pos[1], pos[0], direction);
    Real length = dj::Normalize3(direction);
    if (length < rest_edge_length_[e]) {
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
          tmp[i] *= stiffness_;
        }
      }
    } else {
      Real factor = stiffness_ * rest_edge_length_[e] / length;
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
      matrix[0][0] += stiffness_;
      matrix[1][1] += stiffness_;
      matrix[2][2] += stiffness_;
    }
  }
}

void SubspaceMassSpringVolumetricObject::ComputeForce(std::vector<Vec3> &forces, Real dt)
{
  memset(&forces[0][0], 0, sizeof(Vec3) * v_num_);
  const Real kFloor = 0.02;
  const Real kFloorStiffness = 5;
  //  external_force_.clear();
  //  for (int i = 0; i < segment_num_; ++i) {
  //    net_ext_force_
  //  }
  //  memset(&net_ext_force_[0][0], 0, sizeof(Vec3) * segment_num_);
  if (0) {
    Real mag = force.Magnitude();
    Vec3 pos = vert_[v];
    Vec3 center = pos;
    center[0] = 0;
    center[1] = 0;
    force = center - pos;
    force.Normalize();
    force *= mag;
    forces[v] += force * dt;
  }
  if (0) {
    Real mag = force.Magnitude();
    for (int v = kCrossSectionVerextNum * ring; v < kCrossSectionVerextNum * (ring + 1) - 80; ++v) {
      Vec3 pos = vert_[v];
      Vec3 center = pos;
      center[0] = 0;
      center[1] = 0;
      auto force = center - pos;
      force.Normalize();
      force *= mag;
      forces[v] += force * dt;
    }
  }
  if (0) {
    for (int v : force_vert) {
      //      forces[v][0] = force[0] * dt / mass_[v];
      //      forces[v][1] = force[1] * dt / mass_[v];
      //      forces[v][2] = force[2] * dt / mass_[v];

      forces[v][0] = force[0] * dt;
      forces[v][1] = force[1] * dt;
      forces[v][2] = force[2] * dt;
    }
  }
  OMP_FOR
  for (int v_idx = 0; v_idx < v_num_; ++v_idx) {
    Real* single_force = &forces[v_idx][0];
    if (1) {
      single_force[0] += dt * gravity_[0] * mass_[v_idx];
      single_force[1] += dt * gravity_[1] * mass_[v_idx];
      single_force[2] += dt * gravity_[2] * mass_[v_idx];

      //      net_ext_force_[v_idx / kVertexNumPerSegment][0] += dt * gravity_[0] * mass_[v_idx];
      //      net_ext_force_[v_idx / kVertexNumPerSegment][1] += dt * gravity_[1] * mass_[v_idx];
      //      net_ext_force_[v_idx / kVertexNumPerSegment][2] += dt * gravity_[2] * mass_[v_idx];

      if (vert_[v_idx][1] < kFloor) {
        Real force_magnitude = (kFloor - vert_[v_idx][1]) * mass_[v_idx] * kFloorStiffness;
        single_force[1] += force_magnitude;
        //        net_ext_force_[v_idx / kVertexNumPerSegment][1] += force_magnitude * mass_[v_idx];
        //      external_force_.emplace_back(Vec3(0, force_magnitude, 0), v_idx);
      }
    }
  }
  //    vel_[v_idx] *= kDamping_;
  //    force[0] += vel_[v_idx][0];
  //    force[1] += vel_[v_idx][1];
  //    force[2] += vel_[v_idx][2];

  //    for (int num = 0; num < (int) incident_edge_[v_idx].size(); ++num) {
  for (int e = 0; e < e_num_; ++e) {
    //    int e = incident_edge_[v_idx][num];
    int* verts = &edge_[e][0];
    Real* pos[] = {
      &vert_[verts[0]][0],
      &vert_[verts[1]][0],
    };
    Real* single_force[2] = {
      &forces[verts[0]][0],
      &forces[verts[1]][0],
    };
    Real diff[3];
    dj::SubVec3(pos[1], pos[0], diff);
    Real length = dj::Normalize3(diff);
#if 1
    single_force[0][0] -= -stiffness_ * diff[0] * (length - rest_edge_length_[e]) * dt;
    single_force[0][1] -= -stiffness_ * diff[1] * (length - rest_edge_length_[e]) * dt;
    single_force[0][2] -= -stiffness_ * diff[2] * (length - rest_edge_length_[e]) * dt;

    single_force[1][0] += -stiffness_ * diff[0] * (length - rest_edge_length_[e]) * dt;
    single_force[1][1] += -stiffness_ * diff[1] * (length - rest_edge_length_[e]) * dt;
    single_force[1][2] += -stiffness_ * diff[2] * (length - rest_edge_length_[e]) * dt;
#else
    single_force[0][0] -= -stiffness_ * diff[0] * (length - rest_edge_length_[e]) * dt / mass_[verts[0]];
    single_force[0][1] -= -stiffness_ * diff[1] * (length - rest_edge_length_[e]) * dt / mass_[verts[0]];
    single_force[0][2] -= -stiffness_ * diff[2] * (length - rest_edge_length_[e]) * dt / mass_[verts[0]];

    single_force[1][0] += -stiffness_ * diff[0] * (length - rest_edge_length_[e]) * dt / mass_[verts[1]];
    single_force[1][1] += -stiffness_ * diff[1] * (length - rest_edge_length_[e]) * dt / mass_[verts[1]];
    single_force[1][2] += -stiffness_ * diff[2] * (length - rest_edge_length_[e]) * dt / mass_[verts[1]];
#endif
  }
}

void SubspaceMassSpringVolumetricObject::SimulateInLocalFrame(Real dt)
{
  //    RigidBody::Simulate(dt);return;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // rigid motion, get t->t' ->t'', w
  // get external force and torque
//  gravity_[0] = 0;
//  gravity_[1] = 0;
//  gravity_[2] = 0;
  Mat3 rotation_matrix;
  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  Vec3 net_force = gravity_ * total_mass_;
  Vec3 net_torque(Real(0), Real(0), Real(0));
  //    RigidBody::ApplyGroundCollisionForce(dt, net_force, net_torque);
  std::vector<Vec3> ext_force(v_num_ * 3);// Vec3(0, 0, 0));
  std::vector<Vec3> int_force(v_num_ * 3);// Vec3(0, 0, 0));
  memset(&ext_force[0][0], 0, sizeof(Vec3) * v_num_);
  memset(&int_force[0][0], 0, sizeof(Vec3) * v_num_);
  if (1) {
    const Real kFloor = 0.02;
    const Real kFloorStiffness = 900;
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      if (1) {
        ext_force[v][0] += gravity_[0] * mass_[v];
        ext_force[v][1] += gravity_[1] * mass_[v];
        ext_force[v][2] += gravity_[2] * mass_[v];

        if (vert_[v][1] < kFloor) {
          Real force_magnitude = (kFloor - vert_[v][1]) * mass_[v] * kFloorStiffness;
          ext_force[v][1] += force_magnitude;

          Vec3 collision_force(0, force_magnitude, 0);
          net_force += collision_force;
          Vec3 tmp_torque;
          Vec3 r = vert_[v] - center_of_mass_;
          dj::Cross3(&r[0], &collision_force[0], &tmp_torque[0]);
          net_torque += tmp_torque;
        }
      }
    }

    for (int e = 0; e < e_num_; ++e) {
      int* verts = &edge_[e][0];
      Real pos[][3] = {
        rest_pos_[verts[0]][0] + local_u_[verts[0]][0],
        rest_pos_[verts[0]][1] + local_u_[verts[0]][1],
        rest_pos_[verts[0]][2] + local_u_[verts[0]][2],

        rest_pos_[verts[1]][0] + local_u_[verts[1]][0],
        rest_pos_[verts[1]][1] + local_u_[verts[1]][1],
        rest_pos_[verts[1]][2] + local_u_[verts[1]][2],
      };
      Real* single_force[2] = {
        &int_force[verts[0]][0],
        &int_force[verts[1]][0],
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      Real length = dj::Normalize3(diff);

      single_force[0][0] -= -stiffness_ * diff[0] * (length - rest_edge_length_[e]);
      single_force[0][1] -= -stiffness_ * diff[1] * (length - rest_edge_length_[e]);
      single_force[0][2] -= -stiffness_ * diff[2] * (length - rest_edge_length_[e]);

      single_force[1][0] += -stiffness_ * diff[0] * (length - rest_edge_length_[e]);
      single_force[1][1] += -stiffness_ * diff[1] * (length - rest_edge_length_[e]);
      single_force[1][2] += -stiffness_ * diff[2] * (length - rest_edge_length_[e]);
    }
  }
  // Linear momentum
  acceleration_ = net_force / total_mass_;
  translational_vel_ += acceleration_ * dt;
  translational_vel_ *= (1 - translational_damping_);
  //  translational_vel_.Fill(0);
  // Angular momentum
  Mat3 cur_inv_inertial = rotation_matrix.transpose() * inv_inertia_tensor_ * rotation_matrix;
  Vec3 angular_acceleration;
  Eigen::Map<Vector3> domega_wraper(&angular_acceleration[0]);
  Eigen::Map<Vector3> torque_wraper(&net_torque[0]);
  domega_wraper = cur_inv_inertial * torque_wraper;
  //  angular_acceleration[1] = 1.00;
  Vec3 change_in_angular_vel_ = angular_acceleration * dt;
  //  angular_acceleration.Fill(0);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // compute force
#if 1
  MapVector3 acc(&acceleration_[0]);
  if (0) {
    double angle;
    double axis[3];
    quaternion_.GetRotation(&angle, axis);
    P(angle);
    PV(axis, 3);
//    PMAT(rotation_matrix);
  }
  for (int v = 0; v < v_num_; ++v) {
    // inertial force: F= -m * a
    Vector3 fititious_force = -mass_[v] * acc;
    Vector3 r(vert_[v][0] - center_of_mass_[0], vert_[v][1] - center_of_mass_[1], vert_[v][2] - center_of_mass_[2]);
    if (1) {
      // Euler force: F = -m * dw/dt \times r
      Vector3 euler_force;
      dj::Cross3(&angular_acceleration[0], &r[0], &euler_force[0]);
      euler_force = -mass_[v] * euler_force;
      fititious_force += euler_force;
      //          P(euler_force.norm());
    }
    // centerifugal force: F = -m * w \times (w \times r)
    if (1)  {
      double tmp[3];
      dj::Cross3(&angular_vel_[0], &r[0], tmp);
      Vector3 centerifugal_force;
      dj::Cross3(&angular_vel_[0], &tmp[0], &centerifugal_force[0]);
      centerifugal_force *= -mass_[v];
      fititious_force += centerifugal_force;
    }
    // Coriolis force: F = -2m w\times v
    if (1) {
      Vector3 croilis_force;
      dj::Cross3(&angular_vel_[0], &vel_[v][0], &croilis_force[0]);
      fititious_force -= 2 * mass_[v] * croilis_force;
    }

    ext_force[v][0] += fititious_force[0];
    ext_force[v][1] += fititious_force[1];
    ext_force[v][2] += fititious_force[2];
    Eigen::Map<Vector3> map_ext_force(&ext_force[v][0]);
    // Rotate external force to reference frame
    map_ext_force = rotation_matrix * map_ext_force;
    int_force[v] += ext_force[v];
    int_force[v] *= dt;
    int_force[v] += mass_[v] * vel_[v];
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // local dinaymics
  std::vector<Real> force_gradient(e_num_ * 9, 0);
  ComputeForceGradient(force_gradient);
  const Real t2 = dt * dt;
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      result[v3 + 0] = mass_[v] * x[v3 + 0];
      result[v3 + 1] = mass_[v] * x[v3 + 1];
      result[v3 + 2] = mass_[v] * x[v3 + 2];
    }

    for (int e = 0; e < e_num_; ++e) {
      Real* fg = &force_gradient[e * 9];
      int* verts = &edge_[e][0];
      Real* pos[] = {
        x + verts[0] * 3,
        x + verts[1] * 3,
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      Real tmp[3];
      dj::MulMatrix3x3Vec<Real>((Real (*)[3]) fg, diff, tmp);

      result[verts[0] * 3 + 0] += -t2 * tmp[0];
      result[verts[0] * 3 + 1] += -t2 * tmp[1];
      result[verts[0] * 3 + 2] += -t2 * tmp[2];

      result[verts[1] * 3 + 0] -= -t2 * tmp[0];
      result[verts[1] * 3 + 1] -= -t2 * tmp[1];
      result[verts[1] * 3 + 2] -= -t2 * tmp[2];
    }
  };

  cg_solver_->Solve(&int_force[0][0], &vel_[0][0], StiffnessMatrix, 1000, 1e-10);
  for (int v = 0; v < v_num_; ++v) {
        local_u_[v] += vel_[v] * dt;
        vel_[v] *= 0.98;
    //    local_u_[v] = Vec3(0, 0, 0);
  }
#endif
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  center_of_mass_ += translational_vel_ * dt;
  //  angular_vel_ = Vec3(0, 1, 0);
  angular_vel_ += change_in_angular_vel_;
  angular_vel_ *= (1 - rotational_damping_);
  quaternion_ = quaternion_ + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[0], angular_vel_[1], angular_vel_[2]) * quaternion_;
  quaternion_.Normalize();

  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  rotation_matrix.transposeInPlace();
  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
    Vec3 offset;
    Eigen::Map<Vector3> offset_wrapper(&offset[0]);
    Eigen::Map<Vector3> rest_pos_wrapper(&rest_pos_[v][0]);
    Eigen::Map<Vector3> displacement_wrapper(&local_u_[v][0]);
    offset_wrapper = rotation_matrix * (rest_pos_wrapper + displacement_wrapper);
//        offset_wrapper = (rest_pos_wrapper + displacement_wrapper);
    vert_[v] = offset + center_of_mass_;
  }
  P(center_of_mass_);
}


void SubspaceMassSpringVolumetricObject::Simulate(Real dt)
{
  //  RigidBody::Simulate(dt); return;
  //  EXECUTE_TIMES(400);
  //  ImplicitStepOneSegment(dt);
  SimulateInLocalFrame(dt);
  //  ViewSamlePose("/tmp/q.txt");
  //  GenerateSamplePose("/tmp/q.txt", 200);
  //  SubspaceSimulation(dt);
  //  P(e_num_);
  //        FullSimulation(dt);
  //  collider_->HandleCollision(&vert_[0][0], &tmp_vert_buf_[0][0]);
  //  vert_.swap(tmp_vert_buf_);;
  //  for (int v = 0; v < v_num_; ++v) {
  //    vel_[v] += (vert_[v] - tmp_vert_buf_[v]) / dt;
  //  }
}

void SubspaceMassSpringVolumetricObject::ReconstructObjectFromSubspace()
{
  MapVectorX vert_map(&vert_[0][0], v_num_ * 3);
  MapVectorX rest_vert_map(&rest_pos_[0][0], v_num_ * 3);
  vert_map = full_subspace_ * q_;
  vert_map += rest_vert_map;
  return;
  for (int seg = 0; seg < segment_num_; ++seg) {
    Mat3 rotation_matrix;
    segment_rotation_[seg].Quaternion2Matrix(rotation_matrix.data());
    rotation_matrix.transposeInPlace();
    MapVectorX displacement(&tmp_vert_buffer_[seg * kVertexNumPerSegment][0], kVertexNumPerSegment * 3);
    MapVectorX sub_q_(&q_[seg * subspace_coordinate_per_segment_], subspace_coordinate_per_segment_);
    // TODO delete
    if (seg == 0) {
      displacement = subspace_ * sub_q_;
    } else {
      auto tmp_space_ = full_subspace_.block(kVertexNumPerSegment * 3, subspace_coordinate_per_segment_, kVertexNumPerSegment * 3, subspace_coordinate_per_segment_);
      displacement = tmp_space_ * sub_q_;
    }
    for (int v = 0; v < kVertexNumPerSegment; ++v) {
      MapVector3 one_displacement(&tmp_vert_buffer_[seg * kVertexNumPerSegment + v][0]);
      MapVector3 vert(&vert_[seg * kVertexNumPerSegment + v][0]);
      MapVector3 rest_vert_pos(&rest_pos_[seg * kVertexNumPerSegment + v][0]);
      // TODO: add rotation
      //      vert = rotation_matrix * one_vert + rest_vert_pos;
      vert = one_displacement + rest_vert_pos;
      //      P(one_vert, rest_vert_pos, vert, rest_pos_[0]);
      //      vert_[seg * kVertexNumPerSegment + v] += segment_center_[seg];
    }
  }
}

void SubspaceMassSpringVolumetricObject::ConstructFullSubspace()
{
  full_subspace_ = MatrixX::Zero(v_num_ * 3, segment_num_ * subspace_coordinate_per_segment_);
  full_subspace_transpose_ = MatrixX::Zero(segment_num_ * subspace_coordinate_per_segment_, v_num_ * 3);
#if 0
  int col = 0;
  for (int seg = 0; seg < 1; ++seg) {
    for (int i = 0; i < subspace_coordinate_per_segment_; ++i, col++) {
      Real* column = full_subspace_.col(col).data();
      Real* subcolumn = subspace_.col(i).data();
      memcpy(column + seg * kVertexNumPerSegment * 3, subcolumn, sizeof(Real) * kVertexNumPerSegment * 3);
    }
  }

  for (int seg = 1; seg < 2; ++seg) {
    for (int i = 0; i < subspace_coordinate_per_segment_; ++i, col++) {
      Real* column = full_subspace_.col(col).data() + seg * kVertexNumPerSegment * 3;
      Real* subcolumn = subspace_.col(i).data();
      for (int j = 0; j < kSegmentCrossSectionNum; ++j) {
        memcpy(column + (kSegmentCrossSectionNum - 1 - j) * kCrossSectionVerextNum * 3,
               subcolumn + j * kCrossSectionVerextNum * 3,
               sizeof(Real) * kCrossSectionVerextNum * 3);
      }
    }
  }
#else
  for (int seg = 0, col = 0; seg < segment_num_; ++seg) {
    for (int i = 0; i < subspace_coordinate_per_segment_; ++i, col++) {
      Real* column = full_subspace_.col(col).data();
      Real* subcolumn = subspace_.col(i).data();
      memcpy(column + seg * kVertexNumPerSegment * 3, subcolumn, sizeof(Real) * kVertexNumPerSegment * 3);
    }
  }
#endif
  full_subspace_transpose_ = full_subspace_.transpose();
}

void SubspaceMassSpringVolumetricObject::ComputeRestEdgeLength()
{
  for (int e = 0; e < e_num_; ++e) {
    Vec3 diff = vert_[edge_[e][0]] - vert_[edge_[e][1]];
    rest_edge_length_[e] = diff.Magnitude();
  }
}

void SubspaceMassSpringVolumetricObject::InitializeSubspaceData()
{
  vel_q_ = VectorX::Zero(global_subspace_coord_num_);
  q_ = VectorX::Zero(global_subspace_coord_num_);
  segment_mass_.resize(segment_num_);
  net_ext_force_.reserve(segment_num_);
  //  MapVectorX map_vert(&vert_[0][0], v_num_ * 3);
  //  q_ = full_subspace_transpose_ * map_vert;
  //  map_vert = full_subspace_ * q_;
  //  return;
  segment_center_.resize(segment_num_);
  segment_rotation_.resize(segment_num_);
  for (int seg = 0; seg < segment_num_; ++seg) {
    //    Real* pos = &vert_[seg * kVertexNumPerSegment][0];
    //    MapVectorX pos_map(pos, kVertexNumPerSegment * 3);
    //    AffineTransformer<Real>::GetCenter(pos, kVertexNumPerSegment, &segment_center_[seg][0]);
    //    AffineTransformer<Real>::Translate(pos, kVertexNumPerSegment, -segment_center_[seg][0], -segment_center_[seg][1], -segment_center_[seg][2]);
    //    MapVectorX sub_q_(&q_[seg * subspace_coordinate_per_segment_], subspace_coordinate_per_segment_);
    //    sub_q_ = subspace_transpose_ * pos_map;
    //    pos_map = subspace_ * sub_q_;
    //    AffineTransformer<Real>::Translate(pos, kVertexNumPerSegment, segment_center_[seg][0], segment_center_[seg][1], segment_center_[seg][2]);
    segment_rotation_[seg] = Quaternion<Real>(1, 0, 0, 0);
    segment_mass_[seg] = 0;
    for (int v = 0; v < kVertexNumPerSegment; ++v) {
      segment_mass_[seg] += mass_[seg * kVertexNumPerSegment + v];
    }
    //    Real axis[] = {1, 1, 1};
    //    dj::Normalize3(axis);
    //    segment_rotation_[seg] = Quaternion<Real>(45, axis);
  }
}

void SubspaceMassSpringVolumetricObject::UpdatePosition()
{
  Mat3 rotation_matrix;
  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  rotation_matrix.transposeInPlace();
  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
#if 1
    Vec3 offset;
    Eigen::Map<Vector3> offset_wrapper(&offset[0]);
    Eigen::Map<Vector3> rest_pos_wrapper(&rest_pos_[v][0]);
    Eigen::Map<Vector3> displacement_wrapper(&displacement_[v][0]);
    offset_wrapper = rotation_matrix * (rest_pos_wrapper + displacement_wrapper);
    vert_[v] = offset + center_of_mass_;
#else
    vert_[v] = rest_pos_[v] + center_of_mass_;
#endif
  }
}

void SubspaceMassSpringVolumetricObject::Render(int render_mode)
{
  //  glDisable(GL_LIGHTING);
  //    DrawArrow(&vert_[applied_vert][0], force(), true);
  //  glBegin(GL_POINTS);
  //  Vertex3v(&center_of_mass_[0]);
  //  glEnd();
  if (0) {
    glDisable(GL_LIGHTING);
    glPointSize(6);
    glColor3fv(kRed());
    glBegin(GL_POINTS);
    for (int v : force_vert) {
      Vertex3v(&vert_[v][0]);
    }
    //    for (int v = kCrossSectionVerextNum * ring; v < kCrossSectionVerextNum * (ring + 1) - 80; ++v) {
    //      Vertex3v(&vert_[v][0]);
    //    }
    glEnd();
  }


  Super::Render(render_mode);
}

void SubspaceMassSpringVolumetricObject::FullSimulation(Real dt)
{
  //  ASSERT(false, L("Need to change the stiffness matrix and rhs of equation to work!!!"));
  using namespace dj;
  std::vector<Real> force_gradient(e_num_ * 9, 0);
  ComputeForceGradient(force_gradient);
  const Real t2 = dt * dt;
  const Real kDampingFactor = damping_ * dt;
  UNUSED(kDampingFactor);
  /*
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      const Real factor = kDampingFactor / mass_[v];
      result[v3 + 0] += factor * x[v3 + 0];
      result[v3 + 1] += factor * x[v3 + 1];
      result[v3 + 2] += factor * x[v3 + 2];
      Real tmp_result[3] = {0, 0, 0};
      for (int num = 0; num < (int) incident_edge_[v].size(); ++num) {
        int e = incident_edge_[v][num];
        Real* fg = &force_gradient[e * 9];
        int* verts = &edge_[e][0];
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
  */
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    //    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      result[v3 + 0] = mass_[v] * x[v3 + 0];
      result[v3 + 1] = mass_[v] * x[v3 + 1];
      result[v3 + 2] = mass_[v] * x[v3 + 2];
      //          const float factor = kDampingFactor / mass_[v];
      //          result[v3 + 0] += factor * x[v3 + 0];
      //          result[v3 + 1] += factor * x[v3 + 1];
      //          result[v3 + 2] += factor * x[v3 + 2];
    }

    for (int e = 0; e < e_num_; ++e) {
      Real* fg = &force_gradient[e * 9];
      int* verts = &edge_[e][0];
      Real* pos[] = {
        x + verts[0] * 3,
        x + verts[1] * 3,
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      //      Real sign = (verts[0] == v) ? 1 : -1;
      Real tmp[3];
      dj::MulMatrix3x3Vec<Real>((Real (*)[3]) fg, diff, tmp);
      result[verts[0] * 3 + 0] += -t2 * tmp[0];
      result[verts[0] * 3 + 1] += -t2 * tmp[1];
      result[verts[0] * 3 + 2] += -t2 * tmp[2];

      result[verts[1] * 3 + 0] -= -t2 * tmp[0];
      result[verts[1] * 3 + 1] -= -t2 * tmp[1];
      result[verts[1] * 3 + 2] -= -t2 * tmp[2];
    }
  };
  // Assemble right hand side
  //  memcpy(&tmp_vert_buffer_[0][0], ivel_[0][0], sizeof(Vec3) * v_num_);
  ComputeForce(tmp_vert_buffer_, dt);
  for (int v = 0; v < v_num_; ++v) {
    tmp_vert_buffer_[v] += mass_[v] * vel_[v];
  }
  cg_solver_->Solve(&tmp_vert_buffer_[0][0], &vel_[0][0], StiffnessMatrix, 1000, 1e-6);
  OMP_FOR
  for (int v_idx = 0; v_idx < v_num_; ++v_idx) {
    vert_[v_idx] += vel_[v_idx] * dt;
    //    if (vert_[v_idx][1] < 0) {
    //      vel_[v_idx][1] -= vert_[v_idx][1] / dt;
    //      vert_[v_idx][1] = 0;
    //    }
  }
}
