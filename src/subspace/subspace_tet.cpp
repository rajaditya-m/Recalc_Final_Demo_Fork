//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include "affine_transformer.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <unordered_set>
#include "global.h"
#include "print_macro.h"
#include "binary_file_io.h"
#include "vector_lib.h"
#include "sparseMatrix.h"
#include "tet_mesh_simulator_bridge.h"
#include "subspace_tet.h"
#include "vector_io.h"
#include "generateMassMatrix.h"
#include "config_file.h"
#include "basis_io.h"
#include "matrixIO.h"
#include "polarDecomposition.h"


namespace dummy {
int dummy = 0;
}


SubspaceTet::SubspaceTet(const char *filename, double _limit_threshold, AffineTransformer<double> *affine_transformer,bool initialize_fem_module,bool load_interface, int split)
  : Super(filename, _limit_threshold, affine_transformer, initialize_fem_module, load_interface,split)
  , quaternion_(1, 0, 0, 0) {
    sim_mode_ = 1; //1 means joining mode  2 means stiched mode....
    stop = 0;
    hasPreCompPDK_ = false;
  gravity_ = Vec3(global::gravity[0], global::gravity[1], global::gravity[2]);
  //  PVEC(gravity_);
  translation_acc_.setZero();
  translation_vel_.setZero();
  angular_acc_.setZero();
  angular_vel_.setZero();
  rotation_.setIdentity();
  current_basis_ = -1;
  magnitude_ = 0;
  GenerateMassMatrix::computeVertexMasses((VolumetricMesh*) inv_fem_->vega_mesh_, mass_);
  ComputeCenterOfMass();
  ComputeInertialTensor();
  ComputeInterfaceTetData();
  ComputeInterfaceCOMData();
  vert_offset_from_mass_center_.resize(vertex_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    vert_offset_from_mass_center_[v] = MapVec3(rest_pos_ + v * 3) - center_of_mass_;
  }
}


void SubspaceTet::ComputeInterfaceTetData() {
    for(int t=0;t < tet_number; t++ ) {
        int v1 = tet_[4*t+0];
        int v2 = tet_[4*t+1];
        int v3 = tet_[4*t+2];
        int v4 = tet_[4*t+3];
        //Search interface vertices now
        for(auto it = interface_vertices_.begin(); it!= interface_vertices_.end(); it++) {
            int idx1 = it->first;
            int idx2 = it->second;
            if(idx1 == v1 || idx1 == v2 || idx1 == v3 || idx1 == v4) {
                interface_1_tets.push_back(t);
            }
        }

        //
    }
    P(interface_1_tets.size());
    for(auto it = interface_1_tets.begin(); it!= interface_1_tets.end(); it++) {
        int t = *it;
        int v1 = tet_[4*t+0];
        int v2 = tet_[4*t+1];
        int v3 = tet_[4*t+2];
        int v4 = tet_[4*t+3];
        interface_1_verts_.insert(v1);
        interface_1_verts_.insert(v2);
        interface_1_verts_.insert(v3);
        interface_1_verts_.insert(v4);
    }
}

void SubspaceTet::ComputeInterfaceCOMData() {
   interface_1_restCOM_.setZero();
   double total_mass = 0;
   for(auto it = interface_1_verts_.begin(); it != interface_1_verts_.end(); it++) {
      int v = *it;
      interface_1_restCOM_ += mass_[v] * MapVec3(rest_pos_ + v * 3);
      total_mass += mass_[v];
    }
    interface_1_restCOM_ /= total_mass;
    total_mass = 0;
    for(int v = domain_offset_toggle_;v < vertex_num_; v++) {
        domain_2_restCOM_ += mass_[v] * MapVec3(rest_pos_ + v * 3);
        total_mass += mass_[v];
    }
    domain_2_restCOM_ /= total_mass;
}

void SubspaceTet::ComputeCurrentInterfaceCOMData() {
   interface_1_currentCOM_.setZero();
   double total_mass = 0;
   for(auto it = interface_1_verts_.begin(); it != interface_1_verts_.end(); it++) {
      int v = *it;
      interface_1_currentCOM_ += mass_[v] * MapVec3(X + v * 3);
      total_mass += mass_[v];
    }
    interface_1_currentCOM_ /= total_mass;
    total_mass = 0;
    for(int v = domain_offset_toggle_;v < vertex_num_; v++) {
        domain_2_currentCOM_ += mass_[v] * MapVec3(X + v * 3);
        total_mass += mass_[v];
    }
    domain_2_currentCOM_ /= total_mass;
}

void SubspaceTet::ComputeInterfaceRotationMatrix() {
    Mat3 App;
    App.setZero();
    for(auto it = interface_1_verts_.begin(); it != interface_1_verts_.end(); it++) {
       int v = *it;
       Vec3 p1 = MapVec3(X + v * 3) - interface_1_currentCOM_;
       Vec3 q1 = MapVec3(rest_pos_ + v * 3) - interface_1_restCOM_;
       Mat3 outerPdk = p1 * (q1.transpose());
       outerPdk *= mass_[v];
       App += outerPdk;
     }
     Mat3 F = App.transpose();
     double R[9]; // rotation (row-major)
     double S[9]; // symmetric (row-major)
     double tolerance = 1E-6;
     int forceRotation = 1;
     PolarDecomposition::Compute(F.data(), R, S, tolerance, forceRotation);
     Apq.setZero();
     Apq(0,0) = R[0];  Apq(0,1) = R[1];  Apq(0,2) = R[2];
     Apq(1,0) = R[3];  Apq(1,1) = R[4];  Apq(1,2) = R[5];
     Apq(2,0) = R[6];  Apq(2,1) = R[7];  Apq(2,2) = R[8];
     P(Apq);

}

void SubspaceTet::UpdateRigidMotionAndLocalDeformation(double dt) {
  translation_vel_ += translation_acc_ * dt;
  angular_vel_ += angular_acc_ * dt;
  center_of_mass_ += translation_vel_ * dt;
  quaternion_ = quaternion_ + (0.5 * dt) * Quaternion<Real>(0, angular_vel_[0], angular_vel_[1], angular_vel_[2]) * quaternion_;
  quaternion_.Normalize();
  translation_vel_ *= 0.98;
  angular_vel_ *= 0.98;

  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  rotation_ = rotation_matrix;

  // compute global position and update vertex offset
  OMP_FOR
  for (int v = 0; v < vertex_num_; v++) {
    MapVec3 local_u(&(inv_fem_->u_[v][0]));
    local_u = vert_basis_[v] * q_;
    MapVec3 map_x(X + v * 3);
    map_x = rotation_ * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_;
  }
  //  inv_fem_->Update Offset();
}

void SubspaceTet::UpdatePosition() {
    OMP_FOR
    for (int v = 0; v < vertex_num_; v++) {
      Vec3 local_u = vert_basis_[v] * q_;
      MapVec3 map_x(X + v * 3);
      map_x = rotation_ * (local_u + vert_offset_from_mass_center_[v]) + center_of_mass_;
    }
}

double SubspaceTet::ComputeCubatureError(SubspaceTet::Vec &subspace_force, SubspaceTet::Mat &subspace_k) {
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(0);
  MapVec full_force((double*) inv_fem_->internal_force_, vertex_num_ * 3);
  Vec reduced_force = (-1) * (basis_transpose_ * full_force);
  Vec diff = reduced_force - subspace_force;

  double error = 0;
  for (int i = 0; i < diff.size(); ++i) {
    P(i, diff[i], reduced_force[i], subspace_force[i]);
    if (dj::Abs(reduced_force[i]) < 1e-10) {
      diff[i] = 0;
    } else {
      diff[i] = diff[i] / reduced_force[i];
      error += diff[i];
    }
  }

  Mat k_dot_u = Mat::Zero(vertex_num_ * 3, basis_num_);
  for (int b = 0; b < basis_num_; ++b) {
    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
  }
  Mat reduced_k = basis_transpose_ * k_dot_u;
  Mat diff_k = reduced_k - subspace_k;

  P(diff_k.maxCoeff());
  P(diff_k.minCoeff());
  return error / diff.size();
}

void SubspaceTet::SimulateWithRigidMotion(double dt) {
  // rigid motion
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  Vec3 net_force(0, 0, 0);
  Vec3 net_torque(0, 0, 0);
  net_force += total_mass_ * gravity_;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // external force in world space to each vertex
  std::vector<std::pair<int, Vec3> > collision_force;
  collision_force.reserve(vertex_num_ / 100);
  const double kFloor = 0.00;
  const double kFloorStiffness = 3000;
  for (int v = 0; v < vertex_num_; ++v) {
    // Ground collision
    //    if (0)
    if (X[v * 3 + 1] < kFloor) {
      Vec3 force(0, 0, 0);
      force[1] = (kFloor - X[v * 3 + 1]) * mass_[v] * kFloorStiffness;
      collision_force.emplace_back(v, force);
      net_force[1] += force[1];

      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_;
      net_torque += r.cross(force);
    }
  }
  // user interaction force
  {
#if 1
    Vec3 force;
    int v = GetUIForce(&force[0]);
    if (v >= 0) {
      collision_force.emplace_back(v, force);
      net_force += force;
      MapVec3 map_x(X + v * 3);
      Vec3 r = map_x - center_of_mass_;
      net_torque += r.cross(force);
    }
#endif
  }

  translation_acc_ = net_force / total_mass_;
  angular_acc_ = rotation_ * (inverse_intertial_tensor_ * (rotation_.transpose() * net_torque));

  // local deformation
  Vec subspace_force = Vec::Zero(basis_num_);
  Mat reduced_k = Mat::Zero(basis_num_, basis_num_);
#if 1
  profiler.Start("cubature");
  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(cubature_tets_);
  for (CubaturePoint cub : cubature_) {
    int t = cub.first;
    double weight = cub.second;
    double* force = inv_fem_->element_force_ + t * 12;
    double (*k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
    for (int local_v = 0; local_v < 4; ++local_v) {
      int v = tet_[t * 4 + local_v];
      MapVec3 map_force(force + local_v * 3);
      subspace_force -= vert_basis_transpose_[v] * (weight * map_force);
    }

    for (int i = 0; i < 4; ++i) {
      int v_i = tet_[t * 4 + i];
      for (int j = 0; j < 4; ++j) {
        int v_j = tet_[t * 4 + j];
        Mat sub_k(3, 3);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            sub_k(r, c) = k[i * 3 + r][j * 3 + c];
          }
        }
        sub_k *= weight;
        reduced_k += vert_basis_transpose_[v_i] * sub_k * vert_basis_[v_j];
      }
    }
  }
  profiler.End("cubature");
  //  P(ComputeCubatureError(subspace_force, reduced_k));
#else // compute exact reduced force and reduced K
  profiler.StartTimer("full reduced");
  MapVec map_int_force((double*) inv_fem_->internal_force_, vertex_num_ * 3);
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(0);
  Mat k_dot_u = Mat::Zero(vertex_num_ * 3, basis_num_);
  for (int b = 0; b < basis_num_; ++b) {
    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
  }
  reduced_k = basis_transpose_ * k_dot_u;
  subspace_force = -1 * (basis_transpose_ * map_int_force);
  profiler.EndTimer("full reduced");
#endif
  //  exit(0);
  internal_force_scaling_factor_ = 1.00;
  subspace_force *= internal_force_scaling_factor_;
  Mat3 rotation_transpose = rotation_.transpose();
  //  Vec collision_sub_f = Vec::Zero(basis_num_);
  for (std::pair<int, Vec3>& collision : collision_force) {
    int& v = collision.first;
    Vec3& force = collision.second;
    subspace_force += vert_basis_transpose_[v] * (rotation_transpose * force);
    //    collision_sub_f += vert_basis_transpose_[v] * (rotation_transpose * force);
  }
  //  P(collision_sub_f.norm());
  // fititious force
  if (1) {
    Vec3 acc = rotation_transpose * (translation_acc_ - gravity_);
    Vec3 angular_acc = rotation_transpose * angular_acc_;
    subspace_force += inertial_force_sandwich_ * (acc); // gravity also added here
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Coriolis force: F = -2m w\times v
    Mat3 angular_vel_mat = GetSkewSymmetrixMatrix(rotation_transpose * angular_vel_);
    double* data = angular_vel_mat.data();
    Mat coriolis_subspace_force = Mat::Zero(basis_num_, basis_num_);
    for (int i = 0; i < 9; ++i) {
      coriolis_subspace_force += (coriolis_force_sandwich_[i] * data[i]);
    }
    subspace_force += 2 * (coriolis_subspace_force) * vel_q_;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Euler force: F_v = -m_v * dw/dt \times r_v
    Mat3 angular_acc_mat = GetSkewSymmetrixMatrix(angular_acc);
    data = angular_acc_mat.data();
    Mat euler_force = Mat::Zero(basis_num_, basis_num_);
    for (int i = 0; i < 9; ++i) {
      euler_force += coriolis_force_sandwich_[i] * data[i];
    }
    subspace_force += euler_force * q_ + euler_force_sandwich_ * angular_acc;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // centerifugal force: F = -m * w \times (w \times r)
    angular_vel_mat = angular_vel_mat * angular_vel_mat;
    data = angular_vel_mat.data();
    Mat centrifugal_force0 = Mat::Zero(basis_num_, basis_num_);
    Vec centrifugal_force1 = Vec::Zero(basis_num_);
    for (int i = 0; i < 9; ++i) {
      centrifugal_force0 += coriolis_force_sandwich_[i] * data[i];
      centrifugal_force1 += centrifugal_force_sandwich_[i] * data[i];
    }
    subspace_force += centrifugal_force0 * q_ + centrifugal_force1;
  }

  Vec rhs = dt * subspace_force + vel_q_;
  reduced_k *= (dt * dt * internal_force_scaling_factor_);
  for (int i = 0; i < basis_num_; ++i) {
    reduced_k(i, i) += 1;
  }
  //  Mat A = Mat::Identity(basis_num_, basis_num_) + (dt * dt) * reduced_k;
  //  double rhs_dot = rhs.dot(rhs);
  //  ASSERT(rhs_dot == rhs_dot && std::isfinite(rhs_dot));
  vel_q_ = reduced_k.colPivHouseholderQr().solve(rhs);
  //  double dot = vel_q_.dot(vel_q_);
  //  P(vel_q_.norm());
  //  P(vel_q_.dot(vel_q_), rhs_dot);
  //  PMAT(inverse_intertial_tensor_);
  //  PMAT(inertial_tensor_ );
  //  PVEC(translation_vel_);
  //  PVEC(angular_vel_);
  //  PVEC(translation_acc_);
  //  PVEC(angular_acc_);
  //  ASSERT(dot == dot);
  //  vel_q_.setZero();
  q_ += vel_q_ * dt;
  if (0) {
    Mat3 lagragian_matrix = momentum_matrix_ * momentum_matrix_transpose_;
    Vec3 cm_offset = momentum_matrix_ * q_;
    Vec3 lambda = lagragian_matrix.ldlt().solve(cm_offset);
    Vec offset = momentum_matrix_transpose_ * lambda;
    q_ -= offset;
    vel_q_ -= offset / dt;
    Vec3 m = momentum_matrix_ * q_;
    cm_offset /= total_mass_;
    P(m.norm(), offset.norm(), (basis_ * offset).norm(), cm_offset);
    translation_vel_ += cm_offset / dt;
  }
  //  vel_q_ *= 0.98;

  translation_vel_ += translation_acc_ * dt;
  angular_vel_ += angular_acc_ * dt;
  center_of_mass_ += translation_vel_ * dt;
  quaternion_ = quaternion_ + (0.5 * dt) * Quaternion<double>(0, angular_vel_[0], angular_vel_[1], angular_vel_[2]) * quaternion_;
  quaternion_.Normalize();
  //  translation_vel_ *= 0.98;
  //  angular_vel_ *= 0.98;

  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation_matrix;
  quaternion_.Quaternion2Matrix(rotation_matrix.data());
  rotation_ = rotation_matrix;

  MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
  u = basis_ * q_;
  for (int v = 0; v < vertex_num_; ++v) {
    MapVec3(X + v * 3) = rotation_ * (MapVec3(&u[v * 3]) + vert_offset_from_mass_center_[v]) + center_of_mass_;
  }
}

void SubspaceTet::BuildParallelAccelerationStructures() {
   vv_list_1_.clear();
   vv_cubature_1_.clear();
   v_list_1_.clear();
   v_cubature_1_.clear();

   vv_list_2_.clear();
   vv_cubature_2_.clear();
   v_list_2_.clear();
   v_cubature_2_.clear();

   std::map<std::pair<int, int>, int> vv_map;
   std::map<int, int> v_map;

   int sizeCub = cubature_.size();
   for (int i=0;i<cubature_flip_counter_;i++) {
       int t = cubature_[i].first;
       double weight = cubature_[i].second;
       for (int i = 0; i < 4; ++i) {
           int vi = tet_[t * 4 + i];
           //int pi = vert_part_id_[vi];

           if (v_map.count(vi) > 0) {
               int idx = v_map[vi];
               ASSERT(idx < int(v_list_1_.size()));
               ASSERT(idx < int(v_cubature_1_.size()));
               ASSERT(v_list_1_[idx] == vi);
               //          P(pi, v_cubature_[pi].size(), idx, v_list_[pi].size());
               v_cubature_1_[idx].emplace_back(t, weight, i);
           } else {
               v_map[vi] = int(v_list_1_.size());
               v_list_1_.push_back(vi);
               v_cubature_1_.push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i)));
           }

           for (int j = 0; j < 4; ++j) {
               int vj = tet_[t * 4 + j];
               if (j == i) continue;
               if (vi >= vj) continue;
               auto v_pair = std::make_pair(vi, vj);
               if (vv_map.count(v_pair) > 0) {
                   int idx = vv_map[v_pair];
                   //ASSERT(vv_list_[vv_idx][idx] == v_pair);
                   vv_cubature_1_[idx].emplace_back(t, weight, i, j);
               } else {
                   vv_map[v_pair] = int(vv_list_1_.size());
                   vv_list_1_.emplace_back(vi, vj);
                   vv_cubature_1_.push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i, j)));
               }
           }  // each col
       } // each row
   } // all cubature

   vv_map.clear();
   v_map.clear();

   for (int i=cubature_flip_counter_;i<sizeCub;i++) {
       int t = cubature_[i].first;
       double weight = cubature_[i].second;
       for (int i = 0; i < 4; ++i) {
           int vi = tet_[t * 4 + i];
           //int pi = vert_part_id_[vi];

           if (v_map.count(vi) > 0) {
               int idx = v_map[vi];
               ASSERT(idx < int(v_list_2_.size()));
               ASSERT(idx < int(v_cubature_2_.size()));
               ASSERT(v_list_2_[idx] == vi);
               //          P(pi, v_cubature_[pi].size(), idx, v_list_[pi].size());
               v_cubature_2_[idx].emplace_back(t, weight, i);
           } else {
               v_map[vi] = int(v_list_2_.size());
               v_list_2_.push_back(vi);
               v_cubature_2_.push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i)));
           }

           for (int j = 0; j < 4; ++j) {
               int vj = tet_[t * 4 + j];
               if (j == i) continue;
               if (vi >= vj) continue;
               auto v_pair = std::make_pair(vi, vj);
               if (vv_map.count(v_pair) > 0) {
                   int idx = vv_map[v_pair];
                   //ASSERT(vv_list_[vv_idx][idx] == v_pair);
                   vv_cubature_2_[idx].emplace_back(t, weight, i, j);
               } else {
                   vv_map[v_pair] = int(vv_list_2_.size());
                   vv_list_2_.emplace_back(vi, vj);
                   vv_cubature_2_.push_back(std::vector<VVCubature>(1, VVCubature(t, weight, i, j)));
               }
           }  // each col
       } // each row
   } // all cubature

   P(vv_list_1_.size());
   P(vv_list_2_.size());
   P(v_list_1_.size());
   P(v_list_2_.size());
}

void SubspaceTet::SimulateWithCubature(double dt) {
  //  EXECUTE_TIMES(20);

  Vec subspace_force = Vec::Zero(basis_num_);
  Mat reduced_k = Mat::Zero(basis_num_, basis_num_);
  Vec subspace_force_1 = Vec::Zero(basis_num_);
  Mat reduced_k_1 = Mat::Zero(basis_num_, basis_num_);
  Vec subspace_force_2 = Vec::Zero(basis_num_);
  Mat reduced_k_2 = Mat::Zero(basis_num_, basis_num_);

  std::vector<Mat> result_pairs_1;
  result_pairs_1.resize(v_list_1_.size() +vv_list_1_.size() );

  std::vector<Mat> result_pairs_2;
  result_pairs_2.resize(v_list_2_.size() +vv_list_2_.size() );

  profiler.Start("cubature");
  inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(cubature_tets_);

  /*int kounter = 0;
 for (CubaturePoint cub : cubature_)
{
    int t = cub.first;
    double weight = cub.second;
    double* force = inv_fem_->element_force_ + t * 12;
    double (*k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
    for (int local_v = 0; local_v < 4; ++local_v) {
      int v = tet_[t * 4 + local_v];
      MapVec3 map_force(force + local_v * 3);
      if(kounter<cubature_flip_counter_) {
          subspace_force_1 -= vert_basis_transpose_[v] * (weight * map_force);
      }
      else {
          subspace_force_2 -= vert_basis_transpose_[v] * (weight * map_force);
      }
    }
    for (int i = 0; i < 4; ++i) {
      int v_i = tet_[t * 4 + i];
      for (int j = 0; j < 4; ++j) {
        int v_j = tet_[t * 4 + j];
        Mat sub_k(3, 3);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            sub_k(r, c) = k[i * 3 + r][j * 3 + c];
          }
        }
        sub_k *= weight;
        if(kounter<cubature_flip_counter_) {
            reduced_k_1 += vert_basis_transpose_[v_i] * sub_k * vert_basis_[v_j];
        }
        else {
            reduced_k_2 += vert_basis_transpose_[v_i] * sub_k * vert_basis_[v_j];
        }
      }
    }
    kounter++;
  }

  if(sim_mode_==2 && global::stiffyMult!=1) {
      reduced_k = global::stiffyMult * reduced_k_2 + reduced_k_1;
      subspace_force = global::stiffyMult*subspace_force_2 + subspace_force_1 ;
  }
  else {
      reduced_k = reduced_k_2 + reduced_k_1;
      subspace_force = subspace_force_2 + subspace_force_1 ;
  }*/


int offset_1 = v_list_1_.size();
int offset_2 = v_list_2_.size();

profiler.Start("OMP Assemble");
  OMP_FOR
  for (int j = 0; j < int(v_list_1_.size()); ++j) {
      int v = v_list_1_[j];
      Mat3 sub_k = Mat3::Zero();
      // all cubature tet incident on this v
      for (int k = 0; k < int(v_cubature_1_[j].size()); ++k) {
          VVCubature& cubature = v_cubature_1_[j][k];
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
      result_pairs_1[j] = vert_basis_transpose_[v] * (sub_k) * vert_basis_[v];
  }

  OMP_FOR
  for (int i = 0; i < int(vv_list_1_.size()); ++i) {
    int vi = vv_list_1_[i].first;
    int vj = vv_list_1_[i].second;
    // list of  cubature tets that contains this vert-vert pair
    Mat3 sub_k_1 = Mat3::Zero();
    Mat3 sub_k_2 = Mat3::Zero();
    //profiler.Start("Inner Loop");
    for (int k = 0; k < int(vv_cubature_1_[i].size()); ++k) {
      VVCubature& cubature = vv_cubature_1_[i][k];
      const int& t = cubature.cubature_tet;
      const double& weight = cubature.cubature_weight;
      const int& idx0 = cubature.i;
      const int& idx1 = cubature.j;
      double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          sub_k_1(r, c) += weight * tet_k[idx0 * 3 + r][idx1 * 3 + c];
          sub_k_2(r, c) += weight * tet_k[idx1 * 3 + r][idx0 * 3 + c];
        }
      }
    }
    result_pairs_1[offset_1+i] = vert_basis_transpose_[vi]* sub_k_1 * vert_basis_[vj];
    result_pairs_1[offset_1+i] += vert_basis_transpose_[vj]* sub_k_2 * vert_basis_[vi];
  }

#pragma omp parallel
{
    Mat result_private(basis_num_,basis_num_);
    result_private.setZero();
    #pragma omp for nowait //fill result_private in parallel
    for(int i=0; i<result_pairs_1.size(); i++) result_private += result_pairs_1[i];
    #pragma omp critical
    reduced_k_1 += result_private;
}


//Get help in parallelizatio of this one
/*for(int i = 0; i <result_pairs_1.size();i++) {
      reduced_k_1 +=  result_pairs_1[i];
  }*/

  for (int i=0;i<cubature_flip_counter_;i++) {
    int& t = cubature_[i].first;
    double& weight = cubature_[i].second;
    // reduced internal force
    double* force = inv_fem_->element_force_ + t * 12;
    for (int local_v = 0; local_v < 4; ++local_v) {
      int v = tet_[t * 4 + local_v];
      //int p = vert_part_id_[v];
      MapVec3 map_force(force + local_v * 3);
      subspace_force_1 -= vert_basis_transpose_[v] * (weight * map_force);
    }
  }



  OMP_FOR
  //profiler.Start("Loop time");
  for (int j = 0; j < int(v_list_2_.size()); ++j) {
      int v = v_list_2_[j];
      Mat3 sub_k = Mat3::Zero();
      // all cubature tet incident on this v
      //P(j,v,v_cubature_2_[j].size());
      for (int k = 0; k < int(v_cubature_2_[j].size()); ++k) {
          VVCubature& cubature = v_cubature_2_[j][k];
          const int& t = cubature.cubature_tet;
          const double& weight = cubature.cubature_weight;
          const int& idx = cubature.i;
          double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);

          for (int r = 0; r < 3; ++r) {
              for (int c = 0; c < 3; ++c) {
                  sub_k(r, c) += weight * tet_k[idx * 3 + r][idx * 3 + c];
              }
          }
          //P(k,t,weight,idx);
      }
      //profiler.Start("Mul time");
      result_pairs_2[j] = vert_basis_transpose_[v] * (sub_k) * vert_basis_[v];
      //profiler.End("Mul time");
  }
  //profiler.End("Loop time");

  //L("simulate with cubature");

  OMP_FOR
  for (int i = 0; i < int(vv_list_2_.size()); ++i) {
    int vi = vv_list_2_[i].first;
    int vj = vv_list_2_[i].second;
    // list of  cubature tets that contains this vert-vert pair
    Mat3 sub_k_1 = Mat3::Zero();
    Mat3 sub_k_2 = Mat3::Zero();
    for (int k = 0; k < int(vv_cubature_2_[i].size()); ++k) {
      VVCubature& cubature = vv_cubature_2_[i][k];
      const int& t = cubature.cubature_tet;
      const double& weight = cubature.cubature_weight;
      const int& idx0 = cubature.i;
      const int& idx1 = cubature.j;
      double (*tet_k)[12] = (double (*)[12]) (inv_fem_->element_k_ + t * 144);
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          sub_k_1(r, c) += weight * tet_k[idx0 * 3 + r][idx1 * 3 + c];
          sub_k_2(r, c) += weight * tet_k[idx1 * 3 + r][idx0 * 3 + c];
        }
      }
    }
    result_pairs_2[offset_2+i] = vert_basis_transpose_[vi]* sub_k_1 * vert_basis_[vj];
    result_pairs_2[offset_2+i] += vert_basis_transpose_[vj]* sub_k_2 * vert_basis_[vi];
  }

#pragma omp parallel
{
    Mat result_private(basis_num_,basis_num_);
    result_private.setZero();
    #pragma omp for nowait //fill result_private in parallel
    for(int i=0; i<result_pairs_2.size(); i++) result_private += result_pairs_2[i];
    #pragma omp critical
    reduced_k_2 += result_private;
}

//Get help in parallelizatio of this one
/*for(int i = 0; i <result_pairs_2.size();i++) {
      reduced_k_2 +=  result_pairs_2[i];
  }*/

  for (int i=cubature_flip_counter_;i<cubature_.size();i++) {
    int& t = cubature_[i].first;
    double& weight = cubature_[i].second;
    // reduced internal force
    double* force = inv_fem_->element_force_ + t * 12;
    for (int local_v = 0; local_v < 4; ++local_v) {
      int v = tet_[t * 4 + local_v];
      //int p = vert_part_id_[v];
      MapVec3 map_force(force + local_v * 3);
      subspace_force_2 -= vert_basis_transpose_[v] * (weight * map_force);
    }
  }

  profiler.End("OMP Assemble");
  /*std::ofstream fileeout("TEST_K_RED.txt");
  for(int i=0;i<basis_num_;i++) {
      for(int j=0;j<basis_num_;j++) {
          fileeout << reduced_k(i,j) << " ";
      }
      fileeout << "\n";
  }
  fileeout.close();*/

  profiler.Start("PreSolve");

  if(sim_mode_==2 && global::stiffyMult!=1) {
      reduced_k =  reduced_k_2 + global::stiffyMult *reduced_k_1 ;
      subspace_force = subspace_force_2 +  global::stiffyMult *subspace_force_1 ;
  }
  else {
      reduced_k = reduced_k_2 + reduced_k_1;
      subspace_force = subspace_force_2 + subspace_force_1 ;
  }


  if(sim_mode_ == 2) {
      double stiffness = 5.0e7;
      MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
      Vec my_f_reduced = Vec::Zero(vertex_num_ * 3);
      Mat my_sp_reduced = Mat::Zero(basis_num_, basis_num_);
      for(int i=0;i<basis_num_;i++) {
          for(int j=0;j<basis_num_;j++) {
              for(auto it = interface_vertices_.begin();it!=interface_vertices_.end();it++) {
                  int a = it->first;
                  int b = it->second;
                  int i1 = 3*a+0;
                  int i2 = 3*a+1;
                  int i3 = 3*a+2;
                  int j1 = 3*b+0;
                  int j2 = 3*b+1;
                  int j3 = 3*b+2;
                  double e1 = (basis_(i1,i) - basis_(j1,i))*stiffness;//*global::stiffyMult;
                  double e2 = (basis_(i2,i) - basis_(j2,i))*stiffness;//*global::stiffyMult;
                  double e3 = (basis_(i3,i) - basis_(j3,i))*stiffness;//*global::stiffyMult;
                  double e4 = (basis_(j1,i) - basis_(i1,i))*stiffness;//*global::stiffyMult;
                  double e5 = (basis_(j2,i) - basis_(i2,i))*stiffness;//*global::stiffyMult;
                  double e6 = (basis_(j3,i) - basis_(i3,i))*stiffness;//*global::stiffyMult;
                  my_sp_reduced(i,j) += ((e1*basis_(i1,j)) + (e2*basis_(i2,j)) + (e3*basis_(i3,j)) + (e4*basis_(j1,j)) + (e5*basis_(j2,j)) + (e6*basis_(j3,j)));
              }
          }
      }

      for(int i =0;i < basis_num_;i++) {
          for(auto it = interface_vertices_.begin();it!=interface_vertices_.end();it++) {
              int a = it->first;
              int b = it->second;
              int i1 = 3*a+0;
              int i2 = 3*a+1;
              int i3 = 3*a+2;
              int j1 = 3*b+0;
              int j2 = 3*b+1;
              int j3 = 3*b+2;
              double u1 = (u(i1) - u(j1))*stiffness;//*global::stiffyMult;
              double u2 = (u(i2) - u(j2))*stiffness;//*global::stiffyMult;
              double u3 = (u(i3) - u(j3))*stiffness;//*global::stiffyMult;
              //double u4 = (u(j1) - u(i1))*stiffness;
              //double u5 = (u(j2) - u(i2))*stiffness;
              //double u6 = (u(j3) - u(i3))*stiffness;
              my_f_reduced(i) += ((u1*basis_(i1,i)) + (u2*basis_(i2,i)) + (u3*basis_(i3,i)) - (u1*basis_(j1,i)) - (u2*basis_(j2,i)) - (u3*basis_(j3,i)));
          }
      }
      reduced_k += my_sp_reduced;

      subspace_force -= my_f_reduced;
  }
   profiler.End("PreSolve");

  profiler.Start("Solves");
  subspace_force += subspace_gravity_;

  Vec3 ui_force;
  int ui_force_vert = GetUIForce(&ui_force[0]);
  if (ui_force_vert >= 0) {
    subspace_force += vert_basis_transpose_[ui_force_vert] * ui_force *1.0 ;
  }


  internal_force_scaling_factor_ = 0.1;//0.1;
  Vec rhs = dt * internal_force_scaling_factor_ * subspace_force + vel_q_;
  Mat A = Mat::Identity(basis_num_, basis_num_) + (dt * dt) * reduced_k * internal_force_scaling_factor_;

  vel_q_ = A.colPivHouseholderQr().solve(rhs);
  //  P(vel_q_.dot(vel_q_));
  q_ += vel_q_ * dt;
  //  vel_q_ *= 0.98;

  MapVec map_X(X, vertex_num_ * 3);
  MapVec map_rest_pos_(rest_pos_, vertex_num_ * 3);
  MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
  u = basis_ * q_;
  map_X = u + map_rest_pos_;
  MapVec v((double*) inv_fem_->vel_, vertex_num_ * 3);
  v = basis_ * vel_q_;

  memcpy(render_X,X,sizeof(double)*vertex_num_ * 3);

  if(global::sim_state==3) {
      for(auto it = interface_vertices_.begin(); it!= interface_vertices_.end();it++) {
          int i = it->first;
          int j = it->second;
          double xa = (render_X[3*i+0]+render_X[3*j+0])*0.5;
          render_X[3*i+0] = render_X[3*j+0] = xa;
          double ya = (render_X[3*i+1]+render_X[3*j+1])*0.5;
          render_X[3*i+1] = render_X[3*j+1] = ya;
          double za = (render_X[3*i+2]+render_X[3*j+2])*0.5;
          render_X[3*i+2] = render_X[3*j+2] = za;
      }
  }
    profiler.End("Solves");
    profiler.End("cubature");
}

void SubspaceTet::LoadCubature(const char *file_name) {
  if (!conf.Get<int>("use all cubature points")) {
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    int cubature_point_num;
    in >> cubature_point_num;
    //  P(cubature_point_num);
    cubature_.reserve(cubature_point_num);
    cubature_tets_.reserve(cubature_point_num);
    for (int i = 0; i < cubature_point_num; ++i) {
      int tet = -1;
      double weight = -1e20;
      in >> tet;
      in >> weight;
      if (weight < 1e-10) {
        continue;
        KK;
      } else {
        cubature_.push_back(std::make_pair(tet, weight));
        cubature_tets_.push_back(tet);
      }
    }
    in.close();
  } else {
    //  P(cubature_.size());
    (void) file_name;
    L("use all tets as cubature");
    int cubature_point_num = tet_number;
    cubature_.resize(cubature_point_num);
    cubature_tets_.resize(cubature_point_num);
    for (int i = 0; i < cubature_point_num; ++i) {
      cubature_[i].first = i;
      cubature_[i].second = 1;
      cubature_tets_[i] = cubature_[i].first;
    }
  }

  BuildParallelAccelerationStructures();
}

void SubspaceTet::LoadCubature(std::vector<int> &p, std::vector<double> &w){
    int numCP = p.size();
    for(int i =0; i < numCP; i++) {
        if (w[i] < 1e-10) {
          continue;
        } else {
          cubature_.push_back(std::make_pair(p[i], w[i]));
          cubature_tets_.push_back(p[i]);
        }
    }

    BuildParallelAccelerationStructures();

}

void SubspaceTet::UpdateLocalDeformation(double *new_vel_q, double dt) {
  memcpy(vel_q_.data(), new_vel_q, sizeof(double) * basis_num_);
  q_ += vel_q_ * dt;
  vel_q_ *= 0.95;
  MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
  u = basis_ * q_;
  //  for (int v = 0; v < vertex_num_; ++v) {
  //    MapVec3(X + v * 3) = rotation_ * (MapVec3(&u[v * 3]) + vert_offset_from_mass_center_[v]) + center_of_mass_;
  //  }
}

void SubspaceTet::GetReducedTangentStiffnessMatrix(SubspaceTet::Mat &reduced_k) {
  Mat k_dot_u = Mat::Zero(vertex_num_ * 3, basis_num_);
  for (int b = 0; b < basis_num_; ++b) {
    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
  }
  reduced_k = basis_transpose_ * k_dot_u;
}

void SubspaceTet::AddFictitiousForceAndGravity(double *force) {
  MapVec subspace_force(force, basis_num_);
  Mat3 rotation_transpose = rotation_.transpose();
  Vec3 acc = rotation_transpose * (translation_acc_ - gravity_); // gravity also added here
  Vec3 angular_acc = rotation_transpose * angular_acc_;
  subspace_force += inertial_force_sandwich_ * (acc);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Coriolis force: F = -2m w\times v
  Mat3 angular_vel_mat = GetSkewSymmetrixMatrix(rotation_transpose * angular_vel_);
  double* data = angular_vel_mat.data();
  Mat coriolis_subspace_force = Mat::Zero(basis_num_, basis_num_);
  for (int i = 0; i < 9; ++i) {
    coriolis_subspace_force += (coriolis_force_sandwich_[i] * data[i]);
  }
  subspace_force += 2 * (coriolis_subspace_force) * vel_q_;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Euler force: F_v = -m_v * dw/dt \times r_v
  Mat3 angular_acc_mat = GetSkewSymmetrixMatrix(angular_acc);
  data = angular_acc_mat.data();
  Mat euler_force = Mat::Zero(basis_num_, basis_num_);
  for (int i = 0; i < 9; ++i) {
    euler_force += coriolis_force_sandwich_[i] * data[i];
  }
  subspace_force += euler_force * q_ + euler_force_sandwich_ * angular_acc;
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // centerifugal force: F = -m * w \times (w \times r)
  angular_vel_mat = angular_vel_mat * angular_vel_mat;
  data = angular_vel_mat.data();
  Mat centrifugal_force0 = Mat::Zero(basis_num_, basis_num_);
  Vec centrifugal_force1 = Vec::Zero(basis_num_);
  for (int i = 0; i < 9; ++i) {
    centrifugal_force0 += coriolis_force_sandwich_[i] * data[i];
    centrifugal_force1 += centrifugal_force_sandwich_[i] * data[i];
  }
  subspace_force += centrifugal_force0 * q_ + centrifugal_force1;
}

void SubspaceTet::LoadSubspace(const char *basis_file) {
  std::ifstream in(basis_file);
  ASSERT(in.is_open(), P(basis_file));
  int v_num ;
  in >> v_num >> basis_num_;
  ASSERT(v_num == vertex_num_ * 3, P(vertex_num_, v_num));
  basis_ = Mat::Zero(v_num, basis_num_);
  for (int i = 0; i < v_num; ++i) {
    for (int j = 0; j < basis_num_; ++j) {
      in >> basis_(i, j);
    }
  }
  in.close();
  OnSubspaceLoaded();
}

void SubspaceTet::LoadBinarySubspace(const char *basis_file) {
  BinaryFileReader in(basis_file);
  int row_num;
  //  std::vector<double> basis;
  //    ReadBasisInBinary(basis_file, row_num, basis_num_, basis);
  //    ReadBasisInText(file_name, v_num, basis_num, basis);
  in.Read(&row_num, 1);
  ASSERT(row_num == vertex_num_ * 3);
  in.Read(&basis_num_, 1);
  MatCol tmp_basis(row_num, basis_num_);
  in.Read(tmp_basis.data(), row_num * basis_num_);
  basis_ = tmp_basis;
  OnSubspaceLoaded();
}

void SubspaceTet::LoadBinarySubspace(std::vector<double> &vec, int numVerts, int numBasis) {
    MatCol tmp_basis(numVerts,numBasis);
    memcpy(tmp_basis.data(), &vec[0], sizeof(double)*numVerts*numBasis);
    basis_ = tmp_basis;
    if(sim_mode_==1)
        OnSubspaceLoaded();
    else
        OnSubspaceLoaded2();
}

void SubspaceTet::saveOldBasis() {
    old_basis_ = basis_;
    old_basis_transpose_ = basis_transpose_;
    old_vert_basis_ = vert_basis_;
    old_vert_basis_transpose_ = vert_basis_transpose_;
    P("saved");
}

void SubspaceTet::OnSubspaceLoaded() {
  basis_num_ = basis_.cols();
  basis_transpose_ = basis_.transpose();
  q_ = Vec::Zero(basis_num_);
  vel_q_ = Vec::Zero(basis_num_);
  P(basis_num_);
  //VerifySubspaceDiagonality();
  vert_basis_.resize(vertex_num_);
  vert_basis_transpose_.resize(vertex_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    vert_basis_[v] = basis_.middleRows(v * 3, 3);
    vert_basis_transpose_[v] = vert_basis_[v].transpose();
  }

  Vec total_gravity = Vec::Zero(vertex_num_ * 3);
  for (int v = 0; v < vertex_num_; ++v) {
    total_gravity[v * 3 + 0] = global::gravity[0] * mass_[v];
    total_gravity[v * 3 + 1] = global::gravity[1] * mass_[v];
    total_gravity[v * 3 + 2] = global::gravity[2] * mass_[v];
  }
  subspace_gravity_ = basis_transpose_ * total_gravity;
  ComputeMomentumAndTorqueMatrix();
  PrecomputeFastSanwichTransform();
  P(basis_num_);
}

void SubspaceTet::OnSubspaceLoaded2() {

  basis_num_ = basis_.cols();
  basis_transpose_ = basis_.transpose();
  inv_fem_->Reset();
  //MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
  //MapVec v((double*) inv_fem_->vel_, vertex_num_ * 3);
  q_.setZero();
  vel_q_.setZero();
  //q_ = basis_transpose_ * u;
  //vel_q_ = basis_transpose_ * v;
  VerifySubspaceDiagonality();
  vert_basis_.resize(vertex_num_);
  vert_basis_transpose_.resize(vertex_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    vert_basis_[v] = basis_.middleRows(v * 3, 3);
    vert_basis_transpose_[v] = vert_basis_[v].transpose();
  }
  subspace_gravity_.setZero();
}

void SubspaceTet::VerifySubspaceDiagonality() {
  Mat m_dot_u = basis_;
  for (int v = 0; v < vertex_num_; ++v) {
    m_dot_u.row(v * 3 + 0) *= mass_[v];
    m_dot_u.row(v * 3 + 1) *= mass_[v];
    m_dot_u.row(v * 3 + 2) *= mass_[v];
  }
  Mat zero = basis_transpose_ * m_dot_u;
  zero -= Mat::Identity(basis_num_, basis_num_);
  ASSERT(dj::Abs(zero.minCoeff()) < 1e-5, P(zero.minCoeff(), zero.maxCoeff()));
  ASSERT(dj::Abs(zero.maxCoeff()) < 1e-5, P(zero.minCoeff(), zero.maxCoeff()));
  P(zero.minCoeff(), zero.maxCoeff());
}

void SubspaceTet::setSimMode(int x) {
    sim_mode_ = x;
    inv_fem_->setSimMode(x);
}

void SubspaceTet::LoadMass(const char *file_name) {
  std::ifstream in(file_name);
  int v_num;
  ASSERT(in.is_open(), P(file_name));
  in >> v_num;
  ASSERT(v_num == vertex_num_);
  for (int i = 0; i < vertex_num_; ++i) {
    in >> mass_[i];
  }
  in.close();
  ComputeCenterOfMass();
  ComputeInertialTensor();
}

void SubspaceTet::ComputeInertialTensor() {
  inertial_tensor_.setZero();
  for (int v = 0; v < vertex_num_; ++v) {
    Mat3 tmp = Mat3::Zero();
    Vec3 r = MapVec3(rest_pos_ + v * 3) - center_of_mass_;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        tmp(i, j) = -r[i] * r[j];
      }
    }
    double dot = r.dot(r);
    tmp(0, 0) += dot;
    tmp(1, 1) += dot;
    tmp(2, 2) += dot;
    inertial_tensor_ += mass_[v] * tmp;
  }
  inverse_intertial_tensor_ = inertial_tensor_.inverse();
}

void SubspaceTet::ProjectPositionToSubspace() {
  memcpy(tmp_X, X, sizeof(double) * vertex_num_ * 3);
  for (int v = 0; v < vertex_num_; ++v) {
    tmp_X[v * 3 + 0] *= mass_[v];
    tmp_X[v * 3 + 1] *= mass_[v];
    tmp_X[v * 3 + 2] *= mass_[v];
  }
  MapVec map_X(tmp_X, vertex_num_ * 3);
  q_ = basis_transpose_ * map_X;
}

void SubspaceTet::AnimateBasis() {

  MapVec map_X(X, vertex_num_ * 3);
  MapVec map_rest_pos_(rest_pos_, vertex_num_ * 3);
  map_X = magnitude_ * basis_.col(current_basis_) + map_rest_pos_;
  static double step = 0.02;
  magnitude_ += step;
  if (magnitude_ > 1 || magnitude_ < -1) step *= -1;
}

void SubspaceTet::Simulate(double dt) {
  if (current_basis_ != -1) {
    AnimateBasis();
    return;
  }
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  Mat k_dot_u = Mat::Zero(vertex_num_ * 3, basis_num_);
  for (int b = 0; b < basis_num_; ++b) {
    if(sim_mode_==1)
       inv_fem_->tangent_stiffness_matrix_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
    else if(sim_mode_==2)
       inv_fem_->stitch_tangent_stiffness_matrix_1_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
    //else if(sim_mode_ ==4)
   //     inv_fem_->stitch_tangent_stiffness_matrix_2_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
  }
  Mat reduced_k = basis_transpose_ * k_dot_u;
  Mat A = reduced_k;
  A *= (dt * dt);
  A += Mat::Identity(basis_num_, basis_num_);
  //  L("final A");
  //  PrintMat(&A(0, 0), basis_num_, basis_num_);

  MapVec full_force((double*) inv_fem_->rhs_, vertex_num_ * 3);
  {
    if (0) {
      Vec reduced_force = basis_transpose_ * full_force;
      L("int and ext force")
          PrintVec(&reduced_force[0], basis_num_);
      reduced_force *= dt;
      L("int and ext force times dt")
          PrintVec(&reduced_force[0], basis_num_);
    }
    if (0) {
      Vec t = Vec::Zero(vertex_num_ * 3);
      t[selected_vertex_ * 3 + 2] += conf.Get<double>("force");
      Vec reduced_force = basis_transpose_ * t;
      L("ext force");
      PrintVec(&reduced_force[0], basis_num_);
      //    reduced_force *= dt;
    }
  }
  //  full_force[selected_vertex_ * 3 + 2] += pose_conf.Get<double>("force");
  //  full_force[selected_vertex_ * 3 + 1] -= pose_conf.Get<double>("force");
  double ui_force[3];
  int ui_force_vert = GetUIForce(ui_force);
  if (ui_force_vert >= 0) {
    P(dj::Vec3d(ui_force));
    full_force[ui_force_vert * 3 + 0] += 10.0*ui_force[0];
    full_force[ui_force_vert * 3 + 1] += 10.0*ui_force[1];
    full_force[ui_force_vert * 3 + 2] += 10.0*ui_force[2];
  }

 /* if(global::sim_state==3) {
  full_force[15666 * 3 + 0] += -700.0;
  full_force[15666 * 3 + 1] += 263;
  full_force[15666 * 3 + 2] += 45;
   }*/
Vec reduced_force = basis_transpose_ * full_force;
  if(sim_mode_ == 2) {
      double stiffness = 5.0e6;
      MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
      Vec my_f_reduced = Vec::Zero(vertex_num_ * 3);
      Mat my_sp_reduced = Mat::Zero(basis_num_, basis_num_);
      for(int i=0;i<basis_num_;i++) {
          for(int j=0;j<basis_num_;j++) {
              for(auto it = interface_vertices_.begin();it!=interface_vertices_.end();it++) {
                  int a = it->first;
                  int b = it->second;
                  int i1 = 3*a+0;
                  int i2 = 3*a+1;
                  int i3 = 3*a+2;
                  int j1 = 3*b+0;
                  int j2 = 3*b+1;
                  int j3 = 3*b+2;
                  double e1 = (basis_(i1,i) - basis_(j1,i))*stiffness;//*global::stiffyMult;
                  double e2 = (basis_(i2,i) - basis_(j2,i))*stiffness;//*global::stiffyMult;
                  double e3 = (basis_(i3,i) - basis_(j3,i))*stiffness;//*global::stiffyMult;
                  double e4 = (basis_(j1,i) - basis_(i1,i))*stiffness;//*global::stiffyMult;
                  double e5 = (basis_(j2,i) - basis_(i2,i))*stiffness;//*global::stiffyMult;
                  double e6 = (basis_(j3,i) - basis_(i3,i))*stiffness;//*global::stiffyMult;
                  my_sp_reduced(i,j) += ((e1*basis_(i1,j)) + (e2*basis_(i2,j)) + (e3*basis_(i3,j)) + (e4*basis_(j1,j)) + (e5*basis_(j2,j)) + (e6*basis_(j3,j)));
              }
          }
      }

      for(int i =0;i < basis_num_;i++) {
          for(auto it = interface_vertices_.begin();it!=interface_vertices_.end();it++) {
              int a = it->first;
              int b = it->second;
              int i1 = 3*a+0;
              int i2 = 3*a+1;
              int i3 = 3*a+2;
              int j1 = 3*b+0;
              int j2 = 3*b+1;
              int j3 = 3*b+2;
              double u1 = (u(i1) - u(j1))*stiffness;//*global::stiffyMult;
              double u2 = (u(i2) - u(j2))*stiffness;//*global::stiffyMult;
              double u3 = (u(i3) - u(j3))*stiffness;//*global::stiffyMult;
              //double u4 = (u(j1) - u(i1))*stiffness;
              //double u5 = (u(j2) - u(i2))*stiffness;
              //double u6 = (u(j3) - u(i3))*stiffness;
              my_f_reduced(i) += ((u1*basis_(i1,i)) + (u2*basis_(i2,i)) + (u3*basis_(i3,i)) - (u1*basis_(j1,i)) - (u2*basis_(j2,i)) - (u3*basis_(j3,i)));
          }
      }
      reduced_k += my_sp_reduced;

      reduced_force -= my_f_reduced;
  }

  //  P(selected_vertex_, pose_conf.Get<double>("force"));g

  //  L("int and ext force")
  //  PrintVec(&reduced_force[0], basis_num_);
  reduced_force *= dt;
  //  L("int and ext force times dt")
  //  PrintVec(&reduced_force[0], basis_num_);


  //  exit(0);
#if 1
  reduced_force += vel_q_;
  vel_q_ = A.colPivHouseholderQr().solve(reduced_force);
  P(vel_q_.dot(vel_q_));
#else
  if (0) {
    Vec t = dt * dt * (reduced_k * vel_q_);
    L("(D + h * K) qdot");
    PrintVec(&t[0], basis_num_);
  }
  reduced_force -= dt * dt * (reduced_k * vel_q_);
  //  L("final rhs");
  //  L("int and ext force times dt - dt * dt * K * vel_q");
  //  PrintVec(&reduced_force[0], basis_num_);
  Vec d_qvel = A.colPivHouseholderQr().solve(reduced_force);
  vel_q_ += d_qvel;
  P(vel_q_.dot(vel_q_));
#endif
  q_ += vel_q_ * dt;
  MapVec map_X(X, vertex_num_ * 3);
  MapVec map_rest_pos_(rest_pos_, vertex_num_ * 3);
  MapVec u((double*) inv_fem_->u_, vertex_num_ * 3);
  u = basis_ * q_;
  map_X = u + map_rest_pos_;
  MapVec v((double*) inv_fem_->vel_, vertex_num_ * 3);
  v = basis_ * vel_q_;
/*
  ComputeCurrentInterfaceCOMData();
  ComputeInterfaceRotationMatrix();

  for(int v=domain_offset_toggle_;v<vertex_num_;v++) {
      Vec3 vd =  MapVec3(X + v * 3);
      vd = (Apq * (vd-domain_2_currentCOM_)) + domain_2_currentCOM_;
  }*/

  memcpy(render_X,X,sizeof(double)*vertex_num_ * 3);

  if(global::sim_state==3) {
      for(auto it = interface_vertices_.begin(); it!= interface_vertices_.end();it++) {
          int i = it->first;
          int j = it->second;
          double xa = (render_X[3*i+0]+render_X[3*j+0])*0.5;
          render_X[3*i+0] = render_X[3*j+0] = xa;
          double ya = (render_X[3*i+1]+render_X[3*j+1])*0.5;
          render_X[3*i+1] = render_X[3*j+1] = ya;
          double za = (render_X[3*i+2]+render_X[3*j+2])*0.5;
          render_X[3*i+2] = render_X[3*j+2] = za;
      }
  }
}

void SubspaceTet::SimulateWithReduceMassMatrix(double dt) {
  inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(dt);
  Mat k_dot_u = Mat::Zero(vertex_num_ * 3, basis_num_);
  for (int b = 0; b < basis_num_; ++b) {
    inv_fem_->tangent_stiffness_matrix_->MultiplyVector(basis_.col(b).data(), k_dot_u.col(b).data());
  }
  Mat reduced_k = basis_transpose_ * k_dot_u;
  //  k_dot_u = tsm * basis_;
  //  Mat reduced_k = basis_transpose_ * k_dot_u;
  reduced_k *= dt * dt;
  KK;
  reduced_k += reduce_mass_matrix_;
  KK;
  MapVec full_force((double*) inv_fem_->rhs_, vertex_num_ * 3);
  full_force[selected_vertex_ * 3 + 2] += conf.Get<double>("force") * dt;
  Vec reduced_force = basis_transpose_ * full_force;
  reduced_force *= dt;
  reduced_force += reduce_mass_matrix_ * vel_q_;
  Vec rhs_tmp = reduced_force;
  vel_q_ = reduced_k.colPivHouseholderQr().solve(reduced_force);
  Vec residual = reduced_k * vel_q_ - rhs_tmp;
  P(residual.dot(residual));
  P(vel_q_.dot(vel_q_));
  q_ += vel_q_ * dt;
  vel_q_ *= 0.96;
  MapVec map_X(X, vertex_num_ * 3);
  MapVec map_rest_pos_(rest_pos_, vertex_num_ * 3);
  map_X = basis_ * q_ + map_rest_pos_;

}

void SubspaceTet::PrecomputeFastSanwichTransform() {
  inertial_force_sandwich_ = Mat::Zero(basis_num_, 3);
  euler_force_sandwich_ = Mat::Zero(basis_num_, 3);
  for (int idx = 0; idx < 9; ++idx) {
    coriolis_force_sandwich_[idx] = Mat::Zero(basis_num_, basis_num_);
    centrifugal_force_sandwich_[idx] = Vec::Zero(basis_num_);
  }

  for (int v = 0; v < vertex_num_; ++v) {
    inertial_force_sandwich_ += -mass_[v] * vert_basis_transpose_[v]; // (\sum_{all v} -(m_v * U_v))
    Vec3 r = MapVec3(rest_pos_ + v * 3) - center_of_mass_;
    Mat3 skew_symmetrix_matrix = GetSkewSymmetrixMatrix(r);
    // sandwich = \sum_{all v} (U_v^T * m_v * [r0])
    euler_force_sandwich_ += vert_basis_transpose_[v] * (mass_[v] * skew_symmetrix_matrix);
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        int idx = row * 3 + col;
        Mat3 rotation = Mat3::Zero();
        rotation(row, col) = mass_[v];
        // sandwich = \sum_{all v} (-U_v^T * m_v * [w]^2 * r0)
        centrifugal_force_sandwich_[idx] -= vert_basis_transpose_[v] * (rotation * r);
        // sandwich = \sum_{all v} -(U_v^T * m_v * [w] * U_v)
        coriolis_force_sandwich_[idx] -= vert_basis_transpose_[v] * rotation * vert_basis_[v];
      }
    }
  }
}

void SubspaceTet::ComputeMomentumAndTorqueMatrix() {
  momentum_matrix_ = Mat::Zero(3, basis_num_);
  torque_matrix_ = Mat::Zero(3, basis_num_);
  for (int v = 0; v < vertex_num_; ++v) {
    double mass = mass_[v];
    for (int c = 0; c < basis_num_; ++c) {
      for (int r = 0; r < 3; ++r) {
        momentum_matrix_(r, c) += mass * basis_(v * 3 + r, c);
      }
      Vec3 col(basis_(v * 3 + 0, c), basis_(v * 3 + 1, c), basis_(v * 3 + 2, c));
      col *= mass;
      Vec3 r = MapVec3(rest_pos_ + v * 3) - center_of_mass_;
      Vec3 cross = r.cross(col);
      torque_matrix_(0, c) += cross[0];
      torque_matrix_(1, c) += cross[1];
      torque_matrix_(2, c) += cross[2];
    }
  }
  P(momentum_matrix_.norm());
  P(torque_matrix_.norm());
  momentum_matrix_transpose_ = momentum_matrix_.transpose();
  torque_matrix_transpose_ = torque_matrix_.transpose();
}

void SubspaceTet::ComputeCenterOfMass() {
  center_of_mass_.setZero();
  total_mass_ = 0;
  for (int v = 0; v < vertex_num_; ++v) {
    center_of_mass_ += mass_[v] * MapVec3(rest_pos_ + v * 3);
    total_mass_ += mass_[v];
  }
  center_of_mass_ /= total_mass_;
}

void SubspaceTet::LoadBasis(int basis_id) {
  MapVec map_X(X, vertex_num_ * 3);
  MapVec map_rest_pos_(rest_pos_, vertex_num_ * 3);
  map_X = basis_.col(basis_id) + map_rest_pos_;
}

void SubspaceTet::NextBasis() {
  current_basis_ = (current_basis_ + 1) % basis_num_;
  //  LoadBasis(current_basis_);
}

void SubspaceTet::PrevBasis()
{
  current_basis_--;
  if (current_basis_ < 0) current_basis_ = basis_num_ - 1;
}
