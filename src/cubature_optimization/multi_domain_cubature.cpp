#include "multi_domain_cubature.h"
#include <fstream>
#include "global.h"
#include <map>
#include "tet_mesh_simulator_bridge.h"
#include "string_formatter.h"
#include "random.h"
#include "vector_lib.h"
#include "MERSENNETWISTER.h"

MultiDomainCubature::MultiDomainCubature(const char *mesh_file_name)
  : MultiDomainTet(mesh_file_name, 0, NULL) {
  rigid_cubature_ = false;
}

void MultiDomainCubature::SetFolder(string subspace_pose_file, string partition_folder, string output_folder) {
  output_folder_ = output_folder;
  // over write error file
  std::string file_name = dj::Format("%z/cubature/error.txt", GetDataFolder());
  std::ofstream out(file_name);
  out << "cubature error:" << std::endl;
  out.close();
  ReadInterfaceInfo(partition_folder);
  ReadSubspacePoseData(subspace_pose_file);
}

void MultiDomainCubature::SetFolder(int num_of_samples, string partition_folder, string output_folder,
                                    double scale, bool read_eigen_value, bool over_write_error_file) {
  output_folder_ = output_folder;
  // over write error file
  if (over_write_error_file) {
    std::string file_name = dj::Format("%z/cubature/error.txt", GetDataFolder());
    std::ofstream out(file_name);
    out << "cubature error:" << std::endl;
    out.close();
    ReadInterfaceInfo(partition_folder);
  }

  pose_num_ = num_of_samples;
  pose_q_.resize(pose_num_);
  pose_rotation_.resize(pose_num_);
  pose_center_of_mass_.resize(pose_num_);

  std::vector<Vec> eigen_values(part_num_);
  if (read_eigen_value) {
    for (int p = 0; p < part_num_; ++p) {
      std::string eigen_value_file = dj::Format("%z/modal_basis/partition_%d.eigen_value.bin", GetDataFolder(), p);
      std::ifstream eig(eigen_value_file, std::ios::binary);
      ASSERT(eig.is_open(), P(eigen_value_file));
      int basis_num;
      eig.read((char*) &basis_num, sizeof(int));
      ASSERT(basis_num == part_basis_size_[p], P(basis_num, part_basis_size_[p], p));
      eigen_values[p] = Vec::Zero(basis_num);
      eig.read((char*) &eigen_values[p][0], sizeof(double) * basis_num);
      eig.close();
    }
  } else {
    for (int p = 0; p < part_num_; ++p) {
      eigen_values[p] = Vec::Zero(part_basis_size_[p]);
      for (int i = 0; i < eigen_values[p].size(); ++i) {
        eigen_values[p][i] = 1.0;
      }
    }
  }

  std::ofstream pose_file(dj::Format("%z/random_pose.txt", GetDataFolder()));
  pose_file << pose_num_ << "\n";
  pose_file << part_num_ << "\n";
  for (int p = 0; p < part_num_; ++p) {
    pose_file << part_basis_size_[p] << "\n";
  }
  MERSENNETWISTER rand(time(NULL));
  for (int pose = 0; pose < pose_num_; ++pose) {
    pose_rotation_[pose].resize(part_num_);
    pose_q_[pose].resize(part_num_);
    pose_center_of_mass_[pose].resize(part_num_);
    for (int p = 0; p < part_num_; ++p) {
      pose_rotation_[pose][p].setIdentity();
      pose_center_of_mass_[pose][p] = initial_center_of_mass_[p];
    }

    // genreate random q
    for (int attempt = 0; true; ++attempt) {
      for (int p = 0; p < part_num_; ++p) {
        pose_q_[pose][p] = Vec::Zero(part_basis_size_[p]);
        for (int n = 0; n < part_basis_size_[p]; ++n) {
          //          pose_q_[pose][p][n] = 0.0010 * sqrt(eigen_values[p][n] / eigen_values[p][0]) / 4.0 * rand.randNorm();
          pose_q_[pose][p][n] = scale * sqrt(eigen_values[p][n] / eigen_values[p][0]) / 4.0 * rand.randNorm();
          q_[basis_offset_[p] + n] = pose_q_[pose][p][n];
        }
      }
      //      PVEC(q_);
      UpdatePosition();
      bool is_inverted = false;
      int t = 0;
      for (; t < tet_number; ++t) {
        int* verts = tet_ + t * 4;
        MapVec3 pos[4] = {
          MapVec3(X + verts[0] * 3),
          MapVec3(X + verts[1] * 3),
          MapVec3(X + verts[2] * 3),
          MapVec3(X + verts[3] * 3),
        };
        Vec3 v01 = pos[1] - pos[0];
        Vec3 v02 = pos[2] - pos[0];
        Vec3 v03 = pos[3] - pos[0];
        Vec3 cross = v01.cross(v02);
        if (cross.dot(v03) < 0) {
          is_inverted = true;
          break;
        }
      }
      if (attempt % 100 == 0)
        P(pose, attempt, q_.norm(), t);
      if (!is_inverted) {
        break;
      }
      //      break;
    }
    for (int n = 0; n < total_basis_num_; ++n) {
      pose_file << q_[n] << " ";
    }
    pose_file << "\n";
  }
}

void MultiDomainCubature::GenerateDomainCubature(std::vector<int> max_cubature_point, double relative_error) {
  for (int p = 0; p < part_num_; ++p) {
    file_prefix_ = dj::Format("partition_%d", p);
    current_domain_ = p;
    LoadDomainPartitionInfo(file_prefix_);
    L(dj::Format("Optimizing cubature for domain %d", p));
    std::vector<VECTOR*> training_set;
    VECTOR training_force(pose_num_ * part_basis_size_[p]);
    // generate training data
    for (int pose = 0; pose < pose_num_; ++pose) {
      // update local offset
      for (int v : vert_local_id2global_id_[p]) {
        MapVec3((double*) & (inv_fem_->u_[v][0])) = vert_basis_[v] * pose_q_[pose][p];
      }
      // compute force
      inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(local_tet2global_tet_);
      // assenble subspace force
      Vec subspace_force = Vec::Zero(part_basis_size_[p]);
      for (int t : local_tet2global_tet_) {
        int* verts = tet_ + t * 4;
        for (int i = 0; i < 4; ++i) {
          // vega internal force is in negative direction
          subspace_force -= vert_basis_transpose_[verts[i]] * MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
        }
      }
      memcpy(&training_force(pose * part_basis_size_[p]), &subspace_force[0], sizeof(double) * part_basis_size_[p]);
      VECTOR* q = new VECTOR(part_basis_size_[p]);
      memcpy(&(*q)(0), &pose_q_[pose][p][0], sizeof(double) * part_basis_size_[p]);
      training_set.push_back(q);
      // verify eval
      if (0) {
        int basis_size = part_basis_size_[p];
        VECTOR gout(basis_size);
        Vec sub = Vec::Zero(basis_size);
        for (int t = 0; t < int(local_tet2global_tet_.size()); ++t) {
          evalPointForceDensity(t, *q, gout, pose);
          sub += MapVec(&gout(0), basis_size);
        }
        Vec diff = subspace_force - sub;
        ASSERT(diff.norm() < 1e-5, P(diff.norm()));
      }
    }
    int validation_sample_num = int(0.05 * pose_num_);
    if (validation_sample_num == 0) validation_sample_num = 3;
    run(training_set, training_force, relative_error, max_cubature_point[p], 1000, part_basis_size_[p] / 2, pose_num_ - validation_sample_num);
    {
      std::string file_name = dj::Format("%z/cubature/error.txt", GetDataFolder());
      std::ofstream out(file_name, std::ios::app);
      out << std::endl;
      out.close();
    }
    for (VECTOR * q : training_set) {
      delete q;
    }
  }
}

void MultiDomainCubature::GenerateInterfaceCubature(std::vector<int> max_cubature_point, double relative_error, bool is_rigid_interface) {
  // for each interface domain
  rigid_cubature_ = is_rigid_interface;
  for (int e_id = 0; e_id < interface_num_; ++e_id) {
    if (!is_rigid_interface) {
      file_prefix_ = dj::Format("interface_%z", interface_str_[e_id]);
    } else {
      file_prefix_ = dj::Format("rigid_interface_%z", interface_str_[e_id]);
    }
    current_domain_ = e_id + part_num_;
    LoadDomainPartitionInfo(dj::Format("interface_%z", interface_str_[e_id]));
    interface_opt_stage_ = 0;
    L(dj::Format("Optimizing cubature for %z", file_prefix_));

    int basis_size = 0;
    int q_offset[5] = {0};
    for (int i = 1; i <= int(interface_domains_[e_id].size()); ++i) {
      int part_id = interface_domains_[e_id][i - 1];
      basis_size += part_basis_size_[part_id] ;
      q_offset[i] = q_offset[i - 1] + part_basis_size_[part_id];
    }

    std::vector<VECTOR*> training_set;
    VECTOR training_force(pose_num_ * basis_size);
    for (int pose = 0; pose < pose_num_; ++pose) {
      // update position
      for (int v : local_vert2global_vert_) {
        int p = vert_part_id_[v];
        Vec3 local_u = vert_basis_[v] * pose_q_[pose][p];
        local_u = pose_rotation_[pose][p] * (local_u + vert_offset_from_mass_center_[v]) +
                  pose_center_of_mass_[pose][p] - MapVec3(rest_pos_ + v * 3);
        MapVec3((double*) &inv_fem_->u_[v][0]) = local_u;
      }
      // compute force
      inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(local_tet2global_tet_);
      // assemble subspace force
      Vec subspace_force = Vec::Zero(basis_size);
      for (int t : local_tet2global_tet_) {
        int* verts = tet_ + t * 4;
        for (int i = 0; i < 4; ++i) {
          int v = verts[i];
          int p = vert_part_id_[v];
          Vec force;
          if (!is_rigid_interface) {
            force = -vert_basis_transpose_[verts[i]] * pose_rotation_[pose][p].transpose() * MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
          } else {
            force = Vec::Zero(6);
            Eigen::Vector3d r = MapVec3((double*) &inv_fem_->u_[v][0]) + MapVec3(rest_pos_ + v * 3) - pose_center_of_mass_[pose][p];
            MapVec3(force.data() + 0) = -pose_rotation_[pose][p].transpose() * MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
            MapVec3(force.data() + 3) = r.cross(MapVec3(force.data() + 0));
          }

          int idx = 0;
          for (; idx < int(interface_domains_[e_id].size()) && interface_domains_[e_id][idx] != p; ++idx) {
          }
          ASSERT(idx < int(interface_domains_[e_id].size()));
          MapVec(&subspace_force[q_offset[idx]], part_basis_size_[p]) += force;
        }
      }
      memcpy(&training_force(pose * basis_size), &subspace_force[0], sizeof(double) * basis_size);
      VECTOR* q = new VECTOR(basis_size);
      for (int i = 0; i < int(interface_domains_[e_id].size()); ++i) {
        int p = interface_domains_[e_id][i];
        MapVec(&(*q)(q_offset[i]), part_basis_size_[p]) = pose_q_[pose][p];
      }
      training_set.push_back(q);
      // verify eval
      {
        VECTOR gout(basis_size);
        Vec sub = Vec::Zero(basis_size);
        for (int t = 0; t < int(local_tet2global_tet_.size()); ++t) {
          evalPointForceDensity(t, *q, gout, pose);
          sub += MapVec(&gout(0), basis_size);
        }
        Vec diff = subspace_force - sub;
        ASSERT(diff.norm() < 1e-5, P(diff.norm()));
      }
    }
    //    exit(0);
    int validation_sample_num = int(0.05 * pose_num_);
    if (validation_sample_num == 0) validation_sample_num = 3;
    int max_cubature = dj::Min(int(local_tet2global_tet_.size()), max_cubature_point[e_id]);
    this->run(training_set, training_force, relative_error, max_cubature, 1000, basis_size / 2, pose_num_ - validation_sample_num);
    {
      std::string file_name = dj::Format("%z/cubature/error.txt", GetDataFolder());
      std::ofstream out(file_name, std::ios::app);
      out << std::endl;
      out.close();
    }
    for (VECTOR * q : training_set) {
      delete q;
    }
  }
}

void MultiDomainCubature::GenerateAllCubature(std::vector<int> domain_max_cubature,
                                              double domain_relative_error,
                                              std::vector<int> interface_max_cubature,
                                              double interface_relative_error) {
  ASSERT(int(domain_max_cubature.size()) == part_num_);
  ASSERT(int(interface_max_cubature.size()) == interface_num_);
  GenerateDomainCubature(domain_max_cubature, domain_relative_error);
  GenerateInterfaceCubature(interface_max_cubature, interface_relative_error);
}

void MultiDomainCubature::LoadDomainPartitionInfo(std::string file_prefix) {
  {
    std::string file_name = dj::Format("%z/%z.tet_global_id.txt", partition_folder_, file_prefix);
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    int tet_num;
    in >> tet_num;
    local_tet2global_tet_.resize(tet_num);
    for (int i = 0; i < tet_num; ++i) {
      in >> local_tet2global_tet_[i];
      ASSERT(local_tet2global_tet_[i] < tet_number, P(local_tet2global_tet_[i], i));
    }
    in.close();
  }
  {
    std::string file_name = dj::Format("%z/%z.vert_global_id.txt", partition_folder_, file_prefix);
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    int v_num;
    in >> v_num;
    local_vert2global_vert_.resize(v_num);
    for (int i = 0; i < v_num; ++i) {
      in >> local_vert2global_vert_[i];
    }
    in.close();
  }
}

void MultiDomainCubature::ReadInterfaceInfo(string partition_folder) {
  partition_folder_ = partition_folder;
  std::string file_name = dj::Format("%z/interface_domains.txt", partition_folder);
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  int part_num;
  in >> part_num;
  ASSERT(part_num_ == part_num);
  in >> interface_num_;
  interface_domains_.resize(interface_num_);
  interface_str_.resize(interface_num_);
  for (int interface = 0; interface < interface_num_; ++interface) {
    int interface_id;
    int p[4] = { -1, -1, -1, -1};
    in >> interface_id;
    int domain_count = 0;
    for (int j = 0; j < 4; ++j) {
      in >> p[j];
      if (p[j] != -1) domain_count++;
    }
    ASSERT(domain_count >= 2);
    interface_str_[interface] = dj::Format("%d_%d_%d_%d", p[0], p[1], p[2], p[3]);
    interface_domains_[interface] = std::vector<int>(p, p + domain_count);
  }
  in.close();

}

void MultiDomainCubature::ReadSubspacePoseData(string subspace_pose_file) {
  subspace_pose_file_ = subspace_pose_file;
  std::ifstream in(subspace_pose_file);
  ASSERT(in.is_open(), P(subspace_pose_file));
  int v_num;
  in >> v_num;
  ASSERT(v_num == vertex_num_);
  in >> pose_num_;
  int p_num;
  in >> p_num;
  ASSERT(part_num_ == p_num);
  for (int i = 0; i < part_num_; ++i) {
    int bone_id, basis_size;
    in >> bone_id >> basis_size;
    ASSERT(part_basis_size_[i] == basis_size);
  }
  pose_q_.resize(pose_num_);
  pose_rotation_.resize(pose_num_);
  pose_center_of_mass_.resize(pose_num_);

  for (int pose = 0; pose < pose_num_; ++pose) {
    pose_rotation_[pose].resize(part_num_);
    pose_q_[pose].resize(part_num_);
    pose_center_of_mass_[pose].resize(part_num_);
    for (int p = 0; p < part_num_; ++p) {
      int attached_bone;
      in >> attached_bone;
      // part rotation
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          in >> pose_rotation_[pose][p](r, c);
        }
      }
      // part mass center
      Vec3 initial_center_of_mass;
      in >> initial_center_of_mass[0];
      in >> initial_center_of_mass[1];
      in >> initial_center_of_mass[2];
      ASSERT((initial_center_of_mass - initial_center_of_mass_[p]).norm() < 1e-5);
      in >> pose_center_of_mass_[pose][p][0];
      in >> pose_center_of_mass_[pose][p][1];
      in >> pose_center_of_mass_[pose][p][2];

      pose_q_[pose][p] = Vec::Zero(part_basis_size_[p]);
      for (int n = 0; n < part_basis_size_[p]; ++n) {
        in >> pose_q_[pose][p][n];
      }
    }
  }
  in.close();
}

int MultiDomainCubature::numTotalPoints() {
  return int(local_tet2global_tet_.size());
}

void MultiDomainCubature::evalPointForceDensity(int pointId, VECTOR & q, VECTOR & gOut, int poseIdx) {
  int global_tet_id = local_tet2global_tet_[pointId];
  int* verts = tet_ + global_tet_id * 4;
  if (current_domain_ < part_num_) {
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      int p = vert_part_id_[v];
      ASSERT(p == current_domain_);
      MapVec3((double*) & (inv_fem_->u_[v][0])) = vert_basis_[v] * MapVec(&q(0), part_basis_size_[p]);
    }
    std::vector<int> one_tet(1, global_tet_id);
    inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(one_tet);
    Vec subspace_force = Vec::Zero(part_basis_size_[current_domain_]);
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      // vega internal force is in opposite direction
      subspace_force -= vert_basis_transpose_[v] * MapVec3(inv_fem_->element_force_ + global_tet_id * 12 + i * 3);
    }
    memcpy(&gOut(0), &subspace_force[0], part_basis_size_[current_domain_] * sizeof(double));
  } else {
    int e_id = current_domain_ - part_num_;
    int basis_size = 0;
    int q_offset[5] = {0};
    for (int i = 1; i <= int(interface_domains_[e_id].size()); ++i) {
      int part_id = interface_domains_[e_id][i - 1];
      basis_size += part_basis_size_[part_id] ;
      q_offset[i] = q_offset[i - 1] + part_basis_size_[part_id];
    }
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      int p = vert_part_id_[v];
      int idx = 0;
      for (; idx < int(interface_domains_[e_id].size()) && interface_domains_[e_id][idx] != p; ++idx) {
      }
      ASSERT(idx < int(interface_domains_[e_id].size()));
      Vec3 local_u = vert_basis_[v] * MapVec(&q(q_offset[idx]), part_basis_size_[p]);
      //      P(q_offset[idx], idx);
      //      PVEC(MapVec(&q(0), 80));
      //      PVEC(local_u);
      local_u = pose_rotation_[poseIdx][p] * (local_u + vert_offset_from_mass_center_[v]) +
                pose_center_of_mass_[poseIdx][p] - MapVec3(rest_pos_ + v * 3);
      //      PVEC(local_u);
      MapVec3((double*) &inv_fem_->u_[v][0]) = local_u;
    }

    std::vector<int> one_tet(1, global_tet_id);
    inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(one_tet);
    //    memset(&q(0), 0, sizeof(double) * basis_size);
    gOut.clear();
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      int p = vert_part_id_[v];
      int idx = 0;
      for (; idx < int(interface_domains_[e_id].size()) && interface_domains_[e_id][idx] != p; ++idx) {
      }
      ASSERT(idx < int(interface_domains_[e_id].size()));
      Vec subspace_force;
      if (!rigid_cubature_) {
        subspace_force = vert_basis_transpose_[v] *
                         pose_rotation_[poseIdx][p].transpose() *
                         MapVec3(inv_fem_->element_force_ + global_tet_id * 12 + i * 3);
      } else {
        subspace_force = Vec::Zero(6);
        Eigen::Vector3d r = MapVec3((double*) &inv_fem_->u_[v][0]) + MapVec3(rest_pos_ + v * 3) - pose_center_of_mass_[poseIdx][p];
        MapVec3(subspace_force.data() + 0) = pose_rotation_[poseIdx][p].transpose() * MapVec3(inv_fem_->element_force_ + global_tet_id * 12 + i * 3);
        MapVec3(subspace_force.data() + 3) = r.cross(MapVec3(subspace_force.data() + 0));
      }
      // vega internal force is in opposite direction
      MapVec(&gOut(q_offset[idx]), part_basis_size_[p]) -= subspace_force;
    }
  }
}

void MultiDomainCubature::handleCubature(std::vector<int> &selectedPoints, VECTOR & weights, Real relErr) {
  {
    std::string error_file = dj::Format("%z/cubature/error.txt", GetDataFolder());
    std::ofstream out(error_file, std::ios::app);
    ASSERT(out.is_open(), P(error_file));
    out << dj::Format("%z : %d/%d\t %z\n", file_prefix_, selectedPoints.size(), numTotalPoints(), relErr);
    out.close();
  }
  {
    std::string local_cubature_file = dj::Format("%z/cubature/%z.cubature.local_tet.txt", GetDataFolder(), file_prefix_);
    std::string global_cubature_file = dj::Format("%z/cubature/%z.cubature.global_tet.txt", GetDataFolder(), file_prefix_);
    std::ofstream local_out(local_cubature_file);
    std::ofstream global_out(global_cubature_file);
    ASSERT(local_out.is_open(), P(local_cubature_file));
    ASSERT(global_out.is_open(), P(local_cubature_file));
    local_out << selectedPoints.size() << std::endl;
    global_out << selectedPoints.size() << std::endl;
    for (int i = 0; i < int(selectedPoints.size()); ++i) {
      local_out << selectedPoints[i] << " " << weights(i) << std::endl;
      global_out << local_tet2global_tet_[selectedPoints[i]] << " " << weights(i) << std::endl;
    }
    local_out.close();
    global_out.close();
  }
  L(dj::Format("Write cubature optimization result to file %z", file_prefix_));
}

MultiDomainCubature::~MultiDomainCubature() {
}
