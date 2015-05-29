#include <cstdio>
#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_set>
#include <fstream>
#include "multi_domain_basis_generator.h"
#include "basis_generator.h"
#include "tet_mesh_simulator_bridge.h"
#include "print_macro.h"
#include "basis_io.h"


using std::sprintf;

MultiDomainBasisGenerator::MultiDomainBasisGenerator(const char *mesh_file,
                                                     const char *partition_folder,
                                                     AffineTransformer<double> *transformer)
  : MultiDomainTet(mesh_file, 0, transformer) {
  LoadPartitionInfo(partition_folder);
}

void MultiDomainBasisGenerator::GenerateBasis(const char *folder, int linear_basis_num, int final_basis_num) {
  BasisGenerator basis_generator(inv_fem_->vega_mesh_);
  basis_generator.ComputeAllModes(linear_basis_num, final_basis_num);
  sprintf(file_name_, "%s/global_basis.bin", folder);
  WriteBasisInBinary(file_name_, vertex_num_, final_basis_num, &basis_generator.non_linear_modes_[0]);
  P(file_name_);
  DistributeGlobalBasis2Domains(&basis_generator.non_linear_modes_[0], final_basis_num, folder);
}

void MultiDomainBasisGenerator::DistributeGlobalBasis2Domains(double* global_basis_col_major, int basis_num, const char* basis_output_folder) {
  MapMatCol global_basis(global_basis_col_major, vertex_num_ * 3, basis_num);
  //  std::unordered_set<int> constrainTed({
  //    5 + 16 * 1, 6 + 16 * 1, 9 + 16 * 1, 10 + 16 * 1,
  //    5 + 16 * 5, 6 + 16 * 5, 9 + 16 * 5, 10 + 16 * 5,
  //    5 + 16 * 9, 6 + 16 * 9, 9 + 16 * 9, 10 + 16 * 9,
  //  });
  for (int p = 0; p < part_num_; ++p) {
    MatCol basis(vert_num_per_part_[p] * 3, basis_num);
    for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
      int global_v = vert_local_id2global_id_[p][local_v];
      //      if (is_constrainted_[global_v]) {
      //      if (constrainTed.count(global_v) > 0) {
      if (is_constrainted_[global_v]) {
        basis.row(local_v * 3 + 0).setZero();
        basis.row(local_v * 3 + 1).setZero();
        basis.row(local_v * 3 + 2).setZero();
      } else {
        basis.row(local_v * 3 + 0) = global_basis.row(global_v * 3 + 0);
        basis.row(local_v * 3 + 1) = global_basis.row(global_v * 3 + 1);
        basis.row(local_v * 3 + 2) = global_basis.row(global_v * 3 + 2);
      }
    }
    std::vector<double> verts, mass;
    std::vector<int> tets;
    GenerateDomainTetMesh(p, verts, tets, mass);
    BasisGenerator gen(int(verts.size()) / 3, &verts[0],
                       int(tets.size()) / 4, &tets[0],
                       inv_fem_->young_modulus_,
                       inv_fem_->poisson_ratio_,
                       inv_fem_->density_,
                       &mass[0]);
    std::vector<double> tmp_basis(basis.rows() * basis.cols());
    memcpy(&tmp_basis[0], basis.data(), sizeof(double) * tmp_basis.size());
    gen.MassPCA(basis_num, basis_num, tmp_basis, BasisGenerator::kSimulationData);
    // write eigen value
    {
      sprintf(file_name_, "%s/partition_%d.eigen_value.txt", basis_output_folder, p);
      P(file_name_);
      std::ofstream eig_txt(file_name_);
      sprintf(file_name_, "%s/partition_%d.eigen_value.bin", basis_output_folder, p);
      P(file_name_);
      std::ofstream eig_bin(file_name_, std::ios::binary);
      ASSERT(eig_txt.is_open() && eig_bin.is_open());
      ASSERT(int(gen.eigen_values_.size()) == basis_num);
      eig_txt << basis_num << "\n";
      eig_bin.write((char*) &basis_num, sizeof(int));
      for (int i = 0; i < basis_num; ++i) {
        eig_txt << gen.eigen_values_[i] << "\n";
        // P(i, gen.eigen_values_[i]);
      }
      eig_bin.write((char*) &gen.eigen_values_[0], sizeof(double) * basis_num);
      eig_txt.close();
      eig_bin.close();
    }

    sprintf(file_name_, "%s/partition_%d.basis.bin", basis_output_folder, p);
    WriteBasisInBinary(file_name_, vert_num_per_part_[p], basis_num, &tmp_basis[0]);
    P(file_name_);
    sprintf(file_name_, "%s/partition_%d.basis.txt", basis_output_folder, p);
    WriteBasisInText(file_name_, vert_num_per_part_[p], basis_num, &tmp_basis[0]);
    P(file_name_);
  }
}

void MultiDomainBasisGenerator::GenerateDomainTetMesh(int part_id, std::vector<double> &verts, std::vector<int> &tets, std::vector<double>& mass) {
  tets.clear();
  verts.resize(vert_num_per_part_[part_id] * 3);
  mass.resize(vert_num_per_part_[part_id]);
  for (int v = 0; v < vertex_num_; ++v) {
    if (vert_part_id_[v] == part_id) {
      int local_id = local_vert_idx_[v];
      verts[local_id * 3 + 0] = rest_pos_[v * 3 + 0];
      verts[local_id * 3 + 1] = rest_pos_[v * 3 + 1];
      verts[local_id * 3 + 2] = rest_pos_[v * 3 + 2];
      mass[local_id] = mass_[v];
    }
  }
  for (int t = 0; t < tet_number; ++t) {
    int* v = tet_ + t * 4;
    bool inside_domain = true;
    for (int i = 0; i < 4; ++i) {
      if (vert_part_id_[v[i]] != part_id) {
        inside_domain = false;
        break;
      }
    }
    if (inside_domain) {
      tets.push_back(local_vert_idx_[v[0]]);
      tets.push_back(local_vert_idx_[v[1]]);
      tets.push_back(local_vert_idx_[v[2]]);
      tets.push_back(local_vert_idx_[v[3]]);
    }
  }
}

void MultiDomainBasisGenerator::SetFixedVertex(std::vector<int> fixed_verts) {
  std::fill(is_constrainted_.begin(), is_constrainted_.end(), false);
  for (int v : fixed_verts) {
    is_constrainted_[v] = true;
  }
}

void MultiDomainBasisGenerator::SetFixedVertex(int fixed_vert_num_per_domain) {
  typedef std::pair<double, int> Pair;
  struct Comparator {
    bool operator()(const Pair& a, const Pair& b) const {
      return a.first < b.first;
    }
  };

  std::fill(is_constrainted_.begin(), is_constrainted_.end(), 0);
  for (int p = 0; p < part_num_; ++p) {
    std::priority_queue<Pair, std::vector<Pair>, Comparator> queue;
    for (int local_v = 0; local_v < vert_num_per_part_[p]; ++local_v) {
      int v = vert_local_id2global_id_[p][local_v];
      MapVec3 pos(rest_pos_ + v * 3);
      double distance = (pos - initial_center_of_mass_[p]).norm();
      if (int(queue.size()) < fixed_vert_num_per_domain) {
        queue.push(Pair(distance, v));
      } else {
        if (distance < queue.top().first) {
          queue.pop();
          queue.push(Pair(distance, v));
        }
      }
    }

    while (queue.size() != 0) {
      is_constrainted_[queue.top().second] = true;
      queue.pop();
    }
  }
}

MultiDomainBasisGenerator::~MultiDomainBasisGenerator() {
}
