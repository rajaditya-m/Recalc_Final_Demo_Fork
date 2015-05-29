#include "cubica_basis_generator.h"
#include "basis_generator.h"
#include "tet_mesh_simulator_bridge.h"
#include "subspace_tet.h"
#include <Eigen/Dense>
#include "basis_io.h"
#include "print_macro.h"


CubicaBasisGenerator::CubicaBasisGenerator(const char *mesh_file, const char *folder)
  : CubicaTet(mesh_file, folder) {
}

void CubicaBasisGenerator::GenerateBasis(const char *folder, int linear_mode_num, int final_basis_num) {
  //  std::set<int> fixed_vertex(constrainted_vertex_.begin(), constrainted_vertex_.end());
  //  for (int v : constrainted_vertex_) {
  //    if (vert_partition_info_[v * 4 + 2] != -1) {
  //      fixed_vertex.erase(v);
  //      KK;
  //    }
  //  }
  //  P(fixed_vertex.size());
  BasisGenerator basis_gen_(inv_fem_->vega_mesh_);
  //  basis_gen_.SetFixedVertices(fixed_vertex);
  basis_gen_.ComputeAllModes(linear_mode_num, final_basis_num);
  if (0) {
    char file_name[512];
    sprintf(file_name, "%s/global_basis.bin", folder);
    WriteBasisInBinary(file_name, vertex_num_, final_basis_num, &basis_gen_.non_linear_modes_[0]);
    //    exit(0);
  }
  DistributeGlobalBasisToDomains(final_basis_num, &basis_gen_.non_linear_modes_[0], folder);
}

void CubicaBasisGenerator::DistributeGlobalBasisToDomains(int basis_num, double *basis, const char* basis_output_folder) {
  // Distribute global basis to local basis
  MapMatCol global_basis(basis, vertex_num_ * 3, basis_num);
  for (int p = 0; p < part_num_; ++p) {
    BasisGenerator part_basis_gen_(domain_[p]->inv_fem_->vega_mesh_);
    MatCol basis(domain_[p]->vertex_num_ * 3, basis_num);
    for (int v = 0; v < domain_[p]->vertex_num_; ++v) {
      int global_v = vert_local_id2global_id_[p][v];
      if (is_constrainted_[global_v] && vert_partition_info_[global_v * 4 + 2] == -1) {
        basis.row(v * 3 + 0).setZero();
        basis.row(v * 3 + 1).setZero();
        basis.row(v * 3 + 2).setZero();
      } else {
        basis.row(v * 3 + 0) = global_basis.row(global_v * 3 + 0);
        basis.row(v * 3 + 1) = global_basis.row(global_v * 3 + 1);
        basis.row(v * 3 + 2) = global_basis.row(global_v * 3 + 2);
      }
    }
    std::vector<double> tmp_basis(basis.rows() * basis.cols());
    memcpy(&tmp_basis[0], basis.data(), sizeof(double) * tmp_basis.size());
    part_basis_gen_.MassPCA(basis_num, basis_num, tmp_basis, BasisGenerator::kSimulationData);
    //    double* dummy = &tmp_basis[0];
    //    memcpy(basis.data(), &tmp_basis[0], sizeof(double) * tmp_basis.size());
    char file_name[512];
    sprintf(file_name, "%s/partition_%d.basis.bin", basis_output_folder, parts_[p]);
    WriteBasisInBinary(file_name, domain_[p]->vertex_num_, basis_num, &tmp_basis[0]);
    sprintf(file_name, "%s/partition_%d.basis.txt", basis_output_folder, parts_[p]);
    WriteBasisInText(file_name, domain_[p]->vertex_num_, basis_num, &tmp_basis[0]);
  }
}

CubicaBasisGenerator::~CubicaBasisGenerator() {
}

