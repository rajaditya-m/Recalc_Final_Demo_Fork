#include <queue>
#include <utility>
#include <set>
#include "tet.h"

#include "vector_io.h"
#include "vector_lib.h"
#include "tet_mesh_simulator_bridge.h"
#include "print_macro.h"
#include "string_formatter.h"
#include "basis_io.h"
#include "global.h"
#include "basis_generator.h"
#include "single_domain_basis_generator.h"
#include "sparseMatrix.h"

#define ELT(numRows,i,j) (((long)j)*((long)numRows)+((long)i))

#define COLMAJORIDX

SingleDomainBasisGenerator::SingleDomainBasisGenerator(const char *mesh_file, AffineTransformer<double> *transformer)
  : SubspaceTet(mesh_file, 0, transformer) {
  fixed_verts_.clear();
}

void SingleDomainBasisGenerator::ProcessFixedVertex(const char *filename) {
    int numFixedV;
    int* fixedV;
    loadCommaList(filename,&numFixedV,&fixedV);
    fixed_verts_.resize(numFixedV);
    P(numFixedV);
    for(int i = 0 ;i < numFixedV; i++) {
        fixed_verts_[i] = fixedV[i] - 1;
    }
}

void SingleDomainBasisGenerator::preLoad(const char *basis_prefix) {
    basis_generator = new BasisGenerator(inv_fem_->vega_mesh_, &mass_[0]);
    if (fixed_verts_.size() != 0) {
      basis_generator->SetFixedVertices(std::set<int>(fixed_verts_.begin(), fixed_verts_.end()));
    }
    basis_generator->SetInterfaceVertices(interface_vertices_,interface_vertex_split_);
    basis_generator->preLoad(basis_prefix);
    basis_generator->setPrefix(basis_prefix);

}

void SingleDomainBasisGenerator::GenerateBasis(const char *basis_prefix, int linear_basis_num, int final_basis_num) {
  basis_generator = new BasisGenerator(inv_fem_->vega_mesh_, &mass_[0]);
  int numRows = vertex_num_ * 3;
  int numRows1 = domain_offset_toggle_ * 3;
  int numRows2 = numRows - numRows1;
  if (fixed_verts_.size() != 0) {
    basis_generator->SetFixedVertices(std::set<int>(fixed_verts_.begin(), fixed_verts_.end()));
  }
  basis_generator->setPrefix(basis_prefix);
  basis_generator->SetInterfaceVertices(interface_vertices_,interface_vertex_split_);
  basis_generator->ComputeAllModes(linear_basis_num, final_basis_num);

  sprintf(file_name_, "%s.basis.bin", basis_prefix);
  WriteBasisInBinary(file_name_, vertex_num_, final_basis_num, &(basis_generator->non_linear_modes_[0]));

  sprintf(file_name_, "%s.basis.txt", basis_prefix);
  WriteBasisInText(file_name_, vertex_num_, final_basis_num, &(basis_generator->non_linear_modes_[0]));

  double *basis_split1 = new double[numRows1 * final_basis_num];
  for(int i = 0; i < numRows1; i++) {
      for(int j = 0; j< final_basis_num;j++) {
          basis_split1[ELT(numRows1,i,j)]  = basis_generator->non_linear_modes_[ELT(numRows,i,j)];
      }
  }
  sprintf(file_name_, "%s.basis_1.bin", basis_prefix);
  WriteBasisInBinary(file_name_, domain_offset_toggle_, final_basis_num, basis_split1);

  double *basis_split2 = new double[numRows2 * final_basis_num];
  for(int i = 0; i < numRows2; i++) {
      for(int j = 0; j< final_basis_num;j++) {
          basis_split2[ELT(numRows2,i,j)]  = basis_generator->non_linear_modes_[ELT(numRows,i+numRows1,j)];
      }
  }
  sprintf(file_name_, "%s.basis_2.bin", basis_prefix);
  WriteBasisInBinary(file_name_, vertex_num_-domain_offset_toggle_, final_basis_num, basis_split2);
  sprintf(file_name_, "%s.basis_2.txt", basis_prefix);
  WriteBasisInText(file_name_,vertex_num_-domain_offset_toggle_, final_basis_num, basis_split2);

  sprintf(file_name_, "%s.nonlin_weights.bin", basis_prefix);
  dj::Write1DVectorBinary(file_name_, basis_generator->eigen_values_);

  sprintf(file_name_, "%s.pure_eigen_vals.bin", basis_prefix);
  P(file_name_);
  dj::Write1DVectorBinary(file_name_, basis_generator->pure_eigen_values_);

  sprintf(file_name_, "%s.pure_eigen_vecs.bin", basis_prefix);
  P(file_name_);
  dj::Write1DVectorBinary(file_name_, basis_generator->pure_eigen_vectors_);

  sprintf(file_name_, "%s.lin_freqs.bin", basis_prefix);
  P(file_name_);
  dj::Write1DVectorBinary(file_name_, basis_generator->frequencies_);

  /*
  sprintf(file_name_, "%s.eigen_value.txt", basis_prefix);
  P(file_name_);
  dj::Write1DVectorText(file_name_, basis_generator->eigen_values_);
  */
}

void SingleDomainBasisGenerator::SetFixedVertex(std::function<bool (int, double *)> IsFixed) {
  fixed_verts_.clear();
  for (int v = 0; v < vertex_num_; ++v) {
    if (IsFixed(v, X + v * 3)) {
      fixed_verts_.push_back(v);
      //      P(v, dj::Vec3d(X + v * 3));
    }
  }
  //  P(fixed_verts_.size());
}

void SingleDomainBasisGenerator::SetFixedVertex(std::vector<int> fixed_verts) {
  fixed_verts_ = fixed_verts;
}

void SingleDomainBasisGenerator::SetFixedVertex(int fixed_vert_num) {
  //  P(dj::Vec3d(&center_of_mass_[0]));
  typedef std::pair<double, int> Pair;
  struct Comparator {
    bool operator()(const Pair& a, const Pair& b) const {
      return a.first < b.first;
    }
  };

  std::fill(is_constrainted_.begin(), is_constrainted_.end(), 0);
  std::priority_queue<Pair, std::vector<Pair>, Comparator> queue;
  for (int v = 0; v < vertex_num_; ++v) {
    MapVec3 pos(rest_pos_ + v * 3);
    double distance = (pos - center_of_mass_).norm();
    if (int(queue.size()) < fixed_vert_num) {
      queue.push(Pair(distance, v));
    } else {
      if (distance < queue.top().first) {
        queue.pop();
        queue.push(Pair(distance, v));
      }
    }
  }

  fixed_verts_.clear();
  while (queue.size() != 0) {
    is_constrainted_[queue.top().second] = true;
    fixed_verts_.push_back(queue.top().second);
    queue.pop();
  }
}

SingleDomainBasisGenerator::~SingleDomainBasisGenerator() {
}


void SingleDomainBasisGenerator::RegenerateStitchBasis(int modeID) {
    basis_generator->RegenerateAllModes(modeID);
}

void SingleDomainBasisGenerator::setStitchedStiffnessMatrix(SparseMatrix *spmat) {
    basis_generator->setStitchStiffnessMatrix(spmat);
}
