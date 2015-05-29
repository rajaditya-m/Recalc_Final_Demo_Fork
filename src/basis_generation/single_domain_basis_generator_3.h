#ifndef SINGLEDOMAINBASISGENERATOR_H
#define SINGLEDOMAINBASISGENERATOR_H
#include <vector>
#include <functional>
#include "subspace_tet.h"
#include "basis_generator.h"

class SingleDomainBasisGenerator : public SubspaceTet
{
public:
  SingleDomainBasisGenerator(const char* mesh_file,
                            AffineTransformer<double>* transformer = NULL);
  void GenerateBasis(const char *basis_prefix, int linear_basis_num, int final_basis_num);
  void RegenerateStitchBasis(int id);
  void preLoad(const char* prefix);
  void setStitchedStiffnessMatrix(SparseMatrix* spmat);
  void SetFixedVertex(std::function<bool (int, double *)> IsFixed);
  void SetFixedVertex(std::vector<int> fixed_verts);
  void ProcessFixedVertex(const char *filename);
  /// select the $fixed_vert_num_per_domain$ number vertex that are closest to
  /// the center of mass for each partition as fixed vertices
  void SetFixedVertex(int fixed_vert_num);
  virtual ~SingleDomainBasisGenerator();
  std::vector<int> fixed_verts_;
  char file_name_[1024];

  BasisGenerator *basis_generator;
};

#endif // SINGLEDOMAINBASISGENERATOR_H
