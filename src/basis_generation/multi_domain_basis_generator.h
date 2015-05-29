#ifndef MULTIDOMAINBASISGENERATOR_H
#define MULTIDOMAINBASISGENERATOR_H
#include <vector>
#include "multi_domain_tet.h"

class MultiDomainBasisGenerator : public MultiDomainTet {
public:
  MultiDomainBasisGenerator(const char* mesh_file,
                            const char* partition_folder,
                            AffineTransformer<double>* transformer = NULL);
  void GenerateBasis(const char *folder, int linear_basis_num, int final_basis_num);
  void DistributeGlobalBasis2Domains(double* global_basis_col_major, int basis_num, const char *basis_output_folder);
  void GenerateDomainTetMesh(int part_id, std::vector<double> &verts, std::vector<int>& tets, std::vector<double> &mass);
  void SetFixedVertex(std::vector<int> fixed_verts);
  /// select the $fixed_vert_num_per_domain$ number vertex that are closest to
  /// the center of mass for each partition as fixed vertices
  void SetFixedVertex(int fixed_vert_num_per_domain);
  virtual ~MultiDomainBasisGenerator();
  char file_name_[1024];
};

#endif // MULTIDOMAINBASISGENERATOR_H
