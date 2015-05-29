#ifndef CUBICABASISGENERATOR_H
#define CUBICABASISGENERATOR_H
#include "cubica_tet.h"
class BasisGenerator;
class CubicaBasisGenerator : public CubicaTet {
public:
  CubicaBasisGenerator(const char* mesh_file, const char *folder);
  void GenerateBasis(const char* folder, int linear_mode_num, int final_basis_num_);
  void DistributeGlobalBasisToDomains(int basis_num, double* basis, const char *basis_output_folder);
  virtual ~CubicaBasisGenerator();
};

#endif // CUBICABASISGENERATOR_H
