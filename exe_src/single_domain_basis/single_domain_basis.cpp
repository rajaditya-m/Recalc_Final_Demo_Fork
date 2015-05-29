#include <string>
#include "print_macro.h"


#include <fstream>
#include "string_formatter.h"
#include "config_file.h"
#include "global.h"
#include "single_domain_basis_generator.h"
#include "vector_lib.h"


int main(void) {
  SingleDomainBasisGenerator generator(GetMeshFile());
  std::string basis_prefix = dj::Format("%z/modal_basis/genBasis",GetDataFolder());
  std::string mass_file = dj::Format("%z/vertexmass.M",GetDataFolder());
  generator.LoadMass(mass_file.c_str());
  std::string fixed_vertex_file = dj::Format("%z/fixed_verts.bou",GetModelFolder());
  generator.ProcessFixedVertex(fixed_vertex_file.c_str());
  generator.GenerateBasis(basis_prefix.c_str(),30,60);
  return 0;
}

