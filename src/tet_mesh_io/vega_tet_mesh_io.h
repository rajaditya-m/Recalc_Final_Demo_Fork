#ifndef VEGATETMESHIO_H
#define VEGATETMESHIO_H
#include "tetrahedral_mesh_io.h"
#include "singleton.h"

class VegaTetMeshIO : public TetrahedralMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);

private:
  VegaTetMeshIO() {}
  virtual ~VegaTetMeshIO() {}
  DECLARE_SINGLETON_CLASS(VegaTetMeshIO);
};

#endif // VEGATETMESHIO_H
