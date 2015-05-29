#ifndef CUBICA_TET_MESH_IO_H
#define CUBICA_TET_MESH_IO_H
#include "tetrahedral_mesh_io.h"
#include "singleton.h"

class CubicaTetMeshIO : public TetrahedralMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);

protected:
  CubicaTetMeshIO();
  virtual ~CubicaTetMeshIO() {}
  DECLARE_SINGLETON_CLASS(CubicaTetMeshIO)
};

#endif // CUBICA_TET_MESH_IO_H
