#ifndef TET_GEN_MESH_IO_H
#define TET_GEN_MESH_IO_H
#include <vector>
#include "tetrahedral_mesh_io.h"
#include "singleton.h"
class TetGenMeshIO : public TetrahedralMeshIO
{
public:
  void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
  void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);

protected:
  TetGenMeshIO() {}
  virtual ~TetGenMeshIO() {}
  DECLARE_SINGLETON_CLASS(TetGenMeshIO)
};


#endif // TET_GEN_MESH_IO_H
