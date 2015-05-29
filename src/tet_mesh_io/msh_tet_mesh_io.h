#ifndef MSH_TET_MESH_IO_H
#define MSH_TET_MESH_IO_H
#include "tetrahedral_mesh_io.h"
#include "singleton.h"

class MshTetMeshIO : public TetrahedralMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);

private:
  MshTetMeshIO();
  virtual ~MshTetMeshIO() {}
  DECLARE_SINGLETON_CLASS(MshTetMeshIO);
};

#endif // MSH_TET_MESH_IO_H
