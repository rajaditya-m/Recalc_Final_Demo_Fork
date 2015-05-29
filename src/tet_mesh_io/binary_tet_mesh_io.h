#ifndef BINARYTETMESHIO_H
#define BINARYTETMESHIO_H
#include <vector>
#include "tetrahedral_mesh_io.h"
#include "singleton.h"

class BinaryTetMeshIO : public TetrahedralMeshIO
{
public:
  void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
  void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);

protected:
  BinaryTetMeshIO() {}
  virtual ~BinaryTetMeshIO() {}
  DECLARE_SINGLETON_CLASS(BinaryTetMeshIO)
};

#endif // BINARYTETMESHIO_H
