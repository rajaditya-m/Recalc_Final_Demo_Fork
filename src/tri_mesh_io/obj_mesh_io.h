#ifndef OBJMESHIO_H
#define OBJMESHIO_H
#include "triangular_mesh_io.h"
#include "singleton.h"

class ObjMeshIO : public TriangularMeshIO {
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tri);
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tri);

private:
  ObjMeshIO() {}
  virtual ~ObjMeshIO() {}
  DECLARE_SINGLETON_CLASS(ObjMeshIO)
};

#endif // OBJMESHIO_H
