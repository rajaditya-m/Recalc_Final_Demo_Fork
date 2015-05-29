#ifndef OFFMESHIO_H
#define OFFMESHIO_H
#include "triangular_mesh_io.h"
#include "singleton.h"

class OffMeshIO : public TriangularMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tri);
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tri);

private:
  OffMeshIO() {}
  virtual ~OffMeshIO() {}
  DECLARE_SINGLETON_CLASS(OffMeshIO)
};

#endif // OFFMESHIO_H
