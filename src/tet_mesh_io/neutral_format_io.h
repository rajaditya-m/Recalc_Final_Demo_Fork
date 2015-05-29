#ifndef NEUTRALFORMATIO_H
#define NEUTRALFORMATIO_H
#include "tetrahedral_mesh_io.h"

class NeutralFormatIO : public TetrahedralMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets);
private:
  NeutralFormatIO() {}
  virtual ~NeutralFormatIO() {}
};

#endif // NEUTRALFORMATIO_H
