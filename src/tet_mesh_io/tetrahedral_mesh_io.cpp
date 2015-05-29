#include "tetrahedral_mesh_io.h"
void TetrahedralMeshIO::Read(const char *file_name, std::vector<float> &verts, std::vector<int> &tets)
{
  std::vector<double> double_verts;
  Read(file_name, double_verts, tets);
  verts.resize(double_verts.size());
  for (int i = 0; i < (int) double_verts.size(); ++i) {
    verts[i] = (float) double_verts[i];
  }
}


void TetrahedralMeshIO::Write(const char *file_name, std::vector<float> &verts, std::vector<int> &tets)
{
  std::vector<double> double_verts(verts.size());
  for (int i = 0; i < (int) verts.size(); ++i) {
    double_verts[i] = verts[i];
  }
  Write(file_name, double_verts, tets);
}

