#include "cubica_tet_mesh_io.h"
#include <fstream>
#include <vector>
#include "vector_io.h"
#include "binary_file_io.h"

CubicaTetMeshIO::CubicaTetMeshIO() {
}


void CubicaTetMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tets) {
  BinaryFileReader in(file_name);
  int contrained_vert_num, uncontrained_vert_num;
  in.Read(&uncontrained_vert_num, 1);
  in.Read(&contrained_vert_num, 1);
  int total_vert_num = uncontrained_vert_num + contrained_vert_num;
  verts.resize(total_vert_num * 3);
  in.Read(&verts[0], total_vert_num * 3);
  int tet_num;
  in.Read(&tet_num, 1);
  tets.resize(tet_num * 4);
  in.Read(&tets[0], tet_num * 4);
}


void CubicaTetMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tets) {
  ASSERT(false, P("NOT IMPLEMENTED"));
  (void) file_name;
  (void) verts;
  (void) tets;
}
