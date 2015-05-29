#include "msh_tet_mesh_io.h"
#include <fstream>
#include "print_macro.h"

void MshTetMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  std::ifstream in(file_name);
  ASSERT(in.is_open());
  std::string line;
  in >> line;
  int v_num;
  in >> v_num;
  verts.resize(v_num * 3);
  for (int v = 0; v < v_num; ++v) {
    int junk;
    in >> junk  >> verts[v * 3 + 0] >> verts[v * 3 + 1] >> verts[v * 3 + 2];
  }
  in >> line >> line;
  int tet_num;
  in >> tet_num;
  tets.resize(tet_num * 4);
  for (int t = 0; t < tet_num; ++t) {
    int junk;
    in >> junk; in >> junk; in >> junk; in >> junk; in >> junk;
    in >> tets[t * 4 + 0];
    in >> tets[t * 4 + 1];
    in >> tets[t * 4 + 2];
    in >> tets[t * 4 + 3];
    tets[t * 4 + 0]--;
    tets[t * 4 + 1]--;
    tets[t * 4 + 2]--;
    tets[t * 4 + 3]--;
    ASSERT(!(tets[t * 4 + 0] == 0 && tets[t * 4 + 1] == 0 && tets[t * 4 + 2] == 0 && tets[t * 4 + 3] == 0));
  }
  in.close();
}

void MshTetMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tets) {
  ASSERT(false, P("NOT IMPLEMENTED"));
  (void) file_name;
  (void) verts;
  (void) tets;
}

MshTetMeshIO::MshTetMeshIO()
{
}
