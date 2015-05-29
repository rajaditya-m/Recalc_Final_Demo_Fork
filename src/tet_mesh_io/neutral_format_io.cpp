#include "neutral_format_io.h"
#include <fstream>
#include "print_macro.h"

void NeutralFormatIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  int v_num, tet_num;
  in >> v_num;
  verts.resize(v_num * 3);
  for (int i = 0; i < v_num * 3; ++i) {
    in >> verts[i];
  }
  in >> tet_num;
  tets.resize(tet_num * 4);
  for (int i = 0; i < tet_num * 4; ++i) {
    in >> tets[i];
    tets[i]--; // one indexed to zero indexed
  }
  in.close();
}

void NeutralFormatIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  out << verts.size() / 3 << "\n";
  for (int i = 0; i < int(verts.size()) / 3; ++i) {
    out << verts[i * 3 + 0] << " ";
    out << verts[i * 3 + 1] << " ";
    out << verts[i * 3 + 2] << "\n";
  }
  out << tets.size() / 4 << "\n";
  for (int t = 0; t < int(tets.size()) / 4; ++t) {
    out << tets[t * 4 + 0] + 1 << " ";
    out << tets[t * 4 + 1] + 1 << " ";
    out << tets[t * 4 + 2] + 1 << " ";
    out << tets[t * 4 + 4] + 1 << "\n";
  }
  out.close();
}
