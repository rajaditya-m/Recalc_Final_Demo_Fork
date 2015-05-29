#include "vega_tet_mesh_io.h"
#include "config_file.h"
#include <vector>
#include <sstream>
#include <fstream>

void VegaTetMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  ConfigFile veg(file_name);
  std::vector<std::string> vert_str = *veg.Get<std::vector<std::string>* >("VERTICES");
  std::vector<std::string> tet_str = *veg.Get<std::vector<std::string>* >("ELEMENTS");
  int one_indexed = 0;
  {
    ASSERT(vert_str.size() > 0);
    std::stringstream in(vert_str[0]);
    int v_num;
    in >> v_num;
    verts.resize(v_num * 3);
    ASSERT((int) vert_str.size() >= v_num + 1);
    for (int v = 0; v < v_num; ++v) {
      std::stringstream in(vert_str[v + 1]);
      int junk;
      in >> junk;
      if (v == 0 && junk == 1) one_indexed = 1;
      in >> verts[v * 3 + 0];
      in >> verts[v * 3 + 1];
      in >> verts[v * 3 + 2];
    }
  }

  {
    ASSERT(tet_str.size() > 1);
    if (tet_str[0] != "TET") {
      std::cerr << CURRENT_LINE << " => only support TET files." << std::endl;
      exit(0);
    }
    std::stringstream in(tet_str[1]);
    int tet_num;
    in >> tet_num;
    ASSERT((int) tet_str.size() >= tet_num + 2);
    tets.resize(tet_num * 4);
    for (int t = 0; t < tet_num; ++t) {
      std::stringstream in(tet_str[t + 2]);
      int junk;
      in >> junk;
      in >> tets[t * 4 + 0];
      in >> tets[t * 4 + 1];
      in >> tets[t * 4 + 2];
      in >> tets[t * 4 + 3];
      if (one_indexed) {
        tets[t * 4 + 0]--;
        tets[t * 4 + 1]--;
        tets[t * 4 + 2]--;
        tets[t * 4 + 3]--;
      }
    }
  }
}

void VegaTetMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  ASSERT(verts.size() % 3 == 0);
  ASSERT(tets.size() % 4 == 0);
  char name_with_extension[512];
  sprintf(name_with_extension, "%s.veg", file_name);
  std::ofstream out(name_with_extension);
  ASSERT(out.is_open(), P(name_with_extension));
  out << "*VERTICES" << std::endl;
  out << verts.size() / 3 << " 3 0 0" << std::endl;
  for (int v = 0; v < int(verts.size()) / 3; ++v) {
    out << v + 1 << " "
        << verts[v * 3 + 0] << " "
        << verts[v * 3 + 1] << " "
        << verts[v * 3 + 2] << std::endl;
  }
  out << std::endl;
  out << "*ELEMENTS" << std::endl;
  out << "TET" << std::endl;
  out << tets.size() / 4 << " 4 0" << std::endl;
  for (int t = 0; t < int(tets.size()) / 4; ++t) {
    out << t + 1 << " "
        << tets[t * 4 + 0] + 1 << " "
        << tets[t * 4 + 1] + 1 << " "
        << tets[t * 4 + 2] + 1 << " "
        << tets[t * 4 + 3] + 1 << std::endl;
  }
  out << std::endl;
  out << "*MATERIAL defaultMaterial" << std::endl;
  out << "ENU, 1000, 1.5e5, 0.45" << std::endl;
  out << std::endl;
  out << "*REGION" << std::endl;
  out << "allElements, defaultMaterial" << std::endl;
  out.close();
}

