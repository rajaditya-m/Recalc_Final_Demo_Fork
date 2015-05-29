#include <fstream>
#include <iostream>
#include "binary_tet_mesh_io.h"
#include "print_macro.h"
#include "vector_io.h"

void BinaryTetMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  std::ifstream in(file_name, std::ios::binary);
  ASSERT(in.is_open(), P(file_name));
  int v_num, tet_num;
  in.read((char*) &v_num, sizeof(int));
  in.read((char*) &tet_num, sizeof(int));
  verts.resize(v_num * 3);
  tets.reserve(tet_num * 4);
  in.read((char*) &verts[0], sizeof(double) * v_num * 3);
  in.read((char*) &tets[0], sizeof(double) * tet_num * 3);
  in.close();
}

void BinaryTetMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  std::ofstream out(file_name, std::ios::binary);
  ASSERT(out.is_open(), P(file_name));
  ASSERT(verts.size() % 3 == 0);
  ASSERT(tets.size() % 4 == 0);
  int v_num = int(verts.size()) / 3;
  int tet_num = int(tets.size()) / 4;
  out.write((char*) &v_num, sizeof(int));
  out.write((char*) &tet_num, sizeof(int));
  out.write((char*) &verts[0], sizeof(double) * v_num * 3);
  out.write((char*) &tets[0], sizeof(int) * tet_num * 4);
  out.close();
}
