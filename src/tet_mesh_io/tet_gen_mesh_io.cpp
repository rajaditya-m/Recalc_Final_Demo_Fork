#include <iostream>
#include <fstream>
#include "text_file_io.h"
#include "tet_gen_mesh_io.h"
#include "print_macro.h"
#include "vector_io.h"

void TetGenMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  char filename[1024];
  sprintf(filename, "%s.node", file_name);
  FILE *fp = fopen(filename, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }
  ASSERT(verts.size() % 3 == 0);
  ASSERT(tets.size() % 4 == 0);
  int vertex_num = int(verts.size()) / 3;
  int tet_number = int(tets.size()) / 4;
  fprintf(fp, "%d %d %d %d\n", vertex_num, 3, 0, 0);
  for (int i = 0; i < vertex_num; i++)
    fprintf(fp, "%d %f %f %f\n", i + 1, verts[i * 3 + 0], verts[i * 3 + 1], verts[i * 3 + 2]);
  fclose(fp);

  sprintf(filename, "%s.ele", file_name);
  fp = fopen(filename, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }
  fprintf(fp, "%d %d %d\n", tet_number, 4, 0);

  for (int i = 0; i < tet_number; i++)
    fprintf(fp, "%d %d %d %d %d\n", i + 1, tets[i * 4 + 0] + 1, tets[i * 4 + 1] + 1, tets[i * 4 + 2] + 1, tets[i * 4 + 3] + 1);
  fclose(fp);
}

void TetGenMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tets)
{
  int temp_value;
  int bound;
  std::string prefix(file_name);
  int vertex_num;
  int one_indexed = 0;
  //------------------------------------------------------------------------------
  // node file
  {
    TextFileReader node_file((prefix + ".node").c_str());
    node_file >> vertex_num >> temp_value >> temp_value >> bound;
    verts.resize(vertex_num * 3);
    if (bound == 0) {
      for (int i = 0; i < vertex_num; ++i) {
        double tmp[3];
        node_file >> temp_value >> tmp[0] >> tmp[1] >> tmp[2];
        if (i == 0) one_indexed = (temp_value == 0) ? 0 : 1;
        verts[i * 3 + 0] = tmp[0];
        verts[i * 3 + 1] = tmp[1];
        verts[i * 3 + 2] = tmp[2];
      }
    } else {
      for (int i = 0; i < vertex_num; ++i) {
        double tmp[3];
        node_file >> temp_value >> tmp[0] >> tmp[1] >> tmp[2];
        if (i == 0) one_indexed = (temp_value == 0) ? 0 : 1;
        node_file >> temp_value;
        verts[i * 3 + 0] = tmp[0];
        verts[i * 3 + 1] = tmp[1];
        verts[i * 3 + 2] = tmp[2];
      }
    }
  }
  int tet_number;
  // tets
  {
    std::ifstream ele_file((prefix + ".ele").c_str());
    ASSERT(ele_file.is_open());
    ele_file >> tet_number >> temp_value >> bound;
    tets.resize(tet_number * 4);
    if (bound == 0) {
      for (int i = 0; i < tet_number; ++i) {
        ele_file >> temp_value >> tets[i * 4 + 0] >> tets[i * 4 + 1] >> tets[i * 4 + 2] >> tets[i * 4 + 3];
//        int junk;
//        ele_file >> junk;
      }
    } else {
      for (int i = 0; i < tet_number; ++i) {
        ele_file >> temp_value >> tets[i * 4 + 0] >> tets[i * 4 + 1] >> tets[i * 4 + 2] >> tets[i * 4 + 3] >> temp_value;
//        int junk;
//        ele_file >> junk;
      }
    }
    ele_file.close();
  }
  if (one_indexed) {
    for (int i = 0; i < tet_number; i++) {
      tets[i * 4 + 0] -= 1;
      tets[i * 4 + 1] -= 1;
      tets[i * 4 + 2] -= 1;
      tets[i * 4 + 3] -= 1;
    }
  }
}
