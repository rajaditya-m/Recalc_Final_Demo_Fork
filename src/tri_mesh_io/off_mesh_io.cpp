#include "off_mesh_io.h"
#include <fstream>
#include <sstream>
#include "print_macro.h"

bool StartWithPound(const std::string& str) {
  auto first_pound = str.find("#");
  if (first_pound != std::string::npos) {
    for (unsigned int i = 0; i < first_pound; ++i) {
      char c = str[i];
      if (c != ' ' && c != '\t') {
        return false;
      }
    }
    return true;
  } else {
    return false;
  }
}

void OffMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tri)
{
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));

  std::string line;
  while (true) {
    getline(in, line);
    if (line.size() == 0 || StartWithPound(line)) {
      continue;
    }
    if (line.find("OFF") != std::string::npos || line.find("off") != std::string::npos) {
      getline(in, line);
    }
    break;
  }

  int vertex_num, face_num;
  while (true) {
    if (line.size() > 0 && !StartWithPound(line)) {
      std::istringstream istring(line);
      istring >> vertex_num >> face_num;
      getline(in, line);
      break;
    }
    getline(in, line);
  }
  verts.resize(vertex_num * 3);
  tri.resize(face_num);

  for (int i = 0; i < vertex_num; ++i) {
    if (line.size() == 0 || StartWithPound(line)) {
      --i;
    } else {
      std::istringstream istring(line);
      istring >> verts[i * 3 + 0];
      istring >> verts[i * 3 + 1];
      istring >> verts[i * 3 + 2];
    }
    getline(in, line);
  }

  for (int i = 0; i < face_num; ++i) {
    if (line.size() == 0 || StartWithPound(line)) {
      --i;
    } else {
      std::istringstream istring(line);
      int vertex_count = 0;
      istring >> vertex_count;
      ASSERT(vertex_count == 3, L("Only support triangular mesh."));
      istring >> tri[i * 3 + 0];
      istring >> tri[i * 3 + 1];
      istring >> tri[i * 3 + 2];
    }
    getline(in, line);
  }
  in.close();
}

void OffMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tri)
{
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  ASSERT(verts.size() % 3 == 0);
  ASSERT(tri.size() % 3 == 0);
  out << "OFF\n";
  int v_num = int(verts.size() / 3);
  int tri_num = int(tri.size() / 3);
  out << v_num << " ";
  out << tri_num << " 0\n";

  for (int v = 0; v < v_num; ++v) {
    out << verts[v * 3 + 0] << " ";
    out << verts[v * 3 + 1] << " ";
    out << verts[v * 3 + 2] << "\n";
  }
  for (int t = 0; t < tri_num; ++t) {
    out << "3 ";
    out << tri[t * 3 + 0] << " ";
    out << tri[t * 3 + 1] << " ";
    out << tri[t * 3 + 2] << "\n";
  }
  out.close();
}
