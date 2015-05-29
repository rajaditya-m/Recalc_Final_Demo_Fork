#include "obj_mesh_io.h"
#include <fstream>
#include "print_macro.h"

void ObjMeshIO::Read(const char *file_name, std::vector<double> &verts, std::vector<int> &tri)
{
  /// borrowed from https://github.com/eestrada/wavefront-obj-reader/blob/master/objread.hpp
  std::ifstream objfile(file_name);
  ASSERT(objfile.is_open(), P(file_name));

  std::string line, tmp;
  while (!objfile.eof()) {
    getline(objfile, line);
    if (line.empty()) continue; //if line is empty, restart loop

    std::istringstream strm(line);
    strm >> tmp;
    if (tmp == "v") {
      double x, y, z;
      strm >> x >> y >> z;
      verts.push_back(x);
      verts.push_back(y);
      verts.push_back(z);
    } else if (tmp == "vt") {
      // ignore for now
    } else if (tmp == "vn") {
      // ignore for now
    } else if (tmp == "f") {
      // All faces are assumed to be triangles
      // I may change this later to be more flexible.
      std::string v0, v1, v2;
      int pi, ti;
      strm >> v0 >> v1 >> v2;

      // Obj indices are not 0-based, so compensate by decrementing the value.
      //std::sscanf( v0.c_str(), "%d/%d", &(t.pos[0]), &(t.tex[0]) );
      std::sscanf( v0.c_str(), "%d/%d", &pi, &ti );
      --pi;
      --ti;
      tri.push_back(pi);

      std::sscanf( v1.c_str(), "%d/%d", &pi, &ti );
      --pi;
      --ti;
      tri.push_back(pi);

      std::sscanf( v2.c_str(), "%d/%d", &pi, &ti );
      --pi;
      --ti;
      tri.push_back(pi);
    }
  }
  objfile.close();
}

void ObjMeshIO::Write(const char *file_name, std::vector<double> &verts, std::vector<int> &tri)
{
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  ASSERT(verts.size() % 3 == 0);
  ASSERT(tri.size() % 3 == 0);
  int v_num = int(verts.size() / 3);
  int tri_num = int(tri.size() / 3);
  for (int v = 0; v < v_num; ++v) {
    out << "v "
        << verts[v * 3 + 0] << " "
        << verts[v * 3 + 1] << " "
        << verts[v * 3 + 2] << "\n";
  }
  for (int t = 0; t < tri_num; ++t) {
    out << "f "
        << tri[t * 3 + 0] + 1 << " "
        << tri[t * 3 + 1] + 1 << " "
        << tri[t * 3 + 2] + 1 << "\n";
  }
  out.close();
}
