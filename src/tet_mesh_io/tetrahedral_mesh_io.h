#ifndef TETRAHEDRAL_MESH_IO_H
#define TETRAHEDRAL_MESH_IO_H
#include <string>
#include <vector>
#include "vector_lib.h"
class TetrahedralMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tets) = 0;
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tets) = 0;

  void Read(const char* file_name, std::vector<float>& verts, std::vector<int>& tets);
  void Write(const char* file_name, std::vector<float>& verts, std::vector<int>& tets);

  template <class T>
  void Read(std::string file_name, std::vector<T>& verts, std::vector<int>& tets) {
    Read(file_name.c_str(), verts, tets);
  }

  template <class T>
  void Write(std::string file_name, std::vector<T>& verts, std::vector<int>& tets) {
    Write(file_name.c_str(), verts, tets);
  }

  template <class T>
  void OrientTetrahedron(std::vector<T>& verts, std::vector<int>& tets)
  {
    for (int i = 0; i < (int) tets.size() / 4; ++i) {
      T* v[] = {
        &verts[tets[i * 4 + 0] * 3],
        &verts[tets[i * 4 + 1] * 3],
        &verts[tets[i * 4 + 2] * 3],
        &verts[tets[i * 4 + 3] * 3],
      };
      T vec[3][3];
      dj::SubVec3(v[1], v[0], vec[0]);
      dj::SubVec3(v[2], v[0], vec[1]);
      dj::SubVec3(v[3], v[0], vec[2]);
      T cross[3];
      dj::Cross3(vec[0], vec[1], cross);
      T dot = dj::Dot3(cross, vec[2]);
      if (dot < 0) {
        dj::Swap(tets[i * 4 + 3], tets[i * 4 + 2]);
      }
    }
  }

  template <class T>
  void Write(const char* file_name, T* verts, int vertex_num, int* tets, int tet_num)
  {
    std::vector<T> v(verts, verts + vertex_num * 3);
    std::vector<int> t(tets, tets + tet_num * 4);
    Write(file_name, v, t);
  }

protected:
  virtual ~TetrahedralMeshIO() {}
};

#endif // TETRAHEDRAL_MESH_IO_H
