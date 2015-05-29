#pragma once
#include <vector>
class TriangularMeshIO
{
public:
  virtual void Read(const char* file_name, std::vector<double>& verts, std::vector<int>& tri) = 0;
  virtual void Write(const char* file_name, std::vector<double>& verts, std::vector<int>& tri) = 0;

  void Read(const char* file_name, std::vector<float>& verts, std::vector<int>& tri);
  void Write(const char* file_name, std::vector<float>& verts, std::vector<int>& tri);

  template <class T>
  void Write(const char* file_name, T* verts, int vertex_num, int* tri, int tri_num)
  {
    std::vector<T> v(verts, verts + vertex_num * 3);
    std::vector<int> t(tri, tri + tri_num * 3);
    Write(file_name, v, t);
  }

protected:
  TriangularMeshIO() {}
  virtual ~TriangularMeshIO() {}
};
