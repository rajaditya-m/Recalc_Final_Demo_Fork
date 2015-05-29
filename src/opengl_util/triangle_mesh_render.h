#ifndef TRIANGLEMESHRENDER_H
#define TRIANGLEMESHRENDER_H
#include <vector>
#include <string>

template <class Real>
class TriangleMeshRender
{
  typedef unsigned int uint;
public:
  TriangleMeshRender() {}
  TriangleMeshRender(int vertex_num, Real* verts, Real* normal, int triangle_num, int *triangles, const char* shader_file_name);
  TriangleMeshRender(std::vector<Real>& verts, std::vector<Real>& normal, std::vector<int>& triangles, const char* shader_file_name);
  TriangleMeshRender(int vertex_num, int triangle_num, int *triangles, const char* shader_file_name, Real* colorData = NULL);
  TriangleMeshRender(int vertex_num, int triangle_num, int *triangles, int shader);
  void Construct(int vertex_num, int triangle_num, int* triangles, int shader_id, Real* colorData = NULL);
  void SetShaderSource(const char* shader_name_prefix);
  void set_shader(int shader);
  void UpdateTriangles(int triangle_num, int* triangles);
  void Render();
  void set_vertex_num(int vertex_num);
  void UpdateVertexVBO(Real* pos, Real* normal);
  void set_shader_vert_pos_name(const char* name);
  void set_shader_vert_normal_name(const char* name);
  virtual ~TriangleMeshRender();
private:
  std::string shader_vert_pos_name_;
  std::string shader_vert_normal_name_;
  std::string shader_vert_color_name_;
  int pos_attrib_id_;
  int normal_attrib_id_;
  int color_attrib_id_;

  uint vertex_pos_vbo_;
  uint vertex_normal_vbo_;
  uint tri_element_vbo_;
  uint vertex_color_vbo_;

  int triangle_num_;
  int vertex_num_;
  int* triangles_;
  uint shader_;
  Real* vertex_pos_;
  Real* vertex_normal_;
  Real* color_data_;
};

#endif // TRIANGLEMESHRENDER_H
