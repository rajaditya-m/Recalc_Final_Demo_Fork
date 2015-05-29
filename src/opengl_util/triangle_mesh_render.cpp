#include "triangle_mesh_render.h"
#include "opengl_helper.h"
#include "print_macro.h"

template <class Real>
TriangleMeshRender<Real>::TriangleMeshRender(int vertex_num, Real *verts, Real *normal,
                                             int triangle_num, int *triangles,
                                             const char *shader_file_name)
{
  shader_ = SetupGLSL(shader_file_name);
  Construct(vertex_num, triangle_num, triangles, shader_);
  UpdateVertexVBO(verts, normal);
}

template <class Real>
TriangleMeshRender<Real>::TriangleMeshRender(std::vector<Real> &verts, std::vector<Real> &normal, std::vector<int> &triangles,
                                             const char *shader_file_name)
{
  shader_ = SetupGLSL(shader_file_name);
  Construct(int(verts.size()) / 3, int(triangles.size()) / 3, &triangles[0], shader_);
  UpdateVertexVBO(&verts[0], &normal[0]);
}

template <class Real>
TriangleMeshRender<Real>::TriangleMeshRender(int vertex_num, int triangle_num,
                                             int *triangles,
                                             const char* shader_file_name, Real* colorData) {

  shader_ = SetupGLSL(shader_file_name);
  Construct(vertex_num, triangle_num, triangles, shader_,colorData);
}

template <class Real>
TriangleMeshRender<Real>::TriangleMeshRender(int vertex_num, int triangle_num, int *triangles, int shader) {
  Construct(vertex_num, triangle_num, triangles, shader);
}

template <class Real>
void TriangleMeshRender<Real>::Construct(int vertex_num, int triangle_num, int *triangles, int shader_id,Real* colorData) {
  this->vertex_num_ = vertex_num;
  this->triangle_num_ = triangle_num;
  this->triangles_ = triangles;
  color_data_ = colorData;
  glGenBuffers(1, &vertex_pos_vbo_);
  glGenBuffers(1, &vertex_normal_vbo_);
  glGenBuffers(1, &vertex_color_vbo_);
  glGenBuffers(1, &tri_element_vbo_);

  glBindBuffer(GL_ARRAY_BUFFER, vertex_color_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Real) * triangle_num_ * 3, color_data_, GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tri_element_vbo_);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * triangle_num_ * 3, &triangles[0], GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  shader_vert_pos_name_ = "position";
  shader_vert_normal_name_ = "normal";
  shader_vert_color_name_ = "color";
  set_shader(shader_id);
//  ASSERT(glGetError() == GL_NO_ERROR);
}

template <class Real>
void TriangleMeshRender<Real>::SetShaderSource(const char *shader_name_prefix) {
  shader_ = SetupGLSL(shader_name_prefix);
  set_shader(shader_);
}

template <class Real>
void TriangleMeshRender<Real>::set_shader(int shader) {
  this->shader_ = shader;
  glUseProgram(shader_);
  pos_attrib_id_ = glGetAttribLocation(shader_, shader_vert_pos_name_.c_str());
  normal_attrib_id_ = glGetAttribLocation(shader_, shader_vert_normal_name_.c_str());
  color_attrib_id_ = glGetAttribLocation(shader_, shader_vert_color_name_.c_str());
  glUseProgram(0);
}

template <class Real>
void TriangleMeshRender<Real>::UpdateTriangles(int triangle_num, int *triangles) {
  this->triangle_num_ = triangle_num;
  this->triangles_ = triangles;
  glDeleteBuffers(1, &tri_element_vbo_);
  glGenBuffers(1, &tri_element_vbo_);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tri_element_vbo_);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * triangle_num_ * 3, &triangles[0], GL_STATIC_DRAW);
//  ASSERT(glGetError() == GL_NO_ERROR);
}

template <class Real>
void TriangleMeshRender<Real>::Render() {
#if 1
  const GLenum kPosType = (typeid(Real) == typeid(float)) ? GL_FLOAT : GL_DOUBLE;
  glUseProgram(shader_);
  {
    // position
    glEnableVertexAttribArray(pos_attrib_id_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_pos_vbo_);
    glVertexAttribPointer(pos_attrib_id_, 3, kPosType, GL_FALSE, sizeof(Real) * 3, (char*) NULL + 0);
    // normal
    glEnableVertexAttribArray(normal_attrib_id_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_normal_vbo_);
    glVertexAttribPointer(normal_attrib_id_, 3, kPosType, GL_TRUE, sizeof(Real) * 3, (char*) NULL + 0);
    // colors
    glEnableVertexAttribArray(color_attrib_id_);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_color_vbo_);
    glVertexAttribPointer(color_attrib_id_, 3, kPosType, GL_FALSE, sizeof(Real) * 3, (char*) NULL + 0);

    // triangles
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tri_element_vbo_);
    glDrawElements(GL_TRIANGLES, triangle_num_ * 3, GL_UNSIGNED_INT, (char*) NULL + 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);


  }
  glUseProgram(0);
#else
  glEnable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  for (int t = 0; t < triangle_num_; ++t) {
    Normal(&vertex_normal_[triangles_[t * 3 + 0] * 3]);
    Vertex3v(&vertex_pos_[triangles_[t * 3 + 0] * 3]);
    Normal(&vertex_normal_[triangles_[t * 3 + 1] * 3]);
    Vertex3v(&vertex_pos_[triangles_[t * 3 + 1] * 3]);
    Normal(&vertex_normal_[triangles_[t * 3 + 2] * 3]);
    Vertex3v(&vertex_pos_[triangles_[t * 3 + 2] * 3]);
  }
  glEnd();
#endif
}

template <class Real>
void TriangleMeshRender<Real>::set_vertex_num(int vertex_num) {
  this->vertex_num_ = vertex_num;
}

template <class Real>
void TriangleMeshRender<Real>::UpdateVertexVBO(Real *pos, Real *normal) {
  vertex_pos_ = pos;
  vertex_normal_ = normal;
  glBindBuffer(GL_ARRAY_BUFFER, vertex_pos_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Real) * vertex_num_ * 3, pos, GL_DYNAMIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER, vertex_normal_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Real)* vertex_num_ * 3, normal, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

template <class Real>
void TriangleMeshRender<Real>::set_shader_vert_pos_name(const char *name) {
  this->shader_vert_pos_name_ = std::string(name);
}

template <class Real>
void TriangleMeshRender<Real>::set_shader_vert_normal_name(const char *name) {
  this->shader_vert_normal_name_ = std::string(name);
}

template <class Real>
TriangleMeshRender<Real>::~TriangleMeshRender() {
  glDeleteBuffers(1, &vertex_normal_vbo_);
  glDeleteBuffers(1, &vertex_pos_vbo_);
  glDeleteBuffers(1, &tri_element_vbo_);
}

template class TriangleMeshRender<float>;
template class TriangleMeshRender<double>;
