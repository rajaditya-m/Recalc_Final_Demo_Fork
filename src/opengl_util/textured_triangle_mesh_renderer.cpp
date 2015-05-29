#include "textured_triangle_mesh_renderer.h"
#include "opengl_helper.h"
#include "print_macro.h"

template <class Real>
const TexturedTriangleMeshRenderer<Real> &TexturedTriangleMeshRenderer<Real>::operator=(const TexturedTriangleMeshRenderer<Real> &other)
{
  this->shader_vert_pos_name_ = other.shader_vert_pos_name_;
  this->shader_vert_normal_name_ = other.shader_vert_normal_name_;
  this->shader_vert_tex_coord_name_ = other.shader_vert_tex_coord_name_;

  this->textured_initialized_ = other.textured_initialized_;
  this->texture_names_ = other.texture_names_;
  this->texture_ids_ = other.texture_ids_;
  this->texture_location_in_shader_ = other.texture_location_in_shader_;

  this->pos_attrib_id_ = other.pos_attrib_id_;
  this->normal_attrib_id_ = other.normal_attrib_id_;
  this->tex_attrib_id_ = other.tex_attrib_id_;

  this->vertex_pos_vbo_ = other.vertex_pos_vbo_;
  this->vertex_normal_vbo_ = other.vertex_normal_vbo_;
  this->vertex_tex_vbo_ = other.vertex_tex_vbo_;

  this->triangle_num_ = other.triangle_num_;
  this->shader_ = other.shader_;
  this->vertex_pos_ = other.vertex_pos_;
  this->vertex_tex_coord_ = other.vertex_tex_coord_;
  this->vertex_normal_ = other.vertex_normal_;
  this->own_buffer_ = true;
  other.own_buffer_ = false;
  return *this;
}

template <class Real>
TexturedTriangleMeshRenderer<Real>::TexturedTriangleMeshRenderer() {
  ASSERT(false, L("cannot be used"));
}

template <class Real>
TexturedTriangleMeshRenderer<Real>::TexturedTriangleMeshRenderer(const TexturedTriangleMeshRenderer<Real> &other) {
  this->operator=(other);
}

template <class Real>
TexturedTriangleMeshRenderer<Real>::TexturedTriangleMeshRenderer(const char* shader_file_name) {
  shader_ = SetupGLSL(shader_file_name);
  Construct(shader_);
}

template <class Real>
TexturedTriangleMeshRenderer<Real>::TexturedTriangleMeshRenderer(int shader) {
  Construct(shader);
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::Construct(int shader_id) {
  glGenBuffers(1, &vertex_pos_vbo_);
  glGenBuffers(1, &vertex_normal_vbo_);
  glGenBuffers(1, &vertex_tex_vbo_);
  shader_vert_pos_name_ = "position";
  shader_vert_normal_name_ = "normal";
  shader_vert_tex_coord_name_ = "tex";
  set_shader(shader_id);
  textured_initialized_ = false;
  own_buffer_ = true;
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::SetShaderSource(const char *shader_name_prefix) {
  shader_ = SetupGLSL(shader_name_prefix);
  set_shader(shader_);
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::set_shader(int shader) {
  this->shader_ = shader;
  glUseProgram(shader_);
  normal_attrib_id_ = glGetAttribLocation(shader_, shader_vert_normal_name_.c_str());
  pos_attrib_id_ = glGetAttribLocation(shader_, shader_vert_pos_name_.c_str());
  tex_attrib_id_ = glGetAttribLocation(shader_, shader_vert_tex_coord_name_.c_str());
  glUseProgram(0);
}

#include "global.h"
#include "profiler.h"
template <class Real>
void TexturedTriangleMeshRenderer<Real>::Render() {
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
    // texture coordinate
    if (textured_initialized_) {
      for (int i = 0; i < int(texture_names_.size()); ++i) {
        glActiveTexture(GL_TEXTURE0 + texture_ids_[i]);
        glBindTexture(GL_TEXTURE_2D, texture_ids_[i]);
        texture_location_in_shader_[i] =  glGetUniformLocation(shader_, texture_names_[i].c_str()); // get the location of the biased_MVP matrix
        glUniform1i(texture_location_in_shader_[i], texture_ids_[i]);
      }
      glEnableVertexAttribArray(tex_attrib_id_);
      glBindBuffer(GL_ARRAY_BUFFER, vertex_tex_vbo_);
      glVertexAttribPointer(tex_attrib_id_, 2, kPosType, GL_FALSE, sizeof(Real) * 2, (char*) NULL + 0);
    }
    // triangles
    glDrawArrays(GL_TRIANGLES, 0, triangle_num_ * 3);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    if (textured_initialized_) {
      for (int i = 0; i < int(texture_names_.size()); ++i) {
        glActiveTexture(GL_TEXTURE0 + texture_ids_[i]);
        glBindTexture(GL_TEXTURE_2D, 0);
      }
    }
  }
  glUseProgram(0);
#else
  glEnable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  for (int t = 0; t < triangle_num_; ++t) {
    Normal(&vertex_normal_[t * 9 + 0]);
    Vertex3v(&vertex_pos_[t * 9 + 0]);
    Normal(&vertex_normal_[t * 9 + 3]);
    Vertex3v(&vertex_pos_[t * 9 + 3]);
    Normal(&vertex_normal_[t * 9 + 6]);
    Vertex3v(&vertex_pos_[t * 9 + 6]);
  }
  glEnd();
#endif
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::SetTextures(std::vector<std::string> names, std::vector<unsigned int> texture_ids) {
  ASSERT(names.size() == texture_ids.size());
  this->texture_names_ = names;
  texture_ids_ = texture_ids;
  texture_location_in_shader_.resize(names.size());
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::UpdatePosAndNormal(Real *pos, Real *normal, int triangle_num) {
  this->triangle_num_ = triangle_num;
  vertex_pos_ = pos;
  glBindBuffer(GL_ARRAY_BUFFER, vertex_pos_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Real) * triangle_num_ * 3 * 3, pos, GL_DYNAMIC_DRAW);

  vertex_normal_ = normal;
  glBindBuffer(GL_ARRAY_BUFFER, vertex_normal_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Real) * triangle_num_ * 3 * 3, normal, GL_DYNAMIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::UpdateTexCoord(Real *tex_coord, int triangle_num) {
  textured_initialized_ = true;
  this->triangle_num_ = triangle_num;
  this->vertex_tex_coord_ = tex_coord;
  glBindBuffer(GL_ARRAY_BUFFER, vertex_tex_vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(Real)* triangle_num_ * 3 * 2, tex_coord, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

template <class Real>
int TexturedTriangleMeshRenderer<Real>::shader()
{
  return shader_;
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::set_shader_vert_pos_name(const char *name) {
  this->shader_vert_pos_name_ = std::string(name);
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::set_shader_vert_tex_coord_name(const char *name) {
  this->shader_vert_tex_coord_name_ = std::string(name);
}

template <class Real>
void TexturedTriangleMeshRenderer<Real>::set_shader_vert_normal_name(const char *name) {
  this->shader_vert_normal_name_ = std::string(name);
}

template <class Real>
TexturedTriangleMeshRenderer<Real>::~TexturedTriangleMeshRenderer() {
  if (own_buffer_) {
    glDeleteBuffers(1, &vertex_normal_vbo_);
    glDeleteBuffers(1, &vertex_pos_vbo_);
    glDeleteBuffers(1, &vertex_tex_vbo_);
  }
}

template class TexturedTriangleMeshRenderer<float>;
template class TexturedTriangleMeshRenderer<double>;
