#ifndef TEXTUREDTRIANGLEMESHRENDERER_H
#define TEXTUREDTRIANGLEMESHRENDERER_H
#include <vector>
#include <string>

template <class Real>
class TexturedTriangleMeshRenderer
{
  typedef unsigned int uint;
public:
  TexturedTriangleMeshRenderer();
  TexturedTriangleMeshRenderer(const TexturedTriangleMeshRenderer<Real> &other);
  TexturedTriangleMeshRenderer(const char* shader_file_name);
  TexturedTriangleMeshRenderer(int shader);
  const TexturedTriangleMeshRenderer<Real>& operator=(const TexturedTriangleMeshRenderer& other);
  void Construct(int shader_id);
  void SetShaderSource(const char* shader_name_prefix);
  void set_shader(int shader);
  void Render();
  virtual void SetTextures(std::vector<std::string> names, std::vector<unsigned int> texture_ids);
  void UpdatePosAndNormal(Real* pos, Real* normal, int triangle_num);
  void UpdateTexCoord(Real* tex_coord, int triangle_num);
  int shader();
  void set_shader_vert_pos_name(const char* name);
  void set_shader_vert_tex_coord_name(const char* name);
  void set_shader_vert_normal_name(const char* name);
  virtual ~TexturedTriangleMeshRenderer();
protected:
  mutable bool own_buffer_;
  std::string shader_vert_pos_name_;
  std::string shader_vert_normal_name_;
  std::string shader_vert_tex_coord_name_;

  bool textured_initialized_;
  std::vector<std::string> texture_names_;
  std::vector<unsigned int> texture_ids_;
  std::vector<int> texture_location_in_shader_;

  int pos_attrib_id_;
  int normal_attrib_id_;
  int tex_attrib_id_;

  uint vertex_pos_vbo_;
  uint vertex_normal_vbo_;
  uint vertex_tex_vbo_;

  int triangle_num_;
  uint shader_;
  Real* vertex_pos_;
  Real* vertex_tex_coord_;
  Real* vertex_normal_;
};

#endif // TEXTUREDTRIANGLEMESHRENDERER_H
