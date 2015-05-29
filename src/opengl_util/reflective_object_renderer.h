#ifndef REFLECTIVEOBJECTRENDERER_H
#define REFLECTIVEOBJECTRENDERER_H
#include "enviroment_cube_map.h"
#include "textured_triangle_mesh_renderer.h"

class OpenGLMaterial;
template <class Real>
class ReflectiveObjectRenderer
  : public EnvironmentCubeMap<Real>,
    public TexturedTriangleMeshRenderer<Real> {
  using EnvironmentCubeMap<Real>::cube_map_handle_;
  using TexturedTriangleMeshRenderer<Real>::triangle_num_;
  using TexturedTriangleMeshRenderer<Real>::shader_;
  using TexturedTriangleMeshRenderer<Real>::shader_vert_pos_name_;
  using TexturedTriangleMeshRenderer<Real>::shader_vert_normal_name_;
  using TexturedTriangleMeshRenderer<Real>::shader_vert_tex_coord_name_;
  using TexturedTriangleMeshRenderer<Real>::textured_initialized_;
  using TexturedTriangleMeshRenderer<Real>::texture_names_;
  using TexturedTriangleMeshRenderer<Real>::texture_ids_;
  using TexturedTriangleMeshRenderer<Real>::texture_location_in_shader_;
  using TexturedTriangleMeshRenderer<Real>::pos_attrib_id_;
  using TexturedTriangleMeshRenderer<Real>::normal_attrib_id_;
  using TexturedTriangleMeshRenderer<Real>::tex_attrib_id_;
  using TexturedTriangleMeshRenderer<Real>::vertex_pos_vbo_;
  using TexturedTriangleMeshRenderer<Real>::vertex_normal_vbo_;
  using TexturedTriangleMeshRenderer<Real>::vertex_tex_vbo_;
public:
  ReflectiveObjectRenderer(const char* cube_map_file, const char* shader_source_file,
                           int material_type = -1);
  virtual void SetTextures(std::vector<std::string> names, std::vector<unsigned int> texture_ids);
  void Render();

  ~ReflectiveObjectRenderer();
protected:
  OpenGLMaterial* material_;
};

#endif // REFLECTIVEOBJECTRENDERER_H
