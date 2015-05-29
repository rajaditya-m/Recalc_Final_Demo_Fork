#include "reflective_object_renderer.h"
#include "opengl_header.h"
#include "opengl_material.h"
#include "print_macro.h"


template <class Real>
ReflectiveObjectRenderer<Real>::ReflectiveObjectRenderer(const char* cube_map_file,
                                                         const char* shader_source_file,
                                                         int material_type)
  : EnvironmentCubeMap<Real>(cube_map_file)
  , TexturedTriangleMeshRenderer<Real>(shader_source_file) {
  int loc = glGetAttribLocation(shader_, "has_texture");
  if (loc >= 0) glUniform1i(loc, 0);
  material_ = new OpenGLMaterial();
  if (material_type == -1) {
    *material_ = OpenGLMaterial::CreateMaterial(OpenGLMaterial::kJade);
  } else {
    *material_ = OpenGLMaterial::CreateMaterial(material_type);
  }
  material_->Use();
}

template <class Real>
void ReflectiveObjectRenderer<Real>::SetTextures(std::vector<std::string> names, std::vector<unsigned int> texture_ids) {
  int loc = glGetAttribLocation(shader_, "has_texture");
  glUniform1i(loc, 1);
  TexturedTriangleMeshRenderer<Real>::SetTextures(names, texture_ids);
}

template <class Real>
void ReflectiveObjectRenderer<Real>::Render() {
  glUseProgram(shader_);
  glActiveTexture(GL_TEXTURE0 + cube_map_handle_);
  glBindTexture(GL_TEXTURE_CUBE_MAP, cube_map_handle_);
  int uniloc = glGetUniformLocation(shader_, "cubeMap");
  if (uniloc >= 0)  glUniform1i(uniloc, cube_map_handle_);
  TexturedTriangleMeshRenderer<Real>::Render();
  glActiveTexture(GL_TEXTURE0 + cube_map_handle_);
  glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
  glUseProgram(0);
}

template <class Real>
ReflectiveObjectRenderer<Real>::~ReflectiveObjectRenderer() {
  delete material_;
}

template class ReflectiveObjectRenderer<float>;
template class ReflectiveObjectRenderer<double>;
