#ifndef OPENGLSOFTSHADOWRENDER_H
#define OPENGLSOFTSHADOWRENDER_H
#include <functional>
//#include "opengl_header.h"


class OpenGLSoftShadowRender {
  typedef unsigned int uint;
public:
  OpenGLSoftShadowRender(int screen_with, int screen_height,
                         const float* light_pos_in_world,
                         std::function<void (void)> render_ocluding_objects_);
  void set_light_position(const float *light_pos);
  void ResizeScreen(int width, int height);
  void BeginShadowMap();
  void EndShadowMap();
  void CreateShadowMap();
  void RenderSoftShadow();
  ~OpenGLSoftShadowRender();
private:
  void InitTexture();
  int screen_width_;
  int screen_height_;
  uint texture_id_;
  uint depth_fbo_;
  uint depth_texture_;
  uint shadow_program_;
  float light_position_[3];
  std::function<void (void)> render_ocluding_objects_;
};

#endif // OPENGLSOFTSHADOWRENDER_H
