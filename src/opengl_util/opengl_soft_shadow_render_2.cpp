#include <iostream>
#include <QImageReader>
#include "opengl_soft_shadow_render.h"
#include "opengl_helper.h"
#include "opengl_header.h"
#include "global.h"
#include "print_macro.h"


OpenGLSoftShadowRender::OpenGLSoftShadowRender(int screen_with, int screen_height,
                                               const float *light_pos,
                                               std::function<void (void)> render_ocluding_objects) {
  this->screen_width_ = screen_with;
  this->screen_height_ = screen_height;
  light_position_[0] = light_pos[0];
  light_position_[1] = light_pos[1];
  light_position_[2] = light_pos[2];
  shadow_program_	= SetupGLSL(DATA_DIRECTORY "../shader/shadow");
  InitTexture();
  //Init depth texture and FBO
  glGenFramebuffers(1, &depth_fbo_);
  glBindFramebuffer(GL_FRAMEBUFFER, depth_fbo_);
  // Depth texture. Slower than a depth buffer, but you can sample it later in your shader
  glGenTextures(1, &depth_texture_);
  glActiveTexture(GL_TEXTURE0 + depth_texture_);
  glBindTexture(GL_TEXTURE_2D, depth_texture_);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, screen_width_, screen_height_, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
  //    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_texture_, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_texture_, 0);
  glUseProgram(shadow_program_);
  GLuint uniloc = glGetUniformLocation(shadow_program_, "shadow_texture");
  glUniform1i(uniloc, depth_texture_);
  glUseProgram(0);

  glDrawBuffer(GL_NONE);
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    std::cerr << "OpenGLSoftShadowRender() => create shadow map buffer failed." << std::endl;
    exit(0);
  }

  this->render_ocluding_objects_ = render_ocluding_objects;
}

void OpenGLSoftShadowRender::set_light_position(const float* light_pos) {
  light_position_[0] = light_pos[0];
  light_position_[1] = light_pos[1];
  light_position_[2] = light_pos[2];
}

void OpenGLSoftShadowRender::ResizeScreen(int width, int height) {
  screen_width_ = width;
  screen_height_ = height;
}

void OpenGLSoftShadowRender::BeginShadowMap() {
  glBindFramebuffer(GL_FRAMEBUFFER, depth_fbo_);
  // Render on the whole framebuffer, complete from the lower left corner to the upper right
  glViewport(0, 0, screen_width_, screen_height_);
  // glEnable(GL_CULL_FACE);
  // glCullFace(GL_BACK); // Cull back-facing triangles -> draw only front-facing triangles
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-2, 2, -2, 2, 0, 20);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(light_position_[0], light_position_[1], light_position_[2], 0, 0, 0 , 0, 1, 0);
  //Use fixed program
  glUseProgram(0);

}

#include <QImage>
#include <QImageWriter>
void OpenGLSoftShadowRender::EndShadowMap() {
  if (0) {
    QImage* image;
    QImageWriter* image_writer;
    image = new QImage(screen_width_, screen_height_, QImage::Format_RGB888);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, image->bits());
    //    image_writer = new QImageWriter("/tmp/depth.png", "png");
    image_writer = new QImageWriter("C:/tmp/depth.png", "png");
    image_writer->setQuality(100);
    image_writer->write(image->mirrored(false, true));
    delete image;
    delete image_writer;
    exit(0);
  }
  float bias[16] = {0.5, 0.0, 0.0, 0.0,
                    0.0, 0.5, 0.0, 0.0,
                    0.0, 0.0, 0.5, 0.0,
                    0.5, 0.5, 0.5, 1.0
                   };

  // Grab modelview and transformation matrices
  float	modelView[16];
  float	projection[16];
  float	biased_MVP[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
  glGetFloatv(GL_PROJECTION_MATRIX, projection);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glLoadMatrixf(bias);
  // concatating all matrice into one.
  glMultMatrixf(projection);
  glMultMatrixf(modelView);

  glGetFloatv(GL_MODELVIEW_MATRIX, biased_MVP);

  glUseProgram(shadow_program_);
  GLuint m = glGetUniformLocation(shadow_program_, "biased_MVP"); // get the location of the biased_MVP matrix
  glUniformMatrix4fv(m, 1, GL_FALSE, biased_MVP);

  glUseProgram(0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
}

void OpenGLSoftShadowRender::CreateShadowMap() {
  glBindFramebuffer(GL_FRAMEBUFFER, depth_fbo_);
  // Render on the whole framebuffer, complete from the lower left corner to the upper right
  glViewport(0, 0, screen_width_, screen_height_);
  // glEnable(GL_CULL_FACE);
  // glCullFace(GL_BACK); // Cull back-facing triangles -> draw only front-facing triangles
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-2, 2, -2, 2, 0, 20);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(light_position_[0], light_position_[1], light_position_[2], 0, 0, 0 , 0, 1, 0);
  //Use fixed program
  glUseProgram(0);
  DrawSphere(0.8, 20, 20);
  //  render_ocluding_objects_();

  /*	if(filename)
  {
    float* pixels=new float[1024*1024*3];
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pixels);

    for(int i=0; i<1024; i++)
    for(int j=0; j<1024; j++)
      if(pixels[(i*1024+j)*3+0])
        printf("value i %d, %d: %f\n", i, j, pixels[(i*1024+j)*3+0]);

    BMP_Write(filename, pixels, 1024, 1024);
    delete[] pixels;
    getchar();
  }*/


  //Also we need to set up the projection matrix for shadow texture
  // This is matrix transform every coordinate x,y,z
  // Moving from unit cube [-1,1] to [0,1]
  float bias[16] = {0.5, 0.0, 0.0, 0.0,
                    0.0, 0.5, 0.0, 0.0,
                    0.0, 0.0, 0.5, 0.0,
                    0.5, 0.5, 0.5, 1.0
                   };

  // Grab modelview and transformation matrices
  float	modelView[16];
  float	projection[16];
  float	biased_MVP[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
  glGetFloatv(GL_PROJECTION_MATRIX, projection);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glLoadMatrixf(bias);
  // concatating all matrice into one.
  glMultMatrixf(projection);
  glMultMatrixf(modelView);

  glGetFloatv(GL_MODELVIEW_MATRIX, biased_MVP);

  glUseProgram(shadow_program_);
  GLuint m = glGetUniformLocation(shadow_program_, "biased_MVP"); // get the location of the biased_MVP matrix
  glUniformMatrix4fv(m, 1, GL_FALSE, biased_MVP);

  glUseProgram(0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void OpenGLSoftShadowRender::RenderSoftShadow() {
  glDisable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glUseProgram(shadow_program_);
  GLuint uniloc = glGetUniformLocation(shadow_program_, "light_position");
  glUniform3fv(uniloc, 1, light_position_);

  //    glUseProgram(0);
  const float kHeight = -0.6f;
  const float kTextureScale = 1.0f;
  //  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0 + texture_id_);
  glBindTexture(GL_TEXTURE_2D, texture_id_);
  glActiveTexture(GL_TEXTURE0 + depth_texture_);
  glBindTexture(GL_TEXTURE_2D, depth_texture_);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(-5, kHeight, -5);
  glTexCoord2f(kTextureScale, 0.0);
  glVertex3f(+5, kHeight, -5);
  glTexCoord2f(kTextureScale, kTextureScale);
  glVertex3f(+5, kHeight, +5);
  glTexCoord2f(0.0, kTextureScale);
  glVertex3f(-5, kHeight, +5);
  glEnd();
  //  glDisable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0 + texture_id_);
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0 + depth_texture_);
  glBindTexture(GL_TEXTURE_2D, 0);
  glUseProgram(0);
}

void OpenGLSoftShadowRender::InitTexture() {
  QImageReader reader;
  //  QString qName(DATA_DIRECTORY "texture/floor.jpg");
  //  QString qName(DATA_DIRECTORY "texture/tile.jpg");
  QString qName(DATA_DIRECTORY "texture/checkerIm.png");
  //  QString qName(DATA_DIRECTORY "texture/basketball-floor-textureqq.jpg");
  //  QString qName(DATA_DIRECTORY "texture/irregular_stone_floor.jpg");
  //  QString qName(DATA_DIRECTORY "texture/1.jpg");
  //  QString qName(DATA_DIRECTORY "texture/apple.png");
  reader.setFileName(qName);
  ASSERT(reader.canRead());
  QImage img = reader.read().convertToFormat(QImage::Format_RGBA8888).mirrored(false, true);

  texture_id_ = CreateTexture(img.bits(), img.width(), img.height());

  glUseProgram(shadow_program_);
  int uniloc = glGetUniformLocation(shadow_program_, "texture");
  //  ASSERT(uniloc >= 0, P(uniloc));
  glUniform1i(uniloc, texture_id_);
  glUseProgram(0);
}
