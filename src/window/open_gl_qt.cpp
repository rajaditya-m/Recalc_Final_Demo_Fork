//#include <GL/glew.h>
//#include <QGLFunctions>
#include <QMenu>
#include <QTimer>
#include <QDebug>
#include <QFileDialog>
#include <QString>
#include <cmath>
#include <fstream>
#include <QByteArray>
#include <QIODevice>
#include <QFont>

#include <QLabel>
#include <QByteArray>
#include <QImageWriter>
#include <QString>
#include <QDebug>
#include <QMouseEvent>
#include <sstream>

#include "opengl_helper.h"
#include "opengl_video_recorder.h"
#include "camera.h"
#include "open_gl_qt.h"
#include "input_handler.h"
#include "opengl_helper.h"
#include "qt_gl_3d_scene_navigator.h"
#include "scene.h"
#include "rainbow_color.h"
#include "fps_calculator.h"
#include "global.h"
#include "config_file.h"
#include "print_macro.h"
#include "fps_calculator.h"
#include "timer.h"
#include "string_formatter.h"
#include "opengl_soft_shadow_render.h"
#include "triangle_mesh_render.h"
#include "obj_mesh_io.h"
#include "../demo/hammock_support.h"

#ifndef CURRENT_DIRECTORY
#define CURRENT_DIRECTORY STRINGIZE_TOKEN(DATA_DIRECTORY)
#endif

HammockSupport* hammock_support = nullptr;
OpenGLSoftShadowRender* soft_shadow = nullptr;
TriangleMeshRender<double>* tri_render = nullptr;
FPSCalculator fps_calculator;
const std::string kSceneFileName(std::string(global::kPWD) + "../conf/camera.conf");

unsigned int idle_count = 0;

enum MouseState {
  kDefault,
  kZoom,
  kRotate	,
  kTranslate,
};


OpenGLQt::OpenGLQt(QWidget *parent)
  : QGLWidget(parent)
  , camera_(new dj::Camera(kSceneFileName.c_str())) {
  global::gl = this;
  QtGL3DSceneNavigator* handler = Singleton<QtGL3DSceneNavigator>::Instance();
  video_recorder_ = nullptr;
  handler->set_gl_widget(this);
  handler->set_camera(camera());
  InstallHandler(handler);

  timer_ = new QTimer(this);
  connect(timer_, SIGNAL(timeout()), this, SLOT(idle()));
  timer_->start(1);
}

OpenGLQt::~OpenGLQt() {
  delete timer_;
  if (video_recorder_ != nullptr) {
    delete video_recorder_;
  }
}

void OpenGLQt::ScreenShot(const char *file_name, const char *image_type) {
  QImage* image;
  QImageWriter* image_writer;
  // handle retina display
#ifdef __APPLE__
  image = new QImage(width() * 2, height() * 2, QImage::Format_RGB888);
  glReadPixels(0, 0, width() * 2, height() * 2, GL_RGB, GL_UNSIGNED_BYTE, image->bits());
  image->scaled(width(), height());
#else
  image = new QImage(width() , height(), QImage::Format_RGB888);
  glReadPixels(0, 0, width(), height(), GL_RGB, GL_UNSIGNED_BYTE, image->bits());
#endif
  image_writer = new QImageWriter(file_name, image_type);
  image_writer->setQuality(100);
  image_writer->write(image->mirrored(false, true));
  P(file_name, image_type);
  delete image;
  delete image_writer;
}

void OpenGLQt::saveScene() {
  camera_->Save();
  emit updateStatusMessage(QString("Scene file save"));
}

#include "triangular_mesh_io.h"
void OpenGLQt::initializeGL(void) {
  //#ifdef USE_FREE_GLUT
  //#endif
#ifndef __APPLE__
  int arg = 0;
  glutInit(&arg, NULL);
#endif
  //  makeCurrent();
  ASSERT(glewInit() == GLEW_OK);
  glShadeModel(GL_SMOOTH);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  //  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClearDepth(1.0f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glEnable(GL_NORMALIZE);
  P(global::time_step);
  auto dummy = []() {
  };
  float light_pos[3] = {-1.2f, 2.0f, -0.5f};
  ASSERT(glGetError() == GL_NO_ERROR);
#ifndef __APPLE__
  soft_shadow = new OpenGLSoftShadowRender(width(), height(), light_pos, dummy);
#endif
  if (conf.Get<int>("draw hammock support")) {
    hammock_support = new HammockSupport;
  }
}

void OpenGLQt::resizeGL(int w, int h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // (Field of view, aspect ratio, near plance, far plane)
  gluPerspective(camera_->field_of_view(),  float(w) / h, camera_->near_plane(), camera_->far_plane());
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glEnable(GL_DEPTH_TEST);
  updateGL();
}

void OpenGLQt::paintGL() {
  //    glClearColor(1, 1, 1, 1);
  //  soft_shadow->set_light_position(camera()->eye_pos());
  //  soft_shadow->CreateShadowMap();
  if (soft_shadow) {
    soft_shadow->BeginShadowMap();
    for (InputHandler * handler : input_handlers_) {
      handler->Render();
    }
    if (hammock_support) {
      hammock_support->Render();
    }
    soft_shadow->EndShadowMap();
  }
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //    DrawGradientBackGround();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // (Field of view, aspect ratio, near plance, far plane)
  gluPerspective(camera_->field_of_view(),  float(width()) / height(), camera_->near_plane(), camera_->far_plane());
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (1) {
    float light_diffuse[] = {1.0f, 1.0f, 1.0, 0};
    float light_specular[] = {1.0f, 1.0f, 1.0, 0};
    float light_position[] = {0, 0, 0, 1}; // point light at camera position
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    //    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, );
    glEnable(GL_LIGHT0);
  }
  // (eye_pos, focus_point, view_up)
  gluLookAt(camera_->eye_pos()[0], camera_->eye_pos()[1], camera_->eye_pos()[2],
            camera_->focus()[0], camera_->focus()[1], camera_->focus()[2],
            camera_->up()[0], camera_->up()[1], camera_->up()[2]);

  if (soft_shadow) {
    soft_shadow->RenderSoftShadow();
  } else {
    const float kHeight = -0.41f;
    glPushMatrix();
    glTranslatef(0, kHeight, 0);
    glRotatef(-90, 1, 0, 0);
        DrawCheckBoard(10.0f, 10.0f, 20, 20, kBlue(), kWhite());
    glPopMatrix();
  }
  if (0) {
    float light_diffuse[] = {1.0f, 1.0f, 1.0, 1};
    float light_position[] = {0, 1, 0.00010f, 1};
    float light_specular[] = {1.0f, 1.0f, 1.0, 0};
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position);
    glPushMatrix();
    glTranslatef(light_position[0], light_position[1], light_position[2]);
    DrawSphere(0.08f);
    glPopMatrix();
    glEnable(GL_LIGHT1);
  }
  //  tri_render->Render();
  //  DrawSphere(0.8, 20, 20);
  for (InputHandler * handler : input_handlers_) {
    handler->Render();
  }

  if (hammock_support) {
    hammock_support->Render();
  }
#if 0
  glUseProgram(shader_);
  glBegin(GL_TRIANGLES);
  glNormal3f(0, 0, 1);
  glVertex3f(0, 0, 0);
  glVertex3f(0, 1, 0);
  glVertex3f(1, 0, 0);
  glEnd();
  glUseProgram(0);
#endif
  //  return;
  //------------------------------------------------------------------------------
  // Render texts on screen
#if 1
  QFont sansFont("Helvetica [Cronyx]", 14);
  sansFont.setFamily("sans serif");
  int font_height = QFontMetrics(sansFont).height();
  glColor3fv(kBlack());
  int base = height() - 5;
  QString status;
   switch(global::sim_state)
   {
     case 1: status = QString("Split Mode"); break;
     case 2: status = QString("Recomputing Basis..."); break;
     case 3: status = QString("Join Mode"); break;
     case 4: status = QString("Recomputing Basis..."); break;
     case 5: status = QString("Full Stitched Mode"); break;
   }

   QString texts[] = {
     QString("FPS: %1").arg(fps_calculator.fps()),
     QString("Status: ") + status,
     QString("Stiffness Multiplier: %1x").arg(global::stiffyMult),
   };
  for (unsigned i = 0; i < sizeof(texts) / sizeof(QString); ++i) {
    renderText(5, base - font_height * i, texts[i], sansFont);
  }
#endif
}

void OpenGLQt::keyPressEvent(QKeyEvent *event) {
  int priority = -1;
  for (InputHandler * handler : input_handlers_) {
    if (handler->CanHandle(priority)) {
      priority = handler->HandleKeyPress(event);
    }
  }
  updateGL();
}

void OpenGLQt::mousePressEvent(QMouseEvent *e) {
  QMouseEvent event = adjustPosForRetinaDisplay(e);
  int priority = -1;
  for (InputHandler * handler : input_handlers_) {
    if (handler->CanHandle(priority)) {
      priority = handler->HandleMousePress(&event);
    }
  }
}

void OpenGLQt::mouseReleaseEvent(QMouseEvent *e) {
  QMouseEvent event = adjustPosForRetinaDisplay(e);
  int priority = -1;
  for (InputHandler * handler : input_handlers_) {
    if (handler->CanHandle(priority)) {
      priority = handler->HandleMouseRelease(&event);
    }
  }
}

void OpenGLQt::mouseMoveEvent(QMouseEvent *e) {
  QMouseEvent event = adjustPosForRetinaDisplay(e);
  int priority = -1;
  for (InputHandler * handler : input_handlers_) {
    if (handler->CanHandle(priority)) {
      priority = handler->HandleMouseMove(&event);
    }
  }
  //  QPoint pos = event.pos();
  //  emit updateStatusMessage(QString("Cursor position: [%1, %2]").arg(pos.x()).arg(pos.y()));
  //  updateGL();
}

void OpenGLQt::wheelEvent(QWheelEvent *e) {
  if (e->modifiers() != Qt::ShiftModifier) return;
  if (e->delta() > 0) {
    camera_->Zoom(-camera_->zoom_step());
  } else {
    camera_->Zoom(+camera_->zoom_step());
  }
  updateGL();
}

QMouseEvent OpenGLQt::adjustPosForRetinaDisplay(QMouseEvent *e) {
  QPoint pos = e->pos();
  GLint view_port[4]; // viewport dimensions+pos
  glGetIntegerv(GL_VIEWPORT, view_port);
  pos *= view_port[2] / width();
  QMouseEvent tmp_event(e->type(), pos, e->button(), e->buttons(), e->modifiers());
  return tmp_event;
}

bool OpenGLQt::InstallHandler(InputHandler *handler) {
  for (int i = 0; i < (int) input_handlers_.size(); ++i) {
    if (input_handlers_[i] == handler) {
      return false;
    }
  }
  input_handlers_.push_back(handler);
  auto Comparator = [](const InputHandler * a, const InputHandler * b) -> bool {
    return (*a) > (*b);
  };
  std::sort(input_handlers_.begin(), input_handlers_.end(), Comparator);
  return true;
}

bool OpenGLQt::UninstallHandler(InputHandler *handler) {
  auto iter = std::find(input_handlers_.begin(), input_handlers_.end(), handler);
  if (iter != input_handlers_.end()) {
    input_handlers_.erase(iter);
    return true;
  } else {
    return false;
  }
}

void OpenGLQt::EnableVideoRecording(const char *mp4_file_prefix, int fps, int frame_rate) {
#if defined(_WIN32) || defined(_WIN64)
  const char* default_video_name = "D:/Dropbox/2.video/opengl_video.mp4";
  //  const char* default_video_name = "C:/tmp/opengl_video.mp4";
#else
  const char default_video_name[] = "~/Dropbox/2.video/opengl_video.mp4";
#endif
  if (mp4_file_prefix == nullptr) mp4_file_prefix = default_video_name;
  video_recorder_ = new OpenGLVideoRecorder(mp4_file_prefix, width(), height(), fps, frame_rate);
}

void OpenGLQt::RecordNextFrame() {
  ASSERT(video_recorder_ != nullptr, L("video recorded is not enabled. Call EnableVideoRecording() first"));
  updateGL();
  video_recorder_->RecordNextFrame();
}

void OpenGLQt::FinishRecordingVideo() {
  video_recorder_->FinishRecording();
  std::string msg = dj::Format("video is recorded to %z", video_recorder_->mp4_file_name());
  L(msg);
  delete video_recorder_;
  video_recorder_ = nullptr;
}

void OpenGLQt::idle() {
  if (global::simulate) {
    for (InputHandler * handler : input_handlers_) {
      handler->Idle();
    }
    if (global::pause_per_frame) {
      global::simulate = false;
    }
  }
  fps_calculator.CalculateFPS();
  updateGL();
}

void OpenGLQt::keyReleaseEvent(QKeyEvent *event) {
  for (InputHandler * handler : input_handlers_) {
    handler->HandleKeyRelease(event);
  }
}
