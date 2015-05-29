#ifndef OPENGL_QT_H
#define OPENGL_QT_H
#include "opengl_helper.h"
#include <QGLWidget>
#include <vector>
#include <memory>

namespace dj {
class Camera;
}
class InputHandler;
class OpenGLVideoRecorder;

class OpenGLQt : public QGLWidget {
  Q_OBJECT
public:
  explicit OpenGLQt(QWidget *parent = 0);
  ~OpenGLQt();

  dj::Camera* camera() { return camera_.get(); }
  void ScreenShot(const char* file_name, const char* image_type);
  QMouseEvent adjustPosForRetinaDisplay(QMouseEvent* e);
  bool InstallHandler(InputHandler* handler);
  bool UninstallHandler(InputHandler* handler);

  void EnableVideoRecording(const char* mp4_file_prefix = nullptr, int fps = 30, int frame_rate = 2000);
  void RecordNextFrame();
  void FinishRecordingVideo();

private:
  void keyReleaseEvent(QKeyEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *e);
  void keyPressEvent(QKeyEvent *event);

signals:
  void updateStatusMessage(QString);
  void wheelUp();
  void wheelDown();

public slots:
  void idle();
  void saveScene();

protected:
  void initializeGL();
  void resizeGL(int w, int h);
  void paintGL();

private:
  OpenGLVideoRecorder* video_recorder_;
  std::auto_ptr<dj::Camera> camera_;
  std::vector<InputHandler*> input_handlers_;
  QTimer* timer_;
};

#endif
