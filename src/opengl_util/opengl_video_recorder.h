#ifndef OPENGLVIDEORECORDER_H
#define OPENGLVIDEORECORDER_H
#include <stdio.h>
/**
 * @brief require ffmepg available from command line
 */
class OpenGLVideoRecorder
{
public:
  OpenGLVideoRecorder(const char* mp4_file_name,
                      int width, int height,
                      int fps = 30,
                      int frame_rate = 2000);
  void RecordNextFrame();
  void FinishRecording();
  const char* mp4_file_name();
  ~OpenGLVideoRecorder();
private:
  int fps_;
  int frame_rate_;
  int width_, height_;
  unsigned char* frame_buffer_;
  FILE* ffmpeg_;
  char mp4_file_name_[1024];
};

#endif // OPENGLVIDEORECORDER_H
