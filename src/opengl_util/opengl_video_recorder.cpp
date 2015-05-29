#include <string.h>
#include "opengl_video_recorder.h"
#include "opengl_header.h"

#if defined(_WIN32) || defined(_WIN64)
#define POPEN _popen
#define PCLOSE _pclose
const char* kOpenOption = "wb";
#else
#define POPEN popen
#define PCLOSE pclose
const char* kOpenOption = "w";
#endif
//#include "print_macro.h"
/*
 referenced from:
 http://blog.mmacklin.com/2013/06/11/real-time-video-capture-with-ffmpeg/comment-page-1/#comment-50218
#include <stdio.h>

// start ffmpeg telling it to expect raw rgba 720p-60hz frames
// -i - tells it to read frames from stdin
const char* cmd = "ffmpeg -r 60 -f rawvideo -pix_fmt rgba -s 1280x720 -i - "
                  "-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip output.mp4";

// open pipe to ffmpeg's stdin in binary write mode
FILE* ffmpeg = _popen(cmd, "wb");

int* buffer = new int[width*height];
glutSwapBuffers();
glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer);

fwrite(buffer, sizeof(int)*width*height, 1, ffmpeg);
_pclose(ffmpeg);
 */

OpenGLVideoRecorder::OpenGLVideoRecorder(const char *mp4_file_name, int width, int height, int fps, int frame_rate) {
  width_ = width;
  height_ = height;
  fps_ = fps;
  frame_rate_ = frame_rate;
  strcpy(mp4_file_name_, mp4_file_name);
  char cmd[1024];
//  sprintf(cmd, "ffmpeg -r %d -f rawvideo -pix_fmt rgba -s %dx%d -i - -threads 0 -pix_fmt yuv420p -b:v %dk -vf vflip  -c:v libx264 -preset fast -y %s",
  sprintf(cmd, "ffmpeg -r %d -f rawvideo -pix_fmt rgba -s %dx%d -i - -threads 0 -b:v %dk -vf vflip  -c:v libx264 -preset fast -y %s",
          fps_,
          width_, height_,
          frame_rate_,
          mp4_file_name_);
  // P(std::string(cmd));
  //  const char* cmd = "ffmpeg -r 30 -f rawvideo -pix_fmt rgba -s 1280x724 -i - "
  //                    "-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip C:/tmp/out.mp4";
  ffmpeg_ = POPEN(cmd, kOpenOption);
  if (ffmpeg_ == NULL) {
    fprintf(stderr, "OpenGLVideoRecorder() => failed to open ffmpeg command");
    exit(0);
  }
  // handle retina display
#ifdef __APPLE__
  frame_buffer_ = new unsigned char[width_ * height_ * 4 * 4];
#else
  frame_buffer_ = new unsigned char[width_ * height_ * 4];
#endif
}

void OpenGLVideoRecorder::RecordNextFrame() {
  // handle retina display
#ifdef __APPLE__
  glReadPixels(0, 0, width_ * 2, height_ * 2, GL_RGBA, GL_UNSIGNED_BYTE, frame_buffer_);
  for (int r = 0; r < height_; ++r) {
    for (int c = 0; c < width_; ++c) {
      unsigned int pixels[4][4] = {
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 0) + 0],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 0) + 1],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 0) + 2],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 0) + 3],

        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 1) + 0],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 1) + 1],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 1) + 2],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 0) + 2 * c + 1) + 3],

        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 0) + 0],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 0) + 1],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 0) + 2],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 0) + 3],

        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 1) + 0],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 1) + 1],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 1) + 2],
        frame_buffer_[4 * ((width_ * 2) * (2 * r + 1) + 2 * c + 1) + 3],
      };
      for (int i = 0; i < 4; ++i) {
        unsigned char avg = (unsigned char) ((pixels[0][i] + pixels[1][i] + pixels[2][i] + pixels[3][i]) / 4);
        frame_buffer_[4 * (r * width_ + c) + i] = avg;
      }
    }
  }
#else
  glReadPixels(0, 0, width_, height_, GL_RGBA, GL_UNSIGNED_BYTE, frame_buffer_);
#endif
  fwrite(frame_buffer_, sizeof(unsigned char) * width_ * height_ * 4, 1, ffmpeg_);
}

OpenGLVideoRecorder::~OpenGLVideoRecorder() {
  FinishRecording();
}

void OpenGLVideoRecorder::FinishRecording() {
  PCLOSE(ffmpeg_);
  ffmpeg_ = NULL;
  delete []frame_buffer_;
}

const char *OpenGLVideoRecorder::mp4_file_name() {
  return mp4_file_name_;
}
