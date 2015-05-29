#include "main_window.h"
#include <QDebug>
#include <QApplication>
#include "print_macro.h"
#include <iterator>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include "config_file.h"
#include "global.h"
#include "timer.h"
#include "profiler.h"
#include "config_file.h"
#include "sample_pose_simulator.h"
#include "open_gl_qt.h"
#include "vector_lib.h"

void Init(int argc, char *argv[]) {
  global::gl->InstallHandler(SamplePoseSimulator::Instance());
  SamplePoseSimulator::Instance()->Init(argc, argv);
}

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  w.hide();
  Init(argc, argv); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}


