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
#include "open_gl_qt.h"
#include "pose_viewer_simulator.h"

void Init() {
  global::gl->InstallHandler(PoseViwerSimulator::Instance());
}

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  w.hide();
    #ifdef WINDOW_TITLE
#if defined(_WIN32) || defined(_WIN64)
  QString title(STRINGIZE_TOKEN(WINDOW_TITLE));
#else
  QString title(WINDOW_TITLE);
#endif
  w.setWindowTitle(title);
#endif
  Init(); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}


