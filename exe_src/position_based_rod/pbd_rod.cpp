#include <Eigen/Dense>
#include "main_window.h"
#include <QDebug>
#include <QApplication>
#include <iterator>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <sstream>

#include "global.h"
#include "vector_lib.h"
#include "print_macro.h"
#include "open_gl_qt.h"
#include "pbd_rod_simulator.h"


void Init() {
  global::gl->InstallHandler(PBDRodSimulator::Instance());
}

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
    #ifdef WINDOW_TITLE
#if defined(_WIN32) || defined(_WIN64)
  QString title(STRINGIZE_TOKEN(WINDOW_TITLE));
#else
  QString title(WINDOW_TITLE);
#endif
  w.setWindowTitle(title);
#endif
  Init();
  w.show();
  return a.exec();
}

