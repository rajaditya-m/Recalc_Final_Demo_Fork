#include "main_window.h"
//#include <armadillo>
#include <QDebug>
#include <QApplication>
#include "vector_lib.h"
#include "svd.h"
#include "MY_MATH.h"
#include "print_macro.h"
#include <iterator>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include "config_file.h"
#include "global.h"
#include "timer.h"
#include "open_gl_qt.h"
#include "profiler.h"
#include "VECTOR_TOOLS.h"
#include "tet_mesh_simulator.h"

void Init() {
  global::gl->InstallHandler(TetMeshSimulator::Instance());
}


int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  w.hide();
  Init(); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}

