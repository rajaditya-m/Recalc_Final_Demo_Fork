#include "subspace_simulator.h"
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

void Test() {
  for(int i = 0; i < 10; ++i) {
    std::cout << "\r[";
    for (int j = 0; j < i + 1; ++j) {
      std::cout << "=";
    }
    for (int j = i + 1; j < 10; ++j) {
      std::cout << " ";
    }
    std::cout << "]" << std::flush;
  }
  exit(0);
}

void Init() {
//  Test();
  global::gl->InstallHandler(SubspaceSimulator::Instance());
}


int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.setWindowTitle(QString("Subspace Simulator"));
  w.show();
  w.hide();
  Init(); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}

