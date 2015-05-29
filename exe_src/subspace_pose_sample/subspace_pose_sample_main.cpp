#include "main_window.h"
#include <fstream>
#include <QDebug>
#include <QApplication>
#include <iterator>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <sstream>
#include <Eigen/Dense>

#include "global.h"
#include "vector_lib.h"
#include "print_macro.h"
#include "open_gl_qt.h"
#include "qt_object_selector.h"
#include "subspace_pose_sample_simulator.h"
#include "tet.h"

void Init() {
//  Test(); exit(0);
  global::gl->InstallHandler(SubspacePoseSampleSimulator::Instance());
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

