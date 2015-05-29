#include "subspace_simulator.h"
#include "main_window.h"
#include <QDebug>
#include <QApplication>
#include "print_macro.h"
#include "global.h"
#include "open_gl_qt.h"
#include "cubica_simulator.h"
#include <Eigen/Dense>

void Test() {
  Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
  auto a = mat.row(0);
  a[0] = 1;
  a[1] = 2;
  a[2] = 3;
  PMATCOL(mat);
  exit(0);
}

void Init() {
//  Test();
  global::gl->InstallHandler(CubicaSimulator::Instance());
}


int main(int argc, char *argv[]) {
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
  Init();
  w.show();
  return a.exec();
}

