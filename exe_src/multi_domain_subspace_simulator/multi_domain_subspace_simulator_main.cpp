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
#include "multi_domain_subspace_simulator.h"
#include "qt_object_selector.h"
#include "tet.h"
#include "eig3.h"
#include "string_formatter.h"

void Test() {
#if 0
  test::Profiler<Timer, float> profiler;
  profiler.StartTimer("all");
  profiler.StartTimer("Haha");
  profiler.StartTimer("x");
  double sum = 0;
  for (int i = 0; i < 10000000; ++i) {
    sum += sqrt(i);
  }
  profiler.EndTimer("x");
  profiler.StartTimer("2");
  sum = 0;
  for (int i = 0; i < 20000000; ++i) {
    sum += sqrt(i);
  }
  profiler.EndTimer("2");
  profiler.EndTimer("1");
  profiler.StartTimer("3");
  sum = 0;
  for (int i = 0; i < 30000000; ++i) {
    sum += sqrt(i);
  }
  profiler.StartTimer("4");
  sum = 0;
  for (int i = 0; i < 10000000; ++i) {
    sum += sqrt(i);
  }
  profiler.EndTimer("4");
  profiler.EndTimer("3");
  profiler.EndTimer("all");
#endif
}

#include "profiler.h"
void Init() {
//  Test();
  global::gl->InstallHandler(MultiDomainSimulator::Instance());
  global::gl->InstallHandler(new QtObjectSelector<double>(global::current_body));
}

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  w.hide();
  w.setWindowTitle(QString("multi_domain_subspace_simulator"));
  Init(); // Make sure OpenGL environment is setup before init

  //  int shader = 1;
  //  glUseProgram(shader);
  //  int loc = glGetUniformLocation(shader, "test");
  //  ASSERT(loc >= 0, P(loc));
  //  glUniform1f(loc, 0.8f);
  //  glUseProgram(0);
  w.show();
  return a.exec();
}

