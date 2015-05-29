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
#include "mixed_multi_domain_simulator.h"
#include "qt_object_selector.h"
#include "tet.h"
#include "eig3.h"
#include "string_formatter.h"

void Test() {
  //  P(dj::Format("a\tb"));
  //  P(formatter::Format("is %zsdfk", 'c'));
  //  P(formatter::Format("the %ssdf %z", "abc", 934));
  //  P(formatter::Format("the %05dsdf", 1));
  //  P(formatter::Format("the %.3fsdf", 2.33333333));
  exit(0);
}

void Init() {
//    Test();
  global::gl->InstallHandler(MixedMultiDomainSimulator::Instance());
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
  Init(); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}

