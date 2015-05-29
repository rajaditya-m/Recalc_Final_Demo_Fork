#include <Eigen/Dense>
#include "main_window.h"
#include <QDebug>
#include <QApplication>
#include <iterator>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include "global.h"
#include "vector_lib.h"
#include "print_macro.h"
#include <sstream>
#include "volumetric_rod_simulator.h"
#include "open_gl_qt.h"

void EigenTest() {
  using namespace Eigen;
  using namespace std;
  typedef Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor> RMatXd;
  typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> Mat2d;
  typedef Eigen::Matrix<double, 2, 1> RVec2d;
//  double data[6] = {1, 2, 3, 4, 5, 6};
  P(sizeof(Eigen::Vector2f) / sizeof(float));
  P(sizeof(Eigen::Matrix3f) / sizeof(float));
  P(sizeof(Eigen::Matrix3d) / sizeof(double));
  P(sizeof(Eigen::Vector2d) / sizeof(double));
  Vector3d x(1, 2, 3);
  P(x);
  double* p = &x[0];
  P(dj::Vec3d(p));
}

void Init() {
  global::gl->InstallHandler(VolumetricRodSimulator::Instance());
}

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  Init();
  w.show();
  return a.exec();
}

