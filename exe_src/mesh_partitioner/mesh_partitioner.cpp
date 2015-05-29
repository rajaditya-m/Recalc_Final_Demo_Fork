#include <vector>
#include <string>
#include "main_window.h"
#include <QDebug>
#include <QApplication>
#include <QFileDialog>
#include <QString>
#include "interactive_mesh_partitioner.h"
#include "affine_transformer.h"
#include "global.h"
#include "open_gl_qt.h"
#include "tet.h"
#include "config_file.h"
#include "string_formatter.h"

void Init() {
  //  global::current_body = new Tet(dj::Format("%s/%s", GetModelFolder(), pose_conf.Get<std::string>("body_mesh_file")).c_str(), 1e-15);
  std::vector<std::string> cmd;
//  cmd.push_back("centerize");
//  cmd.push_back("rotate_x 60");
//  cmd.push_back("rotate_z 180");
//  cmd.push_back("range 1.0");
//  cmd.push_back("min_y 0");
  AffineTransformer<double> transformer(cmd);
  if (0) {
    QString fileName = QFileDialog::getOpenFileName(0, "Select tetrahedral mesh",
                                                    DATA_DIRECTORY,
                                                    "TetGen Mesh (*.ele)");
    int last_dot = fileName.lastIndexOf(QString("."));
    fileName = fileName.left(last_dot);
    global::current_body = new Tet(fileName.toStdString().c_str(), 0, &transformer);
  } else {
    global::current_body = new Tet(GetMeshFile(), 0, &transformer, false);
  }
  //  global::current_body = new Tet(DATA_DIRECTORY "octopus/octopus", 0, NULL);
  global::gl->InstallHandler(new InteractiveMeshPartitioner(global::current_body));
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

