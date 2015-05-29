#include "main_window.h"
#include <QDebug>
#include <QFileDialog>
#include <QApplication>
#include "tet.h"
#include "tetrahedral_mesh.h"
#include "affine_transformer.h"
#include "triangular_mesh.h"
#include "tet_gen_mesh_io.h"
#include "obj_mesh_io.h"
#include "off_mesh_io.h"
#include "print_macro.h"
#include "mesh_viewer.h"
#include "global.h"
#include "open_gl_qt.h"
#include "vega_tet_mesh_io.h"

template <int file_type>
void SetUpViewer(std::vector<std::string>& mesh_files) {
  typedef typename MeshChooser<file_type>::Mesh Mesh;
  typedef typename MeshChooser<file_type>::MeshIO MeshIO;
#ifdef _WIN32
  const char* kOutputPath = "C:/tmp/render";
#elif defined(__APPLE__)
  const char* kOutputPath = "/Volumes/ram/";
#else
  const char* kOutputPath = "/tmp/log";
#endif
  Mesh* mesh_ = new Mesh(MeshIO::Instance(), mesh_files[0].c_str());
  MeshViewer<file_type>* viewer = new MeshViewer<file_type>(mesh_, mesh_files);
  global::gl->InstallHandler(viewer);
  viewer->set_output_path(kOutputPath);
}

void Init() {
  const char* kFileFormats =
    "Obj Mesh (*.obj);;"
    "TetGen Mesh (*.ele);;"
    "OFF Mesh (*.off);;"
    "Vega Mesh (*.veg)";
  QStringList selected_files = QFileDialog::getOpenFileNames(0, "Select mesh files", DATA_DIRECTORY, kFileFormats);
  selected_files.clear();
  selected_files.append(QString("/home/dj/Dropbox/subspace_simulator/data/centipede_220k/centipede_220k.ele"));
  if (selected_files.size() == 0) {
    L("no files selected.");
    exit(0);
  }
  std::vector<std::string> mesh_files;
  for (int i = 0; i < selected_files.size(); ++i) {
    mesh_files.push_back(selected_files[i].toStdString());
  }
  QString& first_file = selected_files[0];
  if (first_file.endsWith(".ele")) {
    // strip off .ele subfix
    for (int i = 0; i < selected_files.size(); ++i) {
      auto last_dot = mesh_files[i].find_last_of(".");
      mesh_files[i] = mesh_files[i].substr(0, last_dot);
    }
    if (0) {
      std::vector<std::string> cmd;
//      cmd.push_back("centerize");
//      cmd.push_back("rotate_x 60");
//      cmd.push_back("rotate_z 180");
//      cmd.push_back("range 1.0");
      cmd.push_back("min_y 0");
      AffineTransformer<double> transform(cmd);
      Tet tet(mesh_files[0].c_str(), 0, &transform);
      std::vector<double> verts(tet.X, tet.X + tet.vertex_num_ * 3);
      std::vector<int> tets(tet.tet_, tet.tet_ + tet.tet_number * 4);
      TetGenMeshIO::Instance()->Write(mesh_files[0].c_str(), verts, tets);
      VegaTetMeshIO::Instance()->Write(mesh_files[0].c_str(), verts, tets);
      P(tet.vertex_num_, tet.tet_number, tet.surface_vertices_.size(), tet.surface_vertices_.size());
            exit(0);
    }
    SetUpViewer<MeshFileType::kTetGen>(mesh_files);
  } else if (first_file.endsWith(".obj")) {
    if (0) {
      std::vector<std::string> cmd;
      cmd.push_back("centerize");
      cmd.push_back("range 1");
      cmd.push_back("min_y 0");
      AffineTransformer<double> transform(cmd);
      TriangularMesh tri(ObjMeshIO::Instance(), mesh_files[0].c_str(), &transform);
      std::vector<double> verts(&tri.vert_[0][0], &tri.vert_.back()[0] + 3);
      std::vector<int> tris(&tri.tri_[0][0], &tri.tri_.back()[0] + 3);
      ObjMeshIO::Instance()->Write(mesh_files[0].c_str(), verts, tris);
      P(tri.v_num_, tri.tri_num_, tri.tri_.size());
      exit(0);
    }

    SetUpViewer<MeshFileType::kObj>(mesh_files);
  } else if (first_file.endsWith(".off")) {
    SetUpViewer<MeshFileType::kOff>(mesh_files);
  } else if (first_file.endsWith("*.veg")) {
    SetUpViewer<MeshFileType::kVega>(mesh_files);
  } else {
    ASSERT(false, L("Invalid file type selected."));
  }
}

int main(int argc, char * argv[]) {
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  w.hide();
  w.setWindowTitle(QString("mesh_viewer"));
  Init(); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}

