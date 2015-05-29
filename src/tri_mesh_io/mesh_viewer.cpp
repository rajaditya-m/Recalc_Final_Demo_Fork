#include "mesh_viewer.h"
#include "tet_gen_mesh_io.h"
#include "obj_mesh_io.h"
#include "off_mesh_io.h"
#include "triangular_mesh.h"
#include "vega_tet_mesh_io.h"
#include "tetrahedral_mesh.h"
#include "tet.h"
#include "global.h"
#include "rainbow_color.h"
#include "open_gl_qt.h"


template <int file_type>
MeshViewer<file_type>::MeshViewer(Mesh *mesh, std::vector<std::string> &mesh_files)
  : InputHandler(100)
  , mesh_(mesh)
  , mesh_files_(mesh_files) {
  play_animation_ = false;
  frame_rate_ = 30.0;
  current_mesh_idx_ = 0;
  output_path_ = ".";
  ASSERT(int(mesh_files.size()) > 0);
  LoadMesh(mesh_files[0]);
}

template <int file_type>
void MeshViewer<file_type>::DumpRenderResultToPng(const char *folder) {
  if (folder == NULL) folder = output_path_.c_str();
  dumping_render_result_ = true;
  char file[512];
  for (int i = 0; i < int(mesh_files_.size()); ++i) {
    LoadMesh(mesh_files_[i]);
    //    mesh_->Render();
    global::gl->updateGL();
    //    this->Render();
    sprintf(file, "%s/%05d.png", folder, i);
    global::gl->ScreenShot(file, "png");
  }
  sprintf(file, "finish dumping rendering result to folder %s.", folder);
  L(file);
  LoadMesh(mesh_files_[current_mesh_idx_]);
  dumping_render_result_ = false;
}

template <int file_type>
void MeshViewer<file_type>::set_frame_rate(double frame_rate) {
  this->frame_rate_ = frame_rate;
}

template <int file_type>
void MeshViewer<file_type>::set_output_path(std::string path) {
  output_path_ = path;
}

template <int file_type>
void MeshViewer<file_type>::LoadMesh(std::string mesh_file) {
  std::vector<Real> pos;
  std::vector<int> eles;
  MeshIO::Instance()->Read(mesh_file.c_str(), pos, eles);
  mesh_->LoadPosition(pos);
}

template <int file_type>
int MeshViewer<file_type>::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      mesh_->NextRenderMode();
      break;
    case Qt::Key_D: {
      DumpRenderResultToPng(NULL);
      break;
    }
    case Qt::Key_P:
      play_animation_ = !play_animation_;
      prev_time_ = WorldTime::Instance()->NanoSecond();
      break;
    case Qt::Key_K:
      current_mesh_idx_ = (current_mesh_idx_ + 1) % int(mesh_files_.size());
      LoadMesh(mesh_files_[current_mesh_idx_]);
      break;
    case Qt::Key_J:
      current_mesh_idx_--;
      if (current_mesh_idx_ < 0) current_mesh_idx_ = int(mesh_files_.size()) - 1;
      LoadMesh(mesh_files_[current_mesh_idx_]);
      break;
    default:
      return kNotHandled;
      break;
  }
  return kNotHandled;
}

template <int file_type>
int MeshViewer<file_type>::Render() {
//      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
//    glBegin(GL_QUADS);
//    glNormal3f(0, 0, 1);
//    glVertex3f(0, 0, 0);
//    glVertex3f(0, 1, 0);
//    glVertex3f(1, 1, 0);
//    glVertex3f(1, 0, 0);
//    glEnd();
  mesh_->Render();
  if (dumping_render_result_) {
    return kNotHandled;
  } else {
    QFont sansFont("Helvetica [Cronyx]", 14);
    sansFont.setFamily("sans serif");
    int font_height = QFontMetrics(sansFont).height();
    glColor3fv(kGreen());
    global::gl->renderText(5, font_height + 5,
                           QString("file: %1").arg(mesh_files_[current_mesh_idx_].c_str()),
                           sansFont);
    return kNotHandled;
  }
}

template <int file_type>
int MeshViewer<file_type>::Idle() {
  if (play_animation_) {
    uint64_t current_time = WorldTime::Instance()->NanoSecond();
    double time_elapsed = (current_time - prev_time_) * 1e-9;
    if (time_elapsed > (1.0 / frame_rate_)) {
      current_mesh_idx_ = (current_mesh_idx_ + 1) % int(mesh_files_.size());
      LoadMesh(mesh_files_[current_mesh_idx_]);
      prev_time_ = current_time;
    }
  }
  return kNotHandled;
}

template class MeshViewer<MeshFileType::kTetGen>;
template class MeshViewer<MeshFileType::kVega>;
template class MeshViewer<MeshFileType::kObj>;
template class MeshViewer<MeshFileType::kOff>;
