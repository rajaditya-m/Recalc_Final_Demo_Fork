#include <fstream>
#include "pose_viewer_simulator.h"
#include "tet.h"
#include "pose_sampler.h"
#include "rainbow_color.h"
#include "tet_mesh_simulator.h"
#include "global.h"
#include "string_formatter.h"
#include "open_gl_qt.h"

namespace pose_viwer {
using namespace global;
bool view_local_u = 0;
PoseSampler* pose_sampler;
int pose_id = 0;
int pose_num = 100;
std::vector<std::vector<double> > local_u;

void LoadLocalU() {
  if (!view_local_u) return;
  char file_name[512];
  dj::Format("%s/pose/local_deformation.txt", GetDataFolder());
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  int v_num;
  in >> v_num >> pose_num;
  local_u.resize(pose_num);
  for (int p = 0; p < pose_num; ++p) {
    local_u[p].resize(current_body->vertex_num_ * 3);
    for (int v = 0; v < current_body->vertex_num_; ++v) {
      in >> local_u[p][v * 3 + 0];
      in >> local_u[p][v * 3 + 1];
      in >> local_u[p][v * 3 + 2];
    }
  }
  in.close();
}

void ApplyLocalU(int pose_idx) {
  if (!view_local_u) return;
  for (int i = 0; i < current_body->vertex_num_; ++i) {
    current_body->tmp_X[i * 3 + 0] = current_body->rest_pos_[i * 3 + 0] + local_u[pose_idx][i * 3 + 0];
    current_body->tmp_X[i * 3 + 1] = current_body->rest_pos_[i * 3 + 1] + local_u[pose_idx][i * 3 + 1];
    current_body->tmp_X[i * 3 + 2] = current_body->rest_pos_[i * 3 + 2] + local_u[pose_idx][i * 3 + 2];
  }

  //  dj::Swap(current_body->tmp_X, current_body->X);
  //  current_body->ApplySkeletonTransformationToVertex();
  //  dj::Swap(current_body->tmp_X, current_body->X);
}

std::string GetFileName() {
  return dj::Format("%s/pose/pose_%d.txt", GetDataFolder(), pose_id);
//  char file_name[512];
//  //  sprintf(file_name, DATA_DIRECTORY "armadillo/pose_eu1e6_mass1_possoin_49/pose_%d.txt", pose_id);
//  //  sprintf(file_name, DATA_DIRECTORY "armadillo/pose_eu1e6_rho_100_possion_45/pose_%d.txt", pose_id);
//  sprintf(file_name, DATA_DIRECTORY "armadillo/pose/pose_%d.txt", pose_id);
//  return std::string(file_name);
}

std::string GetNextFileName() {
  pose_id = (pose_id + 1) % pose_num;
  std::string name = GetFileName();
  return name;
}

std::string GetPrevFileName() {
  if (pose_id == 0) pose_id = pose_num - 1;
  else pose_id--;
  std::string name = GetFileName();
  return name;
}

void Init() {
  TetMeshSimulator::Instance();
  int joints[] = {2, 3, 4, 6, 7, 8, 14, 15, 17, 18, 12};
  int num = sizeof(joints) / sizeof(int);
  std::vector<int> moving_joints(&joints[0], &joints[0] + num);
  pose_sampler = new PoseSampler(current_body);
  pose_sampler->LoadPoseResult(GetFileName().c_str());
  LoadLocalU();
  ApplyLocalU(pose_id);
}

}// namespace pose_viwer

using namespace pose_viwer;

int PoseViwerSimulator::HandleKeyPress(QKeyEvent *e) {
  using namespace global;
  switch (e->key()) {
    case Qt::Key_K:
      pose_sampler->LoadPoseResult(GetNextFileName().c_str());
      ApplyLocalU(pose_id);
      break;
    case Qt::Key_J:
      pose_sampler->LoadPoseResult(GetPrevFileName().c_str());
      ApplyLocalU(pose_id);
      break;
    case Qt::Key_R:
      current_body->NextRenderMode();
    default:
      return -1;
  }
  return 0;
}

int PoseViwerSimulator::Render() {
  current_body->Render(Tet::kDefaultRendering, current_body->X);
  if (view_local_u) {
    dj::Swap(current_body->X, current_body->tmp_X);
    glPushMatrix();
    glTranslatef(0.8f, 0, 0);
    current_body->Render(Tet::kDefaultRendering, current_body->X);
    glPopMatrix();
    dj::Swap(current_body->X, current_body->tmp_X);
  }
  QFont sansFont("Helvetica [Cronyx]", 14);
  sansFont.setFamily("sans serif");
  int font_height = QFontMetrics(sansFont).height();
  glColor3fv(kGreen());
  global::gl->renderText(5, font_height + 5,
                         QString("pose id: %1 / %2").arg(pose_id).arg(pose_num - 1),
                         sansFont);
  return kNotHandled;
}


PoseViwerSimulator::PoseViwerSimulator() : InputHandler(2) {
  Init();
}

PoseViwerSimulator::~PoseViwerSimulator() {
}
