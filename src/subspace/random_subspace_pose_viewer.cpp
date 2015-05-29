#include <fstream>
#include "random_subspace_pose_viewer.h"
#include "global.h"
#include "string_formatter.h"
#include "multi_domain_tet.h"

namespace random_subspace_pose_viewer {
MultiDomainTet* multi_tet = NULL;
int pose_num = 100;
int current_pose = 0;
int render_mode = 0;
std::vector<MultiDomainTet::Vec> pose_q;

void Init() {
  multi_tet = new MultiDomainTet(GetMeshFile());
  multi_tet->LoadPartitionInfo(GetPartitionFolder());
  auto GetSubspaceFileName = [](int part_id) -> const char* {
    static std::string name;
    name = dj::Format("%s/modal_basis/partition_%d.basis.bin", GetDataFolder(), part_id);
    return name.c_str();
  };
  multi_tet->LoadSubspace(GetSubspaceFileName, MultiDomainTet::kBinary);
  multi_tet->AssembleGlobalBasis();
  // read randomly sampled q
  {
    std::ifstream in(dj::Format("%z/random_pose.txt"), std::ios::binary);
    int part_num;
    in >> pose_num >> part_num;
    ASSERT(part_num == multi_tet->part_num_);
    for (int p = 0; p < part_num; ++p) {
      int basis_size;
      in >> basis_size;
      ASSERT(basis_size == multi_tet->part_basis_size_[p]);
    }
    pose_q.resize(pose_num);
    for (int i = 0; i < pose_num; ++i) {
      pose_q[i] = MultiDomainTet::Vec::Zero(multi_tet->total_basis_num_);
      for (int n = 0; n < multi_tet->total_basis_num_; ++n) {
        in >> pose_q[i][n];
      }
    }
  }
}

void NextPose() {
  current_pose = (current_pose + 1) % pose_num;
  multi_tet->q_ = pose_q[current_pose];
  multi_tet->UpdatePosition();
}
void PrevPose() {
  current_pose--;
  if (current_pose < 0) current_pose = pose_num - 1;
  multi_tet->q_ = pose_q[current_pose];
  multi_tet->UpdatePosition();
}
}

using namespace random_subspace_pose_viewer;
int RandomSubspacePoseViewer::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      render_mode = (render_mode + 1) / MultiDomainTet::kRenderModeNum;
      break;
    case Qt::Key_J:
      PrevPose();
      break;
    case Qt::Key_K:
      NextPose();
      break;
    default:
      break;
  }
  return kNotHandled;
}

int RandomSubspacePoseViewer::Render() {
  multi_tet->Render(MultiDomainTet::kRenderMode[render_mode], multi_tet->X);
  return kNotHandled;
}

RandomSubspacePoseViewer::RandomSubspacePoseViewer() {
  Init();
}
