#include "subspace_pose_viewer.h"
#include <fstream>
#include "multi_domain_tet.h"
#include "tet.h"
#include "string_formatter.h"
#include "skeleton_node.h"
#include "rainbow_color.h"
#include "skeleton.h"
#include "global.h"
#include "open_gl_qt.h"
#include "pose_sampler.h"
#include "tet_mesh_simulator.h"
#include "affine_transformer.h"
#include "qt_object_selector.h"
#include "config_file.h"

namespace subspace_pose_viewer {
MultiDomainTet* multi_tet;
PoseSampler* pose_sampler;
int pose_id = 0;
int pose_num = 100;

using namespace global;
void CreateBody() {
  L("Start creating bodies");
  std::vector<std::string>* transform_cmd = conf.Get<std::vector<std::string>* >("affine_transform");
  AffineTransformer<double> transform(*transform_cmd);
//  body::body_mesh_file[0] = global::data_directory + body::body_mesh_file[0];
  multi_tet = new MultiDomainTet(GetMeshFile(), 0, &transform);
  multi_tet->AttachSkeleton(global::skeleton);
  L("finish creating body");
}

void CreateSkeleton() {
  L("Create Skeleton");
  skeleton = new Skeleton(global::skeleton_file.c_str(), true);
}

std::string GetFileName(int id = -1) {
  char file_name[512];
  //  sprintf(file_name, DATA_DIRECTORY "armadillo/pose_eu1e6_mass1_possoin_49/pose_%d.txt", pose_id);
  //  sprintf(file_name, DATA_DIRECTORY "armadillo/pose_eu1e6_rho_100_possion_45/pose_%d.txt", pose_id);
  if (id == -1) {
    sprintf(file_name, "%s/pose/pose_%d.txt", GetDataFolder(), pose_id);
  } else {
    sprintf(file_name, "%s/pose/pose_%d.txt", GetDataFolder(), id);
  }
  P(pose_id, file_name);
  return std::string(file_name);
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


void InitializeSubspaceSim() {
  multi_tet->RotateBasis();
  //  multi_tet->ComputeGlobalPositionFromSubspace();
  multi_tet->vel_q_ = MultiDomainTet::Vec::Zero(multi_tet->vel_q_.size());
  multi_tet->ProjectPositionToSubspace();
  multi_tet->ComputeGlobalPositionFromSubspace();
  multi_tet->UpdateVegaVertexOffset();
}

void LoadSubspacePose(int id) {
  char file[512];
  sprintf(file, "%s/subspace_pose/pose_%d.txt", GetDataFolder(), id);
  std::ifstream in(file);
  ASSERT(in.is_open(), P(file));
  int v_num, basis_num;
  in >> v_num >> basis_num;
  ASSERT(v_num == multi_tet->vertex_num_);
  ASSERT(basis_num == multi_tet->q_.size());
  for (int i = 0; i < basis_num; ++i) {
    in >> multi_tet->q_[i];
  }
  in.close();
  multi_tet->RotateBasis();
  //  multi_tet->ComputeGlobalPositionFromSubspace();
  multi_tet->vel_q_ = MultiDomainTet::Vec::Zero(multi_tet->vel_q_.size());
  multi_tet->ComputeGlobalPositionFromSubspace();
}

void ExportSubspacePoseResult() {
  L("exporting pose result");
  char file_name[512];
  sprintf(file_name, "%s/subspace_pose/pose_all.txt", GetDataFolder());
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  out << multi_tet->vertex_num_ << std::endl;
  out << pose_num << std::endl;
  out << multi_tet->part_num_ << std::endl;
  for (int p = 0; p < multi_tet->part_num_; ++p) {
    out << multi_tet->domain_attached_bone_[p] << " " << multi_tet->part_basis_size_[p] << std::endl;
  }
  for (int i = 0; i < pose_num; ++i) {
    // read subspace pose result
    std::string name = GetFileName(i);
    pose_sampler->LoadPoseResult(name.c_str());
    sprintf(file_name, "%s/subspace_pose/pose_%d.txt", GetDataFolder(), i);
    P(file_name);
    std::ifstream in(file_name);
    ASSERT(in.is_open(), P(file_name));
    int v_num, basis_num;
    in >> v_num >> basis_num;
    ASSERT(v_num == multi_tet->vertex_num_);
    ASSERT(basis_num == multi_tet->q_.size());
    for (int j = 0; j < basis_num; ++j) {
      in >> multi_tet->q_[j];
    }
    in.close();
    multi_tet->ComputePartAffineTransformamtion();
    if (0) {
      // output obj file
      multi_tet->RotateBasis();
      multi_tet->ComputeGlobalPositionFromSubspace();
      sprintf(file_name, "%s/subspace_pose/obj/pose_%d.obj", GetDataFolder(), i);
      multi_tet->Save_OBJ_File(file_name);
    }
    //    P(multi_tet->part_rotation_[0]);exit(0);
    // export result for current pose
    for (int p = 0; p < multi_tet->part_num_; ++p) {
      int part_id = multi_tet->domain_attached_bone_[p];
      out << part_id << std::endl;
      // part rotation
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          out << multi_tet->part_rotation_[p](r, c) << " ";
        }
        out << std::endl;
      }
      // part mass center
      out << multi_tet->initial_center_of_mass_[p][0] << " ";
      out << multi_tet->initial_center_of_mass_[p][1] << " ";
      out << multi_tet->initial_center_of_mass_[p][2] << std::endl;
      out << multi_tet->center_of_mass_[p][0] << " ";
      out << multi_tet->center_of_mass_[p][1] << " ";
      out << multi_tet->center_of_mass_[p][2] << std::endl;
      ASSERT(multi_tet->basis_offset_[p + 1] - multi_tet->basis_offset_[p] == multi_tet->part_basis_size_[p]);
      for (int k = multi_tet->basis_offset_[p]; k < multi_tet->basis_offset_[p + 1]; ++k) {
        out << multi_tet->q_[k] << std::endl;
      }
    }
  }
  out.close();
  L("finished exporting subspace pose result");
  exit(0);
}

void Init() {
  CreateSkeleton();
  CreateBody();
  global::current_body = multi_tet;
  //  multi_tet->LoadMeshPartition(DATA_DIRECTORY "armadillo/armadillo_150k.partition.txt");
  multi_tet->LoadPartitionInfo(GetPartitionFolder());
  //        DATA_DIRECTORY "armadillo/partition/vertex_partition_info.txt",
  //                               DATA_DIRECTORY "armadillo/partition/vertex_partition_info.txt");
  //  multi_tet->LoadMass(DATA_DIRECTORY "armadillo/vertex_mass.txt");
  //  multi_tet->LoadLocalMass(DATA_DIRECTORY "armadillo//vertexMass");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/Basis2/Local_Modal_Basis");
  auto GetSubspaceFileName = [](int part_id) -> const char* {
    static char name[512];
    //    sprintf(name, DATA_DIRECTORY "armadillo/Basis2/Local_Basis/%d/ascii_Basis_%d_1.txt", part_id, part_id);
    //    sprintf(name, DATA_DIRECTORY "armadillo/Basis2/Local_Basis/%d/ascii_Basis_%d_2.txt", part_id, part_id);
    sprintf(name, "%s/modal_basis/partition_%d.basis.bin", GetDataFolder(), part_id);
    return name;
  };
  multi_tet->LoadSubspace(GetSubspaceFileName, MultiDomainTet::kBinary);
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/Local_Modal_Basis");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/Local_Modal_Basis");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/Basis2/Local_Pose_Basis");
  multi_tet->AssembleGlobalBasis();
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/ascMat");
  //  multi_tet->SaveConstrainedVertex(DATA_DIRECTORY "armadillo/armadillo_150k.fixed_vertex.txt");
  global::gl->InstallHandler(new QtObjectSelector<Real>(multi_tet));
  pose_sampler = new PoseSampler(multi_tet);
  pose_sampler->LoadPoseResult(GetFileName(pose_id).c_str());
  //  InitializeSubspaceSim();
  LoadSubspacePose(pose_id);
      ExportSubspacePoseResult();
}

} // namespace subspace_pose_viewer

using namespace subspace_pose_viewer;

int SubspacePoseViewer::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      multi_tet->NextRenderMode();
      break;
    case Qt::Key_X: {
      char file_name[512];
      sprintf(file_name, "%s/body%d.obj", GetModelFolder(), 0);
      L("Save body " + std::string(file_name));
      global::current_body->Save_OBJ_File(file_name);
      break;
    }
    case Qt::Key_J: {
      std::string name = GetPrevFileName();
      pose_sampler->LoadPoseResult(name.c_str());
      LoadSubspacePose(pose_id);
      break;
    }
    case Qt::Key_K: {
      std::string name = GetNextFileName();
      pose_sampler->LoadPoseResult(name.c_str());
      LoadSubspacePose(pose_id);
      break;
    }
    default:
      return -1;
  }
  return kPriority_;
}

int SubspacePoseViewer::Render() {
  //  int t = 88813;
  //  int* verts = multi_tet->tet_ + t * 4;
  //  P(verts[0]);
  //  P(verts[1]);
  //  P(verts[2]);
  //  P(verts[3]);
  ////  P(dj::Vec4i(verts));
  //  {
  //    double* pos = multi_tet->X;
  //    P(dj::Vec3d(pos + verts[0] * 3));
  //    P(dj::Vec3d(pos + verts[1] * 3));
  //    P(dj::Vec3d(pos + verts[2] * 3));
  //    P(dj::Vec3d(pos + verts[3] * 3));

  //    pos = multi_tet->rest_pos_;
  //    P(dj::Vec3d(pos + verts[0] * 3));
  //    P(dj::Vec3d(pos + verts[1] * 3));
  //    P(dj::Vec3d(pos + verts[2] * 3));
  //    P(dj::Vec3d(pos + verts[3] * 3));
  //    exit(0);
  //  }
  multi_tet->Render(MultiDomainTet::kDefaultRendering, multi_tet->X);
  //  glPushMatrix();
  //  glTranslatef(0.2, 0, 0);
  //  multi_tet->Render(MultiDomainTet::kDefaultRendering, multi_tet->rest_pos_);
  //  glPopMatrix();
  //  exit(0);
  QFont sansFont("Helvetica [Cronyx]", 14);
  sansFont.setFamily("sans serif");
  int font_height = QFontMetrics(sansFont).height();
  glColor3fv(kGreen());
  global::gl->renderText(5, font_height + 5,
                         QString("pose: %1/%2").arg(pose_id).arg(pose_num - 1),
                         sansFont);
  return kNotHandled;
}


int SubspacePoseViewer::Idle() {
  return -1;
}

SubspacePoseViewer::SubspacePoseViewer()
  : InputHandler(2) {
  Init();
}

SubspacePoseViewer::~SubspacePoseViewer() {
  delete multi_tet;
}
