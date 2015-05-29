#include "subspace_pose_sample_simulator.h"
#include <fstream>
#include "multi_domain_tet.h"
#include "tet.h"
#include "skeleton_node.h"
#include "skeleton.h"
#include "global.h"
#include "string_formatter.h"
#include "string_formatter.h"
#include "open_gl_qt.h"
#include "pose_sampler.h"
#include "tet_mesh_simulator.h"
#include "affine_transformer.h"
#include "qt_object_selector.h"
#include "config_file.h"

namespace subspace_pose_sample_simulator {
MultiDomainTet* multi_tet;
PoseSampler* pose_sampler;
int pose_id = 0;
int pose_num = 100;
int min_pose_idx = -1;
int max_pose_idx = -1;

using namespace global;
void CreateBody() {
  L("Start creating bodies");
  std::vector<std::string>* transform_cmd = conf.Get<std::vector<std::string>* >("affine_transform");
  AffineTransformer<double> transform(*transform_cmd);
//  body::body_mesh_file[0] = global::data_directory + body::body_mesh_file[0];
  std::string mesh_file = dj::Format("%s/%z", GetModelFolder(), conf.Get<std::string>("body_mesh_file"));
  multi_tet = new MultiDomainTet(mesh_file.c_str(), 0, &transform);
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
    sprintf(file_name,  "%s/pose/pose_%d.txt", GetDataFolder(), pose_id);
    //            DATA_DIRECTORY "armadillo/pose/pose_%d.txt", pose_id);
  } else {
    sprintf(file_name,  "%s/pose/pose_%d.txt", GetDataFolder(), id);
  }
  return std::string(file_name);
}

std::string GetNextFileName() {
  std::string name = GetFileName();
  pose_id = (pose_id + 1) % pose_num;
  return name;
}

std::string GetPrevFileName() {
  std::string name = GetFileName();
  if (pose_id == 0) pose_id = pose_num - 1;
  else pose_id--;
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

void Init() {
  CreateSkeleton();
  CreateBody();
  global::current_body = multi_tet;
  //  multi_tet->LoadMeshPartition(DATA_DIRECTORY "armadillo/armadillo_150k.partition.txt");
  multi_tet->LoadPartitionInfo(GetPartitionFolder());
  //  multi_tet->LoadMass(DATA_DIRECTORY "armadillo/vertex_mass.txt");
  //  multi_tet->LoadLocalMass(DATA_DIRECTORY "armadillo//vertexMass");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/Basis2/Local_Modal_Basis");
  auto GetSubspaceFileName = [](int part_id) -> const char* {
    static char name[512];
    // basis from e=1e6, density = 1100, poisson = 0.49
    //    sprintf(name, DATA_DIRECTORY "armadillo/Basis2/Local_Basis/%d/ascii_Basis_%d_1.txt", part_id, part_id);
    // basis from e=1.5e5, density = 1000, poisson = 0.45
    sprintf(name, "%s/modal_basis/partition_%d.basis.bin", GetDataFolder(), part_id);
    return name;
  };
  multi_tet->LoadSubspace(GetSubspaceFileName, MultiDomainTet::kBinary);
  //  multi_tet->LoadSubspace(GetSubspaceFileName);
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/Local_Modal_Basis");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/Local_Modal_Basis");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/Basis2/Local_Pose_Basis");
  multi_tet->AssembleGlobalBasis();
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/ascMat");
  //  multi_tet->SaveConstrainedVertex(DATA_DIRECTORY "armadillo/armadillo_150k.fixed_vertex.txt");
  global::gl->InstallHandler(new QtObjectSelector<Real>(multi_tet));
  pose_sampler = new PoseSampler(multi_tet);
  min_pose_idx = conf.Get<int>("min_pose_idx");
  max_pose_idx = conf.Get<int>("max_pose_idx");;
  pose_id = min_pose_idx;
  pose_sampler->LoadPoseResult(GetFileName(min_pose_idx).c_str());
  InitializeSubspaceSim();
}

void SaveSubspacePoseResult() {
  std::stringstream stream;
  stream << dj::Format("%s/subspace_pose/pose_%d.txt", GetDataFolder(), pose_id);
//            DATA_DIRECTORY << "armadillo/subspace_pose/pose_" << pose_id << ".txt";
  std::ofstream out(stream.str().c_str());
  ASSERT(out.is_open(), P(stream.str()));
  out << multi_tet->vertex_num_ << std::endl;
  out << multi_tet->q_.size() << std::endl;
  for (int i = 0; i < multi_tet->q_.size(); ++i) {
    out << multi_tet->q_[i] << std::endl;
  }
  out.close();
  L(stream.str() + " saved.");
}

} // namespace subspace_pose_sample_simulator

using namespace subspace_pose_sample_simulator;

int SubspacePoseSampleSimulator::HandleKeyPress(QKeyEvent *e) {
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
    default:
      return -1;
  }
  return kPriority_;
}

int SubspacePoseSampleSimulator::Render() {
  multi_tet->Render(MultiDomainTet::kDefaultRendering, multi_tet->X);
  return kNotHandled;
}


int SubspacePoseSampleSimulator::Idle() {
  //  multi_tet->FullSimulation(global::time_step);
  //  multi_tet->Simulate(global::time_step);
  //  return -1;
  int current_pose_idx = pose_id;
  if (current_pose_idx < min_pose_idx || current_pose_idx > max_pose_idx) {
    return -1;
  }
  double sum = 0;
  for (int i = 0; i < 10; ++i) {
    multi_tet->Simulate(global::time_step);
    sum += (multi_tet->vel_q_.dot(multi_tet->vel_q_)) / multi_tet->vel_q_.size();
  }
  // max velocity
  sum /= 10;
  P(pose_id, sum, min_pose_idx, max_pose_idx);
  if (sum < 1e-10) {
    L("finished pose idx = ", current_pose_idx);
    //    std::stringstream stream;
    //    stream << DATA_DIRECTORY << "armadillo/subspace_pose/pose_" << current_pose_idx << ".txt";
    SaveSubspacePoseResult();
    //    pose_sampler->SavePoseResult(stream.str().c_str());
    current_pose_idx++;
    if (current_pose_idx <= max_pose_idx) {
      //      pose_sampler->SetPose(current_pose_idx);
      //            arch::current_body->ApplySkeletonTransformationToVertex();
      //            current_body->vega_inv_fem_->LoadPosition(current_body->X);
      pose_sampler->LoadPoseResult(GetFileName(current_pose_idx).c_str());
      InitializeSubspaceSim();
      pose_id++;
    } else {
//      std::stringstream stream;
//      stream << "finish poses from " << min_pose_idx << " to " << max_pose_idx;
      L(dj::Format("finish poses from %d to %d", min_pose_idx, max_pose_idx));
      exit(0);
    }
  }
  return -1;
}

SubspacePoseSampleSimulator::SubspacePoseSampleSimulator()
  : InputHandler(2) {
  Init();
}

SubspacePoseSampleSimulator::~SubspacePoseSampleSimulator() {
  delete multi_tet;
}
