#include "sample_pose_simulator.h"
#include <fstream>
#include "skeleton_node.h"
#include "tet_mesh_simulator_bridge.h"
#include "tet.h"
#include "skeleton.h"
#include "pose_sampler.h"
#include "tet_mesh_simulator.h"
#include "global.h"
#include "config_file.h"
#include "string_formatter.h"

namespace sample_pose

{
using namespace global;
PoseSampler* pose_sampler;
int pose_id = -1;
int min_pose_idx, max_pose_idx;

const char* GetFileName(int idx) {
  static char file_name[512] ;
  //  sprintf(file_name, DATA_DIRECTORY "armadillo/pose/pose_%d.txt", idx);
  sprintf(file_name, "%s/pose/pose_%d.txt", GetDataFolder(), idx);
  return file_name;
}

void before() {
  skeleton->InitializeBoneRotation();
}

void after() {
  skeleton->UpdatePosition();
  skeleton->AssembleSkeletonTransformation();
  arch::current_body->ApplySkeletonTransformationToVertex();
  global::simulate = false;
}

void NextPose() {
  pose_id = (pose_id + 1) % pose_sampler->PoseNum();
  pose_sampler->SetPose(pose_sampler->poses_[pose_id]);
  //  pose_sampler->SetTargetPose(pose_sampler->poses_[pose_id], 1000);
  L("target pose set as ", pose_id);
  after();

  if (0) {
    std::ofstream out("/Users/dj/deformed.txt");
    out << current_body->vertex_num_ << " " << current_body->constrainted_vertex_.size() << std::endl;
    for (int v : current_body->constrainted_vertex_) {
      out << v << std::endl;
    }
    for (int v = 0; v < current_body->vertex_num_; ++v) {
      out << current_body->X[v * 3 + 0] << std::endl;
      out << current_body->X[v * 3 + 1] << std::endl;
      out << current_body->X[v * 3 + 2] << std::endl;
    }
    out.close();
    exit(0);
  }
  //  pose_sampler->SavePoseResult("/tmp/pose");
  //    pose_sampler->LoadPoseResult(DATA_DIRECTORY "armadillo/pose/pose_90.txt");
  //    pose_sampler->LoadPoseResult("/home/xwu/Dropbox/skeleton/data/armadillo/pose_eu1e6_rho_100_possion_45/pose_79.txt");
  //    std::vector<std::string> in;
  //    in.push_back(DATA_DIRECTORY "armadillo/pose/pose_90.txt");
  //    pose_sampler->ParseLocalDeformations(in, "", true);
  current_body->inv_fem_->LoadPosition(current_body->X);
}

void PrevPose() {
  pose_id--;
  if (pose_id < 0) {
    pose_id = pose_sampler->PoseNum() - 1;
  }
  pose_sampler->SetPose(pose_sampler->poses_[pose_id]);
  //  pose_sampler->SetTargetPose(pose_sampler->poses_[pose_id], 1000);
  L("target pose set as ", pose_id);
  after();
  current_body->inv_fem_->LoadPosition(current_body->X);
}

void Init() {
  TetMeshSimulator::Instance();
  ASSERT(current_body->constrainted_vertex_.size() > 0);
  int joints[] = {2, 3, 4, 6, 7, 8, 14, 15, 17, 18, 12};
  int num = sizeof(joints) / sizeof(int);
  std::vector<int> moving_joints(&joints[0], &joints[0] + num);
  pose_sampler = new PoseSampler(current_body);
  //  const char file_name[] = DATA_DIRECTORY "armadillo_pose_samples.txt";
  std::string file_name = dj::Format("%s/skeleton_pose.txt", GetDataFolder());
  pose_sampler->LoadSkeletonPoses(file_name.c_str());
  L(file_name, " loaded");

  //  return;
  min_pose_idx = conf.Get<int>("min_pose_idx");
  max_pose_idx = conf.Get<int>("max_pose_idx");
  //  pose_sampler->SetPose(min_pose_idx);
  {
    pose_sampler->LoadPoseResult(GetFileName(min_pose_idx));
    current_body->inv_fem_->Reset();
    current_body->inv_fem_->LoadPosition(current_body->X);
    pose_id = min_pose_idx;
  }

  //      pose_sampler->LoadPoseResult(DATA_DIRECTORY "armadillo/pose/pose_0.txt");
  if (0) {
    std::vector<std::string> input;
    for (int i = 0; i < 100; ++i) {
      char file_name[512] ;
      sprintf(file_name, "%s/pose/pose_%d.txt", GetDataFolder(), i);
      P(file_name);
      input.push_back(std::string(file_name));
    }
    pose_sampler->ParseLocalDeformations(input,
                                         dj::Format("%s/pose/local_deformation.txt", GetDataFolder()).c_str(),
                                         false);
    L("finish export data to local deformation");
    exit(0);
  }
}

} // namespace sample_pose


using namespace sample_pose;
int SamplePoseSimulator::HandleKeyPress(QKeyEvent *e) {
  using namespace global;
  switch (e->key()) {
    case Qt::Key_I: break;
    case Qt::Key_R: current_body->NextRenderMode(); break;
    //    case Qt::Key_K:
    //      NextPose();
    //      break;
    //    case Qt::Key_J:
    //      PrevPose();
    //      break;
    //    case Qt::Key_S: {
    //      const char file_name[] = DATA_DIRECTORY "armadillo_pose_samples.txt";
    //      pose_sampler->SaveSkeletonPoses(file_name);
    //      L(file_name, " saved");
    //      break;
    //    }
    //    case Qt::Key_L: {
    //      const char file_name[] = DATA_DIRECTORY "armadillo_pose_samples.txt";
    //      pose_sampler->LoadSkeletonPoses(file_name);
    //      L(file_name, " loaded");
    //      break;
    //    }
    default:
      return -1;
  }
  return 0;
}

int SamplePoseSimulator::Render() {
  current_body->Render(Tet::kDefaultRendering, current_body->X);
  return kNotHandled;
}


SamplePoseSimulator::SamplePoseSimulator() : InputHandler(2) {
  sample_pose::Init();
}

void SamplePoseSimulator::Init(int argc, char *argv[]) {
  if (argc == 1) return;
  if (argc != 3) {
    L("Must have two arguements.");
    return;
  }
  min_pose_idx = atoi(argv[1]);
  max_pose_idx = atoi(argv[2]);
  P(min_pose_idx, max_pose_idx);
  if (1) {
    pose_sampler->LoadPoseResult(GetFileName(min_pose_idx));
    current_body->inv_fem_->Reset();
    current_body->inv_fem_->LoadPosition(current_body->X);
    pose_id = min_pose_idx;
  }
}

int SamplePoseSimulator::Idle() {
  int current_pose_idx = pose_id;
  if (current_pose_idx < min_pose_idx || current_pose_idx > max_pose_idx) {
    current_body->inv_fem_->Simulate(global::time_step);
    return -1;
  }
  double max_velocity = 0;
  for (int i = 0; i < 10; ++i) {
    max_velocity += current_body->inv_fem_->Simulate(global::time_step);
  }
  // max velocity
  max_velocity /= 10;
  P(pose_id, max_velocity, min_pose_idx, max_pose_idx);
  if (max_velocity < 1e-1) {
    L("finished pose idx = ", current_pose_idx);
    std::stringstream stream;
//    stream << DATA_DIRECTORY << "armadillo/pose/pose_" << current_pose_idx << ".txt";
    std::string file_name = dj::Format("%s/pose/pose_%d.txt", GetDataFolder(), current_pose_idx);
    pose_sampler->SavePoseResult(file_name.c_str());
    L(file_name + " saved.");
    current_pose_idx++;
    if (current_pose_idx <= max_pose_idx) {
      //      pose_sampler->SetPose(current_pose_idx);
      //            arch::current_body->ApplySkeletonTransformationToVertex();
      //            current_body->vega_inv_fem_->LoadPosition(current_body->X);
      pose_sampler->LoadPoseResult(GetFileName(current_pose_idx));
      current_body->inv_fem_->Reset();
      current_body->inv_fem_->LoadPosition(current_body->X);
      pose_id++;
    } else {
      std::stringstream stream;
      stream << "finish poses from " << min_pose_idx << " to " << max_pose_idx;
      L(stream.str());
      exit(0);
    }
  }
  return -1;
}


SamplePoseSimulator::~SamplePoseSimulator() {
}
