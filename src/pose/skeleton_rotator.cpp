#include "skeleton_rotator.h"
#include <fstream>
#include <sstream>
#include "config_file.h"
#include "global.h"
#include "string_formatter.h"
#include "tet.h"
#include "skeleton_node.h"
//#include "simulator.h"
#include "opengl_helper.h"
#include "open_gl_qt.h"
#include "pose_sampler.h"
#include "tet_mesh_simulator_bridge.h"
#include "skeleton.h"
using namespace global;
//using namespace simulator;
static PoseSampler* pose_sampler;
int pose_id = -1;
int min_pose_idx, max_pose_idx;


void before() {
  //  global::gl->resetSimulation();
  skeleton->InitializeBoneRotation();
  //  arch::current_body->ReinitializeOffsetToSkeleton();
  //  arch::current_cloth->ReinitializeOffsetToSkeleton();
}

void after() {
  //  skeleton->UpdatePosition();
  //  return;
  //  skeleton->InitializeBoneRotation();
  //  current_body->ReinitializeOffsetToSkeleton();
  //  //  current_cloth->ReinitializeOffsetToSkeleton();
  skeleton->UpdatePosition();
  //  skeleton->AssembleSkeletonTransformation();
  //  current_body->ApplySkeletonTransformationToVertex();
  //      tets0->TetLimiting(500, tets0->X, tets0->X, true);
  //  current_cloth->ApplySkeletonTransformation(true);

  //  simulator::HandleObjectCollision();
  //  current_body->StrainLimiting(100, 50);
  //  for (int i = 0; i < 2; ++i) {
  //    simulator::HandleObjectCollision();
  //    current_cloth->TriangleLimiting(20, current_cloth->vertex_);
  //  }
  //  memset(current_cloth->velocity_, 0, sizeof(double) * current_cloth->vertex_num_ * 3);

  //  cloth->ApplySkeletonTransformation();
  skeleton->AssembleSkeletonTransformation();
  arch::current_body->ApplySkeletonTransformationToVertex();
  //      arch::current_body->StrainLimiting(10, 100);
  //  arch::current_cloth->ApplySkeletonTransformation(true);
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

int SkeletonRotator::HandleKeyPress(QKeyEvent *e) {
  using namespace global;
  int node = 8; //waist
  const double rotate_step = dj::Degree2Radian(5.0f);
  switch (e->key()) {
    case Qt::Key_1: {
      double * angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      before();
      angle[0] -= rotate_step;
      after();
      dj::Vec3d ang(angle);
      ang *= (180.0f / PI);
      P(ang);
    }
    break;
    case Qt::Key_2: {
      before();
      double * angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      angle[0] += rotate_step;
      after();
      dj::Vec3d ang(angle);
      ang *= (180.0f / PI);
      P(ang);
    }
    break;
    case Qt::Key_3: {
      before();
      double* angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      angle[1] -= rotate_step;
      after();
      dj::Vec3d ang(angle);
      ang *= (180.0f / PI);
      P(ang);
    }
    break;
    case Qt::Key_4: {
      before();
      double * angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      angle[1] += rotate_step;
      after();
      angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      dj::Vec3d ang(angle);
      ang *= (180.0f / PI);
      P(ang);
    }
    break;
    case Qt::Key_5: {
      before();
      double * angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      angle[2] -= rotate_step;
      after();
      dj::Vec3d ang(angle);
      ang *= (180.0f / PI);
      P(ang);
    }
    break;
    case Qt::Key_6: {
      before();
      double * angle = current_body->skeleton_->GetNode(node)->rotation_angle();
      angle[2] += rotate_step;
      after();
      dj::Vec3d ang(angle);
      ang *= (180.0f / PI);
      P(ang);
    }
    break;
    case Qt::Key_K:
      NextPose();
      break;
    case Qt::Key_J:
      PrevPose();
      break;
    case Qt::Key_S: {
//      const char file_name[] = DATA_DIRECTORY "armadillo_pose_samples.txt";
      std::string file_name = dj::Format("%s/skeleton_pose.txt", GetDataFolder());
      pose_sampler->SaveSkeletonPoses(file_name.c_str());
      L(file_name, " saved");
      break;
    }
    case Qt::Key_L: {
//      const char file_name[] = DATA_DIRECTORY "armadillo_pose_samples.txt";
      std::string file_name = dj::Format("%s/skeleton_pose.txt", GetDataFolder());
      pose_sampler->LoadSkeletonPoses(file_name.c_str());
      L(file_name, " loaded");
      break;
    }
    case Qt::Key_I :
      Init();
      break;
    default:
      return -1;
  }
  return 0;
}

void SkeletonRotator::Init() {
  delete pose_sampler;
  int joints[] = {2, 3, 4, 6, 7, 8, 14, 15, 17, 18, 12};
  int num = sizeof(joints) / sizeof(int);
  std::vector<int> moving_joints(&joints[0], &joints[0] + num);
  pose_sampler = new PoseSampler(current_body);
  char file_name[512]; // = DATA_DIRECTORY "armadillo_pose_samples.txt";
  sprintf(file_name, "%s/skeleto_pose_txt", GetDataFolder());
  min_pose_idx = conf.Get<int>("min_pose_idx");
  max_pose_idx = conf.Get<int>("max_pose_idx");
  pose_sampler->LoadSkeletonPoses(file_name);
  L(file_name, " loaded");
  pose_sampler->SetPose(min_pose_idx);
  arch::current_body->ApplySkeletonTransformationToVertex();
  current_body->inv_fem_->LoadPosition(current_body->X);

  if (0)
  {
//    char file_name[1024];
    std::vector<std::string> input;
    for (int i = 0; i < 100; ++i) {
//      sprintf(file_name, DATA_DIRECTORY "armadillo/pose/pose_%d.txt", i);
      input.push_back(dj::Format("%s/pose/pose_%d.txt", GetDataFolder(), i));
    }
//    pose_sampler->ParseLocalDeformations(input, DATA_DIRECTORY "armadillo/pose/local_deformation.txt", false);
    pose_sampler->ParseLocalDeformations(input,
//                                         DATA_DIRECTORY "armadillo/pose/local_deformation.txt",
                                         dj::Format("%s/pose/local_deformation.txt", GetDataFolder()).c_str(),
                                         false);
    exit(0);
  }
  //  pose_sampler->SamplePose(100, moving_joints);
}

int SkeletonRotator::Idle() {
  int current_pose_idx = pose_sampler->current_pose_;
  if (current_pose_idx < min_pose_idx || current_pose_idx > max_pose_idx) {
    return -1;
  }
  double sum = 0;
  for (int i = 0; i < 10; ++i) {
    sum += current_body->inv_fem_->Simulate(global::time_step);
  }
  sum /= 10;
  P(pose_sampler->current_pose_, sum, min_pose_idx, max_pose_idx);
  if (sum < 8e-5) {
    L("finished pose idx = ", current_pose_idx);
    std::stringstream stream;
//    stream << DATA_DIRECTORY << "armadillo/pose/pose_" << current_pose_idx << ".txt";
    std::string file_name = dj::Format("%s/pose/pose_%d.txt", current_pose_idx);
    pose_sampler->SavePoseResult(file_name.c_str());
    L(file_name + " saved.");
    current_pose_idx++;
    if (current_pose_idx <= max_pose_idx) {
      pose_sampler->SetPose(current_pose_idx);
      arch::current_body->ApplySkeletonTransformationToVertex();
      current_body->inv_fem_->LoadPosition(current_body->X);
    } else {
      std::stringstream stream;
      stream << "finish poses from " << min_pose_idx << " to " << max_pose_idx;
      L(stream.str());
      exit(0);
    }
  }
  return -1;
}

SkeletonRotator::SkeletonRotator() : InputHandler(10) {
}
