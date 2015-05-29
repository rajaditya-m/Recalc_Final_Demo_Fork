#include "pose_sampler.h"
#include <fstream>
#include "skeleton_node.h"
#include "skeleton.h"
#include "vector_io.h"
#include "random.h"
#include "tet.h"

PoseSampler::PoseSampler(Tet* body) : body_(body), skeleton_(body->skeleton_)
{
  node_num_ = int(skeleton_->nodes_.size());
  rest_pose_ = skeleton_->initial_pose_;
  current_pose_ = -1;
}

void PoseSampler::SetPose(std::vector<double> &pose)
{
  skeleton_->SetPose(pose);
}

void PoseSampler::SetPose(int pose_idx)
{
  ASSERT(pose_idx >= -1 && pose_idx < (int) poses_.size());
  current_pose_ = pose_idx;
  if (pose_idx == -1) {
    skeleton_->Reset();
  } else {
    skeleton_->SetPose(poses_[pose_idx]);
  }
}

void PoseSampler::SamplePose(int number_of_sample, std::vector<int> &moving_joint)
{
  poses_.resize(number_of_sample);
  Random rand;
  for (int i = 0;  i < number_of_sample; ++i) {
    poses_[i] = rest_pose_;
    for (int joint : moving_joint) {
      const double* min_angle = skeleton_->nodes_[joint].min_rotation_angle();
      const double* max_angle = skeleton_->nodes_[joint].max_rotation_angle();
      double* angle = &poses_[i][joint * 3];
      for (int j = 0; j < 3; ++j) {
        angle[j] = rand.GetRandom<double>(max_angle[j], min_angle[j]);
      }
    }
  }
}

void PoseSampler::SaveSkeletonPoses(const char *file_name)
{
  std::ofstream out(file_name);
  ASSERT(out.is_open());
  out << poses_.size() <<  " " << node_num_ << std::endl;
  for (int i = 0; i < (int) poses_.size(); ++i) {
    for (int j = 0; j < node_num_; ++j) {
      out << poses_[i][j * 3 + 0] << " "
          << poses_[i][j * 3 + 1] << " "
          << poses_[i][j * 3 + 2] << std::endl;
    }
  }
  out.close();
}

void PoseSampler::LoadSkeletonPoses(const char *file_name)
{
  std::ifstream in(file_name);
  int pose_num = 0;
  int node_count = 0;
  in >> pose_num >> node_count;
  ASSERT(node_count == node_num_);
  poses_.resize(pose_num);
  for (int i = 0; i < pose_num; ++i) {
    poses_[i] = rest_pose_;
    for (int j = 0; j < node_num_ * 3; ++j) {
      in >> poses_[i][j];
    }
  }
  in.close();
}

int PoseSampler::PoseNum() { return (int) poses_.size(); }

void PoseSampler::SetTargetPose(std::vector<double> &target_pose, int step)
{
  ASSERT(target_pose.size() == rest_pose_.size());
  move_step_.resize(target_pose.size());
  for (int i = 0; i < (int) move_step_.size(); ++i) {
    move_step_[i] = (target_pose[i] - rest_pose_[i]) / step;
  }
  current_step_ = 1;
  total_step_ = step;
}

bool PoseSampler::Move()
{
  if (current_step_ > total_step_) return false;
  std::vector<double> pose = rest_pose_;
  for (int i = 0; i < (int) pose.size(); ++i) {
    pose[i] = rest_pose_[i] + move_step_[i] * current_step_;
  }
  ++current_step_;
  SetPose(pose);
  return true;
}

void PoseSampler::Reset()
{
  skeleton_->Reset();
  body_->Reset();
}

void PoseSampler::SavePoseResult(const char* file_name)
{
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  out << body_->vertex_num_ << " "
      << node_num_ << " "
      << current_pose_ << " "
      << std::endl;
  for (int i = 0; i < node_num_; ++i) {
    out << skeleton_->nodes_[i].rotation_angle()[0] << " ";
    out << skeleton_->nodes_[i].rotation_angle()[1] << " ";
    out << skeleton_->nodes_[i].rotation_angle()[2] << " ";
    out << std::endl;
  }
  for (int v = 0; v < body_->vertex_num_; ++v) {
    out << body_->X[v * 3 + 0] << " ";
    out << body_->X[v * 3 + 1] << " ";
    out << body_->X[v * 3 + 2] << " ";
    out << std::endl;
  }
  out.close();
}

void PoseSampler::LoadPoseResult(const char *file_name)
{
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  int v_num, node_count;
  in >> v_num;
  in >> node_count;
  int pose_idx;
  in >> pose_idx;
  ASSERT(v_num == body_->vertex_num_, P(v_num, body_->vertex_num_));
  ASSERT(node_count == node_num_, P(node_count, node_num_));
  std::vector<double> angles(node_count * 3);
  for (int i = 0; i < node_num_; ++i) {
    in >> angles[i * 3 + 0];
    in >> angles[i * 3 + 1];
    in >> angles[i * 3 + 2];
  }
  skeleton_->SetPose(angles);

  for (int v = 0; v < body_->vertex_num_; ++v) {
    in >> body_->X[v * 3 + 0];
    in >> body_->X[v * 3 + 1];
    in >> body_->X[v * 3 + 2];
  }
  in.close();
}

void PoseSampler::ParseLocalDeformations(std::vector<std::string> &input_file, const char *output_file, bool binary_output)
{
  skeleton_->Reset();
  std::vector<double> rest_transformation = skeleton_->skeleton_transformation_;
  std::vector<std::vector<double> > local_u(input_file.size());
  for (int i = 0; i < (int) input_file.size(); ++i) {
    local_u[i].resize(body_->vertex_num_ * 3);
    LoadPoseResult(input_file[i].c_str());
    std::vector<double>& current_transformation = skeleton_->skeleton_transformation_;
    std::vector<double> rotation = current_transformation;
    for (int j = 0; j < (int) skeleton_->nodes_.size(); ++j) {
      dj::MulMatrixRightTransposed3x3<double>(&rest_transformation[j * 12],
          &current_transformation[j * 12],
          &rotation[j * 12]);
    }
    double* skeleton_transformation = &skeleton_->skeleton_transformation_[0];
    for (int v = 0; v < body_->vertex_num_; ++v) {
      double* pos = body_->X + v * 3;
      int attached_bone = body_->closest_bone_[v];
      double offset[3] = {
        pos[0] - skeleton_transformation[attached_bone * 12 + 9],
        pos[1] - skeleton_transformation[attached_bone * 12 + 10],
        pos[2] - skeleton_transformation[attached_bone * 12 + 11]
      };
      double local_offset[3];
      dj::MulMatrix3x3Vec<double>(
            (double (*)[3]) &rotation[attached_bone * 12],
          offset, local_offset);
      double local_pos[3] = {
        local_offset[0] + rest_transformation[attached_bone * 12 + 9],
        local_offset[1] + rest_transformation[attached_bone * 12 + 10],
        local_offset[2] + rest_transformation[attached_bone * 12 + 11]
      };
      local_u[i][v * 3 + 0] = local_pos[0] - body_->rest_pos_[v * 3 + 0];
      local_u[i][v * 3 + 1] = local_pos[1] - body_->rest_pos_[v * 3 + 1];
      local_u[i][v * 3 + 2] = local_pos[2] - body_->rest_pos_[v * 3 + 2];
//      ASSERT(dj::Abs(local_u[i][v *3 + 0]) < 1e-5);
//      ASSERT(dj::Abs(local_u[i][v *3 + 1]) < 1e-5);
//      ASSERT(dj::Abs(local_u[i][v *3 + 2]) < 1e-5);
//      body_->X[v * 3 + 0] = local_pos[0];
//      body_->X[v * 3 + 1] = local_pos[1];
//      body_->X[v * 3 + 2] = local_pos[2];
    }
  }
//  return;

  int pose_num = int(input_file.size());
  if (binary_output) {
    std::ofstream out(output_file, std::ios::binary);
    ASSERT(out.is_open());
    out.write((char*) &(body_->vertex_num_), sizeof(int) * 1);
    out.write((char*) (&pose_num), sizeof(int) * 1);
    for (int i = 0; i < pose_num; ++i) {
      out.write((char*) &local_u[i][0], sizeof(double) * body_->vertex_num_ * 3);
    }
    out.close();
  } else {
    std::ofstream out(output_file);
    ASSERT(out.is_open());
    out << body_->vertex_num_ << " ";
    out << pose_num << std::endl;
    for (int i = 0; i < pose_num; ++i) {
      for (int v = 0; v < body_->vertex_num_; ++v) {
        out << local_u[i][v * 3 + 0] << " ";
        out << local_u[i][v * 3 + 1] << " ";
        out << local_u[i][v * 3 + 2] << " ";
        out << std::endl;
      }
    }
    out.close();
  }
}
