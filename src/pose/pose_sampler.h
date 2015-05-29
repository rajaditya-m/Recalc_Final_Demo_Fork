#ifndef POSE_SAMPLER_H
#define POSE_SAMPLER_H
#include <vector>
#include <string>

class Skeleton;
class Tet;
class PoseSampler
{
public:
  PoseSampler(Tet* body);
  void SetPose(std::vector<double>& pose);
  void SetPose(int pose_idx);
  void SamplePose(int number_of_sample, std::vector<int>& moving_joint);
  void SaveSkeletonPoses(const char* file_name);
  void LoadSkeletonPoses(const char* file_name);
  int PoseNum();
  void SetTargetPose(std::vector<double>& target_pose, int step);
  bool Move();
  void Reset();

  void SavePoseResult(const char *file_name);
  void LoadPoseResult(const char *file_name);

  void ParseLocalDeformations(std::vector<std::string>& input_file, const char* outout_file, bool binary_output = true);

  int total_step_;
  int current_step_;
  int current_pose_;
  std::vector<double> move_step_;
  std::vector<double> rest_pose_;
  std::vector<std::vector<double> > poses_;
  int node_num_;
  Tet* body_;
  Skeleton* skeleton_;
};

#endif // POSE_SAMPLER_H
