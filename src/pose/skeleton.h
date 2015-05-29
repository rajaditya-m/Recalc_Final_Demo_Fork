#ifndef SKELETON_H_
#define SKELETON_H_
#include "vector_lib.h"
#include <vector>

class PoseSampler;
//extern PoseSampler* pose_sampler;
class SkeletonNode;
class Tet;
class Skeleton
{
  typedef std::pair<int, int> Bone;
public:
  int skeleton_point_number_;
  SkeletonNode* root_;
  std::vector<Bone> bones_;
  std::vector<double> initial_pose_;
  std::vector<double> bone_rotation_transpose_;
  std::vector<double> bone_rotation_;
  std::vector<double> skeleton_points_;
  std::vector<double> skeleton_transformation_;
  std::vector<dj::Vec3d> max_angle;
  std::vector<dj::Vec3d> min_angle;
  std::vector<SkeletonNode> nodes_;
  std::vector<SkeletonNode*> leaf_nodes_containter;
  int node_id2leaf_node_idx[65536];
  int selected_joint_;
  int selected_node_;

  Skeleton(const char *file_name, bool use_initial_posiiton = false);
  int Size()
  {
    return (int) nodes_.size();
  }
  void UpdatePosition();
  SkeletonNode* GetNode(int i)
  {
    return &nodes_[i];
  }

  void Reset();
  void SetPose(std::vector<double> &pose);
  void GetPose(std::vector<double> &pose);

  inline SkeletonNode *GetChildJoint(int bone_idx)
  {
    return GetChildJoint(bones_[bone_idx]);
  }
  void InitializeBoneRotation();

  SkeletonNode *GetChildJoint(const Skeleton::Bone& bone);
  SkeletonNode *GetParentJoint(const Skeleton::Bone& bone);

  inline SkeletonNode *GetParentJoint(int bone_idx)
  {
    return GetParentJoint(bones_[bone_idx]);
  }

  void ReorderBones(void);

  void ExportFakeBVH(const char* file_name);
  void ImportFakeBVH(const char* file_name, bool use_initial_posiiton = false);
  void InitSkeletonFakeBVHData(bool use_initial_posiiton);
  void Scale(double scale);
  void Translate(double dx, double dy, double dz);
  void Translate(double* offset);
  int GetBoneNum();
  double DistanceToBone(const double* pos, int bone_idx);
  void AttachedToBone(const double* pos, int bone_idx, double* initial_offset);
  int AttachedToSkeleton(const double* pos, double* initial_offset);

  void Save(const char* file_name);
  void Load(const char* file_name);

  void InverseKinematics(int moving_joint, double* target_position, const int max_iteration, double step_size);

  void AssembleSkeletonTransformation();
  void Init();
  void LoadAngleConstraints();
  void SaveAngle();

  void SetNewRoot(int root_idx);
  int SelectJoint( const double* start, const double* end);
  void MoveEffector(double dx, double dy, double dz);
  int SelectEffector(const dj::Vec3d& start, const dj::Vec3d& end);
  void HandleKeyPress(unsigned char key);
  void RotateY(double angle);
  void RotateX(double angle);
  void Render();
};

#endif // SKELETON_H
