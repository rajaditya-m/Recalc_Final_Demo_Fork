#include <iostream>
#include <stack>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include "open_gl_qt.h"
#include "opengl_helper.h"
#include "skeleton.h"
#include "skeleton_node.h"
#include "buffer_manager.h"
#include "tet.h"
#include "vector_lib.h"
#include "rainbow_color.h"
#include "global.h"
#include "print_macro.h"
#include "open_gl_qt.h"

#ifndef CURRENT_DIRECTORY
#define CURRENT_DIRECTORY STRINGIZE_TOKEN(DATA_DIRECTORY)
#endif

//PoseSampler* pose_sampler = NULL;
using dj::Mat3d;
using dj::Vec3d;
int skeleton_point_number = 19;
const int skeleton_edge_number = 18;
double skeleton_points[] = {
  0.983108, 0.031926, 0.043085,
  0.982609, 0.033623, -0.047298,
  0.983614, 0.030200, 0.134956,
  0.525866, -0.043298, -0.137139,
  0.065486, 0.039519, -0.186768,
  0.036165, -0.089296, -0.254352,
  0.510611, -0.004893, 0.204899,
  0.075242, 0.084120, 0.211848,
  0.047333, -0.034794, 0.261010,
  1.207314, 0.020588, 0.032412,
  1.447324, 0.009777, 0.026101,
  1.621424, -0.035893, 0.039145,
  1.747674, -0.055367, 0.058087,
  1.497586, -0.017654, -0.153136,
  1.296075, -0.016811, -0.350783,
  1.073056, -0.084916, -0.407875,
  1.468036, 0.009514, 0.226556,
  1.313862, 0.019558, 0.408056,
  1.063088, -0.029135, 0.471262
};


int skeleton_edges[] = {
  1, 2,
  1, 3,
  2, 4,
  4, 5,
  5, 6,
  3, 7,
  7, 8,
  8, 9,
  1, 10,
  10, 11,
  11, 12,
  12, 13,
  11, 14,
  14, 15,
  15, 16,
  11, 17,
  17, 18,
  18, 19
};

static const int kLeafNum = 19;
double end_effector_pos[kLeafNum * 3];
int selected_end_effector = -1;
double target_end_effector_pos[kLeafNum * 3];

double origin_angle[19 * 3];

Skeleton::Skeleton(const char* file_name, bool use_initial_posiiton) {
  if (file_name == NULL) {
    skeleton_point_number_ = skeleton_point_number;
    skeleton_points_ = std::vector<double>(skeleton_points, skeleton_points + skeleton_point_number_ * 3);
    skeleton_transformation_.resize(12 * skeleton_point_number_);
    nodes_.resize(skeleton_point_number_);
    Init();
  } else {
    ImportFakeBVH(file_name, use_initial_posiiton);
    //    ImportFakeBVH(file_name, true);
    //    for (int i = 0; i < int(nodes_.size()); ++i) {
    //      P(i, nodes_[i].name());
    //    }
    //    exit(0);
    //    for (auto n : nodes_) {
    //      P(n.name(), dj::Vec3d(n.world_pos()));
    //    }
  }
  ReorderBones();
  selected_node_ = 0;
  GetPose(initial_pose_);

  InitializeBoneRotation();
  bone_rotation_ = bone_rotation_transpose_;
  for (int i = 0; i < (int) bones_.size(); ++i) {
    dj::Transpose3x3<double>(&bone_rotation_[i * 9]);
  }
}

void Skeleton::UpdatePosition() {
  root_->ComputePosition();
}

void Skeleton::Reset() {
  SetPose(initial_pose_);
}

void Skeleton::SetPose(std::vector<double> &pose) {
  ASSERT(pose.size() == nodes_.size() * 3);
  for (int i = 0; i < (int) nodes_.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      nodes_[i].rotation_angle()[j] = pose[i * 3 + j];
    }
  }
  UpdatePosition();
  AssembleSkeletonTransformation();
}

void Skeleton::GetPose(std::vector<double> &pose) {
  pose.resize(nodes_.size() * 3);
  for (int i = 0; i < (int) nodes_.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      pose[i * 3 + j] = nodes_[i].rotation_angle()[j];
    }
  }
}

double Skeleton::DistanceToBone(const double *pos, int bone_idx) {
  SkeletonNode* joints[2] = {
    &nodes_[bones_[bone_idx].first],
    &nodes_[bones_[bone_idx].second]
  };
  double parameter;
  return dj::Point2SegmentDistance(joints[0]->world_pos(),
                                   joints[1]->world_pos(), pos, parameter);
}

void Skeleton::InitializeBoneRotation() {
  for (int bone_idx = 0; bone_idx < int(bones_.size()); ++bone_idx) {
    SkeletonNode* closest_node = GetChildJoint(bones_[bone_idx]);
    const double* parent_rotation = closest_node->parent()->rotation_matrix();
    const double (*rotation_matrix)[3] = closest_node->rotation_angle_matrix_;
    double* mat = &bone_rotation_transpose_[bone_idx * 9];
    dj::MulMatrix3x3<double>(parent_rotation, rotation_matrix, mat);
//    memcpy(&bone_rotation_[bone_idx * 9], mat, sizeof(double) * 9);
    dj::Transpose3x3<double>(mat);
  }
}

SkeletonNode* Skeleton::GetChildJoint(const Skeleton::Bone &bone) {
  SkeletonNode* joints[2] = {
    &nodes_[bone.first],
    &nodes_[bone.second]
  };
  if (joints[0]->parent() == joints[1]) {
    return joints[0];
  } else {
    return joints[1];
  }
}

SkeletonNode *Skeleton::GetParentJoint(const Skeleton::Bone &bone) {
  SkeletonNode* joints[2] = {
    &nodes_[bone.first],
    &nodes_[bone.second]
  };
  if (joints[0]->parent() != joints[1]) {
    return joints[0];
  } else {
    return joints[1];
  }
}

void Skeleton::ReorderBones() {
  for (int i = 0; i < int(bones_.size()); ++i) {
    SkeletonNode* n[2] = {
      &nodes_[bones_[i].first],
      &nodes_[bones_[i].second]
    };
    if (n[1]->parent() != n[0]) {
      dj::Swap(bones_[i].first, bones_[i].second);
    }
  }
}

void Skeleton::ExportFakeBVH(const char *file_name) {
  std::ofstream out(file_name);
  if (!out.is_open()) {
    std::cerr << "Skeleton::ExportToFakeBVH() => failed to open file " << file_name << std::endl;
    exit(0);
  }
  if (root_ != &nodes_[0]) {
    SetNewRoot(0);
  }
  root_->ComputePosition();
  out << "NODE_NUM " << nodes_.size() << std::endl;
  root_->ExportFakeBVH(out, "");
  out.close();
}

void Skeleton::ImportFakeBVH(const char *file_name, bool use_initial_posiiton) {
  std::ifstream in(file_name);
  if (!in.is_open()) {
    std::cerr << "Skeleton::ImportToFakeBVH() => failed to open file " << file_name << std::endl;
    exit(0);
  }
  std::string token;
  in >> token;
  int node_num;
  in >> node_num;
  nodes_.clear();
  nodes_.reserve(node_num + 30);
  //  P(node_num);
  in >> token;
  if (token == "ROOT")	{
    nodes_.emplace_back();
    root_ = &nodes_.back();
    root_->ImportFakeBVH(in,  NULL, nodes_);
  } else {
    std::cerr << "Skeleton::ImportFakeBVH() => Failed to read ROOT token" << std::endl;
    exit(0);
  }

  skeleton_point_number_ = (int) nodes_.size();
  InitSkeletonFakeBVHData(use_initial_posiiton);
  for (int i = 0; i < (int) nodes_.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT(nodes_[i].min_rotation_angle()[j] <= nodes_[i].max_rotation_angle()[j]);
    }
  }
}

void Skeleton::InitSkeletonFakeBVHData(bool use_initial_posiiton) {
  // Build skeleton hierachy
  std::vector<Vec3d> angles;
  for (int i = 0; i < int(nodes_.size()); ++i) {
    angles.emplace_back(nodes_[i].rotation_angle());
  }
  skeleton_transformation_.resize(nodes_.size() * 12);
  std::queue<SkeletonNode*> queue;
  for (auto node : root_->children()) {
    queue.push(node);
  }
  while (!queue.empty()) {
    SkeletonNode* current_node = queue.front();
    queue.pop();
    const std::vector<SkeletonNode*> children = current_node->children();
    for (unsigned int i = 0; i < children.size(); ++i) {
      queue.push(children[i]);
    }
    SkeletonNode* parent = current_node->parent();
    ASSERT(parent != NULL);
    Mat3d parent_rotation_matrix((double*) parent->rotation_matrix());
    Mat3d parent_rotation_matrix_transpose = parent_rotation_matrix;
    parent_rotation_matrix_transpose.Transpose();
    Vec3d my_pos(current_node->world_pos());
    Vec3d parent_pos(current_node->parent()->world_pos());
    Vec3d world_offset = my_pos - parent_pos;
    // Offset in parent's local frame
    dj::Vec3d local_offset;//= parent_rotation_matrix_transpose * world_offset;
    dj::MulMatrix3x3Vec((double (*)[3]) &parent_rotation_matrix_transpose[0][0], &world_offset[0], &local_offset[0]);
    Vec3d initial_offset(current_node->offset());
    dj::Vec3d projection((double*) local_offset());
    projection[current_node->order()[0]] = 0;
    dj::Vec3d cross_product(0, 0, 0);
    dj::Cross3<double>(initial_offset(), projection(), cross_product());
    double theta = atan2(cross_product[current_node->order()[0]],
                        initial_offset * projection);
    //------------------------------------------------------------------------------
    Mat3d second_rotation_matrix;
    (*(current_node->Rotate_[current_node->order()[0]]))(theta, (double (*)[3]) second_rotation_matrix());
    dj::Cross3(projection(), local_offset(), cross_product());
    Vec3d third_axis = second_rotation_matrix.Col(current_node->order()[1]);
    double sign = dj::Sign(cross_product * third_axis);
    double phi = atan2(sign * cross_product.Magnitude(), projection * local_offset);
    current_node->set_rotation_angle(theta, phi, 0);
    //    current_node->set_rotation_angle(0, 0, 0);
    current_node->AngleChanged();
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // debug
    {
      //      dj::Mat3d ro(current_node->rotation_angle_matrix_);
      //      dj::Mat3d pro(parent->rotation_matrix());
      //      auto diff = pro * ro * initial_offset - world_offset;
      ////      if (!(diff.Norm2() < 1e-6)) {
      //        P(diff, current_node->name());
      //        P(ro);
      //        P(pro);
      //        P(dj::Vec3d(current_node->rotation_angle()));
      ////      }
      //      MY_ASSERT(diff.Norm2() < 1e-6);
      //      KK;
    }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    double reference_frame[3][3] = {
      {world_offset[0], world_offset[1], world_offset[2]},
      {0, 0, 0},
      {0, 0, 0}
    };
    // Compute reference frame in world coordinate as row vectors
    {
      if (children.size() == 0)  {
        // The leaf node can choose any other two axis for reference frame
        reference_frame[2][0] = 0;
        reference_frame[2][1] = 1;
        reference_frame[2][2] = 0;
      } else {
        const double* child_pos = children[0]->world_pos();
        reference_frame[2][0] = child_pos[0] - my_pos[0];
        reference_frame[2][1] = child_pos[1] - my_pos[1];
        reference_frame[2][2] = child_pos[2] - my_pos[2];
      }
      dj::Normalize3(reference_frame[0]);
      dj::Cross3(reference_frame[0], reference_frame[2], reference_frame[1]);
      dj::Normalize3(reference_frame[1]);
      if (dj::Eq(reference_frame[1][0], 0.0, (double) EPSILON)
          && dj::Eq(reference_frame[1][1], 0.0, (double) EPSILON)
          && dj::Eq(reference_frame[1][2], 0.0, (double) EPSILON)) {
        reference_frame[1][0] = 0;
        reference_frame[1][1] = 1;
        reference_frame[1][2] = 0;
      }
      dj::Cross3(reference_frame[0], reference_frame[1], reference_frame[2]);
    }
    dj::Transpose<double, 3, 3>(reference_frame, current_node->reference_frame());
    dj::Transpose<double, 3, 3>(reference_frame, current_node->rotation_matrix());
    Mat3d local_rotation_matrix;
    dj::MulMatrix3x3<double>(parent_rotation_matrix(), current_node->rotation_angle_matrix_,
                            local_rotation_matrix());
    local_rotation_matrix.Transpose();
    dj::MulMatrix3x3<double>(local_rotation_matrix(), current_node->rotation_matrix(),
                            current_node->fixed_rotation_matrix_);

  }


  if (!use_initial_posiiton) {
    for (int i = 0; i < int(nodes_.size()); ++i) {
      nodes_[i].set_rotation_angle(angles[i]());
    }
  }
  root_->ComputePosition();


  bones_.clear();
  for (int i = 0; i < (int) nodes_.size(); ++i) {
    if (nodes_[i].parent()) {
      int parent_idx = nodes_[i].parent() - &nodes_[0];
      if (i < parent_idx) {
        bones_.push_back(Bone(i, parent_idx));
      } else {
        bones_.push_back(Bone(parent_idx, i));
      }
    }
  }
  bone_rotation_transpose_.resize(bones_.size() * 9);
}

void Skeleton::Scale(double scale) {
  std::vector<dj::Vec3d> new_pos_vec;
  new_pos_vec.emplace_back(root_->world_pos());
  for (int i = 0; i < int(nodes_.size()); ++i) {
    if (&nodes_[i] != root_) {
      dj::Vec3d offset(nodes_[i].offset());
      offset *= scale;
      nodes_[i].set_offset(offset());
      dj::Vec3d my_pos(nodes_[i].initial_world_pos());
      dj::Vec3d root_pos(root_->initial_world_pos());
      dj::Vec3d new_pos = root_pos + (my_pos - root_pos) * scale;
      new_pos_vec.push_back(new_pos);
      //      nodes_[i].set_initial_world_pos(new_pos());
      //      nodes_[i].set_world_pos(new_pos());
    }
  }

  for (int i = 0; i < int(nodes_.size()); ++i) {
    nodes_[i].set_initial_world_pos(new_pos_vec[i]());
    nodes_[i].set_world_pos(new_pos_vec[i]());
  }
}

void Skeleton::Translate(double dx, double dy, double dz) {
  for (int i = 0; i < int(nodes_.size()); ++i) {
    dj::Vec3d my_pos(nodes_[i].initial_world_pos());
    my_pos[0] += dx;
    my_pos[1] += dy;
    my_pos[2] += dz;
    nodes_[i].set_initial_world_pos(my_pos());
    nodes_[i].set_world_pos(my_pos());
  }
}

void Skeleton::Translate(double *offset) {
  Translate(offset[0], offset[1], offset[2]);
}

int Skeleton::GetBoneNum() { return (int) bones_.size(); }

void Skeleton::AttachedToBone(const double *pos, int bone_idx, double *initial_offset) {
  // Compute initial offset
  //  if (!(bone_idx >= 0 && bone_idx < int(bones_.size()))) {
  //    P(bone_idx);
  //  }
  //  MY_ASSERT(bone_idx >= 0 && bone_idx < int(bones_.size()))i;
  const double* parent_joint_pos = GetParentJoint(bone_idx)->world_pos();
  double offset[3] = {
    pos[0] - parent_joint_pos[0],
    pos[1] - parent_joint_pos[1],
    pos[2] - parent_joint_pos[2]
  };
  dj::MulMatrix3x3Vec<double>((double (*)[3]) &bone_rotation_transpose_[bone_idx * 9], offset, initial_offset);
}

///   \brief Skeleton::AttachedToSkeleton find closed bone of given position
///   \param pos
///   \param initial_offset
///   \return the index of closest bone
int Skeleton::AttachedToSkeleton(const double *pos, double *initial_offset) {
  // Find the closest bone
  double min_distance = 1e10;
  int closest_bone = -1;
  for (int i = 0; i < int(bones_.size()); ++i) {
    const double* v1 = nodes_[bones_[i].first].world_pos();
    const double* v2 = nodes_[bones_[i].second].world_pos();
    double parameter;
    double distance = dj::Point2SegmentDistance(v1, v2, pos, parameter);
    if (distance < min_distance) {
      min_distance = distance;
      closest_bone = i;
    }
  }
  // Compute initial offset
  AttachedToBone(pos, closest_bone, initial_offset);
  return closest_bone;
}

void Skeleton::Save(const char *file_name) {
  std::ofstream out(file_name);
  if (!out.is_open()) {
    std::cerr << "failed to open file " << file_name << std::endl;
    exit(0);
  }
  for (int i = 0; i < skeleton_point_number_; ++i) {
    double angle_in_degree[3] = {
      dj::Radian2Degree(nodes_[i].rotation_angle()[0]),
      dj::Radian2Degree(nodes_[i].rotation_angle()[1]),
      dj::Radian2Degree(nodes_[i].rotation_angle()[2]),
    };
    out << angle_in_degree[0] << " "
        << angle_in_degree[1] << " "
        << angle_in_degree[2] << " "
        << std::endl;
  }
  out.close();
}

void Skeleton::Load(const char *file_name) {
  std::ifstream in(file_name);
  if (!in.is_open()) {
    std::cerr << "failed to open file " << file_name << std::endl;
    exit(0);
  }
  for (int i = 0; i < skeleton_point_number_; ++i) {
    double angle[3];
    in >> angle[0]
       >> angle[1]
       >> angle[2];
    angle[0] = dj::Degree2Radian(angle[0]);
    angle[1] = dj::Degree2Radian(angle[1]);
    angle[2] = dj::Degree2Radian(angle[2]);
    nodes_[i].set_rotation_angle(angle);
  }
  in.close();
  root_->ComputePosition();
}

void Skeleton::InverseKinematics(int moving_joint, double* target_position,
                                 const int max_iteration, double step_size) {
  //#define ROOT_MOVABLE
  const int kLeafNum = 1;
  root_->ComputePosition();
  const int skeleton_point_number = nodes_[moving_joint].GetDepth();
  const int row_num = kLeafNum * 3;
  const int col_num = skeleton_point_number * 3;
  BufferManager buffer;
  double* theta = buffer.Malloc<double>(col_num);
  double* tmp_theta = buffer.Malloc<double>(col_num);
  double* velocity = buffer.Malloc<double>(row_num);
  double* gradient = buffer.Malloc<double>(col_num);
  double* prev_gradient = buffer.Malloc<double>(col_num);

  std::vector<SkeletonNode*> skeleton_branch;
  std::unordered_map<SkeletonNode*, SkeletonNode*> children_map;
  skeleton_branch.reserve(skeleton_point_number);
  SkeletonNode* current_node = &nodes_[moving_joint];
  SkeletonNode* child = NULL;
  while (current_node) {
    skeleton_branch.push_back(current_node);
    children_map[current_node] = child;
    child = current_node;
    current_node = current_node->parent();
  }
  for (int i = 0; i < skeleton_point_number; ++i) {
    if (skeleton_branch[i]->parent()) {
      const double* angle = skeleton_branch[i]->rotation_angle();
      theta[i * 3 + 0] = angle[0];
      theta[i * 3 + 1] = angle[1];
      theta[i * 3 + 2] = angle[2];
    } else {
      theta[i * 3 + 0] = skeleton_branch[i]->world_pos()[0];
      theta[i * 3 + 1] = skeleton_branch[i]->world_pos()[1];
      theta[i * 3 + 2] = skeleton_branch[i]->world_pos()[2];
    }
  }
  int iteration = 0;
  double obj = 0;

  auto ComputeDistance2Target = [&]() -> double {
    double obj = 0;
    for (int i = 0; i < kLeafNum; ++i) {
      velocity[i * 3 + 0] = -target_position[0] + nodes_[moving_joint].world_pos()[0];
      velocity[i * 3 + 1] = -target_position[1] + nodes_[moving_joint].world_pos()[1];
      velocity[i * 3 + 2] = -target_position[2] + nodes_[moving_joint].world_pos()[2];
      obj += dj::Square(velocity[i * 3 + 0]);
      obj += dj::Square(velocity[i * 3 + 1]);
      obj += dj::Square(velocity[i * 3 + 2]);
    }
    return obj;
  };
  auto UpdatePosition = [&](double * solution) {
    for (int i = 0; i < skeleton_point_number; ++i) {
      if (skeleton_branch[i]->parent()) {
        skeleton_branch[i]->set_rotation_angle(solution + i * 3);
      } else {
        skeleton_branch[i]->set_world_pos(solution + i * 3);
      }
    }

    root_->ComputePosition();
  };

  auto AccumulateGradient = [&](int root_idx, SkeletonNode * n) {
    int start = root_idx * 3;
    gradient[start + 0] += dj::Dot3(n->world_pos_derivative_[0], velocity);
    gradient[start + 1] += dj::Dot3(n->world_pos_derivative_[1], velocity);
    gradient[start + 2] += dj::Dot3(n->world_pos_derivative_[2], velocity);
  };

  auto ComputeGradient = [&]() {
    memset(gradient, 0, sizeof(double) * col_num);
    for (int j = 0; j < skeleton_point_number; ++j) {
      SkeletonNode* rotating_node = skeleton_branch[j];
      if (rotating_node->parent()) {
        // Non root node changes angle
        rotating_node->ComputeDerivative(true);
        std::queue<SkeletonNode*> queue;
        if (children_map[rotating_node] == NULL) {
          AccumulateGradient(j, rotating_node);
        } else {
          queue.push(children_map[rotating_node]);
        }
        while (!queue.empty()) {
          SkeletonNode* current_node = queue.front();
          queue.pop();
          current_node->ComputeDerivative(false);
          if (children_map[current_node] == NULL) {
            AccumulateGradient(j, current_node);
          } else {
            queue.push(children_map[current_node]);
          }
        }
      } else {
        // Root node changes position
#ifdef ROOT_MOVABLE
        int start = j * 3;
        for (int i = 0; i < kLeafNum; ++i) {
          int leaf_idx = i * 3;
          gradient[start + 0] += velocity[leaf_idx + 0];
          gradient[start + 1] += velocity[leaf_idx + 1];
          gradient[start + 2] += velocity[leaf_idx + 2];
        }
#endif
      }
    }
  };

  auto UpdateSolutionVector = [&](const double * direction, const double step_size, double * x) {
    for (int i = 0; i < col_num; ++i) {
      x[i] -= step_size * direction[i];
    }
  };

  double average_step_size = 0;
  int count = 0;
  auto BacktrackingLineSearch = [&](const double * gradient, const double * direction, double * x) {
    double gradient_dot_direction = 0;
    for (int i = 0; i < col_num; ++i) {
      gradient_dot_direction -= gradient[i] * direction[i];
    }
    double beta = 0.9;
    double initial_obj = ComputeDistance2Target();
    while (true) {
      memcpy(tmp_theta, x, sizeof(double) * col_num);
      UpdateSolutionVector(direction, step_size, tmp_theta);
      UpdatePosition(tmp_theta);
      double obj = ComputeDistance2Target();
      //      P(step_size);
      //      P(obj, initial_obj);// + alpha * step_size * gradient_dot_direction);
      if (obj <= initial_obj) { // + alpha * step_size * gradient_dot_direction) {
        break;
      } else {
        if (step_size < 0.00001) {
          break;
        }
        step_size *= beta;
      }
    }
    average_step_size += step_size;
    count++;
    memcpy(x, tmp_theta, sizeof(double) * col_num);
  };
  memset(prev_gradient, 0, sizeof(double) * col_num);
  while (iteration < max_iteration) {
    // Compute velocity that moves leaf nodes to target position
    iteration++;
    obj = ComputeDistance2Target();
    //      P(iteration, obj);
    // Compute gradient
    ComputeGradient();
    // Update theta for each node
    const double kAlpha = 0.6;
    double gradient_norm = 0;
    for (int i = 0; i < col_num; ++i) {
      prev_gradient[i] = (1 - kAlpha) * prev_gradient[i] + kAlpha * gradient[i];
      gradient_norm += dj::Square(prev_gradient[i]);
    }
    BacktrackingLineSearch(prev_gradient, prev_gradient, theta);
    obj /= col_num;
    if (obj < 1e-6 || gradient_norm < 1e-6) {
      //        P(obj);
      break;
    }
  }
  //    if (count > 0) P(average_step_size / count);
  if (iteration >= max_iteration) {
    //      std::cerr << CURRENT_LINE << "=>" << "maximum iteration reached" << std::endl;
  } else {
    //      P(iteration);
  }
}

void Skeleton::AssembleSkeletonTransformation() {
  for (int i = 0; i < int(bones_.size()); ++i) {
    SkeletonNode* current_joint = GetChildJoint(bones_[i]);
    SkeletonNode* parent = current_joint->parent();
    double* parent_rotation = (double*) parent->rotation_matrix();
    double* local_ratation = (double*) current_joint->rotation_angle_matrix_;
    dj::MulMatrix3x3<double>(parent_rotation, local_ratation, &skeleton_transformation_[i * 12]);
    skeleton_transformation_[i * 12 +  9] = parent->world_pos()[0];
    skeleton_transformation_[i * 12 + 10] = parent->world_pos()[1];
    skeleton_transformation_[i * 12 + 11] = parent->world_pos()[2];
  }
}

void Skeleton::Init() {
  max_angle.resize(skeleton_point_number_);
  min_angle.resize(skeleton_point_number_);
  selected_joint_ = -1;

  // Initialize tree structure
  root_ = &nodes_[0];
  {
    nodes_[1 - 1].set_parent(NULL);
    nodes_[2 - 1].set_parent(&nodes_[1 - 1]);
    nodes_[3 - 1].set_parent(&nodes_[1 - 1]);
    nodes_[4 - 1].set_parent(&nodes_[2 - 1]);
    nodes_[5 - 1].set_parent(&nodes_[4 - 1]);
    nodes_[6 - 1].set_parent(&nodes_[5 - 1]);
    nodes_[7 - 1].set_parent(&nodes_[3 - 1]);
    nodes_[8 - 1].set_parent(&nodes_[7 - 1]);
    nodes_[9 - 1].set_parent(&nodes_[8 - 1]);
    nodes_[10 - 1].set_parent(&nodes_[1 - 1]);
    nodes_[11 - 1].set_parent(&nodes_[10 - 1]);
    nodes_[12 - 1].set_parent(&nodes_[11 - 1]);
    nodes_[13 - 1].set_parent(&nodes_[12 - 1]);
    nodes_[14 - 1].set_parent(&nodes_[11 - 1]);
    nodes_[15 - 1].set_parent(&nodes_[14 - 1]);
    nodes_[16 - 1].set_parent(&nodes_[15 - 1]);
    nodes_[17 - 1].set_parent(&nodes_[11 - 1]);
    nodes_[17 - 1].set_parent(&nodes_[11 - 1]);
    nodes_[18 - 1].set_parent(&nodes_[17 - 1]);
    nodes_[19 - 1].set_parent(&nodes_[18 - 1]);
  }
  for (int i = 0; i < skeleton_point_number_; ++i) {
    nodes_[i].set_id(i);
    nodes_[i].set_world_pos(&skeleton_points_[i * 3]);
    nodes_[i].set_initial_world_pos(&skeleton_points_[i * 3]);
    if (nodes_[i].parent()) {
      nodes_[i].parent()->AddChild(&nodes_[i]);
    }
  }
  for (int i = 0; i < skeleton_point_number_; ++i) {
    if (nodes_[i].parent()) {
      double length = dj::Distance3(nodes_[i].world_pos(), nodes_[i].parent()->world_pos());
      nodes_[i].set_length(length);
    }
  }

  // Specify rotation order
  {
    int i;
    i = 1; nodes_[i - 1].set_offset(0, 0, 0);
    i = 2; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 3; nodes_[i - 1].set_offset(-nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 4; nodes_[i - 1].set_offset(0, 0, -nodes_[i - 1].length());
    nodes_[i - 1].set_order(0, 1, 2);
    i = 5; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 6; nodes_[i - 1].set_offset(0, 0, -nodes_[i - 1].length());
    nodes_[i - 1].set_order(0, 1, 2);
    i = 7; nodes_[i - 1].set_offset(0, 0, -nodes_[i - 1].length());
    nodes_[i - 1].set_order(0, 1, 2);
    i = 8; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 9; nodes_[i - 1].set_offset(0, 0, -nodes_[i - 1].length());
    nodes_[i - 1].set_order(0, 1, 2);
    i = 10; nodes_[i - 1].set_offset(0, nodes_[i - 1].length(), 0);
    nodes_[i - 1].set_order(0, 2, 1);
    i = 11; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 12; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 13; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 14; nodes_[i - 1].set_offset(0, nodes_[i - 1].length(), 0);
    nodes_[i - 1].set_order(0, 2, 1);
    i = 15; nodes_[i - 1].set_offset(0, 0, -nodes_[i - 1].length());
    nodes_[i - 1].set_order(0, 1, 2);
    i = 16; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
    i = 17; nodes_[i - 1].set_offset(0, -nodes_[i - 1].length(), 0);
    nodes_[i - 1].set_order(0, 2, 1);
    i = 18; nodes_[i - 1].set_offset(0, 0, -nodes_[i - 1].length());
    nodes_[i - 1].set_order(0, 1, 2);
    i = 19; nodes_[i - 1].set_offset(nodes_[i - 1].length(), 0, 0);
    nodes_[i - 1].set_order(1, 2, 0);
  }

  // Build skeleton hierachy
  std::queue<SkeletonNode*> queue;
  queue.push(&nodes_[2 - 1]);
  queue.push(&nodes_[3 - 1]);
  queue.push(&nodes_[10 - 1]);
  while (!queue.empty()) {
    SkeletonNode* current_node = queue.front();
    queue.pop();
    const std::vector<SkeletonNode*> children = current_node->children();
    for (unsigned int i = 0; i < children.size(); ++i) {
      queue.push(children[i]);
    }
    SkeletonNode* parent = current_node->parent();
    ASSERT(parent != NULL);
    Mat3d parent_rotation_matrix((double*) parent->rotation_matrix());
    Mat3d parent_rotation_matrix_transpose = parent_rotation_matrix;
    parent_rotation_matrix_transpose.Transpose();
    Vec3d my_pos(current_node->world_pos());
    Vec3d parent_pos(current_node->parent()->world_pos());
    Vec3d world_offset = my_pos - parent_pos;
    // Offset in parent's local frame
    dj::Vec3d local_offset;//= parent_rotation_matrix_transpose * world_offset;
    dj::MulMatrix3x3Vec((double (*)[3]) &parent_rotation_matrix_transpose[0][0], &world_offset[0], &local_offset[0]);
    Vec3d initial_offset(current_node->offset());
    dj::Vec3d projection((double*) local_offset());
    projection[current_node->order()[0]] = 0;
    dj::Vec3d cross_product(0, 0, 0);
    dj::Cross3<double>(initial_offset(), projection(), cross_product());
    double theta = atan2(cross_product[current_node->order()[0]],
                        initial_offset * projection);
    Mat3d second_rotation_matrix;
    (*(current_node->Rotate_[current_node->order()[0]]))(theta, (double (*)[3]) second_rotation_matrix());
    dj::Cross3(projection(), local_offset(), cross_product());
    Vec3d third_axis = second_rotation_matrix.Col(current_node->order()[1]);
    double sign = dj::Sign(cross_product * third_axis);
    double phi = atan2(sign * cross_product.Magnitude(), projection * local_offset);
    current_node->set_rotation_angle(theta, phi, 0);
    //    current_node->set_rotation_angle(0, 0, 0);
    current_node->AngleChanged();
    double reference_frame[3][3] = {
      {world_offset[0], world_offset[1], world_offset[2]},
      {0, 0, 0},
      {0, 0, 0}
    };
    // Compute reference frame in world coordinate as row vectors
    {
      if (children.size() == 0)  {
        // The leaf node can choose any other two axis for reference frame
        reference_frame[2][0] = 0;
        reference_frame[2][1] = 1;
        reference_frame[2][2] = 0;
      } else {
        const double* child_pos = children[0]->world_pos();
        reference_frame[2][0] = child_pos[0] - my_pos[0];
        reference_frame[2][1] = child_pos[1] - my_pos[1];
        reference_frame[2][2] = child_pos[2] - my_pos[2];
      }
      dj::Normalize3(reference_frame[0]);
      dj::Cross3(reference_frame[0], reference_frame[2], reference_frame[1]);
      dj::Normalize3(reference_frame[1]);
      if (dj::Eq(reference_frame[1][0], 0.0, (double) EPSILON)
          && dj::Eq(reference_frame[1][1], 0.0, (double) EPSILON)
          && dj::Eq(reference_frame[1][2], 0.0, (double) EPSILON)) {
        reference_frame[1][0] = 0;
        reference_frame[1][1] = 1;
        reference_frame[1][2] = 0;
      }
      dj::Cross3(reference_frame[0], reference_frame[1], reference_frame[2]);
    }
    dj::Transpose<double, 3, 3>(reference_frame, current_node->reference_frame());
    dj::Transpose<double, 3, 3>(reference_frame, current_node->rotation_matrix());
    Mat3d local_rotation_matrix;
    dj::MulMatrix3x3<double>(parent_rotation_matrix(), current_node->rotation_angle_matrix_,
                            local_rotation_matrix());
    local_rotation_matrix.Transpose();
    dj::MulMatrix3x3<double>(local_rotation_matrix(), current_node->rotation_matrix(),
                            current_node->fixed_rotation_matrix_);

  }
  LoadAngleConstraints();
  for (unsigned int i = 0; i < nodes_.size(); ++i) {
    const double* angle = nodes_[i].rotation_angle();
    origin_angle[i * 3 + 0] = angle[0];
    origin_angle[i * 3 + 1] = angle[1];
    origin_angle[i * 3 + 2] = angle[2];
    //    P(i);
    //    P(Vec3d(angle));
    //    P(Vec3d(min_angle[i]));
    //    P(Vec3d(max_angle[i]));
    ASSERT(angle[0] >= min_angle[i][0] && angle[1] >= min_angle[i][1] && angle[2] >= min_angle[i][2]);
    ASSERT(angle[0] <= max_angle[i][0] + EPSILON
           && angle[1] <= max_angle[i][1] + EPSILON
           && angle[2] <= max_angle[i][2] + EPSILON);
    nodes_[i].SetRotationRange(min_angle[i](), max_angle[i]());
  }
  root_->ComputePosition();

  for (int i = 0; i < (int) nodes_.size(); ++i) {
    if (nodes_[i].parent()) {
      int parent_idx = nodes_[i].parent() - &nodes_[0];
      if (i < parent_idx) {
        bones_.push_back(Bone(i, parent_idx));
      } else {
        bones_.push_back(Bone(parent_idx, i));
      }
    }
  }
  bone_rotation_transpose_.resize(bones_.size() * 9);
}

void Skeleton::LoadAngleConstraints() {
  ASSERT(false);
  //  std::ifstream in(CURRENT_DIRECTORY "new_angle_limit.txt");
  //  //  std::ifstream in(CURRENT_DIRECTORY "angles.txt");
  //  if (!in.is_open()) {
  //    std::cerr << CURRENT_LINE << " => failed to open input file" << std::endl;
  //    return;
  //  }
  //  for (int i = 0; i < skeleton_point_number_; ++i) {
  //#if 0
  //    in >> max_angle[i][0];
  //    in >> max_angle[i][1];
  //    in >> max_angle[i][2];

  //    in >> min_angle[i][0];
  //    in >> min_angle[i][1];
  //    in >> min_angle[i][2];
  //#else
  //    in >> max_angle[i][0];
  //    in >> min_angle[i][0];
  //    in >> max_angle[i][1];
  //    in >> min_angle[i][1];
  //    in >> max_angle[i][2];
  //    in >> min_angle[i][2];

  //    max_angle[i][0] = dj::Degree2Radian(max_angle[i][0]);
  //    max_angle[i][1] = dj::Degree2Radian(max_angle[i][1]);
  //    max_angle[i][2] = dj::Degree2Radian(max_angle[i][2]);

  //    min_angle[i][0] = dj::Degree2Radian(min_angle[i][0]);
  //    min_angle[i][1] = dj::Degree2Radian(min_angle[i][1]);
  //    min_angle[i][2] = dj::Degree2Radian(min_angle[i][2]);
  //#endif
  //  }
  //  in.close();
}

void Skeleton::SaveAngle() {
  std::ofstream out("angles.txt");
  if (!out.is_open()) {
    std::cerr << CURRENT_LINE << " => failed to open output file" << std::endl;
    return;
  }
  for (int i = 0; i < skeleton_point_number_; ++i) {
    out << max_angle[i][0] << "\t";
    out << max_angle[i][1] << "\t";
    out << max_angle[i][2] << "\t";
    out << std::endl;
    out << min_angle[i][0] << "\t";
    out << min_angle[i][1] << "\t";
    out << min_angle[i][2] << "\t";
    out << std::endl;
  }
  out.close();
}

void Skeleton::SetNewRoot(int root_idx) {
  //  if (root_ == &nodes_[root_idx]) {
  //    return;
  //  }
  root_->ComputePosition();
  root_ = &nodes_[root_idx];
  std::stack<SkeletonNode*> stack;
  SkeletonNode* current_node = root_;
  stack.push(root_);
  while (current_node->parent()) {
    stack.push(current_node->parent());
    current_node = current_node->parent();
  }
  root_->set_parent(NULL);
  while (stack.size() > 1) {
    current_node = stack.top();
    stack.pop();
    SkeletonNode* new_parent = stack.top();
    new_parent->children().push_back(current_node);

    for (auto iter = current_node->children().begin(); true; ++iter) {
      if (iter == current_node->children().end()) {
        ASSERT(false);
      }
      if (*iter == new_parent) {
        current_node->children().erase(iter);
        break;
      }
    }
    current_node->set_parent(new_parent);

    dj::Vec3i order(new_parent->order());
    std::reverse(order.begin(), order.end());

    dj::Vec3d rotation_angle(new_parent->rotation_angle());
    std::reverse(rotation_angle.begin(), rotation_angle.end());
    rotation_angle *= -1.0;

    // Max angles become min angles
    dj::Vec3d min_angle(new_parent->max_rotation_angle());
    std::reverse(min_angle.begin(), min_angle.end());
    min_angle *= -1.0;

    // Min angles become max angles
    dj::Vec3d max_angle(new_parent->min_rotation_angle());
    std::reverse(max_angle.begin(), max_angle.end());
    max_angle *= -1.0;

    current_node->set_order(order());
    current_node->set_rotation_angle(rotation_angle());
    current_node->SetRotationRange(min_angle(), max_angle());
    current_node->AngleChanged();
    dj::Mat3d rotation_adjustment;
    //    dj::MulMatrix3x3<double>(current_node->rotation_angle_matrix_, new_parent->rotation_matrix(),
    //                             rotation_adjustment());
    dj::MulMatrix3x3<double>(new_parent->rotation_matrix(), current_node->rotation_angle_matrix_,
                            rotation_adjustment());
    rotation_adjustment.Transpose();
    dj::MulMatrix3x3<double>(rotation_adjustment(), current_node->rotation_matrix(),
                            current_node->fixed_rotation_matrix_);
    dj::Vec3d my_pos(current_node->world_pos());
    dj::Vec3d parent_pos(new_parent->world_pos());
    dj::Vec3d offset = my_pos - parent_pos;
    dj::Mat3d parent_rotation(new_parent->rotation_matrix());
    dj::Mat3d offset_transform;//= parent_rotation * dj::Mat3d(current_node->rotation_angle_matrix_);
    dj::MulMatrix3x3<double>(&parent_rotation[0][0], current_node->rotation_angle_matrix_, &offset_transform[0][0]);
    offset_transform.Transpose();
    dj::Vec3d new_offset;//= offset_transform * offset;
    dj::MulMatrix3x3Vec((double (*)[3]) &offset_transform[0][0], &offset[0], &new_offset[0]);
    ASSERT(dj::Abs(new_offset.Magnitude() - offset.Magnitude()) < 1e-2);
    current_node->set_offset(new_offset());
    current_node->set_original_offset(new_offset());
    current_node->set_length(new_parent->length());
  }
  ReorderBones();
}

int Skeleton::SelectJoint( const double* start, const double* end) {
  const double kThreshold = 0.1f;
  int min_idx = -1;
  double closet_distance = 1e10;
  for (int i = 0; i < int(nodes_.size()); ++i) {
    double distance = dj::PointLineDistance<double>(start, end, nodes_[i].world_pos());
    if (distance < kThreshold && distance < closet_distance) {
      closet_distance = distance;
      min_idx = i;
    }
  }
  selected_joint_ = min_idx;
  return selected_joint_;
}


int Skeleton::SelectEffector(const Vec3d &start, const Vec3d &end) {
  const double kThreshold = 0.3;
  double min_distance = 1e10;
  int closet_effector = -1;
  for (int i = 0; i < kLeafNum; ++i) {
    double distance = dj::PointLineDistance<double>(start(), end(), target_end_effector_pos + i * 3);
    if (distance < min_distance && distance < kThreshold) {
      min_distance = distance;
      closet_effector = i;
    }
  }
  selected_end_effector = closet_effector;
  return closet_effector;
}

void Skeleton::HandleKeyPress(unsigned char key) {
  key = tolower(key);
  const double step = 0.01;
#if 0
  switch (key) {
    case '1':
      if (selected_joint >= 0 && selected_joint < skeleton_point_number) {
        SetNewRoot(&nodes[selected_joint]);
      }
      break;
    case 'm':
      selected_end_effector--;
      if (selected_end_effector < 0) {
        selected_end_effector = kLeafNum - 1;
      }
      break;
    case 'n':
      selected_end_effector++;
      if (selected_end_effector >= kLeafNum) {
        selected_end_effector = 0;
      }
      break;
    case 't':
      MoveEffector(-step, 0, 0);
      InverseKinematics();
      break;
    case 'y':
      MoveEffector(step, 0, 0);
      InverseKinematics();
      break;
    case 'u':
      MoveEffector(0, -step, 0);
      InverseKinematics();
      break;
    case 'i':
      MoveEffector(0, +step, 0);
      InverseKinematics();
      break;
    case 'o':
      MoveEffector(0, 0, -step);
      InverseKinematics();
      break;
    case 'p':
      MoveEffector(0, 0, +step);
      InverseKinematics();
      break;
    default:
      break;
  }
  return;
#endif
  int min = 0, max = 18;
  //    static int idx = min;
  const double* angle;
  //  P(Vec3d(nodes[selected_node].rotation_angle()));
  //  P(Vec3d(nodes[selected_node].min_rotation_angle()));
  //  P(Vec3d(nodes[selected_node].max_rotation_angle()));
  //  P(dj::Vec3i(nodes[selected_node].order()));
  switch (key) {
    case '1':
      if (selected_joint_ >= 0 && selected_joint_ < skeleton_point_number_) {
        SetNewRoot(selected_joint_);
      }
      break;
    case 'm':
      selected_node_++;
      if (selected_node_ > max) selected_node_ = min;
      break;
    case 'g':
      min_angle[selected_node_][0] = nodes_[selected_node_].rotation_angle()[0];
      break;
    case 'h':
      max_angle[selected_node_][0] = nodes_[selected_node_].rotation_angle()[0];
      break;
    case 'j':
      min_angle[selected_node_][1] = nodes_[selected_node_].rotation_angle()[1];
      break;
    case 'k':
      max_angle[selected_node_][1] = nodes_[selected_node_].rotation_angle()[1];
      break;
    case 'l':
      min_angle[selected_node_][2] = nodes_[selected_node_].rotation_angle()[2];
      break;
    case ';':
      max_angle[selected_node_][2] = nodes_[selected_node_].rotation_angle()[2];
      break;
    case '\'':
      SaveAngle();
      break;
    case 'n':
      selected_node_--;
      if (selected_node_ < min) selected_node_ = max;
      break;
    case 't':
      angle = nodes_[selected_node_].rotation_angle();
      nodes_[selected_node_].set_rotation_angle(angle[0] - step, angle[1], angle[2]);
      break;
    case 'y':
      angle = nodes_[selected_node_].rotation_angle();
      nodes_[selected_node_].set_rotation_angle(angle[0] + step, angle[1], angle[2]);
      break;
    case 'u':
      angle = nodes_[selected_node_].rotation_angle();
      nodes_[selected_node_].set_rotation_angle(angle[0], angle[1] - step, angle[2]);
      break;
    case 'i':
      angle = nodes_[selected_node_].rotation_angle();
      nodes_[selected_node_].set_rotation_angle(angle[0], angle[1] + step, angle[2]);
      break;
    case 'o':
      angle = nodes_[selected_node_].rotation_angle();
      nodes_[selected_node_].set_rotation_angle(angle[0], angle[1], angle[2] - step);
      break;
    case 'p':
      angle = nodes_[selected_node_].rotation_angle();
      nodes_[selected_node_].set_rotation_angle(angle[0], angle[1], angle[2] + step);
      break;
    case 'b':
      for (unsigned int i = 0; i < nodes_.size(); ++i) {
        nodes_[i].set_rotation_angle(origin_angle + i * 3);
      }
    default:
      break;
  }
  //  for (int i = 0; i < skeleton_point_number; i++) {
  //    P(max_angle[i]);
  //    P(min_angle[i]);
  //  }
  //  Vec3d initial_angle(origin_angle + selected_node * 3);
  //  Vec3d current_angle(nodes[selected_node].rotation_angle());
  //  P(selected_node + 1, current_angle);
  //  P(selected_node + 1, initial_angle);
}

void Skeleton::RotateY(double angle) {
  for (int i = 0; i < skeleton_point_number_; i++) {
    double* p = &skeleton_points_[i * 3];
    double ty, tz;
    ty = p[0] * cosf(angle) - p[2] * sinf(angle);
    tz = p[0] * sinf(angle) + p[2] * cosf(angle);
    p[0] = ty;
    p[2] = tz;
  }
}

void Skeleton::RotateX(double angle) {
  for (int i = 0; i < skeleton_point_number_; i++) {
    double* p = &skeleton_points_[i * 3];
    double ty, tz;
    ty = p[1] * cosf(angle) - p[2] * sinf(angle);
    tz = p[1] * sinf(angle) + p[2] * cosf(angle);
    p[1] = ty;
    p[2] = tz;
  }
}

void Skeleton::Render() {
  //  return;
  // Draw skeleton lines
  glDisable(GL_LIGHTING);
  glLineWidth(1.0);
  glPointSize(4);
  //  glColor3fv(kRed());
  const float* color_map[] = {
    kRed(),
    kGreen(),
    kBlue(),
    kYellow(),
    kOrage(),
    kChocolate(),
    kViolet(),
    kIndigo(),
  };
  const int color_num = sizeof(color_map) / sizeof(double*);

  //  glColor4f(1, 0, 0, 0.8);

  glPushAttrib(GL_LINE_BIT);
  glLineWidth(4.0);
  glBegin(GL_LINES);
  for (int b = 0; b < (int) bones_.size(); ++b) {
    int n1 = bones_[b].first;
    int n2 = bones_[b].second;
    glColor3fv(color_map[b % color_num]);
    Vertex3v(nodes_[n1].world_pos());
    Vertex3v(nodes_[n2].world_pos());
  }
  glEnd();
  glPopAttrib();

  QFont sansFont("Helvetica [Cronyx]", 18);
  sansFont.setFamily("sans serif");
  // draw bone text
  if (1)
    for (int b = 0; b < (int) bones_.size(); ++b) {
      int n1 = bones_[b].first;
      int n2 = bones_[b].second;
      glColor3fv(color_map[b % color_num]);
      dj::Vec3d center(nodes_[n1].world_pos()[0] + nodes_[n2].world_pos()[0],
                       nodes_[n1].world_pos()[1] + nodes_[n2].world_pos()[1],
                       nodes_[n1].world_pos()[2] + nodes_[n2].world_pos()[2]);
      center *= 0.5f;
      global::gl->renderText(center[0], center[1], center[2], QString("%1").arg(b), sansFont);
    }

  // Draw joints
  glColor3fv(kBlue());
  glBegin(GL_POINTS);
  for (int i = 0; i < int(nodes_.size()); ++i) {
    Vertex3v(nodes_[i].world_pos());
  }
  glEnd();

  glColor3fv(kBlack());
  QFont jointFont("Helvetica [Cronyx]", 18);
  jointFont.setFamily("sans serif");
  // draw joint text
  if (0)
    for (int i = 0; i < int(nodes_.size()); ++i) {
      const double* pos = nodes_[i].world_pos();
      global::gl->renderText(pos[0], pos[1], pos[2], QString("%1").arg(i), jointFont);
    }

  // Draw root joint
  {
    //  glColor3fv(kYellow());
    //  glPushMatrix();
    //  const double* root_pos = root_->world_pos();
    //  glTranslatef(root_pos[0], root_pos[1], root_pos[2]);
    //  DrawSphere(0.01, 10, 10);
    //  glPopMatrix();
  }
}
