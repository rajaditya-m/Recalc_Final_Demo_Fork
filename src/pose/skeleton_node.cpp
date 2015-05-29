#include "skeleton_node.h"

#include "global.h"



void SkeletonNode::ExportFakeBVH(std::ostream &out, std::string tab)
{
  const std::string channel_map[] = {"Xrotation", "Yrotation", "Zrotation"};
  if (IsRoot())
  {
    out << tab << "ROOT " << name_ << std::endl;
  }
  else if (IsLeaf())
  {
    out << tab << "End " << name_ <<  std::endl;
  }
  else
  {
    out << tab << "JOINT " << name_ << std::endl;
  }

  out << tab << "{" << std::endl;
  std::string new_tab = tab + "\t";
  out << new_tab << "OFFSET " <<  offset_[0] << " " <<  offset_[1] << " " <<  offset_[2] << std::endl;
  // Output channels
  out << new_tab << "CHANNELS 3 ";
  for (int i = 0; i < 3; ++i)
  {
    out << channel_map[order_[i]] << " ";
  }
  out << std::endl;

  dj::Vec3d angle(rotation_angle_);
  dj::Vec3d max_angle(max_rotation_angle_);
  dj::Vec3d min_angle(min_rotation_angle_);
  angle *= (180.0f / global::kPi);
  max_angle *= (180.0f / global::kPi);
  min_angle *= (180.0f / global::kPi);
  out << new_tab << "ANGLE " << angle[0] << " " << angle[1] << " " << angle[2] << std::endl;
  out << new_tab << "POSITION " << initial_world_pos()[0] << " " << initial_world_pos()[1] << " " << initial_world_pos()[2] << std::endl;
  out << new_tab << "MAX_ANGLE " << max_angle[0] << " " << max_angle[1] << " " << max_angle[2] << std::endl;
  out << new_tab << "MIN_ANGLE " << min_angle[0] << " " << min_angle[1] << " " << min_angle[2] << std::endl;
  // Output children
  for (int i = 0; i < int(children_.size()); ++i)
  {
    children_[i]->ExportFakeBVH(out, new_tab);
  }
  out << tab << "}" << std::endl;
}


void ReadAngleFromDegree2Radian(std::istream &in, double *out)
{
  in >> out[0] >> out[1] >> out[2];
  out[0] *= (global::kPi / 180.0f);
  out[1] *= (global::kPi / 180.0f);
  out[2] *= (global::kPi / 180.0f);
}

void SkeletonNode::ImportFakeBVH(std::istream &in, SkeletonNode *parent, std::vector<SkeletonNode> &nodes)
{
  //Read the joint name
  parent_ = parent;
  in >> name_;
  while (true)
  {
    std::string token;
    in >> token;
    if (token == "{")
    {
    }
    else if (token == "}")
    {
      break;
    }
    else if (token == "OFFSET")
    {
      in >> offset_[0] >> offset_[1] >> offset_[2];
    }
    else if (token == "CHANNELS")
    {
      //Read the data number and order.
      int channel_num;
      in >> channel_num;
      for (int i = 0; i < 3; ++i)
      {
        in >> token;
        if (token == "Xrotation") order_[i] = 0;
        else if (token == "Yrotation") order_[i] = 1;
        else if (token == "Zrotation") order_[i] = 2;
        else
        {
          std::cerr << "SkeletonNode::ImportFakeBVH() => Incorrect fake BVH format" << std::endl;
          exit(0);
        }
      }
    }
    else if (token == "JOINT" || token == "End")
    {
      //Create a new child.
      nodes.emplace_back();
      children_.push_back(&nodes.back());
      //Read the file for the child.
      nodes.back().ImportFakeBVH(in, this, nodes);
    }
    else if (token == "MAX_ANGLE")
    {
      ReadAngleFromDegree2Radian(in, max_rotation_angle_);
    }
    else if (token == "MIN_ANGLE")
    {
      ReadAngleFromDegree2Radian(in, min_rotation_angle_);
    }
    else if (token == "POSITION")
    {
      in >> world_pos_[0] >> world_pos_[1] >> world_pos_[2];
      set_initial_world_pos(world_pos_);
    }
    else if (token == "ANGLE")
    {
      ReadAngleFromDegree2Radian(in, rotation_angle_);
    }
  }
  if (parent)
  {
    dj::Normalize3(offset_);
    const double *my_pos = world_pos_;
    const double *parent_pos = parent_->world_pos();
    double length = dj::Distance3(my_pos, parent_pos);
    offset_[0] *= length;
    offset_[1] *= length;
    offset_[2] *= length;
  }
  //  P(name_, dj::Vec3i(order_), dj::Vec3f(rotation_angle()));
  //  P(this);
  //  P(dj::Mat3f(rotation_matrix()));

  //  using namespace dj;
  //  P(name_, Vec3f(world_pos_), Vec3f(offset()));
}
