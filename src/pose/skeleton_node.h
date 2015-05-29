#ifndef SKELETON_NODE_H_
#define SKELETON_NODE_H_
#include <vector>
#include <cmath>
#include <queue>
#include <functional>
#include <utility>
#include "MY_MATH.h"
#include "vector_lib.h"
#include "print_macro.h"
#ifndef PI
#define PI 3.1415926
#endif



template <class FloatType>
struct RotatorBase
{
  virtual void operator()(FloatType angle, FloatType (*matrix)[3]) = 0;
  virtual void Derivative(FloatType angle, FloatType (*matrix)[3]) = 0;
  virtual ~RotatorBase() {}
};

template <class FloatType>
struct RotateX : RotatorBase<FloatType>
{
  void operator()(FloatType angle, FloatType (*matrix)[3])
  {
    FloatType sin = std::sin(angle);
    FloatType cos = std::cos(angle);
    matrix[0][0] = 1; matrix[0][1] = 0;   matrix[0][2] = 0;
    matrix[1][0] = 0; matrix[1][1] = cos; matrix[1][2] = -sin;
    matrix[2][0] = 0; matrix[2][1] = sin; matrix[2][2] = cos;
  }
  void Derivative(FloatType angle, FloatType (*matrix)[3])
  {
    FloatType sin = std::sin(angle);
    FloatType cos = std::cos(angle);
    matrix[0][0] = 0; matrix[0][1] = 0;   matrix[0][2] = 0;
    matrix[1][0] = 0; matrix[1][1] = -sin; matrix[1][2] = -cos;
    matrix[2][0] = 0; matrix[2][1] = cos; matrix[2][2] = -sin;
  }
  ~RotateX() {}
};

template <class FloatType>
struct RotateY : RotatorBase<FloatType>
{
  void operator()(FloatType angle, FloatType (*matrix)[3])
  {
    FloatType sin = std::sin(angle);
    FloatType cos = std::cos(angle);
    matrix[0][0] = cos;  matrix[0][1] = 0; matrix[0][2] = sin;
    matrix[1][0] = 0;    matrix[1][1] = 1; matrix[1][2] = 0;
    matrix[2][0] = -sin; matrix[2][1] = 0; matrix[2][2] = cos;
  }

  void Derivative(FloatType angle, FloatType (*matrix)[3])
  {
    FloatType sin = std::sin(angle);
    FloatType cos = std::cos(angle);
    matrix[0][0] = -sin; matrix[0][1] = 0; matrix[0][2] = cos;
    matrix[1][0] = 0;    matrix[1][1] = 0; matrix[1][2] = 0;
    matrix[2][0] = -cos; matrix[2][1] = 0; matrix[2][2] = -sin;
  }
  ~RotateY() {}
};

template <class FloatType>
struct RotateZ : public RotatorBase<FloatType>
{
  void operator()(FloatType angle, FloatType (*matrix)[3])
  {
    FloatType sin = std::sin(angle);
    FloatType cos = std::cos(angle);
    matrix[0][0] = cos; matrix[0][1] = -sin; matrix[0][2] = 0;
    matrix[1][0] = sin; matrix[1][1] = cos;  matrix[1][2] = 0;
    matrix[2][0] = 0;   matrix[2][1] = 0;    matrix[2][2] = 1;
  }

  void Derivative(FloatType angle, FloatType (*matrix)[3])
  {
    FloatType sin = std::sin(angle);
    FloatType cos = std::cos(angle);
    matrix[0][0] = -sin; matrix[0][1] = -cos; matrix[0][2] = 0;
    matrix[1][0] = cos;  matrix[1][1] = -sin; matrix[1][2] = 0;
    matrix[2][0] = 0;    matrix[2][1] = 0;    matrix[2][2] = 0;
  }
  ~RotateZ() {};
};

const static double kNullValue = double (1e10);
class SkeletonNode
{
  typedef double Real;
  typedef RotatorBase<Real> RotationFunctor;
public:

  SkeletonNode(void)
  {
    Init();
  }

  SkeletonNode(SkeletonNode *parent, Real length)
  {
    Init();
    parent_ = parent;
    length_ = length;
  }

  ~SkeletonNode()
  {
    delete Rotate_[0];
    delete Rotate_[1];
    delete Rotate_[2];
  }

  void set_id(int id)
  {
    id_ = id;
  }

  int id()
  {
    return id_;
  }

  void set_parent(SkeletonNode *parent)
  {
    parent_ = parent;
  }

  SkeletonNode *parent() const
  {
    return parent_;
  }

  bool IsRoot() const
  {
    return parent() == NULL;
  }

  bool IsLeaf() const
  {
    return children_.size() == 0;
  }

  void AddChild(SkeletonNode *child)
  {
    if (child != NULL) children_.push_back(child);
  }

  Real length()
  {
    return length_;
  }
  void set_length(Real length)
  {
    length_ = length;
  }

  std::vector<SkeletonNode *> &children(void)
  {
    return children_;
  }

  const Real *world_pos(void) const
  {
    return world_pos_;
  }

  Real *reference_frame()
  {
    return (Real *) reference_frame_;
  }

  void set_reference_frame(Real *reference_frame)
  {
    memcpy(reference_frame_, reference_frame, sizeof(Real) * 9);
  }

  void set_rotation_matrix(Real *rotatoin_matrix)
  {
    memcpy(rotation_matrix_, rotatoin_matrix, sizeof(Real) * 9);
  }

  Real *rotation_matrix()
  {
    return (Real *) rotation_matrix_;
  }

  void set_world_pos(Real x, Real y, Real z)
  {
    world_pos_[0] = x;
    world_pos_[1] = y;
    world_pos_[2] = z;
  }


  void set_world_pos(Real *pos)
  {
    world_pos_[0] = pos[0];
    world_pos_[1] = pos[1];
    world_pos_[2] = pos[2];
  }

  const int *order()
  {
    return order_;
  }

  void set_order(int first, int second, int third)
  {
    order_[0] = first;
    order_[1] = second;
    order_[2] = third;
  }

  void set_order(int *order)
  {
    order_[0] = order[0];
    order_[1] = order[1];
    order_[2] = order[2];
  }

  void set_offset(Real x, Real y, Real z)
  {
    offset_[0] = x;
    offset_[1] = y;
    offset_[2] = z;
  }

  void set_offset(Real *offset)
  {
    offset_[0] = offset[0];
    offset_[1] = offset[1];
    offset_[2] = offset[2];
  }

  void set_original_offset(Real *offset)
  {
    original_offset_[0] = offset[0];
    original_offset_[1] = offset[1];
    original_offset_[2] = offset[2];
  }

  void load_original_offset()
  {
    offset_[0] = original_offset_[0];
    offset_[1] = original_offset_[1];
    offset_[2] = original_offset_[2];
  }

  Real *offset()
  {
    return offset_;
  }

  std::pair<Real *, Real *> GetRotationRange(void)
  {
    return std::make_pair<Real *, Real *>((Real *) min_rotation_angle_, (Real *) max_rotation_angle_);
  }

  void SetRotationRange(Real *min_angle, Real *max_angle)
  {
    for (int i = 0; i < 3; ++i)
    {
      min_rotation_angle_[i] = min_angle[i];
      max_rotation_angle_[i] = max_angle[i];
    }
  }

  const Real *min_rotation_angle()
  {
    return min_rotation_angle_;
  }

  const Real *max_rotation_angle()
  {
    return max_rotation_angle_;
  }

  void set_rotation_angle(Real *angle)
  {
    if (angle[0] != kNullValue) rotation_angle_[0] = angle[0];
    if (angle[1] != kNullValue) rotation_angle_[1] = angle[1];
    if (angle[2] != kNullValue) rotation_angle_[2] = angle[2];
  }

  void set_rotation_angle(Real alpha, Real beta, Real gamma)
  {
    if (alpha != kNullValue) rotation_angle_[0] = alpha;
    if (beta != kNullValue) rotation_angle_[1] = beta;
    if (gamma != kNullValue) rotation_angle_[2] = gamma;
  }

  Real *rotation_angle()
  {
    return rotation_angle_;
  }

  void ComputePosition(void)
  {
    // Node is not the root
    if (parent_ != NULL)
    {
      UpdateRotationMatrix();
      Real world_offset[3];
      Real local_offset[3];
      dj::MulMatrix3x3Vec<Real>(rotation_angle_matrix_, offset_, local_offset);
      dj::MulMatrix3x3Vec<Real>((Real (*)[3]) parent_->rotation_matrix(), local_offset, world_offset);
      //      dj::Vec3f off(offset_);
      //      if (off.Norm2() > 10) {
      //        P(id_, off);
      //      }
      const Real *parent_pos = parent_->world_pos();
      //      dj::Vec3f old(world_pos_);
      world_pos_[0] = parent_pos[0] + world_offset[0];
      world_pos_[1] = parent_pos[1] + world_offset[1];
      world_pos_[2] = parent_pos[2] + world_offset[2];
      //          dj::Vec3f n(world_pos_);
      //          auto diff = n - old;
      //          if (!(diff.Norm2() < 1e-6)) {
      //            P(diff, name_);
      //            P(old);
      //            P(n);
      //            P(dj::Vec3f(world_offset));
      //            P(dj::Vec3i(order_));
      //            P(this);
      //            P(dj::Mat3f(parent()->rotation_matrix()));
      //            P(dj::Vec3f(rotation_angle()));
      //            P(dj::Mat3f(fixed_rotation_matrix_));
      //            P(dj::Mat3f(rotation_matrix()));
      //            P(dj::Mat3f(rotation_angle_matrix_));
      //          }
      //            KK;
      //          MY_ASSERT(diff.Norm2() < 1e-6);
    }
    for (unsigned int i = 0; i < children_.size(); ++i)
    {
      children_[i]->ComputePosition();
    }
  }

  void ClearDerivative()
  {
    memset(rotation_matrix_derivative_, 0, sizeof(Real) * 3 * 3 * 3);
    memset(world_pos_derivative_, 0, sizeof(Real) * 3 * 3);
  }

  int  GetNumOfChildren()
  {
    return int(children_.size());
  }

  void PushChildrenIntoQueue(std::queue<SkeletonNode *> &queue) const
  {
    for (unsigned i = 0; i < children_.size(); ++i)
    {
      queue.push(children_[i]);
    }
  }

  int GetDepth()
  {
    int depth = 1;
    SkeletonNode *parent = parent_;
    for (; parent; parent = parent->parent(), ++depth)
    {
    }
    return depth;
  }

  void ComputeDerivative(bool is_changing_node)
  {
    ASSERT(parent_ != NULL);
    if (is_changing_node)
    {
      Real tmp_rotation_matrix_derivative[3][3][3];
      Real angle_rotation[3][3][3];
      for (int i = 0; i < 3; ++i)
      {
        (*Rotate_[order_[i]])(rotation_angle_[i], angle_rotation[i]);
        (*Rotate_[order_[i]]).Derivative(rotation_angle_[i], tmp_rotation_matrix_derivative[i]);
      }
      Real tmp_matrix[2][3][3];
      // First angle
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix(), tmp_rotation_matrix_derivative[0], tmp_matrix[0]);
      dj::MulMatrix3x3<Real>(tmp_matrix[0], angle_rotation[1], tmp_matrix[1]);
      dj::MulMatrix3x3<Real>(tmp_matrix[1], angle_rotation[2], tmp_matrix[0]);
      dj::MulMatrix3x3<Real>(tmp_matrix[0], fixed_rotation_matrix_, rotation_matrix_derivative_[0]);
      dj::MulMatrix3x3Vec<Real>(tmp_matrix[0], offset_, world_pos_derivative_[0]);
      // Second angle
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix(), angle_rotation[0], tmp_matrix[0]);
      dj::MulMatrix3x3<Real>(tmp_matrix[0], tmp_rotation_matrix_derivative[1], tmp_matrix[1]);
      dj::MulMatrix3x3<Real>(tmp_matrix[1], angle_rotation[2], tmp_matrix[0]);
      dj::MulMatrix3x3<Real>(tmp_matrix[0], fixed_rotation_matrix_, rotation_matrix_derivative_[1]);
      dj::MulMatrix3x3Vec<Real>(tmp_matrix[0], offset_, world_pos_derivative_[1]);
      // Third angle
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix(), angle_rotation[0], tmp_matrix[0]);
      dj::MulMatrix3x3<Real>(tmp_matrix[0], angle_rotation[1], tmp_matrix[1]);
      dj::MulMatrix3x3<Real>(tmp_matrix[1], tmp_rotation_matrix_derivative[2], tmp_matrix[0]);
      dj::MulMatrix3x3<Real>(tmp_matrix[0], fixed_rotation_matrix_, rotation_matrix_derivative_[2]);
      dj::MulMatrix3x3Vec<Real>(tmp_matrix[0], offset_, world_pos_derivative_[2]);
    }
    else
    {
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix_derivative_[0],
                             local_rotation_matrix_, rotation_matrix_derivative_[0]);
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix_derivative_[1],
                             local_rotation_matrix_, rotation_matrix_derivative_[1]);
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix_derivative_[2],
                             local_rotation_matrix_, rotation_matrix_derivative_[2]);
      Real tmp_matrix[3][3];
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix_derivative_[0], rotation_angle_matrix_, tmp_matrix);
      dj::MulMatrix3x3Vec<Real>(tmp_matrix, offset_, world_pos_derivative_[0]);
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix_derivative_[1], rotation_angle_matrix_, tmp_matrix);
      dj::MulMatrix3x3Vec<Real>(tmp_matrix, offset_, world_pos_derivative_[1]);
      dj::MulMatrix3x3<Real>(parent_->rotation_matrix_derivative_[2], rotation_angle_matrix_, tmp_matrix);
      dj::MulMatrix3x3Vec<Real>(tmp_matrix, offset_, world_pos_derivative_[2]);
      // Position derivative wrt to alpha
      world_pos_derivative_[0][0] += parent_->world_pos_derivative_[0][0];
      world_pos_derivative_[0][1] += parent_->world_pos_derivative_[0][1];
      world_pos_derivative_[0][2] += parent_->world_pos_derivative_[0][2];
      // Position derivative wrt to beta
      world_pos_derivative_[1][0] += parent_->world_pos_derivative_[1][0];
      world_pos_derivative_[1][1] += parent_->world_pos_derivative_[1][1];
      world_pos_derivative_[1][2] += parent_->world_pos_derivative_[1][2];
      // Position derivative wrt to gamma
      world_pos_derivative_[2][0] += parent_->world_pos_derivative_[2][0];
      world_pos_derivative_[2][1] += parent_->world_pos_derivative_[2][1];
      world_pos_derivative_[2][2] += parent_->world_pos_derivative_[2][2];
    }
  }

  void ComputeWorldFrame(Real (*world_frame)[3])
  {
    dj::MulMatrix3x3<Real>(rotation_matrix_, reference_frame_, world_frame);
  }


  void set_initial_world_pos(double x, double y, double z)
  {
    initial_world_pos_[0] = x;
    initial_world_pos_[1] = y;
    initial_world_pos_[2] = z;
  }

  void set_initial_world_pos(const double *pos)
  {
    initial_world_pos_[0] = pos[0];
    initial_world_pos_[1] = pos[1];
    initial_world_pos_[2] = pos[2];
  }

  double *initial_world_pos()
  {
    return initial_world_pos_;
  }

  std::string name()
  {
    return name_;
  }

  void set_name(std::string &new_name)
  {
    name_ = new_name;
  }


  void AngleChanged()
  {
    Real rotation_matricies[3][3][3];
    for (int i = 0; i < 3; ++i)
    {
      //      double tmp = rotation_angle_[i];
      //      double tmp = rotation_angle_[i];
      //      rotation_angle_[i] = Clamp(rotation_angle_[i], min_rotation_angle_[i], max_rotation_angle_[i]);
      //      if (tmp != rotation_angle_[i]) {
      //        P(i, name(), dj::Radian2Degree(tmp));
      //        P(dj::Radian2Degree(min_rotation_angle_[i]));
      //          P(dj::Radian2Degree(max_rotation_angle_[i]));
      //      }
      //      if (tmp != rotation_angle_[i]) {
      //        P(name_, i, rotation_angle_[i], tmp);
      //        P(max_rotation_angle_[i], min_rotation_angle_[i]);
      //        MY_ASSERT(false);
      //      }
      (*Rotate_[order_[i]])(rotation_angle_[i], rotation_matricies[i]);
    }
    Real tmp_matrix[3][3];
    dj::MulMatrix3x3<Real>(rotation_matricies[0], rotation_matricies[1],
                           tmp_matrix);
    dj::MulMatrix3x3<Real>(tmp_matrix, rotation_matricies[2],
                           rotation_angle_matrix_);
    //    dj::MulMatrix3x3<Real>(fixed_rotation_matrix_, rotation_angle_matrix_,
    dj::MulMatrix3x3<Real>(rotation_angle_matrix_, fixed_rotation_matrix_,
                           local_rotation_matrix_);
  }

  void UpdateRotationMatrix(void)
  {
    if (parent_ == NULL) return;
    AngleChanged();
    dj::MulMatrix3x3<Real>(parent_->rotation_matrix(), local_rotation_matrix_,
                           rotation_matrix_);
  }

public:
  Real rotation_matrix_derivative_[3][3][3];
  Real world_pos_derivative_[3][3];
  Real fixed_rotation_matrix_[3][3];
  Real rotation_angle_matrix_[3][3];
  Real local_rotation_matrix_[3][3];
  RotationFunctor *Rotate_[3];
  void ExportFakeBVH(std::ostream &out, std::string tab);
  void ImportFakeBVH(std::istream &in, SkeletonNode *parent, std::vector<SkeletonNode> &nodes);


private:
  void Init(void)
  {
    name_ = "Unknown";
    parent_ = NULL;
    children_.clear();
    length_ = 0;
    world_pos_[0] = world_pos_[1] = world_pos_[2] = 0.0;
    rotation_angle_[0] = rotation_angle_[1] = rotation_angle_[2] = 0;
    min_rotation_angle_[0] = min_rotation_angle_[1] = min_rotation_angle_[2] = -PI;
    max_rotation_angle_[0] = max_rotation_angle_[1] = max_rotation_angle_[2] = PI;
    memset(reference_frame_, 0, sizeof(Real) * 9);
    reference_frame_[0][0] = 1;
    reference_frame_[1][1] = 1;
    reference_frame_[2][2] = 1;
    memcpy(rotation_matrix_, reference_frame_, sizeof(Real) * 9);
    memcpy(fixed_rotation_matrix_, reference_frame_, sizeof(Real) * 9);
    Rotate_[0] = new RotateX<Real>();
    Rotate_[1] = new RotateY<Real>();
    Rotate_[2] = new RotateZ<Real>();
    order_[0] = 2;
    order_[1] = 1;
    order_[2] = 0;
  }


  SkeletonNode *parent_;
  std::vector<SkeletonNode *> children_;
  int id_;
  Real length_;
  Real offset_[3];
  Real original_offset_[3];
  Real world_pos_[3];
  Real initial_world_pos_[3];
  Real reference_frame_[3][3];
  Real rotation_matrix_[3][3];
  int order_[3]; // Decides which rotation is applied first
  Real rotation_angle_[3];
  Real min_rotation_angle_[3];
  Real max_rotation_angle_[3];
  std::string name_;
};





#endif // SKELETON_NODE_H_
