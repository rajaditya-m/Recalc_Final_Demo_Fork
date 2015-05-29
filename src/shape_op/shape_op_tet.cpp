#include "shape_op_tet.h"
#include "Solver.h"
#include "Constraint.h"
#include "Types.h"
#include "global.h"
//#include "print_macro.h"


//  Tet(const char *filename, double _limit_threshold, AffineTransformer<double>* affine_transform = NULL,  bool initialize_fem_module = true);

ShapeOpTet::ShapeOpTet(const char *mesh_file, AffineTransformer<double> *affine_transform, bool initialize_fem_module)
  : Tet(mesh_file, 0, affine_transform, initialize_fem_module) {
  shape_op_solver_ = new ShapeOp::Solver;
  ShapeOp::Matrix3X pos = ShapeOp::Matrix3X::Zero(3, vertex_num_);
  double total_mass = 0;
  for (int v = 0; v < vertex_num_; ++v) {
    total_mass += mass_[v];
    pos(0, v) = X[v * 3 + 0];
    pos(1, v) = X[v * 3 + 1];
    pos(2, v) = X[v * 3 + 2];
  }
  shape_op_solver_->setPoints(pos);
  for (int t = 0; t < tet_number; ++t) {
    std::vector<int> verts(tet_ + t * 4, tet_ + t * 4 + 4);
    const double kWeight = 10.0;
    std::shared_ptr<ShapeOp::TetrahedronStrainConstraint>
    tet_constraint(new ShapeOp::TetrahedronStrainConstraint(verts, kWeight, pos, 1.0, 1.0));
    shape_op_solver_->addConstraint(tet_constraint);
  }
  shape_op_solver_->initialize(true, total_mass / vertex_num_, 0.95, global::time_step);
}

void ShapeOpTet::Simulate(double dt) {
//  std::vector<std::shared_ptr<ShapeOp::Force> > forces;
  shape_op_solver_->setTimeStep(dt);
  ShapeOp::Vector3 force;
  shape_op_solver_->clearForces();
  int force_vert = GetUIForce(&force[0]);
  if (force_vert >= 0 && force_vert < vertex_num_) {
    std::shared_ptr<ShapeOp::Force> ui_force(new ShapeOp::VertexForce(force, force_vert));
    shape_op_solver_->addForces(ui_force);
//    forces.push_back(ui_force);;
  }
  ShapeOp::Vector3 g(global::gravity[0] * mass_[0], global::gravity[1] * mass_[0], global::gravity[2] * mass_[0]);
  std::shared_ptr<ShapeOp::GravityForce> gravity_force(new ShapeOp::GravityForce(g));
  shape_op_solver_->addForces(gravity_force);
  const double kFloor = -0.10;
  const double kFloorStiffness = 250;
  for (int v = 0; v < vertex_num_; ++v) {
    if (shape_op_solver_->getPoints()(1, v) < kFloor) {
      ShapeOp::Vector3 force(0, 0, 0);
      force[1] = (kFloor - shape_op_solver_->getPoints()(1, v)) * mass_[v] * kFloorStiffness;
      {
        std::shared_ptr<ShapeOp::Force> collision_force(new ShapeOp::VertexForce(force, v));
        shape_op_solver_->addForces(collision_force);
//        forces.push_back(collision_force);;
      }
    }
  }
  shape_op_solver_->solve(10);
  for (int v = 0; v < vertex_num_; ++v) {
    X[v * 3 + 0] = shape_op_solver_->getPoints()(0, v);
    X[v * 3 + 1] = shape_op_solver_->getPoints()(1, v);
    X[v * 3 + 2] = shape_op_solver_->getPoints()(2, v);
  }
}

ShapeOpTet::~ShapeOpTet()
{
  delete shape_op_solver_;
}
