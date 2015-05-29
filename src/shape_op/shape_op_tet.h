#ifndef SHAPEOPTET_H
#define SHAPEOPTET_H
#include <Eigen/Dense>
#include "Solver.h"
#include "Types.h"
#include "tet.h"
template <class T> class AffineTransformer;
namespace ShapeOp {
class Solver;
}
class ShapeOpTet : public Tet
{
public:
  ShapeOpTet(const char* mesh_file, AffineTransformer<double>* affine_transform = NULL, bool initialize_fem_module = false);
  void Simulate(double dt);
  ~ShapeOpTet();
  ShapeOp::Solver *shape_op_solver_;
};

#endif // SHAPEOPTET_H
