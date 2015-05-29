#ifndef MASS_SPRING_OBJ_H
#define MASS_SPRING_OBJ_H
#include "simulation_tetrahedral_mesh.h"

class TetCollider;
class TetrahedralMeshIO;
template <class Float> class ConjugateGradientSolver;
template <class T> class AffineTransformer;

class MassSpringVolumetricObject : public SimulationTetrahedralMesh
{
public:
  typedef SimulationTetrahedralMesh Super;
  MassSpringVolumetricObject(TetrahedralMeshIO* mesh_io, const char* file_name, AffineTransformer<Real>* transformer = NULL);
  MassSpringVolumetricObject(TetrahedralMesh* mesh);

  void Construct();
  virtual ~MassSpringVolumetricObject();
  void ImplicitStepSolvePosition(Real dt);
  void ImplicitStepSolveVelocity(Real dt);
  void ExplicitStep(Real dt);
  void Simulate(Real dt, bool handle_collision = false);

  void ComputeRestEdgeLength();

  std::vector<Real> edge_stiffness_;
  std::vector<Vec3> vel_;
  std::vector<Vec3> prev_vert_;
  std::vector<Vec3> tmp_vert_buf_;
  std::vector<Real> rest_edge_length_;
  ConjugateGradientSolver<Real>* cg_solver_;
  TetCollider* collider_;
  Real stiffness_;
  Real damping_;
};

#endif // MASS_SPRING_OBJ_H
