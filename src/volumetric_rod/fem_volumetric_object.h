#ifndef FEM_VOLUMETRIC_OBJECT_H
#define FEM_VOLUMETRIC_OBJECT_H
#include "global.h"
#include "simulation_tetrahedral_mesh.h"
class TetrahedralMeshIO;
template <class Float> class ConjugateGradientSolver;
template <class T> class AffineTransformer;

class FemVolumetricObject : public SimulationTetrahedralMesh
{
public:
  typedef  Real (*Real3)[3];
  typedef SimulationTetrahedralMesh Super;
  typedef dj::Vec<TetrahedralMesh::Vec3, 3, false> Mat3;
  FemVolumetricObject(TetrahedralMeshIO* mesh_io, const char* file_name, AffineTransformer<Real>* transformer = NULL);
  FemVolumetricObject(TetrahedralMesh* mesh);
  void Construct();
  void Simulate(Real dt);
  inline void ComputeTetForce(int tet_idx, Mat3 &force);
  inline void ComputeTetDifferentialForce(int tet_idx, Real3 dx, Mat3 &force);
  void Render(int render_mode);
  void ComputeInvDm();

  std::vector<Real> mu_;
  std::vector<Real> lambda_;
  std::vector<Mat3> inv_dm_;

  std::vector<Vec3> force_;

  std::vector<Vec3> vel_;
  std::vector<Vec3> prev_vert_;
  std::vector<Vec3> tmp_vert_buf_;

  ConjugateGradientSolver<Real>* cg_solver_;
//  const Real kStiffness_;
//  const Real kDamping_;
};

#endif // FEM_VOLUMETRIC_OBJECT_H
