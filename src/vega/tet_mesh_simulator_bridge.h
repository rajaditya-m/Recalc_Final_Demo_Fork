#ifndef TET_MESH_SIMULATOR_BRIDGE_H
#define TET_MESH_SIMULATOR_BRIDGE_H
#include <vector>
class SparseMatrix;
class Tet;
class BufferManager;
template <class Float> class ConjugateGradientSolver;
class TetMesh;
class StVKForceModel;
class IsotropicHyperelasticFEM;
class CorotationalLinearFEM;
class CorotationalLinearFEMForceModel;
class StVKIsotropicMaterial;
//#define COROTATIONAL_FEM
//#define STVK_FEM

#ifdef COROTATIONAL_FEM
typedef CorotationalLinearFEM Fem;
#elif defined(STVK_FEM)
typedef StVKForceModel Fem;
#else
typedef IsotropicHyperelasticFEM Fem;
#endif
class TetMeshSimulatorBridge
{
public:
  typedef double real;
  TetMeshSimulatorBridge(Tet* tet_mesh, std::vector<std::pair<int,int> > &interfaceV, int split);

  ~TetMeshSimulatorBridge();
  real Simulate(real dt);
  void Reset();

  void ComputeInternalForceAndTangentStiffnessMatrix(double dt);
  void ComputePartialInternalForceAndTangentStiffnessMatrix(std::vector<int>& tets);
  void LoadPosition(double *pos);
  void UpdateOffset();
  void SetInterface(std::vector<std::pair<int,int> > &infVerts);
  void makeLinks(int i, int j);
  void makeSpringConnections();
  void setSimMode(int i)   { sim_mode_ = i;     }

  SparseMatrix* GetStitchStiffnessMatrixPointer(int ID);

  int v_num_;
  int tet_num_;
  int (*tets_)[4];
  real young_modulus_;
  real poisson_ratio_;
  real density_;
  real gravity_[3];
  real (*internal_force_)[3];
  real (*rhs_)[3];
  real (*u_)[3];
  real (*vel_)[3];
  real (*rest_pos_)[3];
  real (*tmp_pos_)[3];
  SparseMatrix* tangent_stiffness_matrix_;
  SparseMatrix* stitch_tangent_stiffness_matrix_1_;
  SparseMatrix* stitch_tangent_stiffness_matrix_2_;
  SparseMatrix* mass_matrix_;
  SparseMatrix* spring_mat_;
  ConjugateGradientSolver<real>* cg_solver_;
  Fem* inv_fem_force_model_;
  StVKIsotropicMaterial* material_;
  Tet* tet_mesh_;
  TetMesh* vega_mesh_;
  real* element_force_;
  real* element_k_;
  BufferManager* buf_;

  std::vector<std::pair<int,int> > interface_points_1_;
  std::vector<std::pair<int,int> > interface_points_2_;
  int split_position_;
  int sim_mode_;

};

#endif // TET_MESH_SIMULATOR_BRIDGE_H
