#ifndef SIMULATION_TETRAHEDRAL_MESH_H
#define SIMULATION_TETRAHEDRAL_MESH_H
#include  "tetrahedral_mesh.h"
class TetrahedralMeshIO;
template <class T> class AffineTransformer;
class SimulationTetrahedralMesh: public TetrahedralMesh
{
public:
  typedef TetrahedralMesh Super;
  typedef dj::Vec<TetrahedralMesh::Vec3, 3, false> Mat3;
  SimulationTetrahedralMesh(TetrahedralMeshIO* mesh_io, const char* file_name, AffineTransformer<Real>* transformer = NULL);
  explicit SimulationTetrahedralMesh(TetrahedralMesh* mesh);
  void Construct();
  void ComputeTetrahedraVolume();
  void ComputeLumpedMass();

  std::vector<Real> tet_vol_;
  std::vector<Real> tet_density_;
  std::vector<Real> mass_;
  Real total_mass_;
  Real total_volume_;
};


#endif // SIMULATION_TETRAHEDRAL_MESH_H
