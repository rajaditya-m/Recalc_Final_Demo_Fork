#include "simulation_tetrahedral_mesh.h"
//#include "affine_transformer.h"

SimulationTetrahedralMesh::SimulationTetrahedralMesh(TetrahedralMeshIO* mesh_io, const char* file_name, AffineTransformer<Real>* transformer)
  : Super(mesh_io, file_name, transformer)
{
  Construct();
}

SimulationTetrahedralMesh::SimulationTetrahedralMesh(TetrahedralMesh* mesh)
  : Super(*mesh)
{
  Construct();
}

void SimulationTetrahedralMesh::Construct() {
  tet_vol_.resize(tet_num_);
  tet_density_.resize(tet_num_);
}

void SimulationTetrahedralMesh::ComputeTetrahedraVolume()
{
  total_volume_ = 0;
  for (int t = 0; t < tet_num_; ++t) {
    Vec3* v[4] = {
      &vert_[tet_[t][0]],
      &vert_[tet_[t][1]],
      &vert_[tet_[t][2]],
      &vert_[tet_[t][3]]
    };
    Vec3 vec[3] = {
      *(v[1]) - *(v[0]),
      *(v[2]) - *(v[0]),
      *(v[3]) - *(v[0]),
    };
    Vec3 cross;
    dj::Cross3(&vec[0][0], &vec[1][0], &cross[0]);
    tet_vol_[t] = (Real(1) / 6) * (vec[2] * cross);
    total_volume_ += tet_vol_[t];
  }
}


void SimulationTetrahedralMesh::ComputeLumpedMass()
{
  mass_.resize(v_num_);
  std::fill(mass_.begin(), mass_.end(), Real(0));
  total_mass_ = 0;
  for (int t = 0; t < tet_num_; ++t) {
    Real m = tet_density_[t] * tet_vol_[t];
    mass_[tet_[t][0]] += m / Real(4.0);
    mass_[tet_[t][1]] += m / Real(4.0);
    mass_[tet_[t][2]] += m / Real(4.0);
    mass_[tet_[t][3]] += m / Real(4.0);
  }
  for (int v = 0; v < v_num_; ++v) {
    total_mass_ += mass_[v];
  }
}
