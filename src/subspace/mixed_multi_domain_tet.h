#ifndef MIXEDMULTIDOMAINTET_H
#define MIXEDMULTIDOMAINTET_H
#include <unordered_set>
#include <vector>
#include <unordered_map>
#include <limits.h>
#include <float.h>
#include "multi_domain_tet.h"

//class SparseMatrix;
//class PardisoSolver;
class MixedSparseMatrix;
template <class Float> class BiconjugateGradientStablized;
//#include "affine_transformer.h"
template <class T> class AffineTransformer;
class MixedMultiDomainTet : public MultiDomainTet {
public:
  enum TetType {
    kFullTet = 0xf,
  };
  struct IntPairHash {
    std::size_t operator()(const std::pair<int, int>& p) const {
      return std::hash<int>()((p.first << 16) + p.second);
    }
  };
  struct EmbededEdge {
    int v0, v1;
    double weight;
    EmbededEdge(int v0_ = INT_MIN, int v1_ = INT_MIN, double weight_ = DBL_MAX) {
      v0 = v0_;
      v1 = v1_;
      weight = weight_;
    }
  };

  typedef Eigen::Matrix<int, 2, 1> Vec2i;
  typedef Eigen::Matrix<int, 4, 1> Vec4i;

  typedef std::unordered_map<int, int> Map;
  typedef MultiDomainTet Super;
  MixedMultiDomainTet(const char* mesh_file, AffineTransformer<double>* transformer = NULL);
  void Build_VN();
  void BuildEdgeTetIndex();
  void BuildVertexAdjacencyList(); void RenderTet(int tet_id);
  void AddFictitiousForceAndGravity(const std::vector<Vec3>& acceleration,
                                    const std::vector<Vec3>& angular_acceleration,
                                    MultiDomainTet::Vec &subspace_rhs);
  void MixSimulaitonOneStepWithCubature(double dt);
  void SimulateCoupledFullRigid(double dt);
  void SimulateMixed(double dt);
  void SimulatePBD(double dt);
  void SimulateFull(double dt);
  void Cut(Vec3 pa, Vec3 pb, Vec3 pc);
  void SetFullSimulationDomains(std::unordered_set<int> domains);
  void SimulationRigidMotion(const double dt, Vec & new_rigid_velocity, const std::vector<ExtForce>* ext_force = nullptr);
  void RenderSurface();
  int Intersect(double * ray_start, double * ray_end,
                double * clicked_world_pos,
                double * selected_pos);
  void TetLimiting(double dt, int iteration);
  inline void UpdateEmbededVertexPosition();
  inline void UpdateBasisOffSet(bool using_pbd = true);
  void UpdateSurfaceTriangleTopology();
  void BuildTetSurfaceTriIndex();
  inline void SplitTriangle(int tri_idx, int v0, int v1, int v2,
                            int v10, int v11, int v20, int v21);
  virtual void SetFixedDomains(std::set<int> fixed_domains);
  ~MixedMultiDomainTet();
  virtual void Render();
  virtual void UpdateSurfaceMeshVBO();
  virtual void EnalbeSurfaceRenderingWithVBO(const char* shader_file);
  virtual void EnalbeTexturedSurfaceRendering(const char* shader_file);
  virtual void UpdateRenderData();
  void UpdateTextureCoordinateVBO();

  BiconjugateGradientStablized<double>* bi_cg_;
  std::vector<int> constrained_domain_offset_;
  std::vector<int> constrained_rigid_offset_;
  std::vector<int> domain_size_;

  std::vector<int> domain_offset_;
  std::vector<int> rigid_offset_;

  int subspace_domain_basis_size_;
  int subspace_domain_num_;
  std::vector<char> is_cut_triangle_;
  std::vector<double> triangle_texture_coordinate_;
  std::vector<std::pair<int, int>> triangle_group_info_; // (group_id, idx within group)
  std::vector<std::pair<int, int> > subspace_subspace_interface_;
  std::vector<std::pair<int, int> > full_subspace_interface_;
  std::vector<int> domain_index_;
  std::vector<int> is_subspace_domain_;

  // For each subspace domain, list all adjacent full domains and adajent vertices in the full domain
  std::vector<std::vector<double> > part_vel_;
  std::vector<std::vector<double> > part_q_;

  std::vector<CubaturePoint> subspace_cubature_; // actuall cubature tets and tets on full-subspace tet
  std::vector<std::vector<std::pair<int, std::vector<int> > > > full_vert_on_domain_;
  std::vector<std::vector<std::vector<std::vector<std::pair<int, int> > > > > vert_incident_tet_domain_;
  std::vector<std::vector<std::vector<Mat> > > full_vert_subspace_matrix_block_;
  std::vector<int> full_subspace_interface_tet_;
  std::vector<int> is_reduced_vertex_;
  std::vector<int> subspace_domains_;
  std::vector<int> full_simulation_domains_;
  std::vector<int> full_tet_;
  std::vector<Mat> domain_k_;
  std::vector<Mat> interface_k_;

  std::vector<int> full_verts_;
  std::vector<int> full_vert_idx_;

  int original_vertex_num_;
  int original_edge_num_;
  int original_tet_num_;
  // vertex data
  int* vertex_id_;
  std::vector<Map> adjacent_vertex_;

  // edge data
  double* mid_point_weight_;
  std::vector<std::unordered_set<int> > incident_tet_on_edge_;
  std::unordered_map<std::pair<int, int>, int, IntPairHash> vertex_pair2edge_;

  int embeded_v_num_;
  std::vector<int> embeded_vertex_;
  std::vector<EmbededEdge> embeded_edge_;

  int* original_edge_;
  int* edge_id_;
  // triangle data
  // tet data
  std::vector<Vec4i> surface_tri_on_tet_;
  int* tet_type_;
  int* edge_on_tet_;
  int* tet_id_;
  MixedSparseMatrix* sparse_k_;

  int constrained_dof_;
  int non_constrained_dof_;
  int* global_vertex_id2surface_vertex_id_;
  char* is_surface_vert_;
  int surface_vertex_num_;
private:
  static const int kMaxVertexNum;
  static const int kMaxEdgeNum;
  static const int kMaxTriNum;
  static const int kMaxTetNum;
};

#endif // MIXEDMULTIDOMAINTET_H
