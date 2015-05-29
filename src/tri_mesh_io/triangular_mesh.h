#ifndef TRIANGULARMESH_H
#define TRIANGULARMESH_H
#include <Eigen/Dense>
#include <vector>
#include "global.h"
//#include "triangular_mesh_io.h"

class TriangularMeshIO;
template <class T> class AffineTransformer;

class TriangularMesh
{
public:
  typedef Eigen::Matrix<Real, 3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Vec3i;
  typedef Eigen::Matrix<int, 2, 1> Vec2i;

  enum MeshRenderMode {
    kRenderEdge = 1 << 0,
    kRenderFlatSurface = 1 << 1,
    kRenderSmoothSurface = 1 << 2,
    kRenderVertex = 1 << 3,
  };
  static const int kRenderMode[];
  static const int kRenderModeNum;

  TriangularMesh(TriangularMeshIO* mesh_io,
                 const char* mesh_file,
                 AffineTransformer<Real>* transformer = NULL);
  virtual ~TriangularMesh();
  void Construct();
  void BuildEdges();
  void BuildIncidentTriangleOnVertex();
  void BuildIncidentEdgeOnVertex();
  void ComputeTriNormal();
  void ComputeVertexNormal();
  void Render(int mode = -1);
  void NextRenderMode();
  void PrevRenderMode();
  void LoadPosition(std::vector<Real>& pos);

  int v_num_;
  int e_num_;
  int tri_num_;
  int render_mode_;

  std::vector<std::vector<int> > incident_edge_on_vert_;
  std::vector<std::vector<int> > incident_tri_on_vert_;

  std::vector<Vec3i> tri_;
  std::vector<Vec3> tri_normal_;
  std::vector<Vec3> vert_;
  std::vector<Vec3> vert_normal_;
  std::vector<Vec2i> edge_;
};

#endif // TRIANGULARMESH_H
