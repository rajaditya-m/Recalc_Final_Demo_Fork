#ifndef SIMPLE_TETRAHEDRA_MESH_H
#define SIMPLE_TETRAHEDRA_MESH_H
#include <vector>
#include "global.h"
#include "vector_lib.h"
#include "selectable_object.h"

class TetrahedralMeshIO;
template <class T> class AffineTransformer;
class TetrahedralMesh : public SelectableObject<Real>
{
public:
  enum MeshRenderMode {
    kWireFrame = 1 << 1,
    kSmoothLighting = 1 << 2,
    kFlatLighting = 1 << 3,
    kSurfaceEdge = 1 << 4,
    kTransparent = 1 << 5,
  };
  static const int kRenderModeNum;
  static const int kRenderMode[];

  typedef dj::Vec<Real, 3, false> Vec3;
  TetrahedralMesh(TetrahedralMeshIO *mesh_io,
                  const char* file_name,
                  AffineTransformer<Real>* transformer = NULL);

  TetrahedralMesh(TetrahedralMeshIO *mesh_io,
                  std::vector<const char*>  file_names,
                  std::vector<AffineTransformer<Real>*> transformers);

  void Construct(TetrahedralMeshIO *mesh_io,
                 std::vector<const char*>  file_names,
                 std::vector<AffineTransformer<Real>*> transformers);

  TetrahedralMesh(const TetrahedralMesh& mesh) : SelectableObject<Real>(mesh)
  {
    *this = mesh;
  }

  int VertexNum()
  {
    return v_num_;
  }
  int TetNum()
  {
    return tet_num_;
  }
  Real* VertexArray()
  {
    return &vert_[0][0];
  }
  int* TetIndex()
  {
    return &tet_[0][0];
  }

  TetrahedralMesh& operator=(const TetrahedralMesh& mesh);
  virtual int Select(Real* ray_start, Real* ray_end, Real* clicked_world_pos, Real* selected_pos);
  virtual bool GetVertexPosition(int vertex_idx, Real *pos);
  virtual ~TetrahedralMesh() {}
  void SaveSurface2Obj(const char* file_name);
  void ComputerVertexNormal();
  void ComputeTriangleNormal();
  void BuildTriangles();
  void BuildSurfaceVertex();
  void BuildEdges();
  void Render(int mode = -1);
  void LoadPosition(std::vector<Real>& pos);
  void NextRenderMode();
  void PrevRenderMode();

  int render_mode_;
  // Vertex
  int v_num_;
  std::vector<Vec3> vert_;
  std::vector<Vec3> v_normal_;
  std::vector<std::vector<int> > incident_edge_;
  std::vector<std::vector<int> > incident_tet_;
  std::vector<int> surface_vert_;

  // Tetrahedra
  std::vector<dj::Vec4i> tet_;
  int tet_num_;

  // Triangle
  int tri_num_;
  std::vector<dj::Vec3i> tri_;
  std::vector<Vec3> tri_normal_;

  // Edge
  int e_num_;
  std::vector<dj::Vec2i> edge_;
};

#endif // SIMPLE_TETRAHEDRA_MESH_H
