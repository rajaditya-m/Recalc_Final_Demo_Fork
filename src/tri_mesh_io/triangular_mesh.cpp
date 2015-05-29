#include "triangular_mesh.h"
#include <set>
#include "affine_transformer.h"
#include "triangular_mesh_io.h"
#include "opengl_helper.h"
#include "rainbow_color.h"
#include "vector_lib.h"

TriangularMesh::TriangularMesh(TriangularMeshIO *mesh_io, const char *mesh_file, AffineTransformer<Real> *transformer)
{
  std::vector<Real> vert;
  std::vector<int> tri;
  mesh_io->Read(mesh_file, vert, tri);
  ASSERT(vert.size() % 3 == 0);
  ASSERT(tri.size() % 3 == 0);
  v_num_ = int(vert.size() / 3);
  tri_num_ = int(tri.size() / 3);
  if (transformer) transformer->Transform(&vert[0], v_num_);
  vert_.resize(v_num_);
  for (int v = 0; v < v_num_; ++v) {
    vert_[v][0] = vert[v * 3 + 0];
    vert_[v][1] = vert[v * 3 + 1];
    vert_[v][2] = vert[v * 3 + 2];
  }
  tri_.resize(tri_num_);
  for (int t = 0; t < tri_num_; ++t) {
    tri_[t][0] = tri[t * 3 + 0];
    tri_[t][1] = tri[t * 3 + 1];
    tri_[t][2] = tri[t * 3 + 2];
    ASSERT(tri_[t][0] >= 0);
    ASSERT(tri_[t][1] >= 0);
    ASSERT(tri_[t][2] >= 0);
  }
  Construct();
  P(mesh_file);
  P(v_num_, tri_num_, e_num_);
  // expand mesh in normal direction
  if (0)
  {
    ComputeVertexNormal();
    double factor = 0.0050;
    for (int v = 0; v < v_num_; ++v) {
      dj::Normalize3(&vert_normal_[v][0]);
      vert_[v] += vert_normal_[v] * factor;
    }
    std::vector<double> verts(&vert_[0][0], &(vert_.back()[0]) + 3);
    std::vector<int> tris(&tri_[0][0], &(tri_[tri_num_ - 1][0]) + 3);
    P(tris.size(), tri_.size());
    mesh_io->Write("/Users/dj/out.obj", verts, tris);
    L("exit");
    exit(0);
  }
}

void TriangularMesh::Construct()
{
  BuildEdges();
  BuildIncidentEdgeOnVertex();
  BuildIncidentTriangleOnVertex();
  tri_normal_.resize(tri_num_);
  vert_normal_.resize(v_num_);
}

void TriangularMesh::BuildEdges()
{
  std::set<std::pair<int, int> > edges;
  for (int t = 0; t < tri_num_; ++t) {
    int three_edge[3][2] = {
      tri_[t][0], tri_[t][1],
      tri_[t][0], tri_[t][2],
      tri_[t][1], tri_[t][2],
    };
    for (int i = 0; i < 3; ++i) {
      if (three_edge[i][0] > three_edge[i][1]) {
        std::swap(three_edge[i][0], three_edge[i][1]);
      }
      edges.insert(std::make_pair(three_edge[i][0], three_edge[i][1]));
    }
  }
  for (const std::pair<int, int>& e : edges) {
    edge_.push_back(Vec2i(e.first, e.second));
  }
  e_num_ = int(edge_.size());
}

void TriangularMesh::BuildIncidentTriangleOnVertex()
{
  incident_tri_on_vert_.resize(v_num_);
  for (int v = 0; v < v_num_; ++v) {
    incident_tri_on_vert_[v].clear();
  }
  for (int t = 0; t < tri_num_; ++t) {
    incident_tri_on_vert_[tri_[t][0]].push_back(t);
    incident_tri_on_vert_[tri_[t][1]].push_back(t);
    incident_tri_on_vert_[tri_[t][2]].push_back(t);
  }
}

void TriangularMesh::BuildIncidentEdgeOnVertex()
{
  incident_edge_on_vert_.resize(v_num_);
  for (int v = 0; v < v_num_; ++v) {
    incident_edge_on_vert_[v].clear();
  }
  for (int e = 0; e < e_num_; ++e) {
    incident_edge_on_vert_[edge_[e][0]].push_back(e);
    incident_edge_on_vert_[edge_[e][1]].push_back(e);
  }
}

void TriangularMesh::ComputeTriNormal()
{
  OMP_FOR
  for (int t = 0; t < tri_num_; ++t) {
    dj::ComputeTriangleNormal<Real>(
      &vert_[tri_[t][0]][0],
      &vert_[tri_[t][1]][0],
      &vert_[tri_[t][2]][0],
      &tri_normal_[t][0]
    );
  }
}

void TriangularMesh::ComputeVertexNormal()
{
  ComputeTriNormal();
  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
    vert_normal_[v] = Vec3(0, 0, 0);
    for (int t : incident_tri_on_vert_[v]) {
      vert_normal_[v] += tri_normal_[t];
    }
    dj::Normalize3(&vert_normal_[v][0]);
  }
}

void TriangularMesh::Render(int mode)
{
  if (mode < 0) mode = render_mode_;
  int render_mode = kRenderMode[mode];
  if (render_mode & kRenderVertex) {
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glColor3fv(kRed());
    glPointSize(5);
    glBegin(GL_POINTS);
    for (int v = 0; v < v_num_; ++v) {
      Vertex3v(&vert_[v][0]);
    }
    glEnd();
    glPopAttrib();
  }

  if (render_mode & kRenderEdge) {
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LINE_WIDTH);
    glDisable(GL_LIGHTING);
    glLineWidth(0.2f);
    glDisable(GL_LIGHTING);
    glColor3fv(kBlack());
    glBegin(GL_LINES);
    for (int e = 0; e < e_num_; ++e) {
      Vertex3v(&vert_[edge_[e][0]][0]);
      Vertex3v(&vert_[edge_[e][1]][0]);
    }
    glEnd();
    glPopAttrib();
    glPopAttrib();
  }

  if ((render_mode & kRenderFlatSurface) || (render_mode & kRenderSmoothSurface)) {
    ComputeVertexNormal();
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LIGHTING_BIT);
    glColor3fv(kBlue());
    float diffuse_color[4] = { 1, 1, 1, 0.6f};
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    float dark[4] = {0, 0, 0, 1};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, dark);
    if (render_mode & kRenderFlatSurface) {
      glShadeModel(GL_FLAT);
    } else {
      glShadeModel(GL_SMOOTH);
    }

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < tri_num_; i++) {
      Normal(&vert_normal_[tri_[i][0]][0]);
      Vertex3v(&vert_[tri_[i][0]][0]);

      Normal(&vert_normal_[tri_[i][1]][0]);
      Vertex3v(&vert_[tri_[i][1]][0]);

      Normal(&vert_normal_[tri_[i][2]][0]);
      Vertex3v(&vert_[tri_[i][2]][0]);
    }
    glEnd();

    glPopAttrib();
    glPopAttrib();
    glDisable(GL_CULL_FACE);
  }
}

void TriangularMesh::NextRenderMode()
{
  render_mode_ = (render_mode_ + 1) % kRenderModeNum;
}

void TriangularMesh::PrevRenderMode()
{
  render_mode_--;
  if (render_mode_ < 0) render_mode_ = kRenderModeNum - 1;
}

void TriangularMesh::LoadPosition(std::vector<Real> &pos)
{
  ASSERT(int(pos.size()) == v_num_ * 3);
  for (int v = 0; v < v_num_; ++v) {
    vert_[v][0] = pos[v * 3 + 0];
    vert_[v][1] = pos[v * 3 + 1];
    vert_[v][2] = pos[v * 3 + 2];
  }
}

TriangularMesh::~TriangularMesh() {}

const int TriangularMesh::kRenderMode[] = {
  TriangularMesh::kRenderSmoothSurface,
  TriangularMesh::kRenderSmoothSurface | TriangularMesh::kRenderEdge,
  TriangularMesh::kRenderSmoothSurface | TriangularMesh::kRenderEdge | TriangularMesh::kRenderVertex,
  TriangularMesh::kRenderSmoothSurface | TriangularMesh::kRenderVertex,
  TriangularMesh::kRenderFlatSurface,
  TriangularMesh::kRenderEdge,
  TriangularMesh::kRenderEdge | TriangularMesh::kRenderVertex,
};

const int TriangularMesh::kRenderModeNum = sizeof(TriangularMesh::kRenderMode) / sizeof(int);
