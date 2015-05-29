#include <algorithm>
#include <unordered_set>
#include <map>
#include <set>
#include "tetrahedral_mesh.h"
#include "opengl_helper.h"
#include "rainbow_color.h"
#include "vector_io.h"
#include "affine_transformer.h"
#include "tetrahedral_mesh_io.h"


TetrahedralMesh::TetrahedralMesh(TetrahedralMeshIO *mesh_io, const char *file_name, AffineTransformer<Real> *transformer) {
  std::vector<const char*> file_names(1, file_name);
  std::vector<AffineTransformer<Real>*> transformers(1, transformer);
  Construct(mesh_io, file_names, transformers);
}

TetrahedralMesh::TetrahedralMesh(TetrahedralMeshIO *mesh_io,
                                 std::vector<const char *> file_names,
                                 std::vector<AffineTransformer<Real> *> transformers) {
  Construct(mesh_io, file_names, transformers);
}

void TetrahedralMesh::Construct(TetrahedralMeshIO *mesh_io, std::vector<const char *> file_names, std::vector<AffineTransformer<Real> *> transformers) {
  render_mode_ = 0;
  std::vector<Real> vert;
  std::vector<int> tet;
  v_num_ = 0;
  tet_num_ = 0;
  for (int i = 0; i < (int) file_names.size(); ++i) {
    mesh_io->Read(file_names[i], vert, tet);
    mesh_io->OrientTetrahedron(vert, tet);
    int old_v_num = v_num_;
    int old_tet_num = tet_num_;
    v_num_ += int(vert.size()) / 3;
    tet_num_ += int(tet.size()) / 4;
    vert_.resize(v_num_);
    tet_.resize(tet_num_);
    if (transformers[i]) {
      transformers[i]->Transform(&vert[0], int(vert.size()) / 3);
    }
    for (int i = 0; i < (int) vert.size(); i += 3) {
      vert_[old_v_num + i / 3][0] = vert[i + 0];
      vert_[old_v_num + i / 3][1] = vert[i + 1];
      vert_[old_v_num + i / 3][2] = vert[i + 2];
    }
    for (int i = 0; i < (int) tet.size(); i += 4) {
      tet_[old_tet_num + i / 4][0] = tet[i + 0] + old_v_num;
      tet_[old_tet_num + i / 4][1] = tet[i + 1] + old_v_num;
      tet_[old_tet_num + i / 4][2] = tet[i + 2] + old_v_num;
      tet_[old_tet_num + i / 4][3] = tet[i + 3] + old_v_num;
    }
  }
  v_normal_.resize(v_num_);
  BuildTriangles();
  BuildSurfaceVertex();
  BuildEdges();
  P(v_num_, tri_num_, e_num_, tet_num_);

  //#define SINGLE_EDGE
#ifdef SINGLE_EDGE
  e_num_ = 1;
  edge_.resize(1);
  vert_.resize(2);
  incident_edge_.resize(2);
  incident_edge_[0] = std::vector<int>(1, 0);
  incident_edge_[1] = std::vector<int>(1, 0);
  v_num_ = 2;
  edge_[0][0] = 0;
  edge_[0][1] = 1;
#endif

  //#define SINGLE_TRIANGLE
#ifdef SINGLE_TRIANGLE
  e_num_ = 3;
  edge_.resize(1);
  vert_.resize(3);
  incident_edge_.resize(3);
  incident_edge_[0] = std::vector<int>({0, 2});
  incident_edge_[1] = std::vector<int>({0, 1});
  incident_edge_[2] = std::vector<int>({1, 2});
  v_num_ = 3;
  vert_[0] = Vec3(0, 1, 0);
  vert_[1] = Vec3(-0.5, 0, 0);
  vert_[2] = Vec3(+0.5, 0, 0);
  edge_[0][0] = 0;
  edge_[0][1] = 1;
  edge_[1][0] = 1;
  edge_[1][1] = 2;
  edge_[2][0] = 2;
  edge_[2][1] = 0;
#endif
}

TetrahedralMesh &TetrahedralMesh::operator=(const TetrahedralMesh &mesh) {
  v_num_ = mesh.v_num_;
  vert_ = mesh.vert_;
  v_normal_ = mesh.v_normal_;
  incident_edge_ = mesh.incident_edge_;
  incident_tet_ = mesh.incident_tet_;

  tet_ = mesh.tet_;
  tet_num_ = mesh.tet_num_;

  tri_num_ = mesh.tri_num_;
  tri_ = mesh.tri_;
  tri_normal_ = mesh.tri_normal_;

  e_num_ = mesh.e_num_;
  edge_ = mesh.edge_;
  return *this;
}

int TetrahedralMesh::Select(Real *ray_start, Real *ray_end,
                            Real *clicked_world_pos, Real *selected_pos) {
  (void) clicked_world_pos;
  Real start[3] = {ray_start[0], ray_start[1], ray_start[2]};
  Real end[3] = {ray_end[0], ray_end[1], ray_end[2]};
  const Real kThreshold = 0.001;
  int selected_vertex = -1;
  std::vector<std::pair<Real, int> > intersect_triangles;
  for (unsigned int tri_idx = 0; tri_idx < tri_.size(); ++tri_idx) {
    Real u, v, t;
    bool intersect = dj::SegmentTriangleIntersection<Real>(start, end,
                                                           &vert_[tri_[tri_idx][0]][0],
                                                           &vert_[tri_[tri_idx][1]][0],
                                                           &vert_[tri_[tri_idx][2]][0],
                                                           kThreshold, u, v, t);
    //    Real length;
    //    Real distance = dj::PointLineDistance<real>(start, end, &vert_[v][0], &length);
    if (intersect) {
      intersect_triangles.push_back(std::pair<Real, int>(t, tri_idx));
    }
  }

  if (intersect_triangles.size() > 0) {
    std::sort(intersect_triangles.begin(), intersect_triangles.end());
    int tri_idx =  intersect_triangles.front().second;
    Real min_dist[3] = {
      dj::PointLineDistance<Real>(start, end, &vert_[tri_[tri_idx][0]][0]),
      dj::PointLineDistance<Real>(start, end, &vert_[tri_[tri_idx][1]][0]),
      dj::PointLineDistance<Real>(start, end, &vert_[tri_[tri_idx][2]][0])
    };
    int min_idx = 0;
    if (min_dist[1] < min_dist[0]) {
      min_idx = 1;
    }
    if (min_dist[2] < min_dist[min_idx]) {
      min_idx = 2;
    }
    selected_vertex = tri_[tri_idx][min_idx];
    selected_pos[0] = vert_[selected_vertex][0];
    selected_pos[1] = vert_[selected_vertex][1];
    selected_pos[2] = vert_[selected_vertex][2];
    return selected_vertex;
  } else {
    return -1;
  }
}

bool TetrahedralMesh::GetVertexPosition(int vertex_idx, Real *pos) {
  if (vertex_idx >= 0 && vertex_idx < v_num_) {
    pos[0] = vert_[vertex_idx][0];
    pos[1] = vert_[vertex_idx][1];
    pos[2] = vert_[vertex_idx][2];
    return true;
  } else {
    return false;
  }
}

void TetrahedralMesh::SaveSurface2Obj(const char *file_name) {
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  std::vector<int> global_vert2surface_vert(v_num_, -1);
  int i = 0;
  for (int v : surface_vert_) {
    global_vert2surface_vert[v] = i;
    out << "v "
        << vert_[v][0] << " "
        << vert_[v][1] << " "
        << vert_[v][2] << "\n";
    ++i;
  }
  for (int i = 0; i < int(tri_.size()); i++) {
    int verts[3] = {
      global_vert2surface_vert[tri_[i][0]] + 1,
      global_vert2surface_vert[tri_[i][1]] + 1,
      global_vert2surface_vert[tri_[i][2]] + 1,
    };
    out << "f " << verts[0] << " " << verts[1] << " " << verts[2] << "\n";
  }
  out.close();
}

void TetrahedralMesh::ComputerVertexNormal() {
  ComputeTriangleNormal();
  memset(&v_normal_[0][0], 0, sizeof(Vec3) * v_num_);
  for (int t = 0; t < tri_num_; ++t) {
    v_normal_[tri_[t][0]] += tri_normal_[t];
    v_normal_[tri_[t][1]] += tri_normal_[t];
    v_normal_[tri_[t][2]] += tri_normal_[t];
  }
}

void TetrahedralMesh::ComputeTriangleNormal() {
  OMP_FOR
  for (int t = 0; t < tri_num_; ++t) {
    Vec3* pos[3] = {
      &vert_[tri_[t][0]],
      &vert_[tri_[t][1]],
      &vert_[tri_[t][2]]
    };
    Vec3 tri_vec[2] = {
      *(pos[1]) - *(pos[0]),
      *(pos[2]) - *(pos[0]),
    };
    dj::Cross3(&tri_vec[0][0], &tri_vec[1][0], &tri_normal_[t][0]);
    tri_normal_[t].Normalize();
  }
}

void TetrahedralMesh::BuildTriangles() {
  std::map<dj::Vec3i, int> surf_tri;
  std::map<dj::Vec3i, dj::Vec3i> ordered_tri;
  for (int i = 0; i < tet_num_; i++) {
    const int idx[4][3] = {
      0, 2, 1,
      0, 3, 2,
      0, 1, 3,
      1, 2, 3,
    };
    for (int j = 0; j < 4; ++j) {
      dj::Vec3i tri(tet_[i][idx[j][0]], tet_[i][idx[j][1]], tet_[i][idx[j][2]]);
      dj::Vec3i sorted_tri(tri);
      std::sort(sorted_tri.begin(), sorted_tri.end());
      ++surf_tri[sorted_tri];
      ordered_tri[sorted_tri] = tri;
    }
  }
  tri_.clear();;
  for (auto iter = surf_tri.begin(); iter != surf_tri.end(); ++iter) {
    if (iter->second == 1) {
      tri_.push_back(ordered_tri[iter->first]);
    }
  }
  tri_num_ = (int) tri_.size();
  tri_normal_.resize(tri_num_);
}

void TetrahedralMesh::BuildSurfaceVertex() {
  std::unordered_set<int> surface_vert_set;
  if (tri_num_ <= 0) {
    std::cerr << CURRENT_LINE << " => triangle number is " << tri_num_ << std::endl;
    return;
  }
  for (int t = 0; t < tri_num_; ++t) {
    surface_vert_set.insert(tri_[t][0]);
    surface_vert_set.insert(tri_[t][1]);
    surface_vert_set.insert(tri_[t][2]);
  }
  surface_vert_.clear();
  surface_vert_.insert(surface_vert_.end(), surface_vert_set.begin(), surface_vert_set.end());
}

void TetrahedralMesh::BuildEdges() {
  const int kEdgePerTet = 6;
  std::vector<dj::Vec2i> temp(tet_num_ * kEdgePerTet);
  for (int i = 0; i < tet_num_; i++) {
    const int idx[kEdgePerTet][2] = {
      0, 1,
      0, 2,
      0, 3,
      1, 2,
      1, 3,
      2, 3,
    };
    for (int j = 0; j < kEdgePerTet; ++j) {
      ASSERT(tet_[i][idx[j][0]] >= 0, P(i, j));
      ASSERT(tet_[i][idx[j][1]] >= 0);
      temp[i * kEdgePerTet + j][0] = tet_[i][idx[j][0]];
      temp[i * kEdgePerTet + j][1] = tet_[i][idx[j][1]];
      if (temp[i * kEdgePerTet + j][0] > temp[i * kEdgePerTet  + j][1]) {
        std::swap(temp[i * kEdgePerTet + j][0], temp[i * kEdgePerTet  + j][1]);
      }
    }
  }
  std::sort(temp.begin(), temp.end());
  std::vector<dj::Vec2i>::iterator end = std::unique(temp.begin(), temp.end());
  edge_.clear();
  edge_.insert(edge_.end(), temp.begin(), end);
  e_num_ = (int) edge_.size();

  incident_edge_.resize(v_num_);
  for (int e = 0; e < e_num_; ++e) {
    incident_edge_[edge_[e][0]].push_back(e);
    incident_edge_[edge_[e][1]].push_back(e);
  }
}

//#include <QString>
//#include "open_gl_qt.h"
void TetrahedralMesh::Render(int mode) {
  int render_mode = (mode < 0) ? kRenderMode[render_mode_] : mode;

//  glPointSize(6.0);
//  glBegin(GL_POINTS);
//  int v = 6688;
//  glColor3fv(kYellow());
//  Vertex3v(&vert_[v][0]);
//  glEnd();
  //  int ss = 0;
  //    if (ss >= 0) {
  //      glPointSize(8);
  //      glBegin(GL_POINTS);
  //      glColor3fv(kRed());
  //      Vertex3v(&vert_[ss][0]);
  //      glColor3fv(kGreen());
  //      Vertex3v(&vert_[1][0]);
  //      glColor3fv(kBlue());
  //      Vertex3v(&vert_[2][0]);
  //      glEnd();
  //    }
  //    for (int e = 0; e < e_num_; ++e) {
  //      auto center = Real(0.5) * (vert_[edge_[e][0]] + vert_[edge_[e][1]]);
  //      QString s = QString("%1").arg(e);
  //      global::gl->renderText(center[0], center[1], center[2], s);
  //    }
  //  glColor3fv(kRed());
  //  for (int v = 0; v < v_num_; ++v) {
  //    QString s = QString("%1").arg(v);
  //    global::gl->renderText(vert_[v][0], vert_[v][1], vert_[v][2], s);
  //  }

  //  if (0)
  if ((render_mode & kSmoothLighting) || (render_mode & kFlatLighting)) {
    ComputerVertexNormal();
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LIGHTING_BIT);
    glPushAttrib(GL_POLYGON_BIT);
    float diffuse_color[4] = {1, 1, 1, 0.6f};
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    float dark[4] = {0, 0, 0, 1};
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, dark);
    if (render_mode & kFlatLighting) {
      glShadeModel(GL_FLAT);
    } else {
      glShadeModel(GL_SMOOTH);
    }
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
    glEnable(GL_NORMALIZE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < tri_num_; i++) {
      Normal(v_normal_[tri_[i][0]]());
      Vertex3v(&vert_[tri_[i][0]][0]);
      Normal(v_normal_[tri_[i][1]]());
      Vertex3v(&vert_[tri_[i][1]][0]);
      Normal(v_normal_[tri_[i][2]]());
      Vertex3v(&vert_[tri_[i][2]][0]);
    }
    glEnd();
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
  }

  if (render_mode & kSurfaceEdge) {
    glColor3fv(kBlack());
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(0.3f);
    //    glLineStipple(3, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    for (int i = 0; i < tri_num_; ++i) {
      Vertex3v(vert_[tri_[i][0]]());
      Vertex3v(vert_[tri_[i][1]]());
      Vertex3v(vert_[tri_[i][0]]());
      Vertex3v(vert_[tri_[i][2]]());
      Vertex3v(vert_[tri_[i][1]]());
      Vertex3v(vert_[tri_[i][2]]());
    }
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glPopAttrib();
    glPopAttrib();
  }

  if (render_mode & kWireFrame) {
    glColor3fv(kBlack());
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(0.2f);
    //    glLineStipple(3, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    for (int e = 0; e < e_num_; ++e) {
      Vertex3v(vert_[edge_[e][0]]());
      Vertex3v(vert_[edge_[e][1]]());
    }
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glPopAttrib();
    glPopAttrib();
  }

  if (render_mode & kTransparent) {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    glDrawBuffer(GL_NONE);
    Render(kSmoothLighting);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDrawBuffer(GL_BACK);
    glEnable(GL_LIGHTING);
    glColor3f(0.0, 0.0, 0.0);
    Render(kSmoothLighting);
    glDisable(GL_BLEND);
  }

  {
#if DRAW_NORMAL
    ComputerVertexNormal();
    glColor3fv(kRed());
    for (int i = 0; i < tri_num_; ++i) {
      //      for (int i : {2428, 2433}) {
      Vec3& v0 = vert_[tri_[i][0]];
      Vec3& v1 = vert_[tri_[i][1]];
      Vec3& v2 = vert_[tri_[i][2]];
      Vec3 center = v0 + v1 + v2;
      center *= (Real(1) / 3);
      //      glVertex3dv(center());
      //      Vertex3v(center());
      //      P(center);
      DrawArrow(center(), tri_normal_[i](), true, 0.1, 0.01, 0.01);
    }

    glColor3fv(kRed());
    glDisable(GL_LIGHTING);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < tri_num_; ++i) {
      //      for (int i : {2428, 2433}) {
      Vec3& v0 = vert_[tri_[i][0]];
      Vec3& v1 = vert_[tri_[i][1]];
      Vec3& v2 = vert_[tri_[i][2]];
      Vec3 center = v0 + v1 + v2;
      center *= (Real(1) / 3);
      glVertex3dv(center());
    }
    glEnd();
#endif
  }
}

void TetrahedralMesh::LoadPosition(std::vector<Real> &pos) {
  ASSERT(int(pos.size()) == v_num_ * 3);
  for (int v = 0; v < v_num_; ++v) {
    vert_[v][0] = pos[v * 3 + 0];
    vert_[v][1] = pos[v * 3 + 1];
    vert_[v][2] = pos[v * 3 + 2];
  }
}

void TetrahedralMesh::NextRenderMode() {
  render_mode_ = (render_mode_ + 1) % kRenderModeNum;
}

void TetrahedralMesh::PrevRenderMode() {
  render_mode_--;
  if (render_mode_ < 0) render_mode_ = kRenderModeNum - 1;
}

const int TetrahedralMesh::kRenderMode[] = {
  TetrahedralMesh::kSmoothLighting,
  TetrahedralMesh::kSurfaceEdge | TetrahedralMesh::kSmoothLighting,
  TetrahedralMesh::kWireFrame | TetrahedralMesh::kSmoothLighting,
  TetrahedralMesh::kFlatLighting,
  TetrahedralMesh::kFlatLighting | TetrahedralMesh::kSurfaceEdge,
  TetrahedralMesh::kSurfaceEdge,
  TetrahedralMesh::kWireFrame,
  TetrahedralMesh::kTransparent,
};
const int TetrahedralMesh::kRenderModeNum = sizeof(TetrahedralMesh::kRenderMode) / sizeof(int);
