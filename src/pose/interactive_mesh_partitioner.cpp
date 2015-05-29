#include "interactive_mesh_partitioner.h"
#include "opengl_helper.h"
#include "rainbow_color.h"
#include <fstream>
#include "string_formatter.h"
#include "global.h"

InteractiveMeshPartitioner::InteractiveMeshPartitioner(Tet *tet)
{
  tet_ = tet;
  mouse_state_ = kDefault;
  current_group_num_ = 0;
  render_mode_ = 0;
  vertex_group_ = std::vector<int>(tet_->vertex_num_, kSelected);
  P(kRenderModeNum);
}

int InteractiveMeshPartitioner::HandleMouseRelease(QMouseEvent *e) { (void) e; return -1;}

int InteractiveMeshPartitioner::HandleKeyPress(QKeyEvent *e)
{
  switch (e->key()) {
    // change render mode
    case Qt::Key_E:
      render_mode_ = (render_mode_ + 1) % kRenderModeNum;
      break;
    // put selected vertex into new partition
    case Qt::Key_S:
      OMP_FOR
      for (int v = 0; v < tet_->vertex_num_; ++v) {
        if (vertex_group_[v] == kSelected)  {
          vertex_group_[v] = current_group_num_;
        } else if (vertex_group_[v] == kUnSelected) {
          vertex_group_[v] = kSelected;
        }
      }
      P(current_group_num_);
      current_group_num_++;
      break;
    // set all unselected vertex to be selected
    case Qt::Key_R:
      OMP_FOR
      for (int v = 0; v < tet_->vertex_num_; ++v) {
        if (vertex_group_[v] == kUnSelected)  {
          vertex_group_[v] = kSelected;
        }
      }
      mouse_state_ = 0;
      break;
    case Qt::Key_D:
//      SavePartition(DATA_DIRECTORY "partition.txt");
      SavePartition(dj::Format("%s/vert_partition.txt", GetPartitionFolder()).c_str());
      break;
    default:
      break;
  }
  return -1;
}

int InteractiveMeshPartitioner::HandleMousePress(QMouseEvent *e)
{
  if (e->modifiers() & Qt::ControlModifier && e->button() == Qt::LeftButton) {
    QPoint pos = e->pos();
    double pp[3] = {0, 0, 0};
    double depth = GetPointDepth(pp);
    GetSelectionRay(pos.x(), pos.y(), &ray_start_[mouse_state_][0], &ray_end_[mouse_state_][0]);
    GetPixelWorldPosition(pos.x(), pos.y(), depth, &clicked_points_[mouse_state_][0]);
    //    P(clicked_points_[mouse_state_]);
    mouse_state_ = (mouse_state_ + 1) % kStateNum;
    if (mouse_state_ == kThirdLeftButton) {
      Vec3 normal;
      Vec3 v0 = ray_end_[0] - ray_start_[0];
      Vec3 v1 = ray_end_[1] - ray_start_[1];
      dj::Cross3(&v0[0], &v1[0], &normal[0]);
      normal.Normalize();
      Vec3 v_ref = clicked_points_[2] - clicked_points_[0];
      double sign = v_ref * normal;
      OMP_FOR
      for (int v = 0; v < tet_->vertex_num_; ++v) {
        if (vertex_group_[v] >= 0 || vertex_group_[v] == kUnSelected) continue;
        Vec3 pos(tet_->X + v * 3);
        Vec3 dir = pos - clicked_points_[0];
        double dot = dir * normal;
        if (dot * sign >= 0) {
        } else {
          vertex_group_[v] = kUnSelected;
        }
      }
      mouse_state_ = 0;
    }
  }
  return kNotHandled;
}

int InteractiveMeshPartitioner::Render()
{
  const float* color_map[] = {
    kRed(),
    kGreen(),
    kBlue(),
    kYellow(),
    kOrage(),
    kChocolate(),
    kViolet(),
    kIndigo(),
  };
  const int color_num = sizeof(color_map) / sizeof(double*);
  glColor3fv(kBlack());
  if (mouse_state_ == kDefault) {
  } else if (mouse_state_ == kFirstLeftButton) {
    glBegin(GL_POINTS);
    Vertex3v(&clicked_points_[0][0]);
    glEnd();
    glBegin(GL_LINES);
    Vertex3v(&ray_start_[0][0]);
    Vertex3v(&ray_end_[0][0]);
    glEnd();
  } else if (mouse_state_ == kSecondLeftButton) {
    glBegin(GL_POINTS);
    Vertex3v(&clicked_points_[0][0]);
    Vertex3v(&clicked_points_[1][0]);
    glEnd();
    glBegin(GL_LINES);
    Vertex3v(&ray_start_[0][0]);
    Vertex3v(&ray_end_[0][0]);
    Vertex3v(&ray_start_[1][0]);
    Vertex3v(&ray_end_[1][0]);
    glEnd();
  } else if (mouse_state_ == kThirdLeftButton) {
    glBegin(GL_POINTS);
    Vertex3v(&clicked_points_[0][0]);
    Vertex3v(&clicked_points_[1][0]);
    Vertex3v(&clicked_points_[2][0]);
    glEnd();
    glBegin(GL_LINES);
    Vertex3v(&ray_start_[0][0]);
    Vertex3v(&ray_end_[0][0]);
    Vertex3v(&ray_start_[1][0]);
    Vertex3v(&ray_end_[1][0]);
    Vertex3v(&ray_start_[2][0]);
    Vertex3v(&ray_end_[2][0]);
    glEnd();
  } else {
  }

  int mode = kRenderModes[render_mode_];
  //  P(render_mode_, mode);
  glPointSize(6.0);
  if (mode & kRenderSelectedVertex) {
    glColor3fv(kBlack());
    glBegin(GL_POINTS);
    for (int v = 0; v < tet_->vertex_num_; ++v) {
      if (vertex_group_[v] == kSelected) {
        Vertex3v(tet_->X + v * 3);
      }
    }
    glEnd();
  }

  if (mode & kRenderUnSelectedVertex) {
    glColor3fv(kGrey());
    glBegin(GL_POINTS);
    for (int v = 0; v < tet_->vertex_num_; ++v) {
      if (vertex_group_[v] == kUnSelected) {
        Vertex3v(tet_->X + v * 3);
      }
    }
    glEnd();
  }

  if (mode & kRenderPartitionedVertex) {
    glBegin(GL_POINTS);
    for (int v = 0; v < tet_->vertex_num_; ++v) {
      if (vertex_group_[v] >= 0) {
        glColor3fv(color_map[vertex_group_[v] % color_num]);
        Vertex3v(tet_->X + v * 3);
      }
    }
    glEnd();
  }


  if (mode & kRenderTetMesh) {
    double* X = tet_->X;
    glBegin(GL_TRIANGLES);
    for (int t = 0; t < tet_->tet_number; ++t) {
      int* verts = tet_->tet_ + t * 4;
      if (vertex_group_[verts[0]] == kSelected) {
        if (vertex_group_[verts[1]] != kSelected) continue;
        if (vertex_group_[verts[2]] != kSelected) continue;
        if (vertex_group_[verts[3]] != kSelected) continue;
        glColor3fv(kBlue());
      }

      if (vertex_group_[verts[0]] == kUnSelected) {
        if (vertex_group_[verts[1]] != kUnSelected) continue;
        if (vertex_group_[verts[2]] != kUnSelected) continue;
        if (vertex_group_[verts[3]] != kUnSelected) continue;
        glColor3fv(kGrey());
      }

      if (vertex_group_[verts[0]] >= 0) {
        int group = vertex_group_[verts[0]];
        if (vertex_group_[verts[1]] != group) continue;
        if (vertex_group_[verts[2]] != group) continue;
        if (vertex_group_[verts[3]] != group) continue;
        glColor3fv(color_map[group % color_num]);
      }

      Vertex3v(X + verts[0] * 3);
      Vertex3v(X + verts[1] * 3);
      Vertex3v(X + verts[2] * 3);

      Vertex3v(X + verts[0] * 3);
      Vertex3v(X + verts[1] * 3);
      Vertex3v(X + verts[3] * 3);

      Vertex3v(X + verts[1] * 3);
      Vertex3v(X + verts[2] * 3);
      Vertex3v(X + verts[3] * 3);

      Vertex3v(X + verts[0] * 3);
      Vertex3v(X + verts[2] * 3);
      Vertex3v(X + verts[3] * 3);
    }
    glEnd();
  }
  return -1;
}

void InteractiveMeshPartitioner::SavePartition(const char *file_name)
{
  std::ofstream out(file_name);
  out << tet_->vertex_num_ << std::endl;
  for (int v = 0; v < tet_->vertex_num_; ++v) {
    ASSERT(vertex_group_[v] >= 0);
    out << vertex_group_[v] << std::endl;
  }
  out.close();
}

const int InteractiveMeshPartitioner::kRenderModes[] = {
  InteractiveMeshPartitioner::kRenderSelectedVertex,
  InteractiveMeshPartitioner::kRenderPartitionedVertex,
  InteractiveMeshPartitioner::kRenderUnSelectedVertex,
  InteractiveMeshPartitioner::kRenderTetMesh,
  InteractiveMeshPartitioner::kRenderSelectedVertex | InteractiveMeshPartitioner::kRenderPartitionedVertex,
  InteractiveMeshPartitioner::kRenderSelectedVertex | InteractiveMeshPartitioner::kRenderPartitionedVertex | InteractiveMeshPartitioner::kRenderUnSelectedVertex,
};
const int InteractiveMeshPartitioner::kRenderModeNum = sizeof(InteractiveMeshPartitioner::kRenderModes) / sizeof(int);
