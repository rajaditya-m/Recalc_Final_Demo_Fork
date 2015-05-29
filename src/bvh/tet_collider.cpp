#include <iostream>
#include <string.h>
#include <set>
#include <assert.h>
#include "tet_collider.h"
#include "bvh.h"
//#define DEBUG
#ifdef DEBUG
#include "print_macro.h"
#include "vector_lib.h"
#endif
#include "collision_util.h"



TetCollider::TetCollider(Float *vert, int v_num, int *tet, int tet_num)
{
  Init(vert, v_num, tet, tet_num);
}

TetCollider::TetCollider(TetCollider::Float *vert, int v_num, int *tet, int tet_num, int* edge, int edge_num)
{
  Init(vert, v_num, tet, tet_num, edge, edge_num);
}

void TetCollider::Init(TetCollider::Float *vert, int v_num, int *tet, int tet_num, int *edge, int edge_num)
{
  vert_ = (Float (*)[3]) vert;
  v_num_ = v_num;
  tet_ = (int (*)[4] )tet;
  tet_num_ = tet_num;
  edge_ = (int (*)[2]) edge;
  edge_num_ = edge_num;
  count_.resize(v_num_);
  tmp_buffer_.resize(v_num_ * 3);
  bvh_ = new DeformBVHTree(*this, false);

  vert_tet_collision_threshold_ = 1e-2;
  vert_tet_collision_force_ = 0.010;
  edge_tri_collision_threshold_ = 2e-2;
  edge_tri_collision_force_ = 3e-2;
}

TetCollider::~TetCollider()
{
  delete bvh_;
}

void TetCollider::HandleCollision(TetCollider::Float *pos, TetCollider::Float *new_pos, const int kMaxIteration)
{
  Float* current_pos = &tmp_buffer_[0];
  this->vert_ = (Float (*)[3]) pos;
#ifndef BRUTE_FORCE
  bvh_->Construct();
#endif
  int v = 0, v3 = 0;
  int e = 0, e2 = 0;

  auto VTetCollision = [&](int tet_idx) {
    if (v == tet_[tet_idx][0] || v == tet_[tet_idx][1] || v == tet_[tet_idx][2] || v == tet_[tet_idx][3]) {
      return;
    }
    Float* tet_pos[4] = {
      &current_pos[tet_[tet_idx][0] * 3],
      &current_pos[tet_[tet_idx][1] * 3],
      &current_pos[tet_[tet_idx][2] * 3],
      &current_pos[tet_[tet_idx][3] * 3],
    };
    Float* v_pos = &current_pos[v3 + 0];
    Float* new_v_pos = &new_pos[v3 + 0];
    Float bary[4];
    if (collision_util::IsInsideTet<Float>(tet_pos[0], tet_pos[1], tet_pos[2], tet_pos[3], v_pos, -vert_tet_collision_threshold_, bary)) {
      int min_idx = 0;
      for (int i = 1; i < 4; ++i) {
        if (collision_util::Abs(bary[i]) < collision_util::Abs(bary[min_idx])) {
          min_idx = i;
        }
      }
      bary[min_idx] = -vert_tet_collision_force_;
      Float sum = bary[0] +  bary[1] +  bary[2] +  bary[3];
      sum  = Float(1) / sum;
      bary[0]  = bary[0] * sum;
      bary[1]  = bary[1] * sum;
      bary[2]  = bary[2] * sum;
      bary[3]  = bary[3] * sum;
      Float factor = bary[0] * bary[0] + bary[1] * bary[1] + bary[2] * bary[2] + bary[3] * bary[3];
      Float offset[3];
      Float* new_tet_pos[4] = {
        &new_pos[tet_[tet_idx][0] * 3],
        &new_pos[tet_[tet_idx][1] * 3],
        &new_pos[tet_[tet_idx][2] * 3],
        &new_pos[tet_[tet_idx][3] * 3],
      };
      for (int j = 0; j < 3; ++j) {
        offset[j] = tet_pos[0][j] * bary[0] + tet_pos[1][j] * bary[1] + tet_pos[2][j] * bary[2] + tet_pos[3][j] * bary[3] - v_pos[j];
        new_v_pos[j] += offset[j] / (factor + 1);
        for (int i = 0; i < 4; ++i) {
          if (i != min_idx) {
            new_tet_pos[i][j] -= offset[j] * bary[i] / (factor + 1);
          }
        }
      }
      ++count_[v];
      for (int i = 0; i < 4; ++i) {
        if (i != min_idx) {
          ++count_[tet_[tet_idx][i]];
        }
      }
    }
  };
  (void) VTetCollision;

  auto EdgeTriCollision = [&](int tet_idx) {
    int* edge_verts = &edge_[e][0];
    for (int f = 0; f < 4; ++f) {
      int tri_verts[3] = {
        tet_[tet_idx][f],
        tet_[tet_idx][(f + 1) % 4],
        tet_[tet_idx][(f + 2) % 4],
      };
      if (edge_verts[0] == tri_verts[0] || edge_verts[0] == tri_verts[1] || edge_verts[0] == tri_verts[2] ||
          edge_verts[1] == tri_verts[0] || edge_verts[1] == tri_verts[1] || edge_verts[1] == tri_verts[2]) {
        continue;
      }
      tri_verts[0] *= 3;
      tri_verts[1] *= 3;
      tri_verts[2] *= 3;
      Float t;
      Float bary[3];
      bool has_collision = collision_util::SegmentTriangleIntersection(
                             &current_pos[edge_verts[0] * 3],
                             &current_pos[edge_verts[1] * 3],
                             &current_pos[tri_verts[0]],
                             &current_pos[tri_verts[1]],
                             &current_pos[tri_verts[2]],
                             -edge_tri_collision_threshold_, bary[1], bary[2], t);
      if (has_collision) {
        bary[0] = 1 - bary[1] - bary[2];
        int min_idx = 0;
        if (collision_util::Abs(bary[1]) < collision_util::Abs(bary[min_idx])) {
          min_idx = 1;
        }
        if (collision_util::Abs(bary[2]) < collision_util::Abs(bary[min_idx])) {
          min_idx = 2;
        }
        bary[min_idx] = -edge_tri_collision_force_;
        Float sum = bary[0] + bary[1] + bary[2];
        bary[0] /= sum;
        bary[1] /= sum;
        bary[2] /= sum;
        Float intersection_point[3];
        Float factor_edge = t * t + (1 - t) * (1 - t);
        Float factor_tri = bary[0] * bary[0] + bary[1] * bary[1] + bary[2] + bary[2];
        Float offset[3];
        for (int i = 0; i < 3; ++i) {
          intersection_point[i] = t * current_pos[edge_verts[1] * 3 + i] + (1 - t) * current_pos[edge_verts[0] * 3 + i];
          offset[i] = bary[0] * current_pos[tri_verts[0] + i] + bary[1] * current_pos[tri_verts[1] + i] + bary[2] * current_pos[tri_verts[2] + i]
                      - intersection_point[i];
          new_pos[edge_verts[0] * 3 + i] += (1 - t) / (factor_edge + factor_tri) * offset[i];
          new_pos[edge_verts[1] * 3 + i] += (t) / (factor_edge + factor_tri) * offset[i];
          for (int j = 0; j < 3; ++j) {
            if (j != min_idx) {
              new_pos[tri_verts[j] + i] -= bary[j] / (factor_edge + factor_tri) * offset[i];
            }
          }
        }
        count_[edge_verts[0]]++;
        count_[edge_verts[1]]++;
        count_[tri_verts[(min_idx + 1) % 3] / 3]++;
        count_[tri_verts[(min_idx + 2) % 3] / 3]++;
      }
    }
  };

  int i = 0;
  memcpy(new_pos, pos, sizeof(Float) * 3 * v_num_);
  for (; i < kMaxIteration; ++i) {
    std::swap(current_pos, new_pos);
    bool has_collision = false;
    memset(new_pos, 0, sizeof(Float) * v_num_ * 3);
    std::fill(count_.begin(), count_.end(), 0);
    //------------------------------------------------------------------------------
    // Handle vertex-tet collision
//    for (v = 0, v3 = 0; v < v_num_; ++v, v3 += 3) {
//#ifdef BRUTE_FORCE
//      for (int t = 0; t < tet_num_; ++t) {
//        VTetCollision(t);
//      }
//#else
//      BOX v_box = bvh_->_root->vertBox(v, false);
//      ForOverlappingFaces(bvh_->_root, v_box, VTetCollision);
//#endif
//    }

    //------------------------------------------------------------------------------
    // Handle  edge triangle collision
    for (e = 0, e2 = 0; e < edge_num_; ++e, e2 += 2) {
#ifdef BRUTE_FORCE
      for (int t = 0; t < tet_num_; ++t) {
        EdgeTriCollision(t);
      }
#else
      BOX e_box = bvh_->_root->edgeBox(&vert_[edge_[e][0]][0], &vert_[edge_[e][1]][0]);
      ForOverlappingFaces(bvh_->_root, e_box, EdgeTriCollision);
#endif
    }

    for (v = 0, v3 = 0; v < v_num_; ++v, v3 += 3) {
      if (count_[v] > 0) {
        new_pos[v3 + 0] /= count_[v];
        new_pos[v3 + 0] += current_pos[v3 + 0];
        new_pos[v3 + 1] /= count_[v];
        new_pos[v3 + 1] += current_pos[v3 + 1];
        new_pos[v3 + 2] /= count_[v];
        new_pos[v3 + 2] += current_pos[v3 + 2];
        has_collision = true;
      } else {
        new_pos[v3 + 0] = current_pos[v3 + 0];
        new_pos[v3 + 1] = current_pos[v3 + 1];
        new_pos[v3 + 2] = current_pos[v3 + 2];
      }
    }
    if (!has_collision) break;
    if (i != kMaxIteration - 1) {
      this->vert_ = (Float (*)[3]) new_pos;
      bvh_->refit();
    }
  }

  if (i >= kMaxIteration) {
    std::cerr << "HandleCollision() => collision is not resolved with " << kMaxIteration << " iterations" << std::endl;
  }
  if (current_pos != &tmp_buffer_[0]) {
    assert(new_pos == &tmp_buffer_[0]);
    memcpy(current_pos, new_pos, sizeof(Float) * 3 * v_num_);
  }
}

#if 0
P(dj::Vec4d(bary));
int vert_idx[3] = {
  tet_[tet_idx][idx[3 - min_idx][0]],
  tet_[tet_idx][idx[3 - min_idx][1]],
  tet_[tet_idx][idx[3 - min_idx][2]],
};
Float dist, normal[3];
Float projection_point[3];
dist = collision_util::PointToTriangleDistance(v_pos, &buffer[vert_idx[0] * 3], &buffer[vert_idx[1] * 3], &buffer[vert_idx[2] * 3],
                                               normal, projection_point, &bary[0]);
const Float kThreshold = 0.020;
P(v, min_idx);
P(dj::Vec3d(&buffer[vert_idx[0] * 3]));
P(dj::Vec3d(&buffer[vert_idx[1] * 3]));
P(dj::Vec3d(&buffer[vert_idx[2] * 3]));
P(dj::Vec3i(vert_idx));
P(dj::Vec3d(normal));
Float offset[3] = {
  projection_point[0] + normal[0] * kThreshold - v_pos[0],
  projection_point[1] + normal[1] * kThreshold - v_pos[1],
  projection_point[2] + normal[2] * kThreshold - v_pos[2],
};
Float* new_tet_pos[3] = {
  &new_pos[vert_idx[0] * 3],
  &new_pos[vert_idx[1] * 3],
  &new_pos[vert_idx[2] * 3],
};
Float factor = bary[0] * bary[0] + bary[1] * bary[1] + bary[2] * bary[2];
for (int j = 0; j < 3; ++j)
{
  new_v_pos[j] += offset[j] / (factor + 1);
  P(j, offset[j] / (factor + 1));
  //        new_tet_pos[0][j] -= offset[j] / (factor + 1);
  //        new_tet_pos[1][j] -= offset[j] / (factor + 1);
  //        new_tet_pos[2][j] -= offset[j] / (factor + 1);
  new_tet_pos[0][j] -= offset[j] * bary[0] / (factor + 1);
  new_tet_pos[1][j] -= offset[j] * bary[1] / (factor + 1);
  new_tet_pos[2][j] -= offset[j] * bary[2] / (factor + 1);
  //                new_tet_pos[3][j] -= offset[j] * bary[3] / (factor + 1);
  //        if (offset[j] > kMaxOffset) offset[j] = kMaxOffset;
  //        else if (offset[j] < -kMaxOffset) offset[j] = -kMaxOffset;
}
//      getchar();
++count_[v];
++count_[vert_idx[0]];
++count_[vert_idx[1]];
++count_[vert_idx[2]];

inline Float PointToTriangleDistance(const Float *p, Float *p0, Float *p1, Float *p2, Float* n, Float *pr, Float* bary)
{
  Float& b0 = bary[0];
  Float& b1 = bary[1];
  Float& b2 = bary[2];
  Float weight;
  Float e[3], e0[3], e1[3], e2[3], temp[3];
  (void) e1;
  for (int i = 0; i < 3; i++) {
    e [i] = p [i] - p0[i];
    e0[i] = p1[i] - p0[i];
    e1[i] = p2[i] - p1[i];
    e2[i] = p0[i] - p2[i];
  }

  dj::Cross3(e0, e1, n);
  Float area = dj::Normalize3(n);
  weight = -dj::Dot3(e, n);

  for (int i = 0; i < 3; i++)
    pr[i] = p[i] + weight * n[i];

  Float er0[3], er1[3], er2[3];
  for (int i = 0; i < 3; i++) {
    er0[i] = pr[i] - p0[i];
    er1[i] = pr[i] - p1[i];
    er2[i] = pr[i] - p2[i];
  }

  dj::Cross3(er1, er2, temp);
  b0 = -dj::Dot3(temp, n);
  dj::Cross3(er2, er0, temp);
  b1 = -dj::Dot3(temp, n);
  b2 = area - b0 - b1;

  if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
    Float inv_sum = 1 / (b0 + b1 + b2);
    b0 *= inv_sum;
    b1 *= inv_sum;
    b2 *= inv_sum;
  } else {
    Float *pp0 = NULL, *pp1 = NULL, e[3];
    Float *bb0 = NULL, *bb1 = NULL;
    if (b0 < 0) {
      b0 = 0;
      if (p1 < p2)	{pp0 = p1; pp1 = p2; bb0 = &b1; bb1 = &b2; }
      else		{pp0 = p2; pp1 = p1; bb0 = &b2; bb1 = &b1; }
    } else if (b1 < 0) {
      b1 = 0;
      if (p0 < p2)	{pp0 = p2; pp1 = p0; bb0 = &b2; bb1 = &b0; }
      else		{pp0 = p0; pp1 = p2; bb0 = &b0; bb1 = &b2; }
    } else if (b2 < 0) {
      b2 = 0;
      if (p0 < p1)	{pp0 = p0; pp1 = p1; bb0 = &b0; bb1 = &b1; }
      else		{pp0 = p1; pp1 = p0; bb0 = &b1; bb1 = &b0; }
    }
    ASSERT(pp0 != NULL, P(b0, b1, b2));
    ASSERT(pp1 != NULL, P(b0, b1, b2));
    ASSERT(bb0 != NULL, P(b0, b1, b2));
    ASSERT(bb1 != NULL, P(b0, b1, b2));
    e[0] = pp1[0] - pp0[0];
    e[1] = pp1[1] - pp0[1];
    e[2] = pp1[2] - pp0[2];

    Float e_length = dj::Normalize3(e);

    Float ep[3];
    ep[0] = p[0] - pp0[0];
    ep[1] = p[1] - pp0[1];
    ep[2] = p[2] - pp0[2];


    *bb1 = dj::Dot3(ep, e) / e_length;
    *bb0 = 1 - *bb1;

    if (*bb1 < 0)	{*bb1 = 0; *bb0 = 1;}
    if (*bb0 < 0)	{*bb1 = 1; *bb0 = 0;}
  }

  pr[0] = p0[0] * b0 + p1[0] * b1 + p2[0] * b2;
  pr[1] = p0[1] * b0 + p1[1] * b1 + p2[1] * b2;
  pr[2] = p0[2] * b0 + p1[2] * b1 + p2[2] * b2;

  return sqrt((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
}

//      if (verbose) {
//        P(tri_verts[0]);
//        P(tri_verts[1]);
//        P(tri_verts[2]);
//        P(dj::Vec3d(&buffer[edge_verts[0] * 3]));
//        P(dj::Vec3d(&buffer[edge_verts[1] * 3]));
//        P(dj::Vec3d(&buffer[tri_verts[0]]));
//        P(dj::Vec3d(&buffer[tri_verts[1]]));
//        P(dj::Vec3d(&buffer[tri_verts[2]]));
//        P(has_collision);
//        P(bary[0]);
//        P(bary[1]);
//        P(bary[2]);
//        P(t);
////        exit(0);
//      }

void TestSegmentTriangleIntersection()
{
  Float tri[3][3] = {
    0, 0, 0,
    1, 0, 0,
    0, 3, 0,
  };
  Float e[2][3] = {
    2.3, 2.3, 2.2,
    3.3, 3.8, -2.4,
  };
  Float u, v, t;
  bool col = collision_util::SegmentTriangleIntersection(e[0], e[1], tri[0], tri[1], tri[2], 0.001, u, v, t);
  dj::Vec3d t0(tri[0]) ;
  dj::Vec3d t1(tri[1]) ;
  dj::Vec3d t2(tri[2]) ;
  dj::Vec3d e0(e[0]);
  dj::Vec3d e1(e[1]);
  dj::Vec3d pe = e0 * (1 - t) + e1 * (t);
  dj::Vec3d pt = t0 * (1 - u - v) + t1 * u + t2 * v;
  P(col);
  P(pe, pt);
  P(u, v, t);
}

#endif
