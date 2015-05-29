// Originally from the SELF-CCD library by Min Tang and Dinesh Manocha
// (http://gamma.cs.unc.edu/SELFCD/).
// Modified by Rahul Narain.

/*************************************************************************\

  Copyright 2010 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
   fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <assert.h>

#include "bvh.h"
#include <climits>
#include <utility>
#include "print_macro.h"
using namespace std;

/*
BOX node_box (const Node *node, bool ccd) {
    BOX box;
    box += node->x;
    if (ccd)
        box += node->x0;
    return box;
}

BOX vert_box (const Vert *vert, bool ccd) {
    return node_box(vert->node, ccd);
}

BOX edge_box (const Edge *edge, bool ccd) {
    BOX box;
    box += node_box(edge->n[0], ccd);
    box += node_box(edge->n[1], ccd);
    return box;
}

BOX face_box (const Face *face, bool ccd) {
    BOX box;
    for (int v = 0; v < 3; v++)
        box += vert_box(face->v[v], ccd);
    return box;
}
*/

BOX dilate (const BOX &box, Float d) {
  static Float sqrt2 = sqrt(2);
  BOX dbox = box;
  for (int i = 0; i < 3; i++) {
    dbox._dist[i] -= d;
    dbox._dist[i + 9] += d;
  }
  for (int i = 0; i < 6; i++) {
    dbox._dist[3 + i] -= sqrt2 * d;
    dbox._dist[3 + i + 9] += sqrt2 * d;
  }
  return dbox;
}

bool overlap (const BOX &box0, const BOX &box1, Float thickness) {
  return box0.overlaps(dilate(box1, thickness));
}

Float
DeformBVHTree::refit() {
  getRoot()->refit(_ccd);
  return 0;
}

BOX
DeformBVHTree::box() {
  return getRoot()->_box;
}


void
DeformBVHNode::refit(bool ccd) {
  if (isLeaf()) {
    _box = tetBox(getTet(), ccd);
  } else {
    getLeftChild()->refit(ccd);
    getRightChild()->refit(ccd);

    _box = getLeftChild()->_box + getRightChild()->_box;
  }
}

bool
DeformBVHNode::find(int face) {
  if (isLeaf())
    return getTet() == face;

  if (getLeftChild()->find(face))
    return true;

  if (getRightChild()->find(face))
    return true;

  return false;
}

DeformBVHNode *DeformBVHNode::findTetBox(int tet)
{
  if (isLeaf()) {
    if (getTet() == tet) return this;
    else return NULL;
  }

  DeformBVHNode* node = getLeftChild()->findTetBox(tet);
  if (node != NULL) return node;
  node = getRightChild()->findTetBox(tet);
  return node;
}

inline Float middle_xyz(char xyz, const vec3f p1, const vec3f p2, const vec3f &p3) {
  Float t0, t1;
  t0 = MIN(p1[(unsigned) xyz], p2[(unsigned)xyz]);
  t0 = MIN(t0, p3[(unsigned)xyz]);
  t1 = MAX(p1[(unsigned)xyz], p2[(unsigned)xyz]);
  t1 = MAX(t1, p3[(unsigned)xyz]);
  return (t0 + t1) * 0.5f;
}

inline Float middle_xyz(char xyz, const vec3f p1, const vec3f p2, const vec3f &p3, const vec3f &p4) {
  Float t0, t1;
  t0 = MIN(p1[(unsigned) xyz], p2[(unsigned)xyz]);
  t0 = MIN(t0, p3[(unsigned)xyz]);
  t0 = MIN(t0, p4[(unsigned)xyz]);
  t1 = MAX(p1[(unsigned)xyz], p2[(unsigned)xyz]);
  t1 = MAX(t1, p3[(unsigned)xyz]);
  t1 = MAX(t1, p4[(unsigned)xyz]);
  return (t0 + t1) * 0.5f;
}

class aap {
public:
  char _xyz;
  Float _p;

  FORCEINLINE aap(const BOX &total) {
    Float center[3];
    total.center(center);
    char xyz = 2;

    if (total.width() >= total.height() && total.width() >= total.depth()) {
      xyz = 0;
    } else if (total.height() >= total.width() && total.height() >= total.depth()) {
      xyz = 1;
    }

    _xyz = xyz;
    _p = center[(unsigned)xyz];
  }

  FORCEINLINE bool inside(const vec3f &mid) const {
    return mid[(unsigned)_xyz] > _p;
  }
};

extern Float middle_xyz(char xyz, const vec3f p1, const vec3f p2, const vec3f &p3);

DeformBVHTree::DeformBVHTree(DeformModel &mdl, bool ccd) {
  _mdl = &mdl;
  _ccd = ccd;
  if (_mdl->v_num_ > 0)
    Construct();
  else
    _root = NULL;
}

void
DeformBVHTree::Construct() {
  BOX total;
  unsigned int count;

  int num_vtx = _mdl->v_num_,
      num_tet = _mdl->tet_num_;

  for (int i = 0; i < num_vtx; i++) {
    total.include((vec3f) &_mdl->vert_[i][0]);
//    if (_ccd)
//      total.include(_mdl->prev_vertex_ + i * 3);
  }

  count = num_tet;

  BOX *tet_boxes = new BOX[count];
  std::vector<Float> tmp_tet_centers(count * 3);
  Float (*tet_centers)[3] = (Float (*)[3]) &tmp_tet_centers[0];
  aap  pln(total);

  tet_buffer = new int[count];
  unsigned int left_idx = 0, right_idx = count;
  unsigned int tet_idx = 0;

  for (int i = 0; i < num_tet; i++) {
    tet_idx++;

//    int* verts = _mdl->triangle_ + i * 3;
    int* verts = &_mdl->tet_[i][0];
    vec3f p1 = &_mdl->vert_[verts[0]][0];
    vec3f p2 = &_mdl->vert_[verts[1]][0];
    vec3f p3 = &_mdl->vert_[verts[2]][0];
    vec3f p4 = &_mdl->vert_[verts[3]][0];
//    vec3f pp1 = _mdl->prev_vertex_ + verts[0] * 3;
//    vec3f pp2 = _mdl->prev_vertex_ + verts[1] * 3;
//    vec3f pp3 = _mdl->prev_vertex_ + verts[2] * 3;

    if (_ccd) {
//      tri_centers[tri_idx - 1][0] = (middle_xyz(0, p1, p2, p3) + middle_xyz(0, pp1, pp2, pp3)) * 0.5f;
//      tri_centers[tri_idx - 1][1] = (middle_xyz(1, p1, p2, p3) + middle_xyz(1, pp1, pp2, pp3)) * 0.5f;
//      tri_centers[tri_idx - 1][2] = (middle_xyz(2, p1, p2, p3) + middle_xyz(2, pp1, pp2, pp3)) * 0.5f;
    } else {
      tet_centers[tet_idx - 1][0] = middle_xyz(0, p1, p2, p3, p4);
      tet_centers[tet_idx - 1][1] = middle_xyz(1, p1, p2, p3, p4);
      tet_centers[tet_idx - 1][2] = middle_xyz(2, p1, p2, p3, p4);
    }

    if (pln.inside(tet_centers[tet_idx - 1]))
      tet_buffer[left_idx++] = i;
    else
      tet_buffer[--right_idx] = i;;

    tet_boxes[tet_idx - 1].include(p1);
    tet_boxes[tet_idx - 1].include(p2);
    tet_boxes[tet_idx - 1].include(p3);
    tet_boxes[tet_idx - 1].include(p4);

//    if (_ccd) {
//      tri_boxes[tri_idx - 1].include(pp1);
//      tri_boxes[tri_idx - 1].include(pp2);
//      tri_boxes[tri_idx - 1].include(pp3);
//    }
  }

  _root = new DeformBVHNode();
  _root->_mesh = _mdl;
  _root->_box = total;
  //_root->_count = count;

  if (count == 1) {
    _root->_tet = 0;
    _root->_left = _root->_right = NULL;
  } else {
    if (left_idx == 0 || left_idx == count)
      left_idx = count / 2;
    _root->_left = new DeformBVHNode(_mdl, _root, tet_buffer, left_idx, tet_boxes, tet_centers);
    _root->_right = new DeformBVHNode(_mdl, _root, tet_buffer + left_idx, count - left_idx, tet_boxes, tet_centers);
  }

  delete [] tet_boxes;
}

DeformBVHTree::~DeformBVHTree() {
  if (!_root)
    return;
  delete _root;
  delete [] tet_buffer;
}

//#################################################################
// called by root
DeformBVHNode::DeformBVHNode() {
  _mesh = NULL;
  _tet = 0;
  _left = _right = NULL;
  _parent = NULL;
  //_count = 0;
  _active = true;
}

DeformBVHNode::~DeformBVHNode() {
  if (_left) delete _left;
  if (_right) delete _right;
}

// called by leaf
DeformBVHNode::DeformBVHNode(DeformModel* mesh, DeformBVHNode *parent, int tet, BOX *tet_boxes, Float (*tet_centers)[3]) {
  (void) tet_centers;
  _mesh = mesh;
  _left = _right = NULL;
  _parent = parent;
  _tet = tet;
  _box = tet_boxes[tet];
  //_count = 1;
  _active = true;
}

// called by nodes
DeformBVHNode::DeformBVHNode(DeformModel* mesh, DeformBVHNode *parent, int* lst, unsigned int lst_num, BOX *tri_boxes, Float (*tri_centers)[3]) {
  assert(lst_num > 0);
  _mesh = mesh;
  _left = _right = NULL;
  _parent = parent;
  _tet = 0;
  //_count = lst_num;
  _active = true;

  if (lst_num == 1) {
    _tet = lst[0];
    _box = tri_boxes[lst[0]];
  } else { // try to split them
    for (unsigned int t = 0; t < lst_num; t++) {
      int i = lst[t];
      _box += tri_boxes[i];
    }

    if (lst_num == 2) { // must split it!
      _left = new DeformBVHNode(_mesh, this, lst[0], tri_boxes, tri_centers);
      _right = new DeformBVHNode(_mesh, this, lst[1], tri_boxes, tri_centers);
    } else {
      aap pln(_box);
      unsigned int left_idx = 0, right_idx = lst_num - 1;

      for (unsigned int t = 0; t < lst_num; t++) {
        int i = lst[left_idx];
        if (pln.inside(tri_centers[i]))
          left_idx++;
        else {// swap it
          int tmp = lst[left_idx];
          lst[left_idx] = lst[right_idx];
          lst[right_idx--] = tmp;
        }
      }

      int hal = lst_num / 2;
      if (left_idx == 0 || left_idx == lst_num) {
        _left = new DeformBVHNode(_mesh, this, lst, hal, tri_boxes, tri_centers);
        _right = new DeformBVHNode(_mesh, this, lst + hal, lst_num - hal, tri_boxes, tri_centers);

      } else {
        _left = new DeformBVHNode(_mesh, this, lst, left_idx, tri_boxes, tri_centers);
        _right = new DeformBVHNode(_mesh, this, lst + left_idx, lst_num - left_idx, tri_boxes, tri_centers);
      }

    }
  }
}


void ForOverlappingFaces(DeformBVHNode *node, BOX& box, RBCallBack call) {
  if (!node || !overlap(node->_box, box, 1e-8f)) return;
  if (node->isLeaf()) {
    call(node->getTet());
  } else {
    ForOverlappingFaces(node->getLeftChild(), box, call);
    if (node->getRightChild() != NULL) {
      ForOverlappingFaces(node->getRightChild(), box, call);
    }
  }
}
