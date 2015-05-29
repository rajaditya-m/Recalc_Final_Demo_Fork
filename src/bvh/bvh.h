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
#ifndef BVH_H
#define BVH_H
#include <functional>
#include <float.h>
#include <stdlib.h>
#include "tet_collider.h"

#ifndef MAX
#define MAX(a,b)	((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b)	((a) < (b) ? (a) : (b))
#endif

#ifndef FORCEINLINE
#  ifndef _WIN32
#    define FORCEINLINE inline
#  else
#    define FORCEINLINE __forceinline
#endif
#endif


typedef TetCollider::Float Float;
typedef Float* vec3f;

class kDOP18
{
public:
  FORCEINLINE static void getDistances(const vec3f p,
                                       Float &d3, Float &d4, Float &d5, Float &d6, Float &d7, Float &d8)
  {
    d3 = p[0] + p[1];
    d4 = p[0] + p[2];
    d5 = p[1] + p[2];
    d6 = p[0] - p[1];
    d7 = p[0] - p[2];
    d8 = p[1] - p[2];
  }

  FORCEINLINE static void getDistances(const vec3f p, Float d[])
  {
    d[0] = p[0] + p[1];
    d[1] = p[0] + p[2];
    d[2] = p[1] + p[2];
    d[3] = p[0] - p[1];
    d[4] = p[0] - p[2];
    d[5] = p[1] - p[2];
  }

  FORCEINLINE static Float getDistances(const vec3f p, int i)
  {
    if (i == 0) return p[0] + p[1];
    if (i == 1) return p[0] + p[2];
    if (i == 2) return p[1] + p[2];
    if (i == 3) return p[0] - p[1];
    if (i == 4) return p[0] - p[2];
    if (i == 5) return p[1] - p[2];
    return 0;
  }

public:
  Float _dist[18];

  FORCEINLINE kDOP18()
  {
    empty();
  }

  FORCEINLINE kDOP18(const vec3f v)
  {
    _dist[0] = _dist[9]  = v[0];
    _dist[1] = _dist[10] = v[1];
    _dist[2] = _dist[11] = v[2];

    Float d3, d4, d5, d6, d7, d8;
    getDistances(v, d3, d4, d5, d6, d7, d8);
    _dist[3] = _dist[12] = d3;
    _dist[4] = _dist[13] = d4;
    _dist[5] = _dist[14] = d5;
    _dist[6] = _dist[15] = d6;
    _dist[7] = _dist[16] = d7;
    _dist[8] = _dist[17] = d8;
  }

  FORCEINLINE kDOP18(const vec3f a, const vec3f b)
  {
    _dist[0]  = MIN(a[0], b[0]);
    _dist[9]  = MAX(a[0], b[0]);
    _dist[1]  = MIN(a[1], b[1]);
    _dist[10] = MAX(a[1], b[1]);
    _dist[2]  = MIN(a[2], b[2]);
    _dist[11] = MAX(a[2], b[2]);

    Float ad3, ad4, ad5, ad6, ad7, ad8;
    getDistances(a, ad3, ad4, ad5, ad6, ad7, ad8);
    Float bd3, bd4, bd5, bd6, bd7, bd8;
    getDistances(b, bd3, bd4, bd5, bd6, bd7, bd8);
    _dist[3]  = MIN(ad3, bd3);
    _dist[12] = MAX(ad3, bd3);
    _dist[4]  = MIN(ad4, bd4);
    _dist[13] = MAX(ad4, bd4);
    _dist[5]  = MIN(ad5, bd5);
    _dist[14] = MAX(ad5, bd5);
    _dist[6]  = MIN(ad6, bd6);
    _dist[15] = MAX(ad6, bd6);
    _dist[7]  = MIN(ad7, bd7);
    _dist[16] = MAX(ad7, bd7);
    _dist[8]  = MIN(ad8, bd8);
    _dist[17] = MAX(ad8, bd8);
  }

  FORCEINLINE bool overlaps(const kDOP18& b) const
  {
    for (int i = 0; i < 9; i++) {
      if (_dist[i] > b._dist[i + 9]) return false;
      if (_dist[i + 9] < b._dist[i]) return false;
    }

    return true;
  }

  FORCEINLINE bool overlaps(const kDOP18 &b, kDOP18 &ret) const
  {
    if (!overlaps(b))
      return false;

    for (int i = 0; i < 9; i++) {
      ret._dist[i] = MAX(_dist[i],  b._dist[i]);
      ret._dist[i + 9] = MIN(_dist[i + 9], b._dist[i + 9]);
    }
    return true;
  }

  FORCEINLINE bool inside(const vec3f p) const
  {
    for (int i = 0; i < 3; i++) {
      if (p[i] < _dist[i] || p[i] > _dist[i + 9])
        return false;
    }

    Float d[6];
    getDistances(p, d);
    for (int i = 3; i < 9; i++) {
      if (d[i - 3] < _dist[i] || d[i - 3] > _dist[i + 9])
        return false;
    }

    return true;
  }

  //  FORCEINLINE kDOP18 &operator += (const vec3f p)
  FORCEINLINE kDOP18& include(const vec3f p)
  {
    _dist[0]  = MIN(p[0], _dist[0]);
    _dist[9]  = MAX(p[0], _dist[9]);
    _dist[1]  = MIN(p[1], _dist[1]);
    _dist[10] = MAX(p[1], _dist[10]);
    _dist[2]  = MIN(p[2], _dist[2]);
    _dist[11] = MAX(p[2], _dist[11]);

    Float d3, d4, d5, d6, d7, d8;
    getDistances(p, d3, d4, d5, d6, d7, d8);
    _dist[3]  = MIN(d3, _dist[3]);
    _dist[12] = MAX(d3, _dist[12]);
    _dist[4]  = MIN(d4, _dist[4]);
    _dist[13] = MAX(d4, _dist[13]);
    _dist[5]  = MIN(d5, _dist[5]);
    _dist[14] = MAX(d5, _dist[14]);
    _dist[6]  = MIN(d6, _dist[6]);
    _dist[15] = MAX(d6, _dist[15]);
    _dist[7]  = MIN(d7, _dist[7]);
    _dist[16] = MAX(d7, _dist[16]);
    _dist[8]  = MIN(d8, _dist[8]);
    _dist[17] = MAX(d8, _dist[17]);

    return *this;
  }

  FORCEINLINE kDOP18 &operator += (const kDOP18 &b)
  {
    _dist[0]  = MIN(b._dist[0], _dist[0]);
    _dist[9]  = MAX(b._dist[9], _dist[9]);
    _dist[1]  = MIN(b._dist[1], _dist[1]);
    _dist[10] = MAX(b._dist[10], _dist[10]);
    _dist[2]  = MIN(b._dist[2], _dist[2]);
    _dist[11] = MAX(b._dist[11], _dist[11]);
    _dist[3]  = MIN(b._dist[3], _dist[3]);
    _dist[12] = MAX(b._dist[12], _dist[12]);
    _dist[4]  = MIN(b._dist[4], _dist[4]);
    _dist[13] = MAX(b._dist[13], _dist[13]);
    _dist[5]  = MIN(b._dist[5], _dist[5]);
    _dist[14] = MAX(b._dist[14], _dist[14]);
    _dist[6]  = MIN(b._dist[6], _dist[6]);
    _dist[15] = MAX(b._dist[15], _dist[15]);
    _dist[7]  = MIN(b._dist[7], _dist[7]);
    _dist[16] = MAX(b._dist[16], _dist[16]);
    _dist[8]  = MIN(b._dist[8], _dist[8]);
    _dist[17] = MAX(b._dist[17], _dist[17]);
    return *this;
  }

  FORCEINLINE kDOP18 operator + ( const kDOP18 &v) const
  { kDOP18 rt(*this); return rt += v; }

  FORCEINLINE Float length(int i) const
  {
    return _dist[i + 9] - _dist[i];
  }

  FORCEINLINE Float width()  const { return _dist[9] - _dist[0]; }
  FORCEINLINE Float height() const { return _dist[10] - _dist[1]; }
  FORCEINLINE Float depth()  const { return _dist[11] - _dist[2]; }
  FORCEINLINE Float volume() const { return width() * height() * depth(); }

  FORCEINLINE void center(vec3f c) const
  {
    c[0] = (_dist[0] + _dist[9]) * 0.5;
    c[1] = (_dist[1] + _dist[10]) * 0.5;
    c[2] = (_dist[2] + _dist[11]) * 0.5;
  }

  FORCEINLINE Float center(int i) const
  {
    return (_dist[i + 9] + _dist[i]) * 0.5;
  }

  FORCEINLINE void empty()
  {
    for (int i = 0; i < 9; i++) {
      _dist[i] = FLT_MAX;
      _dist[i + 9] = -FLT_MAX;
    }
  }
};

#define BOX kDOP18


bool overlap (const BOX &box0, const BOX &box1, Float thickness);


class DeformBVHNode;
typedef TetCollider DeformModel;

class DeformBVHNode
{
public:
  DeformModel* _mesh;
  BOX _box;
  int _tet;
  DeformBVHNode *_parent;
  DeformBVHNode *_left;
  DeformBVHNode *_right;
  bool _active;

public:
  DeformBVHNode();
  DeformBVHNode(DeformModel *mesh, DeformBVHNode *, int, BOX *, Float (*)[3]);
  DeformBVHNode(DeformModel *mesh, DeformBVHNode *, int*, unsigned int, BOX *, Float (*)[3]);

  ~DeformBVHNode();

  void refit(bool = false);
  bool find(int);
  DeformBVHNode* findTetBox(int tet);
  FORCEINLINE DeformBVHNode *getLeftChild() { return _left; }
  FORCEINLINE DeformBVHNode *getRightChild() { return _right; }
  FORCEINLINE DeformBVHNode *getParent() { return _parent; }
  FORCEINLINE BOX vertBox(int vert, bool ccd)
  {
    (void) ccd;
    BOX box;
    //    box.include(_mesh->vertex_ + vert * 3);
    box.include(&_mesh->vert_[vert][0]);
    //    if (ccd)
    //      box.include(_mesh->prev_vertex_ + vert * 3);
    return box;
  }

  BOX edgeBox(Float* v0, Float* v1)
  {
    BOX box;
    box.include(v0);
    box.include(v1);
    return box;
  }

//  BOX edgeBox(int edge, bool ccd)
//  {
//    BOX box;
//    //    int* verts = _mesh->edge_ + edge * 2;
//    int* verts = &_mesh->edge_[edge][0];
//    box += vertBox(verts[0], ccd);
//    box += vertBox(verts[1], ccd);
//    return box;
//  }

  FORCEINLINE BOX tetBox(int tet, bool ccd)
  {
    BOX box;
    int* verts = &_mesh->tet_[tet][0];
    for (int v = 0; v < 4; v++)
      box += vertBox(verts[v], ccd);
    return box;
  }

  FORCEINLINE int getTet() { return _tet; }
  FORCEINLINE bool isLeaf() { return _left == NULL; }
  FORCEINLINE bool isRoot() { return _parent == NULL;}

  friend class DeformBVHTree;
};

class DeformBVHTree
{
public:
  DeformModel		*_mdl;
  DeformBVHNode	*_root;
  int* tet_buffer;

  bool _ccd;

public:
  DeformBVHTree(DeformModel &, bool);
  ~DeformBVHTree();

  void Construct();

  Float refit();

  BOX box();
  FORCEINLINE DeformBVHNode *getRoot() { return _root; }

  friend class DeformBVHNode;
};

//typedef void (*CallBack)(int);
typedef std::function<void (int)> RBCallBack;
void ForOverlappingFaces(DeformBVHNode* node, BOX& box, RBCallBack call);

#endif // BVH_H
