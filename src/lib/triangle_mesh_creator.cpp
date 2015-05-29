#include <cmath>
#include "triangle_mesh_creator.h"
#ifndef M_PI
#define M_PI 3.1415926f
#endif

void DoubleVector2FloatVector(std::vector<double>& double_vector, std::vector<float>& float_vector) {
  float_vector.resize(double_vector.size());
  for (int i = 0; i < int(double_vector.size()); ++i) {
    float_vector[i] = double_vector[i];
  }
}

void CreateCylinder(double radius, double height, int slice, int stack,
                    std::vector<double> &verts, std::vector<double> &normal, std::vector<int> &tri) {
  //                    std::vector<double>* tex_coord) {
  typedef double Real;
  Real base = radius;
  Real top = radius;
  verts.resize(slice * stack * 3);
  normal.resize(slice * stack * 3);
  Real kAngle = 2 * M_PI / slice;
  for (int i = 0, idx = 0; i < stack; i++) {
    Real y = (i * height) / (stack - 1) ;
    Real radius = base + ((top - base) * i) / (stack - 1) ;
    for (int j = 0 ; j < slice; j++, idx += 3) {
      verts[idx + 0] = radius * std::cos(j * kAngle);
      verts[idx + 1] = y;
      verts[idx + 2] = radius * std::sin(j * kAngle);
      normal[idx + 0] = std::cos(j * kAngle);
      normal[idx + 1] = 0;
      normal[idx + 2] = std::sin(j * kAngle);
    }
  }

  tri.resize(3 * (slice * 2 * (stack - 1)));
  for (int y = 0, idx = 0; y < stack - 1 ; y++) {
    for (int x = 0 ; x < slice ; x++, idx += 6) {
      tri[idx + 0] = (y + 0) * slice + (x + 0);
      tri[idx + 1] = (y + 0) * slice + (x + 1) % slice;
      tri[idx + 2] = (y + 1) * slice + (x + 1) % slice;

      tri[idx + 3] = (y + 0) * slice + (x + 0);
      tri[idx + 4] = (y + 1) * slice + (x + 1) % slice;
      tri[idx + 5] = (y + 1) * slice + (x + 0) % slice;
    }
  }
}


void CreateCylinder(float radius, float height, int slice, int stack,
                    std::vector<float>& verts, std::vector<float>& normal, std::vector<int>& tri) {
  std::vector<double> verts_double, normal_double;
  CreateCylinder(radius, height, slice, stack, verts_double, normal_double, tri);
  DoubleVector2FloatVector(verts_double, verts);
  DoubleVector2FloatVector(normal_double, normal);
}


void ExpandTriangleMesh(std::vector<int> &tris, std::vector<double> &verts, std::vector<double> &normal,
                        std::vector<double> &expanded_verts, std::vector<double> &expanded_normal) {
  int triangle_num = int(tris.size()) / 3;
  expanded_verts.resize(triangle_num * 3 * 3);
  expanded_normal.resize(triangle_num * 3 * 3);
  for (int t = 0, idx = 0; t < triangle_num; ++t) {
    int* v = &tris[t * 3];
    for (int i = 0; i < 3; ++i, idx += 3) {
      expanded_verts[idx + 0] = verts[v[i] * 3 + 0];
      expanded_verts[idx + 1] = verts[v[i] * 3 + 1];
      expanded_verts[idx + 2] = verts[v[i] * 3 + 2];

      expanded_normal[idx + 0] = normal[v[i] * 3 + 0];
      expanded_normal[idx + 1] = normal[v[i] * 3 + 1];
      expanded_normal[idx + 2] = normal[v[i] * 3 + 2];
    }
  }
}



void CreateSphere(double radius, int slice, int stack, std::vector<double> &verts, std::vector<double> &normal, std::vector<int> &tri, std::vector<double> *tex_coord) {
  (void) tex_coord;
  typedef double Real;
  int down_stack = (stack - 1) / 2;
  int up_stack = stack - down_stack;
  int vertex_num = slice * stack + 2;
  verts.resize(vertex_num * 3);
  normal.resize(vertex_num * 3);
  verts[(slice * stack) * 3 + 0] = 0;
  verts[(slice * stack) * 3 + 1] = radius;
  verts[(slice * stack) * 3 + 2] = 0;
  normal[(slice * stack) * 3 + 0] = 0;
  normal[(slice * stack) * 3 + 1] = 1;
  normal[(slice * stack) * 3 + 2] = 0;

  verts[(slice * stack + 1) * 3 + 0] = 0;
  verts[(slice * stack + 1) * 3 + 1] = -radius;
  verts[(slice * stack + 1) * 3 + 2] = 0;
  normal[(slice * stack + 1) * 3 + 0] = 0;
  normal[(slice * stack + 1) * 3 + 1] = -1;
  normal[(slice * stack + 1) * 3 + 2] = 0;

  Real kAngle = 2 * M_PI / slice;
  for (int i = 0 ; i < up_stack; i++) {
    Real y = i * 1.0f / up_stack * radius ;
    Real r = sqrtf(radius * radius - y * y) ;
    for (int j = 0; j < slice; j++) {
      verts[(i * slice + j) * 3 + 0] = r * cosf(j * kAngle);
      verts[(i * slice + j) * 3 + 1] = y;
      verts[(i * slice + j) * 3 + 2] = r * sinf(j * kAngle);

      normal[(i * slice + j) * 3 + 0] = verts[(i * slice + j) * 3 + 0] / radius;
      normal[(i * slice + j) * 3 + 1] = verts[(i * slice + j) * 3 + 1] / radius;
      normal[(i * slice + j) * 3 + 2] = verts[(i * slice + j) * 3 + 2] / radius;
    }
  }

  for (int i = 0 ; i < down_stack; i++) {
    float y = -(i + 1) * 1.0f / (down_stack + 1) * radius ;
    float r = sqrtf(radius * radius - y * y) ;
    for (int j = 0 ; j < slice ; j++) {
      verts[(up_stack * slice + i * slice + j) * 3 + 0] = r * cosf(kAngle * j) ;
      verts[(up_stack * slice + i * slice + j) * 3 + 1] = y;
      verts[(up_stack * slice + i * slice + j) * 3 + 2] = r * sinf(kAngle * j) ;

      normal[(up_stack * slice + i * slice + j) * 3 + 0] = verts[(up_stack * slice + i * slice + j) * 3 + 0] / radius;
      normal[(up_stack * slice + i * slice + j) * 3 + 1] = verts[(up_stack * slice + i * slice + j) * 3 + 1] / radius;
      normal[(up_stack * slice + i * slice + j) * 3 + 2] = verts[(up_stack * slice + i * slice + j) * 3 + 2] / radius;
    }
  }

  int triangle_num = slice * (stack - 1) * 2 + slice * 2;
  tri.resize(triangle_num * 3);
  int tri_count = 0;
  for (int y = 0 ; y < up_stack - 1; y++) {
    for (int x = 0 ; x < slice ; x++, tri_count += 2) {
      tri[tri_count * 3 + 0] = (y + 0) * slice + (x + 0);
      tri[tri_count * 3 + 1] = (y + 0) * slice + (x + 1) % slice;
      tri[tri_count * 3 + 2] = (y + 1) * slice + (x + 1) % slice;

      tri[tri_count * 3 + 3] = (y + 0) * slice + (x + 0);
      tri[tri_count * 3 + 4] = (y + 1) * slice + (x + 1) % slice;
      tri[tri_count * 3 + 5] = (y + 1) * slice + (x + 0);
    }
  }
  for (int x = 0; x < slice ; x++, tri_count += 2) {
    tri[3 * tri_count + 0] = (x + 0);
    tri[3 * tri_count + 1] = up_stack * slice + (x + 1) % slice;
    tri[3 * tri_count + 2] = (x + 1) % slice;

    tri[3 * tri_count + 3] = (x + 0);
    tri[3 * tri_count + 4] = up_stack * slice + (x + 0);
    tri[3 * tri_count + 5] = up_stack * slice + (x + 1) % slice;
  }

  for (int y = 0 ; y < down_stack - 1 ; y++) {
    for (int x = 0 ; x < slice ; x++, tri_count += 2) {
      tri[3 * tri_count + 0] = up_stack * slice + (y + 0) * slice + (x + 0);
      tri[3 * tri_count + 1] = up_stack * slice + (y + 1) * slice + (x + 1) % slice;
      tri[3 * tri_count + 2] = up_stack * slice + (y + 0) * slice + (x + 1) % slice;

      tri[3 * tri_count + 3] = up_stack * slice + (y + 0) * slice + (x + 0);
      tri[3 * tri_count + 4] = up_stack * slice + (y + 1) * slice + (x + 0);
      tri[3 * tri_count + 5] = up_stack * slice + (y + 1) * slice + (x + 1) % slice;
    }
  }
  for (int i = 0 ; i < slice ; i++, tri_count += 2) {
    tri[tri_count * 3 + 0] = (up_stack - 1) * slice + i;
    tri[tri_count * 3 + 1] = (up_stack - 1) * slice + (i + 1) % slice;
    tri[tri_count * 3 + 2] = slice * stack ;

    tri[tri_count * 3 + 3] = (stack - 1) * slice + (i + 1) % slice ;
    tri[tri_count * 3 + 4] = (stack - 1) * slice + i;
    tri[tri_count * 3 + 5] = slice * stack + 1;
  }
}


void CreateSphere(float radius, int slice, int stack,
                  std::vector<float>& verts, std::vector<float>& normal, std::vector<int>& tri,
                  std::vector<float>* tex_coord) {
  std::vector<double> verts_double, normal_double, tex_coord_double;
  if (tex_coord == nullptr) {
    CreateSphere(radius, slice, stack, verts_double, normal_double, tri, nullptr);
  } else {
    CreateSphere(radius, slice, stack, verts_double, normal_double, tri, &tex_coord_double);
    DoubleVector2FloatVector(tex_coord_double, *tex_coord);
  }
  DoubleVector2FloatVector(verts_double, verts);
  DoubleVector2FloatVector(normal_double, normal);
}
