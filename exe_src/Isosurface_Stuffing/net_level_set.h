#ifndef NET_LEVEL_SET_H
#define NET_LEVEL_SET_H
#include <Eigen/Dense>
#include <cmath>
#include "vectorObj.h"
#include "lib/vector_lib.h"
#pragma once

struct NetLevelSet {
  typedef Eigen::Vector3d Vec3;
  NetLevelSet(double width, double height, double radius,
              int nx, int ny,
              double x0, double y0, double dx, double dy) {
    width_ = width;
    height_ = height;
    radius_ = radius;
    nx_ = nx;
    ny_ = ny;
    x0_ = x0;
    y0_ = y0;
    dx_ = dx;
    dy_ = dy;
  }

  inline double Point2CylinderDistance(double h, double radius, VectorObj& p) const {
    double distance = 1e10;
    if (p[1] > h / 2 || p[1] < -h / 2) {
      p[1] = (p[1] > h / 2) ? dj::Abs(p[1] - h / 2) : dj::Abs(p[1] + h / 2);
      double r = std::sqrt(p[0] * p[0] + p[2] * p[2]);
      r -= radius;
      if (r <= 0) {
        distance = p[1];
      } else {
        distance = std::sqrt(r * r + p[1] * p[1]);
      }
    } else {
      double r = std::sqrt(p[0] * p[0] + p[2] * p[2]);
      r -= radius;
      distance = r;
    }
    return distance;
  }

  inline double Distance2XCylinder(VectorObj& p) const {
    double dist = 1e10;
    for (int i = 0; i < nx_; ++i) {
      VectorObj px(p[0] - (x0_ + i * dx_), p[1], p[2]);
      double d = Point2CylinderDistance(height_, radius_, px);
      dist = dj::Min(dist, d);
    }
    return dist;
  }

  inline double Distance2YCylinder(VectorObj& p) const {
    double dist = 1e10;
    for (int i = 0; i < ny_; ++i) {
      double y = y0_ + i * dy_;
      VectorObj px(p[1] - y, p[0], p[2]);
      double d = Point2CylinderDistance(width_, radius_, px);
      dist = dj::Min(dist, d);
    }
    return dist;
  }

  double operator()(const VectorObj& p) const {
//    return -std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) + width_;
    VectorObj tmp_p(p[0], p[1], p[2]);
//    return -Point2CylinderDistance(15, 10, tmp_p);
//    return Distance2XCylinder(tmp_p);
    return -dj::Min(Distance2XCylinder(tmp_p), Distance2YCylinder(tmp_p));
  }

  double radius_;
  double width_, height_;
  double x0_, y0_;
  double dx_, dy_;
  int nx_, ny_;
};

#endif // NET_LEVEL_SET_H

