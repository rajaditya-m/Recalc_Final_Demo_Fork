#pragma once
#ifndef HAMMOCK_LEVEL_SET
#define HAMMOCK_LEVEL_SET
#include <Eigen/Dense>
#include <cmath>
#include "vector_lib.h"

class HammockLevelSet {
public:
  typedef Eigen::Vector3d Vec3;
  HammockLevelSet(double width, double height, double radius,
                  int nx, int ny) {
    width_ = width;
    height_ = height;
    radius_ = radius;
    nx_ = nx;
    ny_ = ny;
    dx_ = width_ / nx_;
    dy_ = height_ / ny_;
    cylinder_length_ = std::sqrt(dx_ * dx_ + dy_ * dy_) / 2.0;
    cosine_ = dx_ / (2 * cylinder_length_);
    sine_ = dy_ / (2 * cylinder_length_);
    min_x_ = -width_ / 2.0;
    min_y_ = -height_ / 2.0;
  }

  inline Vec3 Rotate(double cosine, double sine, const Vec3& p) const {
    //    double c = cosine, s = sine;
    //    double theta = atan2(sine, cosine) + 3.1415926 / 4;
    //    cosine = cos(theta);
    //    sine = sin(theta);
    return Vec3(cosine * p[0] - sine   * p[1],
                sine   * p[0] + cosine * p[1],
                p[2]);
  }

  inline double Point2CylinderDistance(double h, double radius, Vec3& p) const {
    double distance = 1e10;
    if (p[1] > h / 2 || p[1] < -h / 2) {
      p[1] = (p[1] > h / 2) ? dj::Abs(p[1] - h / 2) : dj::Abs(p[1] + h / 2);
      double r = std::sqrt(p[0] * p[0] + p[2] * p[2] + p[1] * p[1]) ;
      r -= radius;
      distance = r;
      //      if (r <= 0) {
      //        distance = p[1];
      //      } else {
      //        distance = std::sqrt(r * r + p[1] * p[1]);
      //      }
    } else {
      double r = std::sqrt(p[0] * p[0] + p[2] * p[2]);
      r -= radius;
      distance = r;
    }
    return distance;
  }

  inline double operator()(const double* p) const {
    const double kAngle = 3.1415916 / 4;
    double distance = 1e10;
    const double kScale = 0.6;
    const double X0 = -0.0;
    const double Y0 = +0.0;
    for (int y = 0; y < ny_; ++y) {
      for (int x = 0; x < nx_; ++x) {
        double origin[4][2] = {
          {X0 + min_x_ + x * dx_ - dx_ / 4, Y0 + min_y_ + y * dy_ - dy_ / 4},
          {X0 + min_x_ + x * dx_ + dx_ / 4, Y0 + min_y_ + y * dy_ - dy_ / 4},
          {X0 + min_x_ + x * dx_ + dx_ / 4, Y0 + min_y_ + y * dy_ + dy_ / 4},
          {X0 + min_x_ + x * dx_ - dx_ / 4, Y0 + min_y_ + y * dy_ + dy_ / 4},
        };
        //        if (x == 0) {
        //          origin[0][0] += kScale * 0.5 * (dx_ / 4);
        //          origin[0][1] -= kScale * 0.5 * (dy_ / 4);
        //
        //          origin[3][0] += kScale * 0.5 * (dx_ / 4);
        //          origin[3][1] += kScale * 0.5 * (dy_ / 4);
        //        }
        //        if (y == 0) {// && x > 0 && x < nx_ - 1) {
        //          origin[0][0] -= kScale * 0.5 * (dx_ / 4);
        //          origin[0][1] += kScale * 0.5 * (dy_ / 4);
        //
        //          origin[1][0] += kScale * 0.5 * (dx_ / 4);
        //          origin[1][1] += kScale * 0.5 * (dy_ / 4);
        //        }
        //        if (x == nx_ - 1) {
        //          origin[1][0] -= kScale * 0.5 * (dx_ / 4);
        //          origin[1][1] -= kScale * 0.5 * (dy_ / 4);
        //
        //          origin[2][0] -= kScale * 0.5 * (dx_ / 4);
        //          origin[2][1] += kScale * 0.5 * (dy_ / 4);
        //        }
        double theta[2][2] = {
          { -sine_, cosine_},
          { +sine_, cosine_},
        };
        for (int i = 0; i < 4; ++i) {
          Vec3 tmp_p(p);
          tmp_p = Rotate(cos(kAngle), sin(kAngle), tmp_p);
          Vec3 pos(tmp_p[0] - origin[i][0], tmp_p[1] - origin[i][1], tmp_p[2]);
          Vec3 tmp_pos = Rotate(theta[i % 2][0], theta[i % 2][1], pos);
          //          if (y == 0 && x == 0 && i == 0) continue;
          //          if (y == 0 && x == nx_ - 1 && i == 1) continue;
          //          if (x == 0 || y == 0 || x == nx_ - 1) {
          //            distance = dj::Min(Point2CylinderDistance(kScale * cylinder_length_, radius_, tmp_pos), distance);
          //          } else {
          distance = dj::Min(Point2CylinderDistance(cylinder_length_, radius_, tmp_pos), distance);
          //          }
        }
      }
    }
    // two bars
    double kBarWidth = width_ * 1.2;
    {
      Vec3 tmp_p(p);
      tmp_p = Rotate(cos(kAngle), sin(kAngle), tmp_p);
      Vec3 pos(tmp_p[1] - height_ / 2 + dy_ / 2.5, tmp_p[0] + dx_ / 2, tmp_p[2]);
      distance = dj::Min(Point2CylinderDistance(kBarWidth, radius_, pos), distance);
    }

    if (1) {
      Vec3 tmp_p(p);
      tmp_p = Rotate(cos(kAngle), sin(kAngle), tmp_p);
      //      Vec3 pos(p[1] + height_ / 2 + dy_ / 1.8, p[0] + dx_ / 2, p[2]);
      Vec3 pos(tmp_p[1] + height_ / 2 + dy_ / 1.8, tmp_p[0] + dx_ / 2, tmp_p[2]);
      distance = dj::Min(Point2CylinderDistance(kBarWidth, radius_, pos), distance);
    }
    return distance;
  }

  double cosine_, sine_;
  double width_, height_, radius_;
  double min_x_, min_y_;
  double dx_, dy_;
  int nx_, ny_;
  double cylinder_length_;
};

#endif // HAMMOCK_LEVEL_SET

