#include "fem_volumetric_object.h"
#include "conjugate_gradient_solver.h"
#include "opengl_helper.h"
#include "affine_transformer.h"
#include "vector_io.h"
#include "config_file.h"

#define FOR(i, start, end) for (int i = start; i < end; ++i)

FemVolumetricObject::FemVolumetricObject(TetrahedralMeshIO * mesh_io,
                                         const char * file_name,
                                         AffineTransformer<Real>* transformer)
  : Super(mesh_io, file_name, transformer) {
  Construct();
}

FemVolumetricObject::FemVolumetricObject(TetrahedralMesh *mesh)
  : Super(mesh) {
  Construct();
}

void FemVolumetricObject::Construct() {
  cg_solver_ = new ConjugateGradientSolver<Real>(v_num_ * 3);
  tet_vol_.resize(tet_num_);
  inv_dm_.resize(tet_num_);
  const Real kYoungModulus = 1e5;
  const Real kPoissonRatio = 0.4;
  const Real kLameMu = kYoungModulus / (2 * (1 + kPoissonRatio));
  const Real kLameLambda = kYoungModulus * kPoissonRatio / ((1 + kPoissonRatio) * (1 - 2 * kPoissonRatio));
  mu_ = std::vector<Real>(tet_num_, kLameMu);
  lambda_ = std::vector<Real>(tet_num_, kLameLambda);
  tet_density_ = std::vector<Real>(tet_num_, Real(1000));

  force_.resize(v_num_);
  vel_.resize(v_num_);
  prev_vert_.resize(v_num_);
  tmp_vert_buf_.resize(v_num_);
  prev_vert_ = vert_;
  ComputeInvDm();
  ComputeTetrahedraVolume();
  ComputeLumpedMass();
  //  mass_ = std::vector<Real>(v_num_, real(0.01));
}

void FemVolumetricObject::Simulate(Real dt) {
  const double kFixedVertexThreshold = 0.10;
  Mat3 tmp;
  const Real t2 = dt * dt;
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    for (int t = 0; t < tet_num_; ++t) {
      Real* pos[4] = {
        result + tet_[t][0] * 3,
        result + tet_[t][1] * 3,
        result + tet_[t][2] * 3,
        result + tet_[t][3] * 3,
      };
      ComputeTetDifferentialForce(t, (Real3) x, tmp);
      tmp[0] *= -t2;
      tmp[1] *= -t2;
      tmp[2] *= -t2;
      Real sum[3] = {0, 0, 0};
      for (int i = 0; i < 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
          pos[j][i] += tmp[i][j - 1] / mass_[tet_[t][j]];
          sum[i] += tmp[i][j - 1];
        }
        pos[0][i] -= sum[i] / mass_[tet_[t][0]];
      }
    }
    for (int v = 0; v < v_num_; ++v) {
      if (vert_[v][1] < kFixedVertexThreshold) {
        result[v * 3 + 1] = 0;
      }
    }
  };

  memset(&tmp_vert_buf_[0][0], 0, sizeof(Vec3) * v_num_);
  //  for (int v = 0; v < v_num_; ++v) {
  //    if (vert_[v][1] < kFixedVertexThreshold) {
  //      tmp_vert_buf_[v][1] -= -vert_[v][1] / dt;
  //    }
  //  }
  std::vector<Vec3> tmp_buf(v_num_);
  StiffnessMatrix(&tmp_vert_buf_[0][0], &tmp_buf[0][0]);
  // Compute RHS
  for (int t = 0; t < tet_num_; ++t) {
    ComputeTetForce(t, tmp);
    tmp[0] *= dt;
    tmp[1] *= dt;
    tmp[2] *= dt;
    Real sum[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
      for (int j = 1; j <= 3; ++j) {
        Real force = tmp[i][j - 1];
        vel_[tet_[t][j]][i] += force / mass_[tet_[t][j]];
        sum[i] += force;
      }
      vel_[tet_[t][0]][i] -= sum[i] / mass_[tet_[t][0]];
    }
  }
  for (int v = 0; v < v_num_; ++v) {
    vel_[v][0] += global::gravity[0] * dt + tmp_buf[v][0];
    vel_[v][1] += global::gravity[1] * dt + tmp_buf[v][1] ;
    vel_[v][2] += global::gravity[2] * dt + tmp_buf[v][2];
    if (vert_[v][1] < kFixedVertexThreshold) {
      vel_[v][1] = 0;
    }
  }
  {
    int pulled_vert = conf.Get<int>("selected_vertex");
    double force [3] = {
     300,
      0,
      0,
    };
    vel_[pulled_vert][0] += force[0] / mass_[pulled_vert] * dt;
    vel_[pulled_vert][1] += force[1] / mass_[pulled_vert] * dt;
    vel_[pulled_vert][2] += force[2] / mass_[pulled_vert] * dt;
  }

  //  {
  //    int v = 13019;
  //    vel_[v][1] -= 140000 * dt / mass_[v];
  //  }
  //    PrintImplicitMatrix<Real>(std::cout, StiffnessMatrix, v_num_ * 3);
  //    exit(0);
  //  FOR(v, 0, v_num_) {
  //    Vec3 vec(&vel_[v][0]);
  //    P(v, vec);
  //  }
  tmp_vert_buf_ = vel_;
  auto info = cg_solver_->Solve(&tmp_vert_buf_[0][0], &vel_[0][0], StiffnessMatrix, 2000, Real(1e-6));
  P(info.first, info.second);

  for (int v = 0; v < v_num_; ++v) {
    if (vert_[v][1] < kFixedVertexThreshold) {
      //      vel_[v][1] = -vert_[v][1] / dt;
      //      vert_[v] += vel_[v] * dt;
      //      vert_[v][1] = 0;
    } else {
      //      vel_[v] *= 0.999;
      vert_[v] += vel_[v] * dt;
    }
  }
  //  WriteVectorToMatlab<Real>(v_num_ * 3, &vel_[0][0], std::cout);// exit(0);

  //  Real min[3], max[3];
  //  AffineTransformer<Real>::GetRange(this->VertexArray(), v_num_, min, max);
  //  P(dj::Vec3d(min));
  //  P(dj::Vec3d(max));
}

inline void FemVolumetricObject::ComputeTetForce(int tet_idx, Mat3& force) {
  Vec3* v[4] = {
    &vert_[tet_[tet_idx][0]],
    &vert_[tet_[tet_idx][1]],
    &vert_[tet_[tet_idx][2]],
    &vert_[tet_[tet_idx][3]]
  };
  force[0] = *(v[1]) - *(v[0]);
  force[1] = *(v[2]) - *(v[0]);
  force[2] = *(v[3]) - *(v[0]);
  Real tmp[3][3];
  dj::MulMatrix3x3<Real>(&inv_dm_[tet_idx][0][0], &force[0][0], &tmp[0][0]);
  Real tr = (tmp[0][0] + tmp[1][1] + tmp[2][2] - 3) * lambda_[tet_idx];
  Real two_mu = 2 * mu_[tet_idx];
  Real tmp1[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      tmp1[i][j] = mu_[tet_idx] * (tmp[i][j] + tmp[j][i]);
    }
    tmp1[i][i] += -two_mu + tr;
  }
  dj::MulMatrix3x3<Real>(&tmp1[0][0], &inv_dm_[tet_idx][0][0], &force[0][0]);
  force[0] *= -tet_vol_[tet_idx];
  force[1] *= -tet_vol_[tet_idx];
  force[2] *= -tet_vol_[tet_idx];
}

void FemVolumetricObject::ComputeTetDifferentialForce(int tet_idx, Real3 dx, Mat3 &force) {
  //  force[0] = Vec3(0, 0, 0);
  //  force[0] = Vec3(0, 0, 0);
  //  force[1] = Vec3(0, 0, 0);
  //  return;
  const Real* v[4] = {
    &dx[tet_[tet_idx][0]][0],
    &dx[tet_[tet_idx][1]][0],
    &dx[tet_[tet_idx][2]][0],
    &dx[tet_[tet_idx][3]][0]
  };
  dj::SubVec3<Real>(v[1], v[0], &force[0][0]);
  dj::SubVec3<Real>(v[2], v[0], &force[1][0]);
  dj::SubVec3<Real>(v[3], v[0], &force[2][0]);

  Real tmp[3][3];
  dj::MulMatrix3x3<Real>(&inv_dm_[tet_idx][0][0], &force[0][0], &tmp[0][0]);
  Real tr = (tmp[0][0] + tmp[1][1] + tmp[2][2]) * lambda_[tet_idx];
  Real tmp1[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      tmp1[i][j] = mu_[tet_idx] * (tmp[i][j] + tmp[j][i]);
    }
    tmp1[i][i] += tr;
  }
  dj::MulMatrix3x3<Real>(&tmp1[0][0], &inv_dm_[tet_idx][0][0], &force[0][0]);
  force[0] *= -tet_vol_[tet_idx];
  force[1] *= -tet_vol_[tet_idx];
  force[2] *= -tet_vol_[tet_idx];
}

void FemVolumetricObject::Render(int render_mode) {
  Super::Render(render_mode);
  //  glDisable(GL_LIGHTING);
  //  FOR(v, 0, v_num_) {
  //    DrawArrow<Real>(&vert_[v][0], &force_[v][0], true, 1.0);
  //  }
}

void FemVolumetricObject::ComputeInvDm() {
  for (int t = 0; t < tet_num_; ++t) {
    Vec3* v[4] = {
      &vert_[tet_[t][0]],
      &vert_[tet_[t][1]],
      &vert_[tet_[t][2]],
      &vert_[tet_[t][3]]
    };
    inv_dm_[t][0] = *(v[1]) - *(v[0]);
    inv_dm_[t][1] = *(v[2]) - *(v[0]);
    inv_dm_[t][2] = *(v[3]) - *(v[0]);
    Mat3 tmp(*(v[1]) - * (v[0]), *(v[2]) - * (v[0]), *(v[3]) - * (v[0]));
    dj::Inverse3<Real>((Real3) &tmp[0][0], (Real3) &inv_dm_[t][0][0]);
  }
}
