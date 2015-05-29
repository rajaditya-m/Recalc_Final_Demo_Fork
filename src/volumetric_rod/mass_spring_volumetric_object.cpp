#include "global.h"
#include "mass_spring_volumetric_object.h"
#include "conjugate_gradient_solver.h"
#include "tet_collider.h"
#include "vector_io.h"

MassSpringVolumetricObject::MassSpringVolumetricObject(TetrahedralMeshIO *mesh_io,
                                                       const char *file_name,
                                                       AffineTransformer<Real>* transformer)
  : Super(mesh_io, file_name, transformer)
{
  Construct();
}


MassSpringVolumetricObject::MassSpringVolumetricObject(TetrahedralMesh *mesh)
  : Super(mesh)
{
  Construct();
}

void MassSpringVolumetricObject::Construct()
{
  //  stiffness_ = 10000;
  stiffness_ = 100000000;
  damping_ = 0.02;
  cg_solver_ = new ConjugateGradientSolver<Real>(v_num_ * 3);
  vel_.resize(v_num_);
  prev_vert_.resize(v_num_);
  tmp_vert_buf_.resize(v_num_);
  std::fill(tet_density_.begin(), tet_density_.end(), 1200);
  rest_edge_length_.resize(e_num_);
  ComputeRestEdgeLength();
  ComputeTetrahedraVolume();
  ComputeLumpedMass();
//    mass_.resize(v_num_);
//    std::fill(mass_.begin(), mass_.end(), Real(1));
  prev_vert_ = vert_;
  edge_stiffness_.resize(e_num_);
  Real avg_stiffness = 0;
  Real min_stiff = 99999, max_stiff = -99999;
  for (int e = 0; e < e_num_; ++e) {
    edge_stiffness_[e] = stiffness_ * std::pow(rest_edge_length_[e], 3);
    min_stiff = std::min(min_stiff, edge_stiffness_[e]);
    max_stiff = std::max(max_stiff, edge_stiffness_[e]);
    avg_stiffness += edge_stiffness_[e];
  }
//  std::fill_n(&edge_stiffness_[0], e_num_, 150000);
  avg_stiffness /= e_num_;
  Real avg_mass = total_mass_ / v_num_;
  P(avg_stiffness, min_stiff, max_stiff, avg_stiffness / avg_mass);
  P(total_mass_);
  collider_ = new TetCollider(&vert_[0][0], v_num_, &tet_[0][0], tet_num_, &edge_[0][0], e_num_);
  P(total_mass_, total_volume_, total_mass_ / total_volume_);
}

MassSpringVolumetricObject::~MassSpringVolumetricObject()
{
  delete cg_solver_;
  delete collider_;
}

void MassSpringVolumetricObject::ImplicitStepSolvePosition(Real dt)
{
  using namespace dj;
  std::vector<Real> force_gradient(e_num_ * 9, 0);
  // Compute force gradient for each edge
  // OMP_FOR
  for (int e = 0; e < e_num_; ++e) {
    int e9 = e * 9;
    int* v = &edge_[e][0];
    Real* pos[] = {
      &vert_[v[0]][0],
      &vert_[v[1]][0],
    };
    Real direction[3];
    dj::SubVec3(pos[1], pos[0], direction);
    Real length = dj::Normalize3(direction);
    if (length < rest_edge_length_[e]) {
      Real (*matrix)[3] = (Real (*)[3]) &force_gradient[e9];
      matrix[0][0] = direction[0] * direction[0];
      matrix[0][1] = direction[0] * direction[1];
      matrix[0][2] = direction[0] * direction[2];

      matrix[1][0] = direction[1] * direction[0];
      matrix[1][1] = direction[1] * direction[1];
      matrix[1][2] = direction[1] * direction[2];

      matrix[2][0] = direction[2] * direction[0];
      matrix[2][1] = direction[2] * direction[1];
      matrix[2][2] = direction[2] * direction[2];
      {
        Real * tmp = &force_gradient[e9];
        for (int i = 0; i < 9; ++i) {
          tmp[i] *= edge_stiffness_[e];
        }
      }
    } else {
      Real factor = edge_stiffness_[e] * rest_edge_length_[e] / length;
      Real (*matrix)[3] = (Real (*)[3]) &force_gradient[e9];
      matrix[0][0] = direction[0] * direction[0] - 1;
      matrix[0][1] = direction[0] * direction[1];
      matrix[0][2] = direction[0] * direction[2];

      matrix[1][0] = direction[1] * direction[0];
      matrix[1][1] = direction[1] * direction[1] - 1;
      matrix[1][2] = direction[1] * direction[2];

      matrix[2][0] = direction[2] * direction[0];
      matrix[2][1] = direction[2] * direction[1];
      matrix[2][2] = direction[2] * direction[2] - 1;

      {
        Real * tmp = &force_gradient[e9];
        for (int i = 0; i < 9; ++i) {
          tmp[i] *= factor;
        }
      }
      matrix[0][0] += edge_stiffness_[e];
      matrix[1][1] += edge_stiffness_[e];
      matrix[2][2] += edge_stiffness_[e];
    }
  }

  Real t2 = dt * dt;
  //    t2 = dt;
  const Real kDampingFactor = damping_ * dt;
  //  const Real kDampingFactor = 0;
  // result = (I - h^2 * M^-1 * dF/dx) * x
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      const Real  factor = kDampingFactor / mass_[v];
      result[v3 + 0] += factor * x[v3 + 0];
      result[v3 + 1] += factor * x[v3 + 1];
      result[v3 + 2] += factor * x[v3 + 2];
      Real tmp_result[3] = {0, 0, 0};
      for (int num = 0; num < (int) incident_edge_[v].size(); ++num) {
        int e = incident_edge_[v][num];
        Real* fg = &force_gradient[e * 9];
        int* verts = &edge_[e][0];
        Real* pos[] = {
          x + verts[0] * 3,
          x + verts[1] * 3,
        };
        Real diff[3];
        dj::SubVec3(pos[1], pos[0], diff);
        Real sign = (verts[0] == v) ? 1 : -1;
        Real tmp[3];
        dj::MulMatrix3x3Vec<Real>((Real (*)[3]) fg, diff, tmp);
        tmp_result[0] += tmp[0] * sign;
        tmp_result[1] += tmp[1] * sign;
        tmp_result[2] += tmp[2] * sign;
      }
      result[v3 + 0] += -t2 * tmp_result[0] / mass_[v];
      result[v3 + 1] += -t2 * tmp_result[1] / mass_[v];
      result[v3 + 2] += -t2 * tmp_result[2] / mass_[v];
    }
  };
  //  PrintImplicitMatrix<Real>(std::cout, StiffnessMatrix, v_num_ * 3);
  //  exit(0);

  // Assemble right hand side
  // rhs = (I - h^2 * M^-1 * dF/dx) * x_0 + x_0 - x_-1 + h^2 * M^-1 F(x_0) + h^2 * g
  std::vector<Real> rhs(v_num_ * 3, 0);
  //  memcpy(&rhs[0], &vert_[0][0], sizeof(Real) * 3 * v_num_);
  StiffnessMatrix(&vert_[0][0], &rhs[0]); // rhs = (I - h^2 * M^-1 * dF/dx) * x_0
  OMP_FOR
  for (int v_idx = 0; v_idx < v_num_; ++v_idx) {
    int v3 = v_idx * 3;
    Real* force = &rhs[v3];
    force[0] += t2 * global::gravity[0];
    force[1] += t2 * global::gravity[1];
    force[2] += t2 * global::gravity[2];
    //    if (v_idx >= 4) {
    //      force[0] += t2 * 9.8;
    //    } else {
    //      force[0] -= t2 * 9.8;
    //    }

    //    force[0] += vertex_[v3 + 0] - prev_vertex_[v3 + 0];
    //    force[1] += vertex_[v3 + 1] - prev_vertex_[v3 + 1];
    //    force[2] += vertex_[v3 + 2] - prev_vertex_[v3 + 2];
    vel_[v_idx] *= damping_;
    force[0] += vel_[v_idx][0] * dt;
    force[1] += vel_[v_idx][1] * dt;
    force[2] += vel_[v_idx][2] * dt;

    for (int num = 0; num < (int) incident_edge_[v_idx].size(); ++num) {
      int e = incident_edge_[v_idx][num];
      int* verts = &edge_[e][0];
      Real* pos[] = {
        &vert_[verts[0]][0],
        &vert_[verts[1]][0],
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      Real sign = (verts[1] == v_idx) ? 1 : -1;
      Real length = dj::Normalize3(diff);

      force[0] += -edge_stiffness_[e] * sign * diff[0] * (length - rest_edge_length_[e]) * t2 / mass_[v_idx];
      force[1] += -edge_stiffness_[e] * sign * diff[1] * (length - rest_edge_length_[e]) * t2 / mass_[v_idx];
      force[2] += -edge_stiffness_[e] * sign * diff[2] * (length - rest_edge_length_[e]) * t2 / mass_[v_idx];
    }
  }

  prev_vert_ = vert_;
  cg_solver_->Solve(&rhs[0], &vert_[0][0], StiffnessMatrix, 1000, 1e-40);

  //    vert_[0] = prev_vert_[0];
  OMP_FOR
  for (int v_idx = 0; v_idx < v_num_; ++v_idx) {
    if (vert_[v_idx][1] <= 0) {
      vert_[v_idx][1] = 0;
    }
    vel_[v_idx] = (vert_[v_idx] - prev_vert_[v_idx]) / dt;
  }
}

void MassSpringVolumetricObject::ImplicitStepSolveVelocity(Real dt)
{
  using namespace dj;
  std::vector<Real> force_gradient(e_num_ * 9, 0);
  // Compute force gradient for each edge
  // OMP_FOR
  for (int e = 0; e < e_num_; ++e) {
    int e9 = e * 9;
    int* v = &edge_[e][0];
    Real* pos[] = {
      &vert_[v[0]][0],
      &vert_[v[1]][0],
    };
    Real direction[3];
    dj::SubVec3(pos[1], pos[0], direction);
    Real length = dj::Normalize3(direction);
    if (length < rest_edge_length_[e]) {
      Real (*matrix)[3] = (Real (*)[3]) &force_gradient[e9];
      matrix[0][0] = direction[0] * direction[0];
      matrix[0][1] = direction[0] * direction[1];
      matrix[0][2] = direction[0] * direction[2];

      matrix[1][0] = direction[1] * direction[0];
      matrix[1][1] = direction[1] * direction[1];
      matrix[1][2] = direction[1] * direction[2];

      matrix[2][0] = direction[2] * direction[0];
      matrix[2][1] = direction[2] * direction[1];
      matrix[2][2] = direction[2] * direction[2];
      {
        Real * tmp = &force_gradient[e9];
        for (int i = 0; i < 9; ++i) {
          tmp[i] *= edge_stiffness_[e];
        }
      }
    } else {
      Real factor = edge_stiffness_[e] * rest_edge_length_[e] / length;
      Real (*matrix)[3] = (Real (*)[3]) &force_gradient[e9];
      matrix[0][0] = direction[0] * direction[0] - 1;
      matrix[0][1] = direction[0] * direction[1];
      matrix[0][2] = direction[0] * direction[2];

      matrix[1][0] = direction[1] * direction[0];
      matrix[1][1] = direction[1] * direction[1] - 1;
      matrix[1][2] = direction[1] * direction[2];

      matrix[2][0] = direction[2] * direction[0];
      matrix[2][1] = direction[2] * direction[1];
      matrix[2][2] = direction[2] * direction[2] - 1;
      {
        Real * tmp = &force_gradient[e9];
        for (int i = 0; i < 9; ++i) {
          tmp[i] *= factor;
        }
      }
      matrix[0][0] += edge_stiffness_[e];
      matrix[1][1] += edge_stiffness_[e];
      matrix[2][2] += edge_stiffness_[e];
    }
  }

  const Real t2 = dt * dt;
  const Real kDampingFactor = damping_ * dt;
  auto StiffnessMatrix = [&](Real * x, Real * result) {
    memcpy(result, x, sizeof(Real) * v_num_ * 3);
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      const Real factor = kDampingFactor / mass_[v];
      result[v3 + 0] += factor * x[v3 + 0];
      result[v3 + 1] += factor * x[v3 + 1];
      result[v3 + 2] += factor * x[v3 + 2];
      Real tmp_result[3] = {0, 0, 0};
      for (int num = 0; num < (int) incident_edge_[v].size(); ++num) {
        int e = incident_edge_[v][num];
        Real* fg = &force_gradient[e * 9];
        int* verts = &edge_[e][0];
        Real* pos[] = {
          x + verts[0] * 3,
          x + verts[1] * 3,
        };
        Real diff[3];
        dj::SubVec3(pos[1], pos[0], diff);
        Real sign = (verts[0] == v) ? 1 : -1;
        Real tmp[3];
        dj::MulMatrix3x3Vec<Real>((Real (*)[3]) fg, diff, tmp);
        tmp_result[0] += tmp[0] * sign;
        tmp_result[1] += tmp[1] * sign;
        tmp_result[2] += tmp[2] * sign;
      }
      result[v3 + 0] += -t2 * tmp_result[0] / mass_[v];
      result[v3 + 1] += -t2 * tmp_result[1] / mass_[v];
      result[v3 + 2] += -t2 * tmp_result[2] / mass_[v];
    }
  };
  //    PrintImplicitMatrix<Real>(std::cout, StiffnessMatrix, v_num_ * 3);
  //    exit(0);

  // Assemble right hand side
  memcpy(&tmp_vert_buf_[0][0], &vel_[0][0], sizeof(Vec3) * v_num_);
  //    for (int v : force_vert)
//  {
//    int v = 3019;
//    tmp_vert_buf_[v][1] = -100000 * dt / mass_[v];
//  }
  OMP_FOR
  for (int v_idx = 0; v_idx < v_num_; ++v_idx) {
    Real* force = &tmp_vert_buf_[v_idx][0];
    //    force[0] += dt * global::gravity[0];
    //    force[1] += dt * global::gravity[1];
    //    force[2] += dt * global::gravity[2];
    //    vel_[v_idx] *= kDamping_;
    //    force[0] += vel_[v_idx][0];
    //    force[1] += vel_[v_idx][1];
    //    force[2] += vel_[v_idx][2];

    for (int num = 0; num < (int) incident_edge_[v_idx].size(); ++num) {
      int e = incident_edge_[v_idx][num];
      int* verts = &edge_[e][0];
      Real* pos[] = {
        &vert_[verts[0]][0],
        &vert_[verts[1]][0],
      };
      Real diff[3];
      dj::SubVec3(pos[1], pos[0], diff);
      Real sign = (verts[1] == v_idx) ? 1 : -1;
      Real length = dj::Normalize3(diff);

      force[0] += -edge_stiffness_[e] * sign * diff[0] * (length - rest_edge_length_[e]) * dt / mass_[v_idx];
      force[1] += -edge_stiffness_[e] * sign * diff[1] * (length - rest_edge_length_[e]) * dt / mass_[v_idx];
      force[2] += -edge_stiffness_[e] * sign * diff[2] * (length - rest_edge_length_[e]) * dt / mass_[v_idx];
    }
  }
  cg_solver_->Solve(&tmp_vert_buf_[0][0], &vel_[0][0], StiffnessMatrix, 1000, 1e-10);
  //  vel_[0] = Vec3(0, 0, 0);
  OMP_FOR
  for (int v_idx = 0; v_idx < v_num_; ++v_idx) {
    vert_[v_idx] += vel_[v_idx] * dt;
    if (vert_[v_idx][1] < 0) {
      vel_[v_idx][1] -= vert_[v_idx][1] / dt;
      vert_[v_idx][1] = 0;
    }
  }
}

void MassSpringVolumetricObject::ExplicitStep(Real dt)
{
  for (int v = 0; v < v_num_; ++v) {
    vel_[v][0] += global::gravity[0] * dt;
    vel_[v][1] += global::gravity[1] * dt;
    vel_[v][2] += global::gravity[2] * dt;
  }

  for (int e = 0; e < e_num_; ++e) {
    Vec3 v0 = vert_[edge_[e][0]];
    Vec3 v1 = vert_[edge_[e][1]];
    Vec3 direction = v1 - v0;
    Real length = direction.Magnitude();
    direction *= Real(1 / length);
    Vec3 force = (length - rest_edge_length_[e]) * edge_stiffness_[e] * direction / dt;
    vel_[edge_[e][0]] += force;
    vel_[edge_[e][1]] -= force;
  }

  for (int v = 0; v < v_num_; ++v) {
    vel_[v] *= damping_;
    vert_[v] += vel_[v] * dt;
    if (vert_[v][1] < 0) {
      vel_[v][1] -= (vert_[v][1]) / dt;
      vert_[v][1] = 0;
    }
  }
}

void MassSpringVolumetricObject::Simulate(Real dt, bool handle_collision)
{
  //    ExplicitSpringStep(dt);
  //  EXECUTE_TIMES(1);
  //  ImplicitSpringStepSolvePosition(dt);
  ImplicitStepSolveVelocity(dt);
  if (handle_collision) {
    collider_->HandleCollision(&vert_[0][0], &tmp_vert_buf_[0][0]);
    vert_.swap(tmp_vert_buf_);;
    for (int v = 0; v < v_num_; ++v) {
      vel_[v] += (vert_[v] - tmp_vert_buf_[v]) / dt;
    }
  }
}

void MassSpringVolumetricObject::ComputeRestEdgeLength()
{
  for (int e = 0; e < e_num_; ++e) {
    Vec3 diff = vert_[edge_[e][0]] - vert_[edge_[e][1]];
    rest_edge_length_[e] = diff.Magnitude();
  }
}
