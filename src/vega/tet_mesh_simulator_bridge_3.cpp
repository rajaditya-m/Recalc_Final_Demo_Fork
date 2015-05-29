#include "tet_mesh_simulator_bridge.h"
#include "tet.h"
#include "conjugate_gradient_solver.h"
#include "StVKForceModel.h"
#include "buffer_manager.h"
#include "forceModel.h"
#include "StVKIsotropicMaterial.h"
#include "StVKElementABCD.h"
#include "StVKTetABCD.h"
#include "StVKInternalForces.h"
#include "StVKStiffnessMatrix.h"
#include "corotationalLinearFEM.h"
#include "vector_lib.h"
#include "isotropicHyperelasticFEM.h"
#include "tetMesh.h"
#include "global.h"
#include "generateMassMatrix.h"
#include "tetMesh.h"
#include "config_file.h"

TetMeshSimulatorBridge::TetMeshSimulatorBridge(Tet *tet_mesh,std::vector<std::pair<int,int> > &interfaceV) {
   sim_mode_ = 1;
  buf_ = new BufferManager;
  gravity_[0] = global::gravity[0];
  gravity_[1] = global::gravity[1];
  gravity_[2] = global::gravity[2];
  {
    dj::Vec3d vega_gravity(gravity_);
    P(vega_gravity);
  }
  tet_mesh_ = tet_mesh;
  v_num_ = tet_mesh->VertexNum();
  tet_num_ = tet_mesh->TetNum();
  tets_ = (int (*)[4]) tet_mesh->TetIndex();
  u_ = (real (*)[3]) buf_->Malloc<real>(v_num_ * 3);
  vel_ = (real (*)[3]) buf_->Malloc<real>(v_num_ * 3);
  internal_force_ = (real (*)[3]) buf_->Malloc<real>(v_num_ * 3);
  rhs_ = (real (*)[3]) buf_->Malloc<real>(v_num_ * 3);
  rest_pos_ = (real (*)[3]) buf_->Malloc<real>(v_num_ * 3);
  tmp_pos_ = (real (*)[3]) buf_->Malloc<real>(v_num_ * 3);
  for (int i = 0; i < v_num_ * 3; ++i) {
    rest_pos_[i / 3][i % 3] = real(tet_mesh->VertexArray()[i]);
  }
  cg_solver_ = new ConjugateGradientSolver<real>(v_num_ * 3);
  // TODO change parameter
  //    double young_module = 1e6;
  //    double poisson_ratio = 0.49;
  //  double young_module = 1.5e5;
  //  double poisson_ratio = 0.45;
  young_modulus_ = conf.Get<double>("young's modulus");
  poisson_ratio_ = conf.Get<double>("poisson ratio");
  density_ = tet_mesh->density_;
  P(young_modulus_, poisson_ratio_, density_);
  vega_mesh_ = new TetMesh(v_num_, (real*) rest_pos_, tet_num_, (int*) tets_, young_modulus_, poisson_ratio_, density_);
  //  vega_mesh_ = new TetMesh(v_num_, (real*) rest_pos_, tet_num_, (int*) tets_, 1e5, 0.49, 1100);
  material_ = new StVKIsotropicMaterial(vega_mesh_, 1, 500);//, enableCompressionResistance, compressionResistance);
  GenerateMassMatrix::computeMassMatrix(vega_mesh_, &mass_matrix_, true);
  element_force_ = NULL;
  element_k_ = NULL;
  interface_points_.resize(interfaceV.size());
  std::copy(interfaceV.begin(),interfaceV.end(),interface_points_.begin());

#ifdef COROTATIONAL_FEM
  L("using corotational FEM");
  inv_fem_force_model_ = new CorotationalLinearFEM(vega_mesh_);
  inv_fem_force_model_->GetStiffnessMatrixTopology(&tangent_stiffness_matrix_);
#elif defined(STVK_FEM)
  L("using StVK FEM");
  StVKTetABCD* tet_abcd = new StVKTetABCD(vega_mesh_);
  StVKInternalForces* st_vk_force = new StVKInternalForces(vega_mesh_, tet_abcd);
  inv_fem_force_model_ = new StVKForceModel(st_vk_force);
  inv_fem_force_model_->GetTangentStiffnessMatrixTopology(&tangent_stiffness_matrix_);
#else
  L("using invertible FEM");
  // TODO uncomment
  inv_fem_force_model_ = new IsotropicHyperelasticFEM(vega_mesh_, material_, 0.1);
  inv_fem_force_model_->GetStiffnessMatrixTopology(&tangent_stiffness_matrix_);
  inv_fem_force_model_->GetStiffnessMatrixTopology(&stitch_tangent_stiffness_matrix_,interface_points_);
  //Precompute this
  double *u_rest = new double[3 * tet_mesh_->vertex_num_];
  memset(u_rest, 0, sizeof(double) * 3 * tet_mesh_->vertex_num_);
  inv_fem_force_model_->GetTangentStiffnessMatrix((real*) u_rest, stitch_tangent_stiffness_matrix_);
  element_force_ = inv_fem_force_model_->elementForce;
  element_k_ = inv_fem_force_model_->elementK;
#endif
  P(dj::Vec3d(gravity_));
}

void TetMeshSimulatorBridge::makeSpringConnections() {
    for(auto it = interface_points_.begin(); it!= interface_points_.end(); it++) {
        this->makeLinks(it->first, it->second);
    }
}

void TetMeshSimulatorBridge::makeLinks(int i, int j) {

    double k = 5.0e4;

    int colIdx;

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*i+0,3*i+0);
    stitch_tangent_stiffness_matrix_->AddEntry(3*i+0,colIdx,k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*i+1,3*i+1);
    stitch_tangent_stiffness_matrix_->AddEntry(3*i+1,colIdx,k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*i+2,3*i+2);
    stitch_tangent_stiffness_matrix_->AddEntry(3*i+2,colIdx,k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*j+0,3*j+0);
    stitch_tangent_stiffness_matrix_->AddEntry(3*j+0,colIdx,k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*j+1,3*j+1);
    stitch_tangent_stiffness_matrix_->AddEntry(3*j+1,colIdx,k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*j+2,3*j+2);
    stitch_tangent_stiffness_matrix_->AddEntry(3*j+2,colIdx,k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*i+0,3*j+0);
    stitch_tangent_stiffness_matrix_->AddEntry(3*i+0,colIdx,-k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*i+1,3*j+1);
    stitch_tangent_stiffness_matrix_->AddEntry(3*i+1,colIdx,-k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*i+2,3*j+2);
    stitch_tangent_stiffness_matrix_->AddEntry(3*i+2,colIdx,-k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*j+0,3*i+0);
    stitch_tangent_stiffness_matrix_->AddEntry(3*j+0,colIdx,-k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*j+1,3*i+1);
    stitch_tangent_stiffness_matrix_->AddEntry(3*j+1,colIdx,-k);

    colIdx = stitch_tangent_stiffness_matrix_->GetInverseIndex(3*j+2,3*i+2);
    stitch_tangent_stiffness_matrix_->AddEntry(3*j+2,colIdx,-k);

}

TetMeshSimulatorBridge::~TetMeshSimulatorBridge() {
  delete tangent_stiffness_matrix_;
  delete stitch_tangent_stiffness_matrix_;
  delete cg_solver_;
  delete material_;
  delete inv_fem_force_model_;
  delete vega_mesh_;
  delete buf_;
}

TetMeshSimulatorBridge::real TetMeshSimulatorBridge::Simulate(real dt) {
  if (1) {
    std::fill(tet_mesh_->is_constrainted_.begin(), tet_mesh_->is_constrainted_.end(), false);
    for (int v = 0; v < tet_mesh_->vertex_num_; ++v) {
      if (tet_mesh_->X[v * 3 + 2] < 1e-3) {
        tet_mesh_->is_constrainted_[v] = true;
      }
    }
  }

#if 1
  real dt_2 = dt * dt;
  auto K = [&](real * x, real * result) {
    //    memcpy((real*) result, (real*) x, sizeof(real) * 3 * v_num_);
    //    return;

    tangent_stiffness_matrix_->MultiplyVector(x, result);
    OMP_FOR
    for (int v = 0; v < v_num_; ++v) {
      int v3 = v * 3;
      if (tet_mesh_->is_constrainted_[v]) {
        result[v3 + 0] = 0;
        result[v3 + 1] = 0;
        result[v3 + 2] = 0;
      } else {
        result[v3 + 0] = x[v3 + 0] * tet_mesh_->mass_[v] + dt_2 * result[v3 + 0];
        result[v3 + 1] = x[v3 + 1] * tet_mesh_->mass_[v] + dt_2 * result[v3 + 1];
        result[v3 + 2] = x[v3 + 2] * tet_mesh_->mass_[v] + dt_2 * result[v3 + 2];
      }
    }
  };
  (void) K;
#ifdef COROTATIONAL_FEM
  inv_fem_force_model_->ComputeForceAndStiffnessMatrix((real*) u_, (real*) internal_force_, tangent_stiffness_matrix_);
#elif defined(STVK_FEM)
  inv_fem_force_model_->GetForceAndMatrix((real*) u_, (real*) internal_force_, tangent_stiffness_matrix_);
#else
  inv_fem_force_model_->GetForceAndTangentStiffnessMatrix((real*) u_, (real*) internal_force_, tangent_stiffness_matrix_);
#endif

  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
    if (tet_mesh_->is_constrainted_[v]) {
      rhs_[v][0] = 0;
      rhs_[v][1] = 0;
      rhs_[v][2] = 0;
    } else {
      rhs_[v][0] = dt * (tet_mesh_->mass_[v] * gravity_[0] - internal_force_[v][0]) + vel_[v][0] * tet_mesh_->mass_[v];
      rhs_[v][1] = dt * (tet_mesh_->mass_[v] * gravity_[1] - internal_force_[v][1]) + vel_[v][1] * tet_mesh_->mass_[v];
      rhs_[v][2] = dt * (tet_mesh_->mass_[v] * gravity_[2] - internal_force_[v][2]) + vel_[v][2] * tet_mesh_->mass_[v];
    }
  }

  double ui_force[3];
  int ui_force_vert = tet_mesh_->GetUIForce(ui_force);
  if (ui_force_vert >= 0) {
    rhs_[ui_force_vert][0] += ui_force[0] * dt;
    rhs_[ui_force_vert][1] += ui_force[1] * dt;
    rhs_[ui_force_vert][2] += ui_force[2] * dt;
  }

#if 0
  // ground collision
  //  double total_collision_force = 0;
  const double kFloor = 0.10;
  const double kFloorStiffness = 5000;
  for (int v = 0; v < v_num_; ++v) {
    if (tet_mesh_->VertexArray()[v * 3 + 1] < kFloor) {
      double collision_force = (kFloor - tet_mesh_->VertexArray()[v * 3 + 1]) * tet_mesh_->mass_[v]  * kFloorStiffness;
      rhs_[v][1] += collision_force * dt;
      //      total_collision_force += collision_force;
    }
  }
  //  P(total_collision_force);
#endif

  //  memcpy((real*) vel_, (real*) rhs_, sizeof(real) * 3 * v_num_);
//  std::pair<int, double> info = cg_solver_->Solve((real*) rhs_, (real*) vel_, K, 3000, 1e-11);
//  P(info.first, info.second);
#endif

  //  for (int i = 0; i < tet_mesh_->attached_tet_of_joints_.size(); ++i) {
  //    int v = tet_mesh_->attached_tet_of_joints_[i];
  //    vel_[v][0] = tet_mesh_->joint_vel_[i * 3 + 0];
  //    vel_[v][1] = tet_mesh_->joint_vel_[i * 3 + 1];
  //    vel_[v][2] = tet_mesh_->joint_vel_[i * 3 + 2];
  //    P(dj::Vec3d(vel_[v]));
  //  }

  //  for (int i = 0; i < (int) tet_mesh_->constrainted_vertex_.size(); ++i) {
  //    int v = tet_mesh_->constrainted_vertex_[i];
  //    vel_[v][0] = tet_mesh_->constrainted_vert_vel_[i * 3 + 0];
  //    vel_[v][1] = tet_mesh_->constrainted_vert_vel_[i * 3 + 1];
  //    vel_[v][2] = tet_mesh_->constrainted_vert_vel_[i * 3 + 2];
  //  }


  real max_velocity = 0;
  for (int v = 0; v < v_num_; ++v) {
    max_velocity = dj::Max(vel_[v][0] * vel_[v][0] + vel_[v][1] * vel_[v][1] + vel_[v][2] * vel_[v][2], max_velocity);
  }
  max_velocity = std::sqrt(max_velocity);
  //  real diff = 0.0f;
  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
    if (tet_mesh_->is_constrainted_[v]) {
      vel_[v][0] = 0;
      vel_[v][1] = 0;
      vel_[v][2] = 0;
    } else {
      u_[v][0] += vel_[v][0] * dt;
      u_[v][1] += vel_[v][1] * dt;
      u_[v][2] += vel_[v][2] * dt;

      const real kDamping = 1.00;
      vel_[v][0] *= kDamping;
      vel_[v][1] *= kDamping;
      vel_[v][2] *= kDamping;

      tet_mesh_->VertexArray()[v * 3 + 0] = u_[v][0] + rest_pos_[v][0];
      tet_mesh_->VertexArray()[v * 3 + 1] = u_[v][1] + rest_pos_[v][1];
      tet_mesh_->VertexArray()[v * 3 + 2] = u_[v][2] + rest_pos_[v][2];
    }
  }
  return max_velocity;
}

void TetMeshSimulatorBridge::Reset() {
  memset((real*) u_, 0, sizeof(real) * v_num_ * 3);
  memset((real*) vel_, 0, sizeof(real) * v_num_ * 3);
}

SparseMatrix* TetMeshSimulatorBridge::GetStitchStiffnessMatrixPointer() {
    //double *u_rest = new double[3 * tet_mesh_->vertex_num_];
    //memset(u_rest, 0, sizeof(double) * 3 * tet_mesh_->vertex_num_);
    //inv_fem_force_model_->GetForceAndTangentStiffnessMatrix((real*) u_rest, (real*) internal_force_, stitch_tangent_stiffness_matrix_);
    makeSpringConnections();
    //stitch_tangent_stiffness_matrix_->SaveToMatlabFormat("sav.mm");
    return stitch_tangent_stiffness_matrix_;

}



void TetMeshSimulatorBridge::ComputeInternalForceAndTangentStiffnessMatrix(double dt) {
  (void) dt;
#ifdef STVK_FEM
  inv_fem_force_model_->GetForceAndMatrix((real*) u_, (real*) internal_force_, tangent_stiffness_matrix_);
#else
  //  memset(u_, 0, sizeof(double) * 3 * tet_mesh_->vertex_num_);
  if(sim_mode_==1)
    inv_fem_force_model_->GetForceAndTangentStiffnessMatrix((real*) u_, (real*) internal_force_, tangent_stiffness_matrix_);
  else {
      inv_fem_force_model_->GetForceAndTangentStiffnessMatrix((real*) u_, (real*) internal_force_, stitch_tangent_stiffness_matrix_);
      makeSpringConnections();
  }
  //  P(blas::dotd(tet_mesh_->vertex_num_ * 3, (real*) internal_force_, (real*) internal_force_)); exit(0);
#endif

  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
    //    if (tet_mesh_->is_constrainted_[v]) {
    //      rhs_[v][0] = 0;
    //      rhs_[v][1] = 0;
    //      rhs_[v][2] = 0;
    //    } else
    {
      rhs_[v][0] = (tet_mesh_->mass_[v] * gravity_[0] - internal_force_[v][0]);
      rhs_[v][1] = (tet_mesh_->mass_[v] * gravity_[1] - internal_force_[v][1]);
      rhs_[v][2] = (tet_mesh_->mass_[v] * gravity_[2] - internal_force_[v][2]);
    }
  }
}

void TetMeshSimulatorBridge::ComputePartialInternalForceAndTangentStiffnessMatrix(std::vector<int>& tets) {
  //  int computationMode = Fem::COMPUTE_INTERNALFORCES | Fem::COMPUTE_TANGENTSTIFFNESSMATRIX;
  //  inv_fem_force_model_->GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue((real*) u_, NULL, (real*) internal_force_, tangent_stiffness_matrix_, computationMode);

  //  for (int i = 0; i < v_num_ * 3; i++)
  //    inv_fem_force_model_->currentVerticesPosition[i] = inv_fem_force_model_->restVerticesPosition[i] + u_[i / 3][i % 3];

  inv_fem_force_model_->ComputePartialForceAndStiffnessMatrix(tets, (real*) u_);
}

void TetMeshSimulatorBridge::LoadPosition(double * pos) {
  for (int v = 0; v < v_num_; ++v) {
    u_[v][0] = pos[v * 3 + 0] - rest_pos_[v][0];
    u_[v][1] = pos[v * 3 + 1] - rest_pos_[v][1];
    u_[v][2] = pos[v * 3 + 2] - rest_pos_[v][2];
  }
}

void TetMeshSimulatorBridge::UpdateOffset() {
  OMP_FOR
  for (int v = 0; v < v_num_; ++v) {
    u_[v][0] = tet_mesh_->VertexArray()[v * 3 + 0] - rest_pos_[v][0];
    u_[v][1] = tet_mesh_->VertexArray()[v * 3 + 1] - rest_pos_[v][1];
    u_[v][2] = tet_mesh_->VertexArray()[v * 3 + 2] - rest_pos_[v][2];
  }
}

void TetMeshSimulatorBridge::SetInterface(std::vector<std::pair<int,int> > &infVerts) {
    interface_points_ = infVerts;
}
