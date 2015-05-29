#include "armadillo"
#include "basis_generator.h"
#include "RedSVD_2.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include "matrixPCA.h"
#include "tetMesh.h"
#include "linearSolver.h"
#include "CGSolver.h"
#include "performanceCounter.h"
#include "PardisoSolver.h"
#include "sparseMatrix.h"
#include "generateMassMatrix.h"
#include "StVKStiffnessMatrix.h"
#include "StVKHessianTensor.h"
#include "computeStiffnessMatrixNullspace.h"
#include "StVKCubeABCD.h"
#include "print_macro.h"
#include "matrixIO.h"
#include "insertRows.h"
#include "StVKElementABCDLoader.h"
#include "invKSolver.h"
#include "ARPACKSolver.h"
#include "basis_io.h"
#include "binary_file_io.h"
#include "MatrixOps.h"
#include "string_formatter.h"

#define COLMAJORIDX(r,c,nrows) (c*nrows+r)

#include "arpackpp/argsym.h"
#include "global.h"
#include "profiler.h"
class Profiler;
extern Profiler profiler;

#ifndef M_PI
#define M_PI 3.1415926
#endif

#define COLMAJORIDX(r,c,nrows) ((c)*(nrows)+(r))

using std::set;


BasisGenerator::BasisGenerator(int vert_num, double* vert_pos,
                               int tet_num, int* tets,
                               double youngs_modulus,
                               double poisson_ratio,
                               double density,
                               double* mass,
                               std::set<int> *fixed_vertex_set) {
  mesh_ = new TetMesh(vert_num, vert_pos, tet_num, tets, youngs_modulus, poisson_ratio, density);
  Construct(mass, fixed_vertex_set);
  stiffyMultiplier_ = 1;
}

BasisGenerator::BasisGenerator(TetMesh* mesh, double *mass, std::set<int> *fixed_vertex_set) {
  mesh_ = new TetMesh(*mesh);
  Construct(mass, fixed_vertex_set);
  stiffyMultiplier_ = 1;
}

void BasisGenerator::Construct(double *mass, std::set<int> *fixed_vertex_set) {
  linear_mode_num_ = -1;
  non_linear_mode_num_ = -1;
  modal_derivative_num_ = -1;
  rigid_mode_num_ = -1;
  if (fixed_vertex_set) fixed_vertices_ = *fixed_vertex_set;
  if (mass == NULL) {
    GenerateMassMatrix::computeMassMatrix(mesh_, &mass_matrix_, true);
  } else {
    SparseMatrixOutline matrix_outline(mesh_->getNumVertices() * 3);
    for (int i = 0; i < mesh_->getNumVertices() * 3; ++i) {
      matrix_outline.AddEntry(i, i, mass[i / 3]);
    }
    mass_matrix_ = new SparseMatrix(&matrix_outline);
    mass_matrix_eig_.resize(mesh_->getNumVertices() * 3,mesh_->getNumVertices() * 3);
    for (int i = 0; i < mesh_->getNumVertices() * 3; ++i) {
      mass_matrix_eig_.insert(i, i) = mass[i / 3];
    }
  }

  double *ones = new double[mesh_->getNumVertices() * 3];
  massInv_.resize(mesh_->getNumVertices() * 3);
  massFlat_.resize(mesh_->getNumVertices() * 3);
  massSqrt_.resize(mesh_->getNumVertices() * 3);
  for(int i =0 ; i< mesh_->getNumVertices() * 3; i++ ) {
      ones[i] = 1.0;
  }
  mass_matrix_->MultiplyVector(ones,&massInv_[0]);
   mass_matrix_->MultiplyVector(ones,&massFlat_[0]);
   mass_matrix_->MultiplyVector(ones,&massSqrt_[0]);
  for(int i =0 ; i< mesh_->getNumVertices() * 3; i++ ) {
      massInv_[i] = 1.0/massInv_[i];
      massSqrt_[i] = sqrt(massSqrt_[i]);
  }



}

void BasisGenerator::preLoad(const char *basis_prefix){
    //Read the non_linear_modes
    int rbck;
    char file_name[512];
    sprintf(file_name, "%s.basis.bin", basis_prefix);
    BinaryFileReader in(file_name);
    int row_num,basis_num_;
    in.Read(&row_num, 1);
    in.Read(&basis_num_, 1);
    non_linear_modes_.resize(row_num*basis_num_);
    in.Read(&non_linear_modes_[0], row_num * basis_num_);
    non_linear_mode_num_ = basis_num_;
    rbck = row_num;

    //Read the pure_eigen_values
    sprintf(file_name, "%s.pure_eigen_vals.bin", basis_prefix);
    BinaryFileReader in2(file_name);
    in2.Read(&row_num, 1);
    pure_eigen_values_.resize(row_num);
    in2.Read(&pure_eigen_values_[0], row_num);


    //Read the pure eigen vectors
    sprintf(file_name, "%s.pure_eigen_vecs.bin", basis_prefix);
    BinaryFileReader in3(file_name);
    in3.Read(&row_num, 1);
    pure_eigen_vectors_.resize(row_num);
    in3.Read(&pure_eigen_vectors_[0], row_num);

    /*sprintf(file_name,"%s_write_lin_basis.txt",basis_prefix);
    std::ofstream fout(file_name);
    fout << rbck << "\n";
    fout << "30\n";
    for(int i = 0; i<rbck;i++) {
        for(int j=0;j < 30; j++) {
            fout << pure_eigen_vectors_[COLMAJORIDX(i,j,rbck)] << " ";
        }
        fout << "\n";
    }
    fout.close();*/


    //Read the RHS Files
    sprintf(file_name, "%s.rhs_original.bin", basis_prefix);
    BinaryFileReader in4(file_name);
    in4.Read(&row_num, 1);
    in4.Read(&basis_num_, 1);
    Eigen::MatrixXd tmp_rhs(row_num, basis_num_);
    in4.Read(tmp_rhs.data(), row_num * basis_num_);
    rhsOriginal_ = tmp_rhs;
    numColsOriginalRHS_ = basis_num_;

    //Read the hessian tensors also
    sprintf(file_name, "%s.default_hessian.bin", basis_prefix);
    stVKStiffnessHessian = new StVKHessianTensor(file_name,mesh_->getNumVertices(),-1,-1); //Remeber last two are dummies

    //Read the pure linear frequencies
    sprintf(file_name, "%s.lin_freqs.bin", basis_prefix);
    BinaryFileReader in5(file_name);
    in5.Read(&row_num, 1);
    frequencies_.resize(row_num);
    in5.Read(&frequencies_[0], row_num);
    linear_mode_num_ = row_num;
}

BasisGenerator::~BasisGenerator() {
  delete mesh_;
  delete mass_matrix_;
  delete stVKStiffnessHessian;
  delete stitched_stiffness_matrix_;
}

 void BasisGenerator::setStitchStiffnessMatrix(SparseMatrix* spm) {
     stitched_stiffness_matrix_ = new SparseMatrix(*spm);
     //stitched_stiffness_matrix_->SaveToMatlabFormat("BAD.mm");
 }

void BasisGenerator::SetFixedVertices(std::set<int> fixed_verts) {
  fixed_vertices_ = fixed_verts;
  numVertsToRemove_ = fixed_vertices_.size()*3;
  vertsToRemove_.resize(numVertsToRemove_);
  int i = 0;
  for(int fv : fixed_vertices_) {
      vertsToRemove_[3*i+0] = 3*fv+0;
      vertsToRemove_[3*i+1] = 3*fv+1;
      vertsToRemove_[3*i+2] = 3*fv+2;
      i++;
  }
  std::cout <<vertsToRemove_.size() << "\n";
}

void BasisGenerator::SetInterfaceVertices(std::vector<std::pair<int, int> > &iv, int split) {
    P(split);
    interface_vertex_ = iv;
    if(split) {
        for(auto it = interface_vertex_.begin(); it!= interface_vertex_.begin()+split; it++) {
            int v1 = it->first;
            int v2 = it->second;
            for(int d=0;d<3;d++) {
                nonZeroColumns_1_.push_back(3*v1+d);
                columnInformation_1_.push_back(3*v2+d);
            }
        }
        for(auto it = interface_vertex_.begin()+split; it!= interface_vertex_.end(); it++) {
            int v1 = it->first;
            int v2 = it->second;
            for(int d=0;d<3;d++) {
                nonZeroColumns_2_.push_back(3*v1+d);
                columnInformation_2_.push_back(3*v2+d);
            }
        }
    } else {
        for(auto it = interface_vertex_.begin(); it!= interface_vertex_.end(); it++) {
            int v1 = it->first;
            int v2 = it->second;
            for(int d=0;d<3;d++) {
                nonZeroColumns_1_.push_back(3*v1+d);
                columnInformation_1_.push_back(3*v2+d);
            }
        }
    }
}

void BasisGenerator::RemoveSixRigidModes(int numVectors, double *x) {
  int n3 = 3 * mesh_->getNumVertices();

  // remove six rigid modes from rhs
  double * nullspace6 = (double*) malloc (sizeof(double) * n3 * 6);
  double * defoPos6 = (double*) malloc (sizeof(double) * n3);
  for (int i = 0; i < n3 / 3; i++) {
    Vec3d restPos = *((mesh_)->getVertex(i));
    for (int j = 0; j < 3; j++)
      defoPos6[3 * i + j] = restPos[j];
  }

  ComputeStiffnessMatrixNullspace::ComputeNullspace(n3 / 3, defoPos6, nullspace6, 1, 1);
  free(defoPos6);

  for (int i = 0; i < numVectors; i++) {
    ComputeStiffnessMatrixNullspace::RemoveNullspaceComponent(n3 / 3, 6, nullspace6, &x[ELT(n3, 0, i)]);
  }

  free(nullspace6);
}

void BasisGenerator::SetFixedVertexBasisZero(int basis_num, double *basis) {
  for (int v : fixed_vertices_) {
    for (int i = 0; i < basis_num; ++i) {
      basis[i * mesh_->getNumVertices() * 3 + v * 3 + 0] = 0;
      basis[i * mesh_->getNumVertices() * 3 + v * 3 + 1] = 0;
      basis[i * mesh_->getNumVertices() * 3 + v * 3 + 2] = 0;
    }
  }
}

void BasisGenerator::ComputeModalDerivatives() {
  if (linear_mode_num_ < 0) {
    std::cerr << "BasisGenerator::ComputeModalDerivatives() => linear modes are not computed yet!!!" << std::endl;
    return;
  }
  int dummy = 0;
  int* code = &dummy;

  // create stiffness matrix
  StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(mesh_);
  StVKInternalForces * internalForces = new StVKInternalForces(mesh_, precomputedIntegrals);

  // create stiffness matrix
  int n3 = 3 * mesh_->getNumVertices();
  SparseMatrix * stiffnessMatrix;
  StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
  stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
  double * zero = (double*) calloc(n3, sizeof(double));
  stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);
  free(zero);

  // now, the stiffness matrix is computed
  // constrain the degrees of freedom
  int numConstrainedVertices = (int) (fixed_vertices_.size());
  int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
  set<int> :: iterator iter;
  int i = 0;
  for (iter = fixed_vertices_.begin(); iter != fixed_vertices_.end(); iter++) {
    constrainedDOFs[3 * i + 0] = 3 * (*iter) + 1;
    constrainedDOFs[3 * i + 1] = 3 * (*iter) + 2;
    constrainedDOFs[3 * i + 2] = 3 * (*iter) + 3;
    i++;
  }

  int oneIndexed = 1;
  stiffnessMatrix->RemoveRowsColumns(3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  int numRetainedDOFs = stiffnessMatrix->Getn();

  // generate rhs side
  bool computeHessianAtZero = true;
  stVKStiffnessHessian = new StVKHessianTensor(stiffnessMatrixClass);

  int numUsedLinearModes = linear_mode_num_ - rigid_mode_num_;
  //modal_derivative_num_ = numUsedLinearModes * (numUsedLinearModes + 1) / 2;
  modal_derivative_num_ = 321;
  double * rhs = (double*) malloc (sizeof(double) * n3 * modal_derivative_num_);
  if (!rhs) {
    printf("Error: could not allocate space for all modal derivatives.\n");
    *code = 1;
    delete(precomputedIntegrals);
    delete(stiffnessMatrixClass);
    delete(internalForces);
    return;
  }

  if (computeHessianAtZero) {
    // compute hessian at zero
    if (stVKStiffnessHessian->ComputeHessianAtZero(0) != 0) {
      printf("Error: failed to evaluate the Hessian at the origin.\n");
      *code = 1;
      delete(precomputedIntegrals);
      delete(stiffnessMatrixClass);
      delete(internalForces);
      return;
    }
  }

  char file_name_[512];
  sprintf(file_name_, "%s.default_hessian.bin", prefix_);
  stVKStiffnessHessian->SaveHessianAtZeroToFile(file_name_);

  printf("Preparing to compute %d modal derivatives...\n", modal_derivative_num_);

  profiler.Start("Modal Ders naive");
  profiler.Start("Hessian Pdk Naive");

  double * Ulin = &linear_modes_[0];
  if (computeHessianAtZero) {
    printf("Using the high-memory version.\n");
    int pos = 0;

      for (int i = 0; i < 6; i++) {
  //      printf("%d: ", i); fflush(NULL);
        for (int j = i; j < 6; j++) {
  //        printf("%d ", j); fflush(NULL);
          stVKStiffnessHessian->EvaluateHessianQuadraticForm(
            &Ulin[n3 * (rigid_mode_num_ + i)], &Ulin[n3 * (rigid_mode_num_ + j)], &rhs[ELT(n3, 0, pos)]);

          for (int k = 0; k < n3; k++) //multiply by -1
            rhs[ELT(n3, k, pos)] *= -1.0;
          pos++;
        }
      }

        for (int i = 6; i < numUsedLinearModes; i++) {
        //      printf("%d: ", i); fflush(NULL);
              for (int j = i; j < numUsedLinearModes; j++) {
        //        printf("%d ", j); fflush(NULL);
                stVKStiffnessHessian->EvaluateHessianQuadraticForm(
                  &Ulin[n3 * (rigid_mode_num_ + i)], &Ulin[n3 * (rigid_mode_num_ + j)], &rhs[ELT(n3, 0, pos)]);

                for (int k = 0; k < n3; k++) //multiply by -1
                  rhs[ELT(n3, k, pos)] *= -1.0;

                pos++;
              }
        }
        P(pos);
        //modal_derivative_num_ = pos;
//      printf("\n");
  } else {
    printf("Using the low-memory version.\n");
    stVKStiffnessHessian->EvaluateHessianQuadraticFormDirectAll(Ulin, linear_mode_num_, rhs, rigid_mode_num_, 1);

    if (n3 * modal_derivative_num_ < 0) {
      printf("Error: data too large to be indexed with the word size of your machine.\n");
      *code = 2;
      delete(stVKStiffnessHessian);
      delete(stiffnessMatrix);
      free(rhs);
      delete(precomputedIntegrals);
      delete(stiffnessMatrixClass);
      delete(internalForces);
      return;
    }

    /*
    if ((n3 > 200000) || (precomputationState.numDeriv > 1000))
    {
      printf("Warning: size of data %d might be too large to be indexed with the word size of your machine.\n",n3*precomputationState.numDeriv);
    }
    */

    //multiply by -1
    for (int i = 0; i < n3 * modal_derivative_num_; i++)
      rhs[i] *= -1.0;
  }
  printf("Right-hand sides for modal derivatives computed.\n"); fflush(NULL);

  delete(stVKStiffnessHessian);
  delete(precomputedIntegrals);
  delete(stiffnessMatrixClass);
  delete(internalForces);

  // create mass matrix



  if (rigid_mode_num_ < 6) {
    double * buffer0 = (double*) malloc (sizeof(double) * n3);
    for (int i = 0; i < modal_derivative_num_; i++) {
      mass_matrix_->MultiplyVector(&rhs[ELT(n3, 0, i)], buffer0);
      for (int j = 0; j < rigid_mode_num_; j++) {
        // rhs -= <rhs, rigid mode j> * rigid mode j
        // rigid modes are mass-orthonormal
        double dotp = 0.0;
        for (int k = 0; k < n3; k++)
          dotp += buffer0[k] * Ulin[ELT(n3, k, j)];
        for (int k = 0; k < n3; k++)
          rhs[ELT(n3, k, i)] -= dotp * Ulin[ELT(n3, k, j)];
      }
    }
    free(buffer0);
  } else {
    RemoveSixRigidModes(modal_derivative_num_, rhs);
  }

  sprintf(file_name_, "%s.rhs_original.bin", prefix_);
  std::ofstream out(file_name_,std::ios::binary);
  out.write((char*) &n3, sizeof(int));
  out.write((char*) &modal_derivative_num_, sizeof(int));
  out.write((char*) rhs, sizeof(double) * n3 * modal_derivative_num_);
  out.close();
  numColsOriginalRHS_ = modal_derivative_num_;
  rhsOriginal_.resize(n3,modal_derivative_num_);
  memcpy(rhsOriginal_.data(),rhs,sizeof(double)*n3 * modal_derivative_num_);

  // constrain rhs
  double * rhsConstrained = (double*) malloc (sizeof(double) * numRetainedDOFs * modal_derivative_num_);
  for (int i = 0; i < modal_derivative_num_; i++)
    RemoveRows(n3, &rhsConstrained[numRetainedDOFs * i],
               &rhs[n3 * i], 3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  free(rhs);

  profiler.End("Hessian Pdk Naive");

  // make room for (uninflated) derivatives
  double * modalDerivativesConstrained = (double*) malloc (sizeof(double) * modal_derivative_num_ * numRetainedDOFs);

  // solve K * modesTemp = rhs
  //printf("Factoring the %d x %d stiffness matrix...\n", numRetainedDOFs, numRetainedDOFs);
  //SPOOLESSolver * solver = new SPOOLESSolver(stiffnessMatrix);

  LinearSolver * solver;

  profiler.Start("Solve time Naive");

#define PARDISO_SOLVER_IS_AVAILABLE
#ifdef PARDISO_SOLVER_IS_AVAILABLE
  int positiveDefinite = 0;
  int directIterative = 0;
  int numThreads = 7;
  PardisoSolver * pardisoSolver = new PardisoSolver(stiffnessMatrix, numThreads, positiveDefinite, directIterative);
  pardisoSolver->ComputeCholeskyDecomposition(stiffnessMatrix);
  solver = pardisoSolver;
#elif defined(SPOOLES_SOLVER_IS_AVAILABLE)
  int numThreads = wxThread::GetCPUCount();
  if (numThreads > 1)
    solver = new SPOOLESSolverMT(stiffnessMatrix, numThreads);
  else
    solver = new SPOOLESSolver(stiffnessMatrix);
#else
  solver = new CGSolver(stiffnessMatrix);
#endif

  pardisoSolver->SolveLinearSystemMultipleRHS(modalDerivativesConstrained,rhsConstrained,modal_derivative_num_);

  free(rhsConstrained);
  delete(solver);
  delete(stiffnessMatrix);

  modal_derivatives_.resize(modal_derivative_num_ * n3);

  // insert zero rows into the computed derivatives
  for (int i = 0; i < modal_derivative_num_; i++) {
    InsertRows(n3, &modalDerivativesConstrained[numRetainedDOFs * i],
               &(modal_derivatives_[n3 * i]),
               3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
  }

  free(modalDerivativesConstrained);

   profiler.End("Solve time Naive");
   profiler.End("Modal Ders naive");

  // remove rigid modes from modal derivatives
  if (rigid_mode_num_ > 0)
    printf("Removing rigid modes from modal derivatives...\n"); fflush(NULL);

  if (rigid_mode_num_ < 6) {
    double * buffer1 = (double*) malloc (sizeof(double) * n3);
    for (int i = 0; i < modal_derivative_num_; i++) {
      mass_matrix_->MultiplyVector(&(modal_derivatives_[n3 * i]), buffer1);
      for (int j = 0; j < rigid_mode_num_; j++) {
        // rhs -= <rhs, rigid mode j> * rigid mode j
        // rigid modes are mass-orthonormal
        double dotp = 0.0;
        for (int k = 0; k < n3; k++)
          dotp += buffer1[k] * Ulin[ELT(n3, k, j)];
        for (int k = 0; k < n3; k++)
          modal_derivatives_[ELT(n3, k, i)] -= dotp * Ulin[ELT(n3, k, j)];
      }
    }
    free(buffer1);
  } else {
    RemoveSixRigidModes(modal_derivative_num_, &modal_derivatives_[0]);
  }

  //printf("Mass-normalizing modal derivatives...\n");fflush(NULL);

  // mass-normalize modal derivatives
  //for(int i=0; i < precomputationState.numDeriv; i++)
  //massMatrix->NormalizeVector(&((*modalDerivatives)[n3 * i]));

  //  delete(mass_matrix_);

  free(constrainedDOFs);
}

bool BasisGenerator::ComputeNonLinearModes(int numNonLinearModes) {
  if (modal_derivative_num_ < 0) {
    std::cerr << "BasisGenerator::ComputeNonLinearModes() => modal derivatives are not computed yet" << std::endl;
  }
  non_linear_mode_num_ = numNonLinearModes;
  MassPCA(modal_derivative_num_, numNonLinearModes, non_linear_modes_, kModalAnalysis);
  return true;
}


bool BasisGenerator::MassPCA(int input_basis_num, int output_basis_num, std::vector<double>& basis, int data_origin) {
  int dummy = 0;
  int* code = &dummy;
  int n3 = 3 * mesh_-> getNumVertices();
  int numDataVectors = 0;

  double * dataMatrix;
  if (data_origin == kModalAnalysis) {
    // use linear modes and derivatives
    // construct PCA data matrix:
    // mass-normalize modal derivatives
    double * normalizedModalDerivatives = (double*) malloc (sizeof(double) * n3 * input_basis_num);
    memcpy(normalizedModalDerivatives, &modal_derivatives_[0], sizeof(double) * n3 * input_basis_num);
    for (int i = 0; i < input_basis_num; i++)
      mass_matrix_->NormalizeVector(&(normalizedModalDerivatives[n3 * i]));
    int numUsedLinearModes = linear_mode_num_ - rigid_mode_num_;
    //numDataVectors = numUsedLinearModes + numUsedLinearModes * (numUsedLinearModes + 1) / 2;
    numDataVectors = 321+30;
    printf("Number of PCA datamatrix columns: %d.\n", numDataVectors);
    dataMatrix = (double*) malloc (sizeof(double) * n3 * numDataVectors);

    printf("Generating datamatrix for SVD...\n");

    double lambda0 = frequencies_[6] * frequencies_[6];
    double max = 0.008;
    double min = 0.00008;
    double diff = (max - min)/5;
    //double lambda0 = 100;
    // scale linear modes
    double * Ulin = &linear_modes_[0];
    for (int i = 0; i < numUsedLinearModes; i++) {
      double lambda,factor;
      if(i<6) {
          factor = (max) - (i*diff);
      }
      else {
          lambda = frequencies_[rigid_mode_num_ + i] * frequencies_[rigid_mode_num_ + i];
          factor = lambda0 / lambda;
      }
      //std::cout << "LM i:" << i << " Factor:" << factor << "\n";
      for (int vertex = 0; vertex < n3; vertex++)
        dataMatrix[ELT(n3, vertex, i)] = factor * Ulin[ELT(n3, vertex, rigid_mode_num_ + i)];
    }

    // scale modal derivatives
    int pos = 0;
    max = 0.0008;
    min = 0.000008;
    diff = (max - min)/20;
    int ctr = 0;
    for (int i = 0; i < 6; i++) {
      for (int j = i; j < 6; j++) {
        double factor = 1.0;//(max) - (ctr*diff);
        //std::cout << "MD i:" << i << " j:" << j << " Factor:" << factor << "\n";
        for (int vertex = 0; vertex < n3; vertex++)
          dataMatrix[ELT(n3, vertex, numUsedLinearModes + pos)] =
            factor * normalizedModalDerivatives[ELT(n3, vertex, pos)];
        pos++;
        ctr++;
      }
    }
    for (int i = 6; i < numUsedLinearModes; i++) {
      for (int j = i; j < numUsedLinearModes; j++) {
        double lambdai = frequencies_[rigid_mode_num_ + i] * frequencies_[rigid_mode_num_ + i];
        double lambdaj = frequencies_[rigid_mode_num_ + j] * frequencies_[rigid_mode_num_ + j];
        double factor = lambda0 * lambda0 / (lambdai * lambdaj);
        //std::cout << "MD i:" << i << " j:" << j << " Factor:" << factor << "\n";
        for (int vertex = 0; vertex < n3; vertex++)
          dataMatrix[ELT(n3, vertex, numUsedLinearModes + pos)] =
            factor * normalizedModalDerivatives[ELT(n3, vertex, pos)];
        pos++;
      }
    }
    free(normalizedModalDerivatives);
  } else {
    // data from external simulation
    numDataVectors = input_basis_num;
    dataMatrix = (double*) malloc (sizeof(double) * n3 * numDataVectors);
    memcpy(dataMatrix, &basis[0], sizeof(double) * n3 * numDataVectors);
  }

  profiler.Start("Regular PCA");

  // do lumped-mass-PCA on dataMatrix ( n3 x numDataVectors )

  double * ones = (double*) malloc (sizeof(double) * n3);
  for (int i = 0; i < n3; i++)
    ones[i] = 1.0;

  double * LTDiagonal = (double*) malloc (sizeof(double) * n3);
  mass_matrix_->MultiplyVector(ones, LTDiagonal);
  free(ones);

  // sqrt
  for (int i = 0; i < n3; i++)
    LTDiagonal[i] = sqrt(LTDiagonal[i]);

  // number of retained dimensions can't be more than num linear modes + num derivatives
  if (output_basis_num > numDataVectors)
    output_basis_num = numDataVectors;

  // premultiply by LT
  for (int i = 0; i < n3; i++)
    for (int j = 0; j < numDataVectors; j++)
      dataMatrix[ELT(n3, i, j)] *= LTDiagonal[i];

  // do SVD on dataMatrix ( n3 x numDataVectors ), retain uiState.numComputedNonLinearModes modes
  ThresholdingSpecification thresholdingSpecification;
  thresholdingSpecification.tresholdingType = ThresholdingSpecification::numberOfModesBased;
  thresholdingSpecification.rDesired = output_basis_num;

  int outputr;
  eigen_values_.resize(numDataVectors);
  int matrixPCACode = MatrixPCA(&thresholdingSpecification, n3, numDataVectors, dataMatrix, &outputr, NULL, &eigen_values_[0]);
  if ((matrixPCACode != 0) || (outputr != output_basis_num)) {
    printf("Error performing SVD. Code: %d\n", matrixPCACode);
    *code = matrixPCACode;
    free(dataMatrix);
    free(LTDiagonal);
    ASSERT(false);
    return NULL;
  }
  eigen_values_.resize(output_basis_num);
  for(int i = 0 ; i< output_basis_num;i++) {
      std::cout << "[" << i << "]" << eigen_values_[i] << "\n";
  }
  // solve L^T U = V
  for (int i = 0; i < n3; i++)
    for (int j = 0; j < output_basis_num; j++)
      dataMatrix[ELT(n3, i, j)] /= LTDiagonal[i];

  free(LTDiagonal);

  // export data
  basis.clear();
  basis.insert(basis.end(), dataMatrix, dataMatrix + n3 * output_basis_num);

   profiler.End("Regular PCA");

  *code = 0;
  return true;
}


bool BasisGenerator::MassPCA2(int output_basis_num, Eigen::MatrixXd dataMatrix,std::vector<double> &basis) {

    int n3 = dataMatrix.rows();
    Eigen::MatrixXd dataTranspose = dataMatrix.transpose();
    RedSVD::RedSVD<Eigen::MatrixXd> rrr(dataTranspose,output_basis_num);
    Eigen::MatrixXd pcaData = rrr.matrixV();

    massOrthogonalizationFast(pcaData);

    basis.clear();
    basis.insert(basis.end(), pcaData.data(), pcaData.data() + n3 * output_basis_num);
}

/*
bool BasisGenerator::MassPCA2(int output_basis_num, double* dataMatrix, int numDataVectors, std::vector<double>& basis) {
    int dummy = 0;
    int* code = &dummy;
    int n3 = 3 * mesh_-> getNumVertices();
    // do lumped-mass-PCA on dataMatrix ( n3 x numDataVectors )
    double * ones = (double*) malloc (sizeof(double) * n3);
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for (int i = 0; i < n3; i++)
        ones[i] = 1.0;

    double * LTDiagonal = (double*) malloc (sizeof(double) * n3);
    mass_matrix_->MultiplyVector(ones, LTDiagonal);
    free(ones);

    // sqrt
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for (int i = 0; i < n3; i++)
        LTDiagonal[i] = sqrt(LTDiagonal[i]);

    // number of retained dimensions can't be more than num linear modes + num derivatives
    if (output_basis_num > numDataVectors)
        output_basis_num = numDataVectors;

    // premultiply by LT
   // for (int i = 0; i < n3; i++)
     //   for (int j = 0; j < numDataVectors; j++)
       //     dataMatrix[ELT(n3, i, j)] *= LTDiagonal[i];

    // do SVD on dataMatrix ( n3 x numDataVectors ), retain uiState.numComputedNonLinearModes modes
    ThresholdingSpecification thresholdingSpecification;
    thresholdingSpecification.tresholdingType = ThresholdingSpecification::numberOfModesBased;
    thresholdingSpecification.rDesired = output_basis_num;

    if(0) {
    int outputr;
    stitched_cubature_weights_.resize(numDataVectors);
    //profiler.Start("Act PCA");
    int matrixPCACode = MatrixPCA(&thresholdingSpecification, n3, numDataVectors, dataMatrix, &outputr, NULL, &stitched_cubature_weights_[0]);
    if ((matrixPCACode != 0) || (outputr != output_basis_num)) {
        printf("Error performing SVD. Code: %d\n", matrixPCACode);
        *code = matrixPCACode;
        free(dataMatrix);
        free(LTDiagonal);
        ASSERT(false);
        return NULL;
    }
    }

    Eigen::MatrixXd dataInEig(n3,numDataVectors);
    memcpy(dataInEig.data(),dataMatrix,sizeof(double)*numDataVectors*n3);
    Eigen::MatrixXd dataInEigt = dataInEig.transpose();
    RedSVD::RedSVD<Eigen::MatrixXd> rrr(dataInEigt,output_basis_num);

    //RedSVD::RedPCA<Eigen::MatrixXd> rrr(dataInEigt,output_basis_num);
    Eigen::MatrixXd pcaData = rrr.matrixV();
    Eigen::MatrixXd svals = rrr.singularValues();

    for (int i = 0; i < 40; i++) {
            std::cout << svals(i, 0) << "\n";
        }
    std::cout << pcaData.rows() << "x" << pcaData.cols() << "\n";
    std::ofstream outana("S_RED.txt");
    for(int i = 0; i< output_basis_num;i++ ) {
        outana << svals(i,i) << "\n";
    }
    outana.close();
    for(int i = 0 ; i<n3;i++) {
        for(int j=0; j < output_basis_num;j++) {
            //outana << pcaData(i,j) << " ";
            outana << dataMatrix[ELT(n3,i,j)] << " ";
        }
        outana << "\n";
    }
    outana.close();

    //Eigen::MatrixXd pcaData = Eigen::MatrixXd::Random(n3,output_basis_num);

    //profiler.End("Act PCA");
    //stitched_cubature_weights_.resize(output_basis_num);
    //for(int i = 0 ; i< output_basis_num;i++) {
    //    std::cout << "[" << i << "]" << stitched_cubature_weights_[i] << "\n";
   // }
    // solve L^T U = V
  //  for (int i = 0; i < n3; i++)
   //     for (int j = 0; j < output_basis_num; j++) {
            //dataMatrix[ELT(n3, i, j)] /= LTDiagonal[i];
      //      pcaData(i, j) /= LTDiagonal[i];
     //   }

    free(LTDiagonal);

    Eigen::MatrixXf rm =Eigen::MatrixXf::Random(1000,40);
    RedSVD::RedPCA<Eigen::MatrixXf> rrr(rm,40);
    Eigen::MatrixXf sc = rrr.scores();
    for(int i = 0; i<10;i++) {
        std::cout << "REDPCA:" << sc(i,0) << "\n";
    }


    Eigen::MatrixXd massOrthoNonLin;
    massOrthogonalization(pcaData,mass_matrix_eig_,massOrthoNonLin);

    // export data
    basis.clear();
    //basis.insert(basis.end(), dataMatrix, dataMatrix + n3 * output_basis_num);
    //basis.insert(basis.end(), pcaData.data(), pcaData.data() + n3 * output_basis_num);
    basis.insert(basis.end(), massOrthoNonLin.data(), massOrthoNonLin.data() + n3 * output_basis_num);

    *code = 0;
    return true;
}
*/

void BasisGenerator::ComputeAllModes(int linear_mode_num, int non_linear_mode_num) {
  //L("Compute linear modes.");
   profiler.Start("Basic Lin Modes");
  ComputeLinearModes(linear_mode_num);
   profiler.End("Basic Lin Modes");
   ComputeModalDerivatives();
  ComputeNonLinearModes(non_linear_mode_num);
}

void BasisGenerator::SaveLinearModes(const char *file_name_prefix) {
  char file_name[1024];
  sprintf(file_name, "%s.Ulin", file_name_prefix);
  WriteBasisInBinary(file_name, mesh_->getNumVertices(), linear_mode_num_, &linear_modes_[0]);
}

bool BasisGenerator::ComputeLinearModes(int num_of_desired_modes) {

  profiler.Start("Linear Modes");
  linear_mode_num_ = -1;
  // create stiffness matrix
  StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(mesh_);
  StVKInternalForces * internalForces = new StVKInternalForces(mesh_, precomputedIntegrals);

  SparseMatrix * stiffnessMatrix;
  StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
  stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
  double * zero = (double*) calloc(3 * mesh_->getNumVertices(), sizeof(double));
  int numRows = mesh_->getNumVertices() * 3;
  stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);

  free(zero);
  delete(precomputedIntegrals);
  delete(stiffnessMatrixClass);
  delete(internalForces);

  // constrain the degrees of freedom
  int numConstrainedVertices = (int) (fixed_vertices_.size());
  int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
  set<int> :: iterator iter;
  int i = 0;
  for (iter = fixed_vertices_.begin(); iter != fixed_vertices_.end(); iter++) {
    constrainedDOFs[3 * i + 0] = 3 * (*iter) + 1;
    constrainedDOFs[3 * i + 1] = 3 * (*iter) + 2;
    constrainedDOFs[3 * i + 2] = 3 * (*iter) + 3;
    i++;
  }

  int oneIndexed = 1;
  SparseMatrix* tmp_mass_matrix = new SparseMatrix(*mass_matrix_);
  tmp_mass_matrix->RemoveRowsColumns(3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  stiffnessMatrix->RemoveRowsColumns(3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
 int numRetainedDOFs = stiffnessMatrix->Getn();
 int numLinearSolverThreads = 3;
 PerformanceCounter ARPACKCounter;
  //if(0) {
  // call ARPACK
  double * frequenciesTemp = (double*) malloc (sizeof(double) * num_of_desired_modes);
  double * modesTemp = (double*) malloc (sizeof(double) * num_of_desired_modes * numRetainedDOFs);
  printf("Computing linear modes using ARPACK: ...\n");

  double sigma = -1.0;
  ARPACKSolver generalizedEigenvalueProblem;
  int nconv = generalizedEigenvalueProblem.SolveGenEigShInv(stiffnessMatrix, tmp_mass_matrix,
                                                            num_of_desired_modes, frequenciesTemp,
                                                            modesTemp, sigma, numLinearSolverThreads);
 //}

  // THis is portion I am adding
  /*double sigma = -1.0;
  stiffnessMatrix->AddDiagonalMatrix(-sigma);
  PardisoSolver* pardisoSolver = new PardisoSolver(stiffnessMatrix,numLinearSolverThreads);
  pardisoSolver->ComputeCholeskyDecomposition(stiffnessMatrix);

  //Now create class invK
  InvKSolver *invKSolver = new InvKSolver(pardisoSolver);

  //Generate the Eigen Vectors and the Eigen Values
  double *frequenciesTemp = new double[num_of_desired_modes];
  double *modesTemp = new double[num_of_desired_modes*numRetainedDOFs];
  ARSymStdEig<double, InvKSolver> solver2(numRetainedDOFs, num_of_desired_modes,
                                          invKSolver, (void (InvKSolver::*)(double*,double*)) &InvKSolver::ComputeInvK,sigma,"LM",0,1e-10);
  int nconv = solver2.EigenValVectors(modesTemp, frequenciesTemp);*/

  ARPACKCounter.StopCounter();
  double ARPACKTime = ARPACKCounter.GetElapsedTime();
  printf("ARPACK time: %G s.\n", ARPACKTime); fflush(NULL);

  if (nconv < num_of_desired_modes) {
    free(modesTemp);
    free(frequenciesTemp);
    linear_mode_num_ = -3;
    free(constrainedDOFs);
    delete(stiffnessMatrix);
    delete(tmp_mass_matrix);
    std::cerr << "ComputeLinearModes() => failed in computing linear modes.";
    return false;
  }
  int n3 = 3 * mesh_->getNumVertices();
  frequencies_.resize(num_of_desired_modes);
  linear_modes_.resize(num_of_desired_modes * n3);

  for (int i = 0; i < num_of_desired_modes; i++) {
    // insert zero rows into the computed modes
    int oneIndexed = 1;
    InsertRows(n3, &modesTemp[numRetainedDOFs * i], &((linear_modes_)[n3 * i]),
               3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
  }

  //Save the pure eigen values and the pure eigen vectors for the recomputation stuff
   pure_eigen_vectors_ = linear_modes_;
   pure_eigen_values_.resize(num_of_desired_modes);
   memcpy(&pure_eigen_values_[0],frequenciesTemp,sizeof(double)*num_of_desired_modes);
   profiler.End("Linear Modes");

  for (int i = 0; i < num_of_desired_modes; i++) {
    if (frequenciesTemp[i] <= 0)
      (frequencies_)[i] = 0.0;
    else
      (frequencies_)[i] = sqrt((frequenciesTemp)[i]) / (2 * M_PI);
  }

  //Convert this to a  eigen dense matrix to mass orthogonalization
  /*Eigen::MatrixXd linModes(numRows,num_of_desired_modes);
  for(int r = 0; r<numRows; ++r) {
      for(int c = 0; c<num_of_desired_modes; ++c) {
          linModes(r,c) = linear_modes_[c*numRows+r];
      }
  }*/
  //Eigen::MatrixXd massOrthoLinModes_;
  //massOrthogonalization(linModes,mass_matrix_eig_,massOrthoLinModes_);
  //memcpy(&linear_modes_[0],massOrthoLinModes_.data(),sizeof(double)*(num_of_desired_modes)*numRows);
  //memcpy(&linear_modes_[0],massOrthoLinModes_.data(),sizeof(double)*(num_of_desired_modes)*numRows);

  free(modesTemp);
  free(frequenciesTemp);
  free(constrainedDOFs);

  delete(tmp_mass_matrix);
  delete(stiffnessMatrix);

  linear_mode_num_ = num_of_desired_modes;
  if (fixed_vertices_.size() >= 3) {
    rigid_mode_num_ = 0;
  } else if (fixed_vertices_.size() == 2) {
    rigid_mode_num_ = 1;
  } else if (fixed_vertices_.size() == 1) {
    rigid_mode_num_ = 3;
  } else {
    rigid_mode_num_ = 6;
  }
  //  ComputeDot(mesh_->getNumVertices() * 3, linear_mode_num_, &linear_modes_[0]);
  return true;
}

void BasisGenerator::massOrthogonalizationFast(Eigen::MatrixXd &linBasis) {
    int n3 = 3 * mesh_-> getNumVertices();

    // premultiply by LT
    for (int i = 0; i < n3; i++)
        for (int j = 0; j < linBasis.cols(); j++)
            linBasis(i, j) *= massSqrt_[i];

    //Do Gram Schmidt
    RedSVD::gram_schmidt(linBasis);

    //Perform division
    for (int i = 0; i < n3; i++)
         for (int j = 0; j < linBasis.cols(); j++)
             linBasis(i, j) /= massSqrt_[i];

}

void BasisGenerator::massOrthogonalization(Eigen::MatrixXd linBasis, Eigen::SparseMatrix<double> massMatrix, Eigen::MatrixXd &massOrtho) {
    int numRows = linBasis.rows();
    int numCols = linBasis.cols();

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > choleskySolver(massMatrix);
    Eigen::SparseMatrix<double> LT = choleskySolver.matrixU();
    Eigen::MatrixXd LTU = LT*linBasis;

    //Gram Schmidt Process
    Eigen::MatrixXd ISpecial = Eigen::MatrixXd::Identity(numRows,numCols);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(LTU);
    Eigen::HouseholderSequence<Eigen::MatrixXd,Eigen::VectorXd> Q = qr.householderQ();
    Eigen::MatrixXd orthoBasis = Q*ISpecial;

    massOrtho = Eigen::MatrixXd::Zero(numRows,numCols);
    for(int i=0;i<numCols;i++)
    {
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(LT);
        Eigen::VectorXd VCol = orthoBasis.col(i);
        Eigen::VectorXd UCol = solver.solve(VCol);
        massOrtho.col(i) = UCol;
    }
}

double DiffFromFile(const char *file_name0, const char *file_name1) {
  int m0, n0;
  int m1, n1;
  double* b0, *b1;
  ReadMatrixFromDisk(file_name0, &m0, &n0, &b0);
  ReadMatrixFromDisk(file_name1, &m1, &n1, &b1);
  if (m0 != m1 || n1 != n0) {
    return 1e10;
  }
  P(m0, n0);
  //  ASSERT(m0 == m1 && n1 == n0);
  double diff = 0;
  for (int i = 0; i < m0 * n0; ++i) {
    diff += (b0[i] - b1[i]) * (b0[i] - b1[i]);
  }
  return std::sqrt(diff / (m0 * n0));
}


double ComputeDiff(int n3, int r, double *basis0, double *basis1) {
  double diff = 0;
  for (int i = 0; i < n3 * r; ++i) {
    diff += (basis0[i] - basis1[i]) * (basis0[i] - basis1[i]);
  }
  return std::sqrt(diff / (n3 * r));
}

void BasisGenerator::RegenerateAllModes(int modeId){
    //First lets work with the linear modes
    profiler.Start("Lin Mode Recomputation");
    //profiler.Start("Eigen Val Recomp");
    //profiler.Start("Intro");
    int numVertices = mesh_->getNumVertices();
    int numRows = numVertices * 3;
    int numRetainedDOFs = numRows - numVertsToRemove_;
    int c ;
    int r = 30;
    std::vector<int> nonZeroColumns_;
    std::vector<int> columnInformation_;
    double *eigenVectors = new double[r*numRows];
    double *eigenValues = new double[r];
    //if(modeId==2) {
    nonZeroColumns_ = nonZeroColumns_1_;
    columnInformation_ = columnInformation_1_;
    memcpy(eigenVectors,&pure_eigen_vectors_[0],sizeof(double)*r*numRows);
    memcpy(eigenValues,&pure_eigen_values_[0],sizeof(double)*r);
    //} else if(modeId == 4) {
    //    nonZeroColumns_ = nonZeroColumns_2_;
    //    columnInformation_ = columnInformation_2_;
    //    memcpy(eigenVectors,&stitched_eigen_vectors_1_[0],sizeof(double)*r*numRows);
    //    memcpy(eigenValues,&stitched_eigen_values_1_[0],sizeof(double)*r);
    //}
    c = nonZeroColumns_.size();

    if(stiffyMultiplier_!=1) {
        for(int i=6;i<22;i++) {
            eigenValues[i] *= stiffyMultiplier_;
        }
        eigenValues[25] *= stiffyMultiplier_;
        eigenValues[27] *= stiffyMultiplier_;
        eigenValues[28] *= stiffyMultiplier_;
        eigenValues[29] *= stiffyMultiplier_;
    }
   // profiler.End("Intro");

   // profiler.Start("Product Assembly");
    arma::mat product(numRows,c);
    P(stiffyMultiplier_);
    double stiffness = 5.0e7;// * (double)stiffyMultiplier_;
    double stiffnessSqrt_ = sqrt(stiffness);

    OMP_FOR
            for(int i=0;i<c;i++) {
        int rowVal,colVal;
        colVal = nonZeroColumns_[i];
        rowVal = columnInformation_[i];

        for(int ri=0;ri<numRows;ri++) {
            double valrc1 = 0.0, valrc2 = 0.0;
            double r1,r2;
            for(int ci = 0;ci<r;ci++) {
                valrc1 += (eigenVectors[COLMAJORIDX(ri,ci,numRows)]*eigenVectors[COLMAJORIDX(colVal,ci,numRows)]);
                valrc2 += (eigenVectors[COLMAJORIDX(ri,ci,numRows)]*eigenVectors[COLMAJORIDX(rowVal,ci,numRows)]);
            }
            if(ri==colVal)
                r1 = (1.0-valrc1);//*massInv_[colVal];
            else
                r1 = -valrc1;// * massInv_[colVal];
            if(ri==rowVal)
                r2 = (1.0-valrc2);// * massInv_[rowVal];
            else
                r2 = -valrc2;// * massInv_[rowVal];
            product(ri,i) = (r1-r2)*stiffnessSqrt_;
        }
    }
   // profiler.End("Product Assembly");

    //profiler.Start("QR Decomp");
    arma::mat Q1,R1;
    qr_econ(Q1,R1,product);
    //profiler.End("QR Decomp");

    //profiler.Start("UTARA Assembly");
    arma::mat UTARA(r+c,c);
    UTARA.zeros();
    for(int i=0;i<c;i++) {
        int rowVal,colVal;
        rowVal = nonZeroColumns_[i];
        colVal = columnInformation_[i];
        for(int j=0;j<r;j++) {
            //UTARA(j,i) = ((eigenVectors[COLMAJORIDX(rowVal,j,numRows)]*massInv_[rowVal]) - (eigenVectors[COLMAJORIDX(colVal,j,numRows)]*massInv_[colVal] ))*stiffnessSqrt_;
            UTARA(j,i) = ((eigenVectors[COLMAJORIDX(rowVal,j,numRows)]) - (eigenVectors[COLMAJORIDX(colVal,j,numRows)]))*stiffnessSqrt_;
            //UTARA(j,i) = (tempRes[COLMAJORIDX(rowVal,j,numRows)]-tempRes[COLMAJORIDX(colVal,j,numRows)])*stiffnessSqrt_;
        }
    }

    //delete[] tempRes;

    for(int rows = 0;rows <c;rows++) {
        for(int cols = 0;cols<c;cols++) {
            UTARA(r+rows,cols) = R1(rows,cols);
        }
    }

    arma::mat UTARA_Transpose = UTARA.t();
    arma::mat KPRIME = UTARA * UTARA_Transpose;
    for(int i=0;i<r;i++) {
        KPRIME(i,i) += eigenValues[i];
    }



    //profiler.End("UTARA Assembly");

    //profiler.Start("Eigen Solve");
    arma::vec SPRIME;
    arma::mat SVECTOR;
    arma::eig_sym(SPRIME, SVECTOR, KPRIME);
    //std::cout << "Size of KPRIME:" << KPRIME.n_rows << "x" << KPRIME.n_cols << "\n";
    //profiler.End("Eigen Solve");

    //profiler.Start("Final Assembly");
    arma::mat finalLeft(numRows,r+c);
    finalLeft.zeros();
    memcpy(finalLeft.memptr(),eigenVectors,sizeof(double)*numRows*r);
    int offset = numRows*r;
    memcpy(finalLeft.memptr()+offset,Q1.memptr(),sizeof(double)*numRows*c);
    arma::mat finalEV = finalLeft * SVECTOR;


    Eigen::MatrixXd linBasis(numRows,r);
    memcpy(linBasis.data(),finalEV.memptr(),sizeof(double)*numRows*r);
   // profiler.End("Final Assembly");
    //profiler.End("Eigen Val Recomp");

    //profiler.Start("Mass Ortho");
    //Eigen::Matrix massLin
    massOrthogonalizationFast(linBasis);
    //profiler.End("Mass Ortho");



    //if(modeId == 2) {
    stitched_eigen_vectors_1_.resize(numRows*r);
    stitched_eigen_values_1_.resize(r);
    stitched_frequencies_1_.resize(r);
    //}
    //else {
    //    stitched_eigen_vectors_2_.resize(numRows*r);
    //    stitched_eigen_values_2_.resize(r);
    //   stitched_frequencies_2_.resize(r);
    //}
    stitched_linear_modes_.resize(numRows*r);
    stitched_linear_mode_num_ = r;


    // Then place the modes here also to calculate the actual eigen values
    //if(modeId==2) {
    memcpy(&stitched_eigen_values_1_[0],SPRIME.memptr(),sizeof(double)*stitched_linear_mode_num_);
    memcpy(&stitched_eigen_vectors_1_[0],finalEV.memptr(),sizeof(double)*numRows*r);
    //} else {
    //    memcpy(&stitched_eigen_values_2_[0],SPRIME.memptr(),sizeof(double)*stitched_linear_mode_num_);
    //    memcpy(&stitched_eigen_vectors_2_[0],finalEV.memptr(),sizeof(double)*numRows*r);
    //}
    memcpy(&stitched_linear_modes_[0],linBasis.data(),sizeof(double)*r*numRows);


    /*std::ofstream if_file("Lin_Modes_Recomp_Eigen.txt");
    if_file << numRows << "\n";
    if_file << r << "\n";
    for(int i = 0; i < numRows; i++) {
        for(int j = 0; j < r; j++) {
            if_file << finalEV(i,j) << " ";
        }
        if_file << "\n";
    }
    if_file.close();*/



    if(0) {

        //Insert the linear modes hhere to dump...
        std::ofstream if_file("Lin_Modes_Recomp_Eigen.txt");
        for(int i = 0; i < numRows; i++) {
            for(int j = 0; j < r; j++) {
                //if_file << stitched_linear_modes_[ELT(numRows,i,j)] << " ";

                if_file << finalEV(i,j) << " ";
            }
            if_file << "\n";
        }
        if_file.close();

        SparseMatrix *tempSpMat = new SparseMatrix(*stitched_stiffness_matrix_);
        SparseMatrix *tempMassMat = new SparseMatrix(*mass_matrix_);

        tempSpMat->RemoveRowsColumns(numVertsToRemove_,&vertsToRemove_[0],0);
        tempMassMat->RemoveRowsColumns(numVertsToRemove_,&vertsToRemove_[0],0);
        // call ARPACK
        double * frequenciesTemp = (double*) malloc (sizeof(double) * r);
        double * modesTemp = (double*) malloc (sizeof(double) * r * numRetainedDOFs);

        double sigma = -1.0;
        tempSpMat->AddDiagonalMatrix(-sigma);
        PardisoSolver* pardisoSolver2 = new PardisoSolver(tempSpMat,3);
        pardisoSolver2->ComputeCholeskyDecomposition(tempSpMat);

        //Now create class invK
        InvKSolver *invKSolver2 = new InvKSolver(pardisoSolver2);

        //Generate the Eigen Vectors and the Eigen Values
        ARSymStdEig<double, InvKSolver> solver2(numRetainedDOFs, r,
                                                invKSolver2, (void (InvKSolver::*)(double*,double*)) &InvKSolver::ComputeInvK,sigma,"LM",0,1e-10);
        int nconv = solver2.EigenValVectors(modesTemp, frequenciesTemp);

        double * modes = (double*) malloc (sizeof(double) * r * numRows);

        for (int i = 0; i < r; i++) {
            // insert zero rows into the computed modes
            InsertRows(numRows, &modesTemp[numRetainedDOFs * i], &modes[numRows * i],
                    numVertsToRemove_, &vertsToRemove_[0], 0);
        }
        //Insert the linear modes hhere to dump...
        if_file.open("Lin_Modes_Naive_Eig.txt");
        for(int i = 0; i < numRows; i++) {
            for(int j = 0; j < r; j++) {
                if_file << modes[ELT(numRows,i,j)] << " ";
            }
            if_file << "\n";
        }
        if_file.close();

        free(frequenciesTemp);
        free(modesTemp);
        free(modes);
    }



    //Calculate the frequencies
    for(int i = 0; i<stitched_linear_mode_num_;i++) {
        //if(modeId==2) {
        stitched_frequencies_1_[i] = sqrt(stitched_eigen_values_1_[i]) / (2 * M_PI);
        //}
        //else {
        //   stitched_frequencies_2_[i] = sqrt(stitched_eigen_values_2_[i]) / (2 * M_PI);
        //}
    }

    /*for(int i = 0; i<stitched_linear_mode_num_;i++) {
        std::cout << "TR:" << stitched_frequencies_1_[i] << "\n";
    }*/

    for(int i = 0;i<6;i++) {
        stitched_frequencies_1_[i] = stitched_frequencies_1_[6];//*(i+15);
        frequencies_[i] = frequencies_[6];//*(i+15);
    }




    profiler.End("Lin Mode Recomputation");

    //Now we compute the non-linear modes
    //Get the numRows
    //First create the additional derivatives using product
    profiler.Start("Modal Derivative");
    int numNewModesToUse = 5;
    int rhsColsToUse;
    //numColsOriginalRHS_ = 0;
    if(numColsOriginalRHS_) {
        rhsColsToUse = 0.4 * numColsOriginalRHS_;
        //P(rhsColsToUse);
    }
    else {
        rhsColsToUse = 0;
    }
    //int numNewModalDerivs = 21 + ((numNewModesToUse-6)*((numNewModesToUse-6)+1))/2;
    int numNewModalDerivs = ((numNewModesToUse)*((numNewModesToUse)+1))/2;
    double *finalRHS = new double[numRows*(numNewModalDerivs+rhsColsToUse)];
    //int counter = 0;
    //profiler.Start("Hessian Pdk");

    OMP_FOR
    for(int k=0; k < numNewModalDerivs;k++) {
        int n = numNewModesToUse;
        int i = n - floor(sqrt(-8*k + 4*n*(n+1)-7)/2.0 - 0.5) - 1;
        int j = k + i - n*(n+1)/2 + (n-i)*((n-i)+1)/2;
        stVKStiffnessHessian->EvaluateHessianQuadraticForm(&stitched_eigen_vectors_1_[i*numRows],&stitched_eigen_vectors_1_[j*numRows],&finalRHS[numRows*k]);
    }

  /* OMP_FOR
    for(int i = 0 ; i < numNewModesToUse; i++) {
        OMP_FOR
        for(int j = i; j < numNewModesToUse; j++) {
            int n = numNewModesToUse;
            int linIDX = (n*(n+1)/2) - (n-i)*((n-i)+1)/2 + j - i ;
            //stVKStiffnessHessian->EvaluateHessianQuadraticForm(&stitched_linear_modes_[i*numRows],&stitched_linear_modes_[j*numRows],&finalRHS[numRows*counter]);
            //if(modeId==2) {
            stVKStiffnessHessian->EvaluateHessianQuadraticForm(&stitched_eigen_vectors_1_[i*numRows],&stitched_eigen_vectors_1_[j*numRows],&finalRHS[numRows*linIDX]);
            // } else if(modeId==4) {
            //    stVKStiffnessHessian->EvaluateHessianQuadraticForm(&stitched_eigen_vectors_2_[i*numRows],&stitched_eigen_vectors_2_[j*numRows],&finalRHS[numRows*counter]);
            // }
            //P(counter,linIDX);
           // counter++;
        }
    }*/
    //profiler.End("Hessian Pdk");



   // if(numNewModalDerivs!=counter) {
   //     std::cout << counter << " " << numNewModalDerivs << "\n";
   // }

    OMP_FOR
            for(int k = 0; k < (numRows*numNewModalDerivs); ++k) {
        finalRHS[k] *= -1.0;
    }

    if(rhsColsToUse)
        memcpy(&finalRHS[numRows*numNewModalDerivs],rhsOriginal_.data(),sizeof(double)*numRows*rhsColsToUse);
    //int numRetainedDOFs = numRows - numVertsToRemove_;
    double * rhsConstrained = new double[numRetainedDOFs * (numNewModalDerivs+rhsColsToUse) ];
    //OMP_FOR
    for(int i=0; i<(numNewModalDerivs+rhsColsToUse); i++) {
        RemoveRows(numRows, &rhsConstrained[numRetainedDOFs*i],
                &finalRHS[numRows*i], numVertsToRemove_, &vertsToRemove_[0],0 );
    }

    //std::cout << "[INFO] Starting to decompose and solve the modal derivatives.\n";

    //Atthis point there should be a stiffness matrix
    stitched_stiffness_matrix_->RemoveRowsColumns(numVertsToRemove_,&vertsToRemove_[0],0);

    //profiler.Start("Solve time");

    PardisoSolver* modalDerivativeSolver_ = new PardisoSolver(stitched_stiffness_matrix_,7);
    modalDerivativeSolver_->ComputeCholeskyDecomposition(stitched_stiffness_matrix_);
    //Now solve all the rows and cols using our precomputed cholesky
    double *newModalDerivsConstrained = new double[numRetainedDOFs*(numNewModalDerivs+rhsColsToUse)];
    modalDerivativeSolver_->SolveLinearSystemMultipleRHS(newModalDerivsConstrained,rhsConstrained,numNewModalDerivs+rhsColsToUse);

    //profiler.End("Solve time");

    //Insert zeros in the cols
    double *newModalDerivs = new double[numRows*(numNewModalDerivs+rhsColsToUse)];
    //OMP_FOR
    for(int i=0; i<(numNewModalDerivs+rhsColsToUse); i++) {
        InsertRows(numRows, &newModalDerivsConstrained[numRetainedDOFs*i],
                &newModalDerivs[numRows*i],
                numVertsToRemove_, &vertsToRemove_[0], 0);
    }

    profiler.End("Modal Derivative");

    profiler.Start("Non LIn MOdes Act");
   // profiler.Start("Scaling");
    int numColsInSoup = numNewModalDerivs + rhsColsToUse + stitched_linear_mode_num_;
    Eigen::MatrixXd derivSoup(numRows,numColsInSoup);
    memcpy(derivSoup.data(),newModalDerivs,sizeof(double)*(numNewModalDerivs + rhsColsToUse)*numRows);
    offset = (numNewModalDerivs + rhsColsToUse)*numRows;
    memcpy(derivSoup.data()+offset,&stitched_linear_modes_[0],sizeof(double)*numRows*stitched_linear_mode_num_);



    //Frequency Scaling must happen here
    int counter = 0;


    for(int i = 0 ; i < numNewModesToUse; i++) {
        for(int j = i; j < numNewModesToUse; j++) {
            double lambdai,lambdaj,lambda0;
            if(modeId == 2) {
                lambdai = stitched_frequencies_1_[i] * stitched_frequencies_1_[i];
                lambdaj = stitched_frequencies_1_[j] * stitched_frequencies_1_[j];
                lambda0 = stitched_frequencies_1_[0] * stitched_frequencies_1_[0];
            } else if(modeId==4) {
                lambdai = stitched_frequencies_2_[i] * stitched_frequencies_2_[i];
                lambdaj = stitched_frequencies_2_[j] * stitched_frequencies_2_[j];
                lambda0 = stitched_frequencies_2_[0] * stitched_frequencies_2_[0];
            }
            double factor = lambda0 * lambda0 / (lambdai * lambdaj);
            OMP_FOR
                    for(int row = 0 ; row < numRows; ++row) {
                derivSoup(row,counter) *= factor;
            }
            counter++;
        }
    }
    bool stopLoop = false;


    if(rhsColsToUse) {
        int counter2 = 0;
        for(int i = 0 ; i < linear_mode_num_; i++) {
            for(int j = i; j < linear_mode_num_; j++) {
                double lambdai = frequencies_[i] * frequencies_[i];
                double lambdaj = frequencies_[j] * frequencies_[j];
                double lambda0 = frequencies_[0] * frequencies_[0];
                double factor = lambda0 * lambda0 / (lambdai * lambdaj);
                OMP_FOR
                        for(int row = 0 ; row < numRows; ++row) {
                    derivSoup(row,counter) *= factor;
                }
                counter++;
                counter2++;
                if(counter2==rhsColsToUse) {
                    stopLoop = true;
                    break;
                } //Exit inner loop
            }
            if(stopLoop)
                break;
        }
    }
    for(int i = 0; i < r;++i) {
        double lambdai,lambda0,factor;
        if(modeId == 2) {
            lambdai = stitched_frequencies_1_[i] * stitched_frequencies_1_[i];
            lambda0 = stitched_frequencies_1_[0] * stitched_frequencies_1_[0];
            factor = lambda0 / lambdai ;

        } else if(modeId == 4) {
            lambdai = stitched_frequencies_2_[i] * stitched_frequencies_2_[i];
            lambda0 = stitched_frequencies_2_[0] * stitched_frequencies_2_[0];
            factor = lambda0 / lambdai ;
        }
        OMP_FOR
                for(int row = 0 ; row < numRows; ++row) {
            derivSoup(row,counter) *= factor;
        }
        counter++;
    }

    if(counter!=numColsInSoup) {
        std::cout << "[WARNING] Scaling coverage is not complete " << counter << " " << numColsInSoup << "\n";
    }


   // profiler.End("Scaling");

    // P("Here..");

   // profiler.Start("Rand PCA");
    stitched_non_linear_mode_num_ = 90;
    MassPCA2(stitched_non_linear_mode_num_, derivSoup,stitched_non_linear_modes_);
   // profiler.End("Rand PCA");

    profiler.End("Non LIn MOdes Act");


    //std::cout << "[INFO] Finished computing stitched nonlinear modes.\n";
    delete[] rhsConstrained;
    delete[] finalRHS;
    delete[] newModalDerivsConstrained;
    delete[] newModalDerivs;

}


inline void BasisGenerator::generateColsOfIMUUT(int colIdx, double *x, double *eigs, int numRows, int numCols) {
    for(int r=0;r<numRows;r++) {
        double valrc = 0.0;
        for(int c = 0;c<numCols;c++) {
            valrc += (eigs[COLMAJORIDX(r,c,numRows)]*eigs[COLMAJORIDX(colIdx,c,numRows)]);
        }
        if(r==colIdx)
            x[r] = 1.0-valrc;
        else
            x[r] = -valrc;
    }
}

 void BasisGenerator::generateColsOfIMUUT2(int colIdx1, int colIdx2,double *x, double *eigs, int numRows, int numCols,double multiplier) {
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for(int r=0;r<numRows;r++) {
        double valrc1 = 0.0, valrc2 = 0.0;
        double r1,r2;
        for(int c = 0;c<numCols;c++) {
            valrc1 += (eigs[COLMAJORIDX(r,c,numRows)]*eigs[COLMAJORIDX(colIdx1,c,numRows)]);
            valrc2 += (eigs[COLMAJORIDX(r,c,numRows)]*eigs[COLMAJORIDX(colIdx2,c,numRows)]);
        }
        if(r==colIdx1)
            r1 = (1.0-valrc1);//*massInv_[colIdx1];
        else
            r1 = -valrc1;//*massInv_[colIdx1];
        if(r==colIdx2)
            r2 = (1.0-valrc2);//*massInv_[colIdx2];
        else
            r2 = -valrc2;//*massInv_[colIdx2];
        x[r] = (r1-r2)*multiplier;
    }
}
