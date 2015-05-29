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
}

BasisGenerator::BasisGenerator(TetMesh* mesh, double *mass, std::set<int> *fixed_vertex_set) {
  mesh_ = new TetMesh(*mesh);
  Construct(mass, fixed_vertex_set);
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
}

void BasisGenerator::preLoad(const char *basis_prefix){
    //Read the non_linear_modes
    char file_name[512];
    sprintf(file_name, "%s.basis.bin", basis_prefix);
    BinaryFileReader in(file_name);
    int row_num,basis_num_;
    in.Read(&row_num, 1);
    in.Read(&basis_num_, 1);
    non_linear_modes_.resize(row_num*basis_num_);
    in.Read(&non_linear_modes_[0], row_num * basis_num_);
    non_linear_mode_num_ = basis_num_;

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
     stitched_stiffness_matrix_->SaveToMatlabFormat("GT.mm");
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

void BasisGenerator::SetInterfaceVertices(std::vector<std::pair<int, int> > &iv) {
    interface_vertex_ = iv;

    for(auto it = interface_vertex_.begin(); it!= interface_vertex_.end(); it++) {
        int v1 = it->first;
        int v2 = it->second;
        for(int d=0;d<3;d++) {
            nonZeroColumns_.push_back(3*v1+d);
            columnInformation_.push_back(3*v2+d);
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
 // bool computeHessianAtZero = (mesh_->getNumElements() < 5000);
  bool computeHessianAtZero = true; //(mesh_->getNumElements() < 5000);
/*
  if (computeHessianAtZero)
    printf("Hessian at zero will be computed explicitly.\n");
  else
    printf("Hessian at zero will not be computed explicitly.\n");
*/
  stVKStiffnessHessian = new StVKHessianTensor(stiffnessMatrixClass);

  int numUsedLinearModes = linear_mode_num_ - rigid_mode_num_;
  modal_derivative_num_ = numUsedLinearModes * (numUsedLinearModes + 1) / 2;
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

  double * Ulin = &linear_modes_[0];
  if (computeHessianAtZero) {
    printf("Using the high-memory version.\n");
    int pos = 0;
    for (int i = 0; i < numUsedLinearModes; i++) {
//      printf("%d: ", i); fflush(NULL);
      for (int j = i; j < numUsedLinearModes; j++) {
//        printf("%d ", j); fflush(NULL);
        stVKStiffnessHessian->EvaluateHessianQuadraticForm(
          &Ulin[n3 * (rigid_mode_num_ + i)], &Ulin[n3 * (rigid_mode_num_ + j)], &rhs[ELT(n3, 0, pos)]);

        for (int k = 0; k < n3; k++) //multiply by -1
          rhs[ELT(n3, k, pos)] *= -1.0;

        pos++;
      }
//      printf("\n");
    }
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

  sprintf(file_name_, "%s.rhs_original.bin", prefix_);
  std::ofstream out(file_name_,std::ios::binary);
  out.write((char*) &n3, sizeof(int));
  out.write((char*) &modal_derivative_num_, sizeof(int));
  out.write((char*) rhs, sizeof(double) * n3 * modal_derivative_num_);
  out.close();
  numColsOriginalRHS_ = modal_derivative_num_;
  rhsOriginal_.resize(n3,modal_derivative_num_);
  memcpy(rhsOriginal_.data(),rhs,sizeof(double)*n3 * modal_derivative_num_);

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

  // constrain rhs
  double * rhsConstrained = (double*) malloc (sizeof(double) * numRetainedDOFs * modal_derivative_num_);
  for (int i = 0; i < modal_derivative_num_; i++)
    RemoveRows(n3, &rhsConstrained[numRetainedDOFs * i],
               &rhs[n3 * i], 3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  free(rhs);

  // make room for (uninflated) derivatives
  double * modalDerivativesConstrained = (double*) malloc (sizeof(double) * modal_derivative_num_ * numRetainedDOFs);

  // solve K * modesTemp = rhs
  //printf("Factoring the %d x %d stiffness matrix...\n", numRetainedDOFs, numRetainedDOFs);
  //SPOOLESSolver * solver = new SPOOLESSolver(stiffnessMatrix);

  LinearSolver * solver;

#define PARDISO_SOLVER_IS_AVAILABLE
#ifdef PARDISO_SOLVER_IS_AVAILABLE
  int positiveDefinite = 0;
  int directIterative = 0;
  int numThreads = 3;
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

#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < modal_derivative_num_; i++) {
    //printf("Solving for derivative #%d out of %d.\n", i + 1, modal_derivative_num_); fflush(NULL);
    solver->SolveLinearSystem(&modalDerivativesConstrained[ELT(numRetainedDOFs, 0, i)], &rhsConstrained[ELT(numRetainedDOFs, 0, i)]);
  }

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
    numDataVectors = numUsedLinearModes + numUsedLinearModes * (numUsedLinearModes + 1) / 2;
    printf("Number of PCA datamatrix columns: %d.\n", numDataVectors);
    dataMatrix = (double*) malloc (sizeof(double) * n3 * numDataVectors);

    printf("Generating datamatrix for SVD...\n");

    double lambda0 = frequencies_[rigid_mode_num_] * frequencies_[rigid_mode_num_];

    // scale linear modes
    double * Ulin = &linear_modes_[0];
    for (int i = 0; i < numUsedLinearModes; i++) {
      double lambda = frequencies_[rigid_mode_num_ + i] * frequencies_[rigid_mode_num_ + i];
      double factor = lambda0 / lambda;
      for (int vertex = 0; vertex < n3; vertex++)
        dataMatrix[ELT(n3, vertex, i)] = factor * Ulin[ELT(n3, vertex, rigid_mode_num_ + i)];
    }

    // scale modal derivatives
    int pos = 0;
    for (int i = 0; i < numUsedLinearModes; i++) {
      for (int j = i; j < numUsedLinearModes; j++) {
        double lambdai = frequencies_[rigid_mode_num_ + i] * frequencies_[rigid_mode_num_ + i];
        double lambdaj = frequencies_[rigid_mode_num_ + j] * frequencies_[rigid_mode_num_ + j];
        double factor = lambda0 * lambda0 / (lambdai * lambdaj);

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
   profiler.Start("Basic Modal Ders");
   ComputeModalDerivatives();
  profiler.End("Basic Modal Ders");
  profiler.Start("Basic NonLin Modes");
  ComputeNonLinearModes(non_linear_mode_num);
  profiler.End("Basic NonLin Modes");
}

void BasisGenerator::SaveLinearModes(const char *file_name_prefix) {
  char file_name[1024];
  sprintf(file_name, "%s.Ulin", file_name_prefix);
  WriteBasisInBinary(file_name, mesh_->getNumVertices(), linear_mode_num_, &linear_modes_[0]);
}

bool BasisGenerator::ComputeLinearModes(int num_of_desired_modes) {
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
  if(0) {
  // call ARPACK
  double * frequenciesTemp = (double*) malloc (sizeof(double) * num_of_desired_modes);
  double * modesTemp = (double*) malloc (sizeof(double) * num_of_desired_modes * numRetainedDOFs);
  printf("Computing linear modes using ARPACK: ...\n");

  double sigma = -1.0;
  ARPACKSolver generalizedEigenvalueProblem;
  int nconv = generalizedEigenvalueProblem.SolveGenEigShInv(stiffnessMatrix, tmp_mass_matrix,
                                                            num_of_desired_modes, frequenciesTemp,
                                                            modesTemp, sigma, numLinearSolverThreads);
  }

  // THis is portion I am adding
  double sigma = -1.0;
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
  int nconv = solver2.EigenValVectors(modesTemp, frequenciesTemp);

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

  for (int i = 0; i < num_of_desired_modes; i++) {
    if (frequenciesTemp[i] <= 0)
      (frequencies_)[i] = 0.0;
    else
      (frequencies_)[i] = sqrt((frequenciesTemp)[i]) / (2 * M_PI);
  }

  //Convert this to a  eigen dense matrix to mass orthogonalization
  Eigen::MatrixXd linModes(numRows,num_of_desired_modes);
  for(int r = 0; r<numRows; ++r) {
      for(int c = 0; c<num_of_desired_modes; ++c) {
          linModes(r,c) = linear_modes_[c*numRows+r];
      }
  }
  Eigen::MatrixXd massOrthoLinModes_;
  massOrthogonalization(linModes,mass_matrix_eig_,massOrthoLinModes_);
  memcpy(&linear_modes_[0],massOrthoLinModes_.data(),sizeof(double)*(num_of_desired_modes)*numRows);

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
    int n3chk = linBasis.rows();
    ASSERT(n3==n3chk);
    double * ones = (double*) malloc (sizeof(double) * n3);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n3; i++)
        ones[i] = 1.0;

    double * LTDiagonal = (double*) malloc (sizeof(double) * n3);
    mass_matrix_->MultiplyVector(ones, LTDiagonal);
    free(ones);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n3; i++)
        LTDiagonal[i] = sqrt(LTDiagonal[i]);

    // premultiply by LT
    for (int i = 0; i < n3; i++)
        for (int j = 0; j < linBasis.cols(); j++)
            linBasis(i, j) *= LTDiagonal[i];

    //Do Gram Schmidt
    RedSVD::gram_schmidt(linBasis);

    //Perform division
    for (int i = 0; i < n3; i++)
         for (int j = 0; j < linBasis.cols(); j++)
             linBasis(i, j) /= LTDiagonal[i];

    free(LTDiagonal);

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
    int numVertices = mesh_->getNumVertices();
    int numRows = numVertices * 3;

    int c = nonZeroColumns_.size();
    //@TODO::This value is hardcoded for now, later on we will modify this
    int r = 30;

    double *eigenVectors = new double[r*numRows];
    if(modeId == 2) {
        memcpy(eigenVectors,&pure_eigen_vectors_[0],sizeof(double)*r*numRows);
    } else if(modeId == 4) {
        memcpy(eigenVectors,&stitched_eigen_vectors_1_[0],sizeof(double)*r*numRows);
    }
    double *eigenValues = new double[r];
    if(modeId == 2) {
        memcpy(eigenValues,&pure_eigen_values_[0],sizeof(double)*r);
    }
    else if(modeId == 4) {
        memcpy(eigenValues,&stitched_eigen_values_1_[0],sizeof(double)*r);
    }

    profiler.Start("Product Assembly");
    arma::mat product(numRows,c);
    double stiffness = 50000.0;
    double stiffnessSqrt_ = sqrt(stiffness);
    for(int i=0;i<c;i++) {
        int colVal = nonZeroColumns_[i];
        int rowVal = columnInformation_[i];

        double *x = new double[numRows];
        double *y = new double[numRows];
        profiler.Start("Col Gen");
        generateColsOfIMUUT(colVal,x,eigenVectors,numRows,r);
        generateColsOfIMUUT(rowVal,y,eigenVectors,numRows,r);
        profiler.End("Col Gen");

        profiler.Start("Pdk ass");
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
        for(int r_=0;r_<numRows;r_++) {
            product(r_,i) = (x[r_]-y[r_])*stiffnessSqrt_;
        }
        profiler.End("Pdk ass");
    }
     profiler.End("Product Assembly");

     //profiler.Start("QR Decomp");
    arma::mat Q1,R1;
    qr_econ(Q1,R1,product);
    //profiler.End("QR Decomp");

    //profiler.Start("UTARA Assembly");
    arma::mat UTARA(r+c,c);
    UTARA.zeros();
    for(int i=0;i<c;i++) {
        int rowVal = nonZeroColumns_[i];
        int colVal = columnInformation_[i];
        for(int j=0;j<r;j++) {
            UTARA(j,i) = eigenVectors[COLMAJORIDX(rowVal,j,numRows)]*stiffnessSqrt_;
        }
        for(int j=0;j<r;j++) {
            UTARA(j,i) += (-eigenVectors[COLMAJORIDX(colVal,j,numRows)]*stiffnessSqrt_);
        }
    }

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

   // profiler.Start("Eigen Solve");
    arma::vec SPRIME;
    arma::mat SVECTOR;
    arma::eig_sym(SPRIME, SVECTOR, KPRIME);
    //std::cout << "Size of KPRIME:" << KPRIME.n_rows << "x" << KPRIME.n_cols << "\n";
  //  profiler.End("Eigen Solve");

   // profiler.Start("Final Assembly");
    arma::mat finalLeft(numRows,r+c);
    finalLeft.zeros();
    memcpy(finalLeft.memptr(),eigenVectors,sizeof(double)*numRows*r);
    int offset = numRows*r;
    memcpy(finalLeft.memptr()+offset,Q1.memptr(),sizeof(double)*numRows*c);
    arma::mat finalEV = finalLeft * SVECTOR;
    Eigen::MatrixXd linBasis(numRows,r);
    memcpy(linBasis.data(),finalEV.memptr(),sizeof(double)*numRows*r);
   // profiler.End("Final Assembly");

  //  profiler.Start("Mass Ortho");
    massOrthogonalizationFast(linBasis);
   // profiler.End("Mass Ortho");

    if(modeId == 2) {
        stitched_eigen_vectors_1_.resize(numRows*r);
        stitched_eigen_values_1_.resize(r);
        stitched_frequencies_1_.resize(r);
    }
    else if(modeId == 4) {
        stitched_eigen_vectors_2_.resize(numRows*r);
        stitched_eigen_values_2_.resize(r);
        stitched_frequencies_2_.resize(r);
    }

    stitched_linear_mode_num_ = r;
    stitched_linear_modes_.resize(numRows*r);

    //if(modeId==2) {
        memcpy(&stitched_eigen_values_1_[0],SPRIME.memptr(),sizeof(double)*stitched_linear_mode_num_);
        //memcpy(&stitched_eigen_vectors_1_[0],finalEV.memptr(),sizeof(double)*numRows*r);
    //} else {
        //memcpy(&stitched_eigen_values_2_[0],SPRIME.memptr(),sizeof(double)*stitched_linear_mode_num_);
        //memcpy(&stitched_eigen_vectors_2_[0],finalEV.memptr(),sizeof(double)*numRows*r);
    //}
    //memcpy(&stitched_eigen_values_1_[0],SPRIME.memptr(),sizeof(double)*stitched_linear_mode_num_);
    memcpy(&stitched_linear_modes_[0],linBasis.data(),sizeof(double)*r*numRows);

    //Calculate the frequencies
    for(int i = 0; i<stitched_linear_mode_num_;i++) {
        stitched_frequencies_1_[i] = sqrt(stitched_eigen_values_1_[i]) / (2 * M_PI);
    }

    //std::cout << "[INFO] Linear Modes computation finished.\n";
    profiler.End("Lin Mode Recomputation");

    //Now we compute the non-linear modes
    //Get the numRows
    //First create the additional derivatives using product
    int numNewModesToUse = 5;
    int rhsColsToUse;
    if(numColsOriginalRHS_) {
        rhsColsToUse = 0.4 * numColsOriginalRHS_;
    }
    else {
        rhsColsToUse = 0;
    }
    int numNewModalDerivs = (numNewModesToUse*(numNewModesToUse+1))/2;
    double *finalRHS = new double[numRows*(numNewModalDerivs+rhsColsToUse)];
    int counter = 0;
    for(int i = 0 ; i < numNewModesToUse; i++) {
        for(int j = i; j < numNewModesToUse; j++) {
            //stVKStiffnessHessian->EvaluateHessianQuadraticForm(&stitched_linear_modes_[i*numRows],&stitched_linear_modes_[j*numRows],&finalRHS[numRows*counter]);
            stVKStiffnessHessian->EvaluateHessianQuadraticForm(&stitched_eigen_vectors_1_[i*numRows],&stitched_eigen_vectors_1_[j*numRows],&finalRHS[numRows*counter]);
            counter++;
        }
    }



    if(numNewModalDerivs!=counter) {
        std::cout << counter << " " << numNewModalDerivs << "\n";
    }
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for(int k = 0; k < (numRows*counter); ++k) {
        finalRHS[k] *= -1.0;
    }

    if(rhsColsToUse)
        memcpy(&finalRHS[numRows*counter],rhsOriginal_.data(),sizeof(double)*numRows*rhsColsToUse);
    int numRetainedDOFs = numRows - numVertsToRemove_;
    double * rhsConstrained = new double[numRetainedDOFs * (numNewModalDerivs+rhsColsToUse) ];
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for(int i=0; i<(numNewModalDerivs+rhsColsToUse); i++) {
        RemoveRows(numRows, &rhsConstrained[numRetainedDOFs*i],
                &finalRHS[numRows*i], numVertsToRemove_, &vertsToRemove_[0],0 );
    }

    //std::cout << "[INFO] Starting to decompose and solve the modal derivatives.\n";

    //Atthis point there should be a stiffness matrix
    stitched_stiffness_matrix_->RemoveRowsColumns(numVertsToRemove_,&vertsToRemove_[0],0);

    PardisoSolver* modalDerivativeSolver_ = new PardisoSolver(stitched_stiffness_matrix_,7);
    modalDerivativeSolver_->ComputeCholeskyDecomposition(stitched_stiffness_matrix_);

    //Now solve all the rows and cols using our precomputed cholesky
    double *newModalDerivsConstrained = new double[numRetainedDOFs*(numNewModalDerivs+rhsColsToUse)];
    modalDerivativeSolver_->SolveLinearSystemMultipleRHS(newModalDerivsConstrained,rhsConstrained,numNewModalDerivs+rhsColsToUse);

    //Insert zeros in the cols
    double *newModalDerivs = new double[numRows*(numNewModalDerivs+rhsColsToUse)];
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for(int i=0; i<(numNewModalDerivs+rhsColsToUse); i++) {
        InsertRows(numRows, &newModalDerivsConstrained[numRetainedDOFs*i],
                &newModalDerivs[numRows*i],
                numVertsToRemove_, &vertsToRemove_[0], 0);
    }

    int numColsInSoup = numNewModalDerivs + rhsColsToUse + stitched_linear_mode_num_;
    Eigen::MatrixXd derivSoup(numRows,numColsInSoup);
    memcpy(derivSoup.data(),newModalDerivs,sizeof(double)*(numNewModalDerivs + rhsColsToUse)*numRows);
    offset = (numNewModalDerivs + rhsColsToUse)*numRows;
    memcpy(derivSoup.data()+offset,&stitched_linear_modes_[0],sizeof(double)*numRows*stitched_linear_mode_num_);

    //Frequency Scaling must happen here
    counter = 0;
    for(int i = 0 ; i < numNewModesToUse; i++) {
        for(int j = i; j < numNewModesToUse; j++) {
            double lambdai = stitched_frequencies_1_[i] * stitched_frequencies_1_[i];
            double lambdaj = stitched_frequencies_1_[j] * stitched_frequencies_1_[j];
            double lambda0 = stitched_frequencies_1_[0] * stitched_frequencies_1_[0];
            double factor = lambda0 * lambda0 / (lambdai * lambdaj);
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
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
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
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
        double lambdai = stitched_frequencies_1_[i] * stitched_frequencies_1_[i];
        double lambda0 = stitched_frequencies_1_[0] * stitched_frequencies_1_[0];
        double factor = lambda0 / lambdai ;
#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
        for(int row = 0 ; row < numRows; ++row) {
            derivSoup(row,counter) *= factor;
        }
        counter++;
    }

    if(counter!=numColsInSoup) {
        std::cout << "[WARNING] Scaling coverage is not complete " << counter << " " << numColsInSoup << "\n";
    }

    stitched_non_linear_mode_num_ = 90;
    MassPCA2(stitched_non_linear_mode_num_, derivSoup,stitched_non_linear_modes_);

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
