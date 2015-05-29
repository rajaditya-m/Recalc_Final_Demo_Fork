#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>


#include "cubature_base/GreedyCubop.h"
#include "cubature_base/VECTOR.h"
#include "cubature_base/VEC3.h"
#include "cubature_base/MATRIX.h"
#include "cubature_base/MATRIX3.h"
#include "tetMesh.h"
#include "StVKIsotropicMaterial.h"
#include "isotropicHyperelasticFEM.h"

class CubatureObjInterface :
  public GreedyCubop {
public:
  CubatureObjInterface(int da, int db);
  ~CubatureObjInterface(void);

  //Set all the initialization data + (constructor) should init everything
  void setBasisMatrix(const char* basisFile_A, const char* basisFile_B);
  void setRotationMatrixData(std::vector<MATRIX3> &rma, std::vector<MATRIX3> &rmb);
  void setGlobalMeshData(int numNode, double* nodes, int numElems, int* elems);
  void setGlobalDeformations(std::vector<VECTOR> &globaldefo);
  void setGlobalInterfaceTets(std::vector<int> &intTet);
  void setGlobalVertexDomainID(std::vector<int> &did);
  void setTransformMaps(std::vector<std::pair<int, int> > &mapA, std::vector<std::pair<int, int> > &mapB);

  //The three functions we will need to actually implmement
  virtual int numTotalPoints();
  virtual void evalPointForceDensity( int pointId, VECTOR& q, VECTOR& gOut, int poseIdx );
  virtual void handleCubature( std::vector<int>& selectedPoints, VECTOR& weights, Real relErr );

  //Additional functions that we will need to implement
  VECTOR getGlobalForceByTetID(int tetId, int poseId);
  VECTOR getsubspaceInternalForce(VECTOR& q, int poseId);
  MATRIX getVertexElementBasis(MATRIX &basis, int vertexId);

private:
  int domainIdx_A;
  int domainIdx_B;

  MATRIX U_A;
  MATRIX U_B;
  int reducedDim;
  int effectiveReducedDim;

  std::vector<MATRIX3> rotMat_A;
  std::vector<MATRIX3> rotMat_B;

  TetMesh *globalMesh;
  IsotropicHyperelasticFEM* globalFEMSolver;
  StVKIsotropicMaterial* globalMaterialDS;
  int numGlobalNodes;
  int numGlobalElems;
  double density;
  double ym;
  double pr;

  std::vector<VECTOR> globalDeformations;

  std::vector<int> globalInterfaceTetId;
  int numElems;
  std::unordered_map<int, int> global2LocalTetrahedrons;

  std::vector<int> globalVertexDomainId;

  std::unordered_map<int, int> global2LocalVertexMap_A;
  std::unordered_map<int, int> global2LocalVertexMap_B;

};

