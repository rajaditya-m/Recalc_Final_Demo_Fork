#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "cubature_base/GreedyCubop.h"
#include "tetMesh.h"
#include "StVKIsotropicMaterial.h"
#include "isotropicHyperelasticFEM.h"
#include "cubature_base/MATRIX3.h"
#include "cubature_base/MATRIX.h"

//This is my concreted implemantion of the Greedy CubeOP

class CubatureObjDomain : public GreedyCubop {
public:
  //This needs a node and ele file as input
  CubatureObjDomain(const char* fileNameNode, const char* fileNameEle, const char* fileNameBasis);
  ~CubatureObjDomain(void);

  //All the pure virtual functions in the GreedyCubop must be implemented
  int numTotalPoints();
  void evalPointForceDensity( int pointId, VECTOR& q, VECTOR& gOut, int poseIdx );
  virtual void handleCubature( std::vector<int>& selectedPoints, VECTOR& weights, Real relErr );

  //These are some additional helper setters
  void setDomainId(int i)							{ domainId = i;	}

  //These are additional helper API's that we will need
  MATRIX getElementBasisMatrix(int pointId);
  VECTOR getsubspaceInternalForce(VECTOR& q, int poseId);
  //VECTOR getsubspaceInternalForceByCubature(VECTOR &q);
  //MATRIX getsubspaceStiffnessMatrixByCubature(VECTOR &q);
  void readCubatureFromFile(const char* file);

private:
  double* nodes;
  int* elems;
  int numElems;
  int numNodes;

  //VEGA STUFF That we need
  TetMesh* vegaMesh;
  IsotropicHyperelasticFEM* femSolver;
  StVKIsotropicMaterial* materialDS;

  double density;
  double ym;
  double pr;

  //The raison-de-entre the basis transformation matrix
  MATRIX U;
  int rsize;

  int domainId;

  //This contains the cubature values trained
  bool hasCubature;
  std::vector<int> cubaturePoints;
  VECTOR cubatureWeights;
};

