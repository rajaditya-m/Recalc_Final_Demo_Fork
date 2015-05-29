#include "CubatureObjInterface.h"
#include "global.h"

static const char kFolder[] = DATA_DIRECTORY "/armadillo";

CubatureObjInterface::CubatureObjInterface(int da, int db) {
  domainIdx_A = da;
  domainIdx_B = db;

  density = 1.1e3;
  ym = 1.0e6;
  pr = 0.49;
}

CubatureObjInterface::~CubatureObjInterface(void) {
}

void CubatureObjInterface::setBasisMatrix(const char* basisFile_A, const char* basisFile_B) {
  std::string line;
  std::ifstream basisFile(basisFile_A);
  std::getline(basisFile, line);
  std::istringstream iss(line);
  int numRows, numCols;
  iss >> numRows ;
  std::getline(basisFile, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> numCols;
  U_A = MATRIX(numRows, numCols);
  for (int r = 0; r < numRows; r++) {
    std::getline(basisFile, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    for (int c = 0; c < numCols; c++) {
      double val;
      iss >> val;
      U_A(r, c) = val;
    }
  }
  basisFile.close();

  basisFile.open(basisFile_B);
  std::getline(basisFile, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> numRows ;
  std::getline(basisFile, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> numCols;
  U_B = MATRIX(numRows, numCols);
  for (int r = 0; r < numRows; r++) {
    std::getline(basisFile, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    for (int c = 0; c < numCols; c++) {
      double val;
      iss >> val;
      U_B(r, c) = val;
    }
  }
  basisFile.close();

  reducedDim = U_A.cols();
  effectiveReducedDim = 2 * reducedDim;
}

void CubatureObjInterface::setRotationMatrixData(std::vector<MATRIX3> &rma, std::vector<MATRIX3> &rmb) {
  rotMat_A = std::vector<MATRIX3>(rma);
  rotMat_B = std::vector<MATRIX3>(rmb);
}

void CubatureObjInterface::setGlobalMeshData(int numNode, double* nodes, int numElems, int* elems) {
  numGlobalNodes = numNode;
  numGlobalElems = numElems;

  globalMesh = new TetMesh(numGlobalNodes, nodes, numGlobalElems, elems, ym, pr, density);
  globalMaterialDS = new StVKIsotropicMaterial(globalMesh, 1, 500);
  globalFEMSolver = new IsotropicHyperelasticFEM(globalMesh, globalMaterialDS);
}

void CubatureObjInterface::setGlobalDeformations(std::vector<VECTOR> &globaldefo) {
  globalDeformations = std::vector<VECTOR>(globaldefo);
}

void CubatureObjInterface::setGlobalInterfaceTets(std::vector<int> &intTet) {
  numElems = int(intTet.size());

  std::cout << "Number of elems:" <<  numElems << "\n";

  globalInterfaceTetId = std::vector<int>(intTet);

  int counter = 0;
  while (counter < numElems) {
    global2LocalTetrahedrons.insert(std::pair<int, int>(globalInterfaceTetId[counter], counter));
    counter++;
  }

}

void CubatureObjInterface::setGlobalVertexDomainID(std::vector<int> &did) {
  globalVertexDomainId = std::vector<int>(did);
}

void CubatureObjInterface::setTransformMaps(std::vector<std::pair<int, int> > &mapA, std::vector<std::pair<int, int> > &mapB) {
  for (int i = 0; i < int(mapA.size()); i++) {
    global2LocalVertexMap_A.insert(std::pair<int, int>(mapA[i].second, mapA[i].first));
  }
  for (int i = 0; i < int(mapB.size()); i++) {
    global2LocalVertexMap_B.insert(std::pair<int, int>(mapB[i].second, mapB[i].first));
  }
}

int CubatureObjInterface::numTotalPoints() {
  return numElems;
}

void CubatureObjInterface::evalPointForceDensity( int pointId, VECTOR& q, VECTOR& gOut, int poseIdx ) {
  (void) q;
  int globalTetIdx = globalInterfaceTetId[pointId];

  VECTOR forceInGlobalCoords = getGlobalForceByTetID(globalTetIdx, poseIdx);
  VECTOR localDomainForceA(reducedDim);
  VECTOR localDomainForceB(reducedDim);
  localDomainForceA.clear();
  localDomainForceB.clear();

  for (int i = 0; i < 4; i++) {
    int globalVertexId = globalMesh->getVertexIndex(globalTetIdx, i);
    int domainId = globalVertexDomainId[globalVertexId];
    VEC3 vertexForceInGlobalCoords(forceInGlobalCoords(i * 3 + 0), forceInGlobalCoords(i * 3 + 1), forceInGlobalCoords(i * 3 + 2));
    if (domainId == domainIdx_A) {
      auto it = global2LocalVertexMap_A.find(globalVertexId);
      int localVertexId = it->second;
      MATRIX3 rotMat = rotMat_A[poseIdx];
      MATRIX3 rotMatTrans = rotMat.transpose();
      VEC3 tmpRes = rotMatTrans * vertexForceInGlobalCoords;
      MATRIX UAElem = getVertexElementBasis(U_A, localVertexId);
      MATRIX UAElemTranspose = UAElem.transpose();
      VECTOR fullSpaceLocalForce(3);
      fullSpaceLocalForce(0) = tmpRes[0];
      fullSpaceLocalForce(1) = tmpRes[1];
      fullSpaceLocalForce(2) = tmpRes[2];
      VECTOR finalVertexForce = UAElemTranspose * fullSpaceLocalForce;
      localDomainForceA += finalVertexForce;
    } else {
      auto it = global2LocalVertexMap_B.find(globalVertexId);
      int localVertexId = it->second;
      MATRIX3 rotMat = rotMat_B[poseIdx];
      MATRIX3 rotMatTrans = rotMat.transpose();
      VEC3 tmpRes = rotMatTrans * vertexForceInGlobalCoords;
      MATRIX UBElem = getVertexElementBasis(U_B, localVertexId);
      MATRIX UBElemTranspose = UBElem.transpose();
      VECTOR fullSpaceLocalForce(3);
      fullSpaceLocalForce(0) = tmpRes[0];
      fullSpaceLocalForce(1) = tmpRes[1];
      fullSpaceLocalForce(2) = tmpRes[2];
      VECTOR finalVertexForce = UBElemTranspose * fullSpaceLocalForce;
      localDomainForceB += finalVertexForce;
    }
  }

  VECTOR qRes(effectiveReducedDim);
  for (int i = 0; i < reducedDim; i++) {
    qRes(i) = localDomainForceA(i);
    qRes(i + reducedDim) = localDomainForceB(i);
  }

  gOut = qRes;

}

void CubatureObjInterface::handleCubature(std::vector<int>& selectedPoints, VECTOR& weights, Real relErr ) {
  cout << "[INFO] n = " << selectedPoints.size() << " \t relerr = " << (relErr * 100) << "%" << endl;
  char fileName[500];
  sprintf(fileName, "%s/Cubatures/interface_%d_%d.txt", kFolder, domainIdx_A, domainIdx_B);
  std::ofstream fileOut(fileName);
  fileOut << selectedPoints.size() << "\n";
  for (int i = 0; i < int(selectedPoints.size()); i++) {
    fileOut << selectedPoints[i] << " " << weights(i) << "\n";
  }
  fileOut.close();
}

VECTOR CubatureObjInterface::getGlobalForceByTetID(int tetId, int poseId) {
  double* u = globalDeformations[poseId].data();
  double* intForces = new double[globalMesh->getNumVertices() * 3];
  int computationMode = IsotropicHyperelasticFEM::COMPUTE_INTERNALFORCES;
  globalFEMSolver->GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue(u, NULL, intForces, NULL, computationMode);
  int retVal = globalFEMSolver->GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse(tetId, tetId + 1, u, NULL, intForces, NULL, computationMode);
  if (retVal) {
    std::cout << "[INFO] Internal Force computation failed.\n";
  }
  VECTOR intForceVec(12);
  int counter = 0;
  for (int i = 0; i < 4; i++) {
    for (int d = 0; d < 3; d++) {
      intForceVec(counter) = -intForces[globalMesh->getVertexIndex(tetId, i) * 3 + d]; //The negative sign is VV IMP
      counter++;
    }
  }
  delete[] intForces;
  return intForceVec;
}

VECTOR CubatureObjInterface::getsubspaceInternalForce(VECTOR& unused, int poseId) {
  (void) unused;
  VECTOR domainAForces(reducedDim);
  VECTOR domainBForces(reducedDim);

  domainAForces.clear();
  domainBForces.clear();

  for (int t = 0; t < int(globalInterfaceTetId.size()); t++) {
    int globalTetIdx = globalInterfaceTetId[t];
    VECTOR forceInGlobalCoords = getGlobalForceByTetID(globalTetIdx, poseId);
    for (int i = 0; i < 4; i++) {
      int globalVertexId = globalMesh->getVertexIndex(globalTetIdx, i);
      int domainId = globalVertexDomainId[globalVertexId];
      VEC3 vertexForceInGlobalCoords(forceInGlobalCoords(i * 3 + 0), forceInGlobalCoords(i * 3 + 1), forceInGlobalCoords(i * 3 + 2));
      if (domainId == domainIdx_A) {
        auto it = global2LocalVertexMap_A.find(globalVertexId);
        int localVertexId = it->second;
        MATRIX3 rotMat = rotMat_A[poseId];
        MATRIX3 rotMatTrans = rotMat.transpose();
        VEC3 tmpRes = rotMatTrans * vertexForceInGlobalCoords;
        MATRIX UAElem = getVertexElementBasis(U_A, localVertexId);
        MATRIX UAElemTranspose = UAElem.transpose();
        VECTOR fullSpaceLocalForce(3);
        fullSpaceLocalForce(0) = tmpRes[0];
        fullSpaceLocalForce(1) = tmpRes[1];
        fullSpaceLocalForce(2) = tmpRes[2];
        VECTOR finalVertexForce = UAElemTranspose * fullSpaceLocalForce;
        domainAForces += finalVertexForce;
      } else {
        auto it = global2LocalVertexMap_B.find(globalVertexId);
        int localVertexId = it->second;
        MATRIX3 rotMat = rotMat_B[poseId];
        MATRIX3 rotMatTrans = rotMat.transpose();
        VEC3 tmpRes = rotMatTrans * vertexForceInGlobalCoords;
        MATRIX UBElem = getVertexElementBasis(U_B, localVertexId);
        MATRIX UBElemTranspose = UBElem.transpose();
        VECTOR fullSpaceLocalForce(3);
        fullSpaceLocalForce(0) = tmpRes[0];
        fullSpaceLocalForce(1) = tmpRes[1];
        fullSpaceLocalForce(2) = tmpRes[2];
        VECTOR finalVertexForce = UBElemTranspose * fullSpaceLocalForce;
        domainBForces += finalVertexForce;
      }
    }
  }

  VECTOR qRes(effectiveReducedDim);
  for (int i = 0; i < reducedDim; i++) {
    qRes(i) = domainAForces(i);
    qRes(i + reducedDim) = domainBForces(i);
  }

  return qRes;

}

MATRIX CubatureObjInterface::getVertexElementBasis(MATRIX &basis, int vertexId) {
  MATRIX res(3, reducedDim);
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < reducedDim; c++) {
      res(r, c) = basis(3 * vertexId + r, c);
    }
  }
  return res;
}
