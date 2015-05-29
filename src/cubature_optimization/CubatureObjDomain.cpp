#include "CubatureObjDomain.h"
#include "print_macro.h"
#include <sstream>
#include "global.h"

static const char kFolder[] = DATA_DIRECTORY "/armadillo";

CubatureObjDomain::CubatureObjDomain(const char* fileNameNode, const char* fileNameEle, const char* fileNameBasis) {
  //Read the node file
  std::ifstream nodeFile(fileNameNode);
  std::string line;
  std::getline(nodeFile, line);
  {
    std::istringstream iss(line);
    iss >> numNodes ;
  }
  nodes = new double[numNodes * 3];
  int counter = 0 ;
  while (counter < numNodes) {
    std::getline(nodeFile, line);
    std::istringstream iss(line);
    double x, y, z;
    unsigned int idx;
    iss >> idx >> x >> y >> z;
    nodes[counter * 3 + 0] = x;
    nodes[counter * 3 + 1] = y;
    nodes[counter * 3 + 2] = z;
    counter++;
  }
  nodeFile.close();
  std::cout << "[INFO] Successfully processed " << numNodes << " vertices\n";

  //Read the element file
  std::ifstream eleFile(fileNameEle);
  std::getline(eleFile, line);
  {
    std::istringstream iss(line);
    iss >> numElems;
  }
  elems = new int[numElems * 4];
  counter = 0;
  while (counter < numElems) {
    int p1, p2, p3, p4, mat, idx;
    std::getline(eleFile, line);
    std::istringstream iss(line);
    iss >> idx >> p1 >> p2 >> p3 >> p4 >> mat;
    elems[4 * counter + 0] = (p1);
    elems[4 * counter + 1] = (p2);
    elems[4 * counter + 2] = (p3);
    elems[4 * counter + 3] = (p4);
    counter++;
  }
  eleFile.close();
  std::cout << "[INFO] Successfully processed " <<  numElems << " tets.\n";

  //Material properties - Assume its an ENU material
  density = 1.0e3;
  ym = 1.5e5;
  pr = 0.45;

  //Create the VEGA Mesh FILE
  vegaMesh = new TetMesh(numNodes, nodes, numElems, elems, ym, pr, density);
  materialDS = new StVKIsotropicMaterial(vegaMesh, 1, 500);
  femSolver = new IsotropicHyperelasticFEM(vegaMesh, materialDS);

  //Now read the matrix U from the system
  std::ifstream basisFile(fileNameBasis);
  std::getline(basisFile, line);
  int numRows, numCols;
  {
    std::istringstream iss(line);
    iss >> numRows ;
  }
  std::getline(basisFile, line);
  {
    std::istringstream iss(line);
    iss >> numCols;
  }
  U = MATRIX(numRows, numCols);
  for (int r = 0; r < numRows; r++) {
    std::getline(basisFile, line);
    std::istringstream iss(line);
    for (int c = 0; c < numCols; c++) {
      double val;
      iss >> val;
      U(r, c) = val;
    }
  }
  rsize = U.cols();
  basisFile.close();
  std::cout << "[INFO] Successfully processed " <<  U.rows() << "x" << U.cols() << " basis matrix.\n";
  /*
  //Get the precomputed deformation constants also
  for(int t=0;t<numElems;t++)
  {
  	double xi = nodes[elems[4*t+0]*3+0];
  	double xj = nodes[elems[4*t+1]*3+0];
  	double xk = nodes[elems[4*t+2]*3+0];
  	double xl = nodes[elems[4*t+3]*3+0];
  	double yi = nodes[elems[4*t+0]*3+1];
  	double yj = nodes[elems[4*t+1]*3+1];
  	double yk = nodes[elems[4*t+2]*3+1];
  	double yl = nodes[elems[4*t+3]*3+1];
  	double zi = nodes[elems[4*t+0]*3+2];
  	double zj = nodes[elems[4*t+1]*3+2];
  	double zk = nodes[elems[4*t+2]*3+2];
  	double zl = nodes[elems[4*t+3]*3+2];

  	MATRIX3 dm(VEC3(xi-xl,xj-xl,xk-xl),VEC3(yi-yl,yj-yl,yk-yl),VEC3(zi-zl,zj-zl,zk-zl));
  	MATRIX3 bm = dm.inverse();

  	double w = det(dm);
  	w *= 0.167;

  	BmValues.push_back(bm);
  	wValues.push_back(w);
  }
  */
}

CubatureObjDomain::~CubatureObjDomain(void) {
}

int CubatureObjDomain::numTotalPoints() {
  return numElems;
}

void CubatureObjDomain::evalPointForceDensity(int pointId, VECTOR& q, VECTOR& gOut, int poseIdx ) {
  (void) poseIdx;
  VECTOR fullSpaceDeformations = U * q;
  double* u = fullSpaceDeformations.data();
  //std::cout <<  q(12) << "\n";
  double* intForces = new double[vegaMesh->getNumVertices() * 3];
  int computationMode = IsotropicHyperelasticFEM::COMPUTE_INTERNALFORCES;

  femSolver->GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue(u, NULL, intForces, NULL, computationMode);
  int retVal = femSolver->GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse(pointId, pointId + 1, u, NULL, intForces, NULL, computationMode);
  if (retVal) {
    std::cout << "[INFO] Internal Force computation failed.\n";
  }
  VECTOR intForceVec(12);
  int counter = 0;
  for (int i = 0; i < 4; i++) {
    for (int d = 0; d < 3; d++) {
      intForceVec(counter) = -intForces[elems[4 * pointId + i] * 3 + d]; //The negative sign is VV IMP
      counter++;
    }
  }

  MATRIX UElem = getElementBasisMatrix(pointId);
  MATRIX UElemTrans = UElem.transpose();

  gOut = UElemTrans * intForceVec;

  delete[] intForces;

  return ;
}

void CubatureObjDomain::handleCubature(std::vector<int>& selectedPoints, VECTOR& weights, Real relErr ) {
  hasCubature = true;
  cubaturePoints.resize(selectedPoints.size());
  cubatureWeights = VECTOR(int(selectedPoints.size()));
  cout << "[INFO] n = " << selectedPoints.size() << " \t relerr = " << (relErr * 100) << "%" << endl;
  char fileName[500];
  sprintf(fileName, "%s/Cubatures/domain_%d.txt", kFolder, domainId);
  std::ofstream fileOut(fileName);
  ASSERT(fileOut.is_open(), P(fileName));
  fileOut << selectedPoints.size() << "\n";
  for (int i = 0; i < int(selectedPoints.size()); i++) {
    fileOut << selectedPoints[i] << " " << weights(i) << "\n";
    cubaturePoints[i] = selectedPoints[i];
    cubatureWeights(i) = weights(i);
  }
  fileOut.close();
}

MATRIX CubatureObjDomain::getElementBasisMatrix(int pointId) {
  MATRIX result(12, rsize);
  int counter = 0;
  for (int i = 0; i < 4; i++) {
    for (int d = 0; d < 3; d++) {
      for (int k = 0; k < rsize; k++) {
        result(counter, k) = U(elems[4 * pointId + i] * 3 + d, k);
      }
      counter++;
    }
  }
  return result;
}

VECTOR CubatureObjDomain::getsubspaceInternalForce(VECTOR& q, int poseId) {
  (void) poseId;
  //Project this to full space
  VECTOR fullSpaceDeformations = U * q;
  double* u = fullSpaceDeformations.data();
  double* intForces = new double[vegaMesh->getNumVertices() * 3];
  femSolver->ComputeForces(u, intForces);

  VECTOR fullSpaceForce(vegaMesh->getNumVertices() * 3);
  for (int i = 0; i < vegaMesh->getNumVertices() * 3; i++) {
    fullSpaceForce(i) = -intForces[i];
  }

  //Project this back to the subspace
  MATRIX u_transpose = U.transpose();
  VECTOR subspaceForce = (u_transpose) * fullSpaceForce;
  return subspaceForce;
}

/*
VECTOR CubatureObjDomain::getsubspaceInternalForceByCubature(VECTOR &q)
{
	VECTOR result(q.size());

	for(int i=0;i<rsize;i++)
		result(i) = 0.0;

	if(hasCubature)
	{
		for(int i=0;i<cubaturePoints.size();i++)
		{
			VECTOR fiHat;
			evalPointForceDensity(cubaturePoints[i],q,fiHat);
			fiHat *= cubatureWeights(i);
			result += fiHat;
		}
		return result;
	}
	else
	{
		std::cout << "[ERROR] No cubature is defined. Falling back to full space computation.\n";
		result =  getsubspaceInternalForce(q);
	}
	return result;
}


MATRIX CubatureObjDomain::getsubspaceStiffnessMatrixByCubature(VECTOR &q)
{
	MATRIX result(rsize,rsize);
	result.clear();

	if(hasCubature)
	{
		for(int i=0;i<cubaturePoints.size();i++)
		{
			//Evaluate the element stiffness matrix
			double elemK[144];
			femSolver->ComputeTetK(cubaturePoints[i],elemK,0);
			MATRIX Ke(12,12);
			for(int r=0;r<12;r++)
			{
				for(int c=0;c<12;c++)
				{
					Ke(r,c) = elemK[12*c+r];
				}
			}
			MATRIX Ue = getElementBasisMatrix(cubaturePoints[i]);
			MATRIX pdk = (Ue.transpose())*Ke*Ue;
			pdk *= cubatureWeights(i);
			result += pdk;
		}
		return result;
	}
	else
	{
		std::cout << "[ERROR] No cubature is defined. Stiffness Matrix cannot be computed.\n";
		result.clear();
	}
	return result;
}
*/

void CubatureObjDomain::readCubatureFromFile(const char* fileName) {
  int numCubPoints;
  std::ifstream fileIn(fileName);
  std::string line;
  std::getline(fileIn, line);
  {
    std::istringstream iss(line);
    iss >> numCubPoints ;
  }
  nodes = new double[numNodes * 3];
  cubaturePoints.resize(numCubPoints);
  cubatureWeights = VECTOR(numCubPoints);
  int counter = 0 ;
  while (counter < numCubPoints) {
    std::getline(fileIn, line);
    std::istringstream iss(line);
    double w;
    unsigned int idx;
    iss >> idx >> w;
    cubaturePoints[counter] = idx;
    cubatureWeights(counter) = w;
    counter++;
  }
  fileIn.close();
  hasCubature = true;
}
