#include "MatrixOps.h"

double vecDotPdk(double *A, double *B, int len) {
	double res = 0.0;
	for(int i = 0; i < len; i++) {
		res += A[i]*B[i];
	}
	return res;
}

void vecMinus(double *A, double *B, int len) {
	for(int i = 0; i < len; i++) {
		A[i] -= B[i];
	}
}

void vecScalarPdk(double *A, double alpha, int len) {
	for(int i = 0; i < len; i++) {
		A[i] *= alpha;
	}
}

double vec2Norm(double *A, int len) {
	double res = 0.0;
	for(int i = 0; i < len; i++) {
		res += (A[i]*A[i]);
	}
	return sqrt(res);
}

void vecScalarDivide(double* A, double alpha, int len) {
	for(int i = 0; i < len; i++) {
		A[i] /= alpha;
	}
}

void QRDecomposition(int rows, int cols, double *U, double *Q, double *R) {
	int M = rows;
	int N = cols;

	for(int j = 0; j < N; j++ ) {
		//Get the vector vJ
		double *vJ = new double[M];
		memcpy(vJ,U + (M*j),sizeof(double)*M);
		for(int i = 0; i < j; i++) {
			double *qI = new double[M];
			memcpy(qI,Q + (M*i),sizeof(double)*M);
			R[COLMAJORIDX(i,j,N)] =  vecDotPdk(vJ,qI,M);
			vecScalarPdk(qI,R[COLMAJORIDX(i,j,N)],M);
			vecMinus(vJ,qI,M);
			delete[] qI;
		}
		R[COLMAJORIDX(j,j,N)] = vec2Norm(vJ,M);
		if(fabs(R[COLMAJORIDX(j,j,N)])>1.0e-6) {
			vecScalarDivide(vJ,R[COLMAJORIDX(j,j,N)],M);
		}
		memcpy(Q + (M*j),vJ,sizeof(double)*M);
		delete[] vJ;
	}
	
}

void computeRotationInvariantPCA(Eigen::MatrixXd &dataPoints, int targetDims, Eigen::MatrixXd &fullToSub)
{
	Eigen::MatrixXd mean;
	Eigen::MatrixXd U0; //Initial PCA Guess

	//Initial PCA 
	dimensionalityReductionByPCASVD(dataPoints,targetDims,U0,mean);

	//Compute the value of c 
	double c;
	int numPoints = dataPoints.rows();
	int numDims = dataPoints.cols();
	Eigen::MatrixXd U0t = U0.transpose();
	std::vector<double> cValues;
	for(int p=0;p<numPoints;p++)
	{
		Eigen::MatrixXd xi = dataPoints.row(p).transpose();
		Eigen::MatrixXd xit = xi.transpose();
		Eigen::MatrixXd valMat = (xit*xi) - (xit*U0*U0t*xi);
		cValues.push_back(valMat(0,0));
	}
	std::sort(cValues.begin(),cValues.end());
	if(numPoints%2)
	{
		c = cValues[numPoints/2];
	}
	else
	{
		c = cValues[numPoints/2]+cValues[(numPoints+1)/2];
		c /=2.0;
	}
	Eigen::MatrixXd ISpecial = Eigen::MatrixXd::Identity(numPoints,targetDims);
	//Now iterate using the values of c
	Eigen::MatrixXd U = U0;
	for(int i=0;i<10;i++)
	{
		Eigen::MatrixXd Ut = U.transpose();
		Eigen::MatrixXd Cr = Eigen::MatrixXd::Zero(numDims,numDims);
		//Construct the matrix Cr iteratively 
		for(int p=0;p<numPoints;p++)
		{
			Eigen::MatrixXd xi = dataPoints.row(p).transpose();
			Eigen::MatrixXd xit = xi.transpose();
			//Compute the cauchy estimator 
			Eigen::MatrixXd intRes = xi - (U*Ut*xi) ;
			double normIntRes = intRes.norm();
			double wi = (1.0+(normIntRes/(c*c)));
			wi = 1.0/wi;
			Eigen::MatrixXd crwi = wi*(xi*xit);
			Cr += crwi;
		}
		Eigen::MatrixXd UtHalf = Cr*U;
		//std::cout << "R:" << UtHalf.rows() << " C:" << UtHalf.cols() << "\n";
		//Orthogonalize
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(UtHalf);
		Eigen::MatrixXd Q = qr.householderQ();
		Q = Q*ISpecial;
		U = Q;
		std::cout << "[INFO] Rotation Invariant PCA Iteration#" << (i+1) << "\n";
	}
	//std::cout <<  U.rows() << "x" << U.cols() << "\n";
	fullToSub = dataPoints * U;
}



void dimensionalityReductionByPCASVD(Eigen::MatrixXd dataPoints, int targetDims,Eigen::MatrixXd &fullToSub,Eigen::MatrixXd &mean)
{

	//Get statistics data 
	int numPoints = dataPoints.rows();
	int numDims = dataPoints.cols();

	//Compute the mean vector which will be used for reconstruction
	mean = Eigen::VectorXd::Zero(numDims);
	for(int d=0;d<numDims;d++)
	{
		Eigen::VectorXd dthCol = dataPoints.col(d);
		double meanCol = dthCol.sum()/numPoints;
		mean(d) = meanCol;
	}

	for(int c=0;c<dataPoints.cols();c++)
	{
		for(int r = 0; r<dataPoints.rows();r++)
		{
			dataPoints(r,c) -= mean(c);
		}
	}
	//Perform the SVD 
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(dataPoints,Eigen::ComputeThinU|Eigen::ComputeThinV);

	std::cout << "SVD Completed\n";
	Eigen::MatrixXd VM = svd.matrixV();
	Eigen::VectorXd sval = svd.singularValues();
	//Select the k-best ones
	fullToSub = Eigen::MatrixXd::Zero(numDims,targetDims);
	for(int k=0;k<targetDims;k++)
	{
		fullToSub.col(k) = VM.col(k);
	}
}

void massPCA(Eigen::MatrixXd dataPoints, Eigen::SparseMatrix<double> massMatrix, int targetDims,Eigen::MatrixXd &fullToSub)
{
	//First compute the cholesky decompositio of massMatrix 
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > choleskySolver(massMatrix);
	//Eigen::LLT<Eigen::MatrixXd> choleskySolver(massMatrix);
	Eigen::SparseMatrix<double> L = choleskySolver.matrixL();
	Eigen::SparseMatrix<double> LT = L.transpose();
	//std::cout << "[INFO] Cholesky Solve finished.\n";

	//Then keeping in mind OUR row major data format compute z
	//Eigen::MatrixXd smallU = dataPoints.transpose();
	Eigen::MatrixXd smallZ = LT*dataPoints;
	//asciiMatrixWriter("dataForMatlab.txt",smallZ);
	Eigen::MatrixXd fuckingHell = smallZ.transpose();


	//Do standard PCA on this 
	Eigen::MatrixXd fstsp;
	Eigen::MatrixXd mean;
	dimensionalityReductionByPCASVD(fuckingHell,targetDims,fstsp,mean);
	//computeRotationInvariantPCA(smallZ,targetDims,fstsp);
	
	//Now solve a bunch of equations 
	//Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(LT);
	fullToSub = Eigen::MatrixXd::Zero(fstsp.rows(),fstsp.cols());
	for(int i=0;i<targetDims;i++)
	{
		Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(LT);
		Eigen::VectorXd VCol = fstsp.col(i);
		Eigen::VectorXd UCol = solver.solve(VCol);
		fullToSub.col(i) = UCol;
	}
	
}

void massOrthogonalization(Eigen::MatrixXd linBasis, Eigen::SparseMatrix<double> massMatrix, Eigen::MatrixXd &massOrtho) {
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