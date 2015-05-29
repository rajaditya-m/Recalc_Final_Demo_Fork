#ifndef SINGLE_DOMAIN_CUBATURE_H
#define SINGLE_DOMAIN_CUBATURE_H

#include "subspace_tet.h"
#include "cubature_base/GreedyCubop.h"


class SingleDomainCubature : public GreedyCubop, public SubspaceTet
{
public:
    SingleDomainCubature(const char* filename,AffineTransformer<double> *affine_transformer = NULL,bool initialize_fem_module = true,bool load_interface = true, int split=0);
    ~SingleDomainCubature();

    void SetFolder(int num_of_samples, string output_folder, double scale, bool read_eigen_value);

    void GenerateCubature(int max_cubature_point, double relative_error);

    int numTotalPoints();
    void evalPointForceDensity(int pointId, VECTOR & q, VECTOR & gOut, int poseIdx);
    virtual void handleCubature(std::vector<int>& selectedPoints, VECTOR& weights, Real relErr);

    std::vector<int>& getCubaturePoints()   { return cubPoints1_;           }
    std::vector<double>& getCubatureWeights()   { return cubWeights1_;           }

private:
     std::string output_folder_;

     int pose_num_;
     std::vector<Vec> pose_rigid_q_;
     std::vector<Vec> pose_q_;
     std::vector<Mat3> pose_rotation_;

     bool stage1CubatureDetected_;
     std::vector<int> cubPoints1_;
     std::vector<double> cubWeights1_;
};

#endif // SINGLE_DOMAIN_CUBATURE_H
