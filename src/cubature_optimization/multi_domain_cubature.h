#ifndef MULTIDOMAINCUBATURE_H
#define MULTIDOMAINCUBATURE_H
#include <string>
#include "multi_domain_tet.h"
#include "cubature_base/GreedyCubop.h"

class MultiDomainCubature : public MultiDomainTet, public GreedyCubop  {
public:
  typedef Eigen::Vector4d Vec4;
  MultiDomainCubature(const char* mesh_file_name);
  void GenerateDomainCubature(std::vector<int> max_cubature_point, double relative_error);
  void GenerateInterfaceCubature(std::vector<int> max_cubature_point, double relative_error, bool is_rigid_interface = false);
  void GenerateAllCubature(std::vector<int> domain_max_cubature, double domain_relative_error,
                           std::vector<int> interface_max_cubature, double interface_relative_error);
  void SetFolder(std::string subspace_pose_file,
                 std::string partition_folder,
                 std::string output_folder);

  void SetFolder(int num_of_samples, string partition_folder, string output_folder,
                 double scale, bool read_eigen_value = true,
                 bool over_write_error_file = true);


  virtual int numTotalPoints();
  virtual void evalPointForceDensity(int pointId, VECTOR& q, VECTOR& gOut, int poseIdx);
  virtual void handleCubature(std::vector<int>& selectedPoints, VECTOR& weights, Real relErr);
  virtual ~MultiDomainCubature();

private:
  void ReadInterfaceInfo(std::string partition_folder);
  void ReadSubspacePoseData(std::string subspace_pose_file);
  void LoadDomainPartitionInfo(std::string file_prefix);
  std::vector<std::vector<Vec> > pose_rigid_q_;
  std::vector<std::vector<Vec> > pose_q_;
  std::vector<std::vector<Vec3> > pose_center_of_mass_;
  std::vector<std::vector<Mat3> > pose_rotation_;
  std::vector<int> local_tet2global_tet_;
  std::vector<int> local_vert2global_vert_;
  std::vector<std::string> interface_str_;
  std::vector<std::vector<int> > interface_domains_;
  bool rigid_cubature_;
  int pose_num_;
  std::string file_prefix_;
  std::string subspace_pose_file_;
  std::string partition_folder_;
  std::string output_folder_;

  std::vector<int> interface_cubature_point2local_tet_id_;
  std::vector<Vec4> interface_cubature_weight_;
  int interface_num_;
  int interface_opt_stage_;
  int current_domain_;
};

#endif // MULTIDOMAINCUBATURE_H
