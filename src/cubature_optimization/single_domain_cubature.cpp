#include "single_domain_cubature.h"
#include "string_formatter.h"
#include "global.h"
#include "tet_mesh_simulator_bridge.h"
#include <fstream>

SingleDomainCubature::SingleDomainCubature(const char* filename,AffineTransformer<double> *affine_transformer,bool initialize_fem_module,bool load_interface, int split)
    :SubspaceTet(filename,0,affine_transformer,initialize_fem_module,load_interface,split)
{
}

SingleDomainCubature::~SingleDomainCubature()
{
}

void SingleDomainCubature::SetFolder(int num_of_samples, string output_folder, double scale, bool read_eigen_value) {
    output_folder_ = output_folder;

    pose_num_ = num_of_samples;
    pose_q_.resize(pose_num_);
    pose_rotation_.resize(pose_num_);
    // pose_center_of_mass_.resize(pose_num_);

    Vec eigen_values;
    if (read_eigen_value) {
        std::string eigen_value_file = dj::Format("%z/modal_basis/genBasis.nonlin_weights.bin", GetDataFolder());
        std::ifstream eig(eigen_value_file, std::ios::binary);
        ASSERT(eig.is_open(), P(eigen_value_file));
        int basis_num;
        eig.read((char*) &basis_num, sizeof(int));
        ASSERT(basis_num == basis_num_);
        eigen_values = Vec::Zero(basis_num);
        eig.read((char*)&eigen_values[0], sizeof(double) * basis_num);
        eig.close();
    } else {
        eigen_values = Vec::Zero(basis_num_);
        for (int i = 0; i < eigen_values.size(); ++i) {
            eigen_values[i] = 1.0;
        }
    }

    std::ofstream pose_file(dj::Format("%z/random_pose.txt", GetDataFolder()));
    pose_file << pose_num_ << "\n";
    pose_file << basis_num_ << "\n";
    MERSENNETWISTER rand(time(NULL));
    for (int pose = 0; pose < pose_num_; ++pose) {
        pose_rotation_[pose].setIdentity();
        //   pose_center_of_mass_[pose] = initial_center_of_mass_[p];

        // generate random q
        for (int attempt = 0; true; ++attempt) {
            pose_q_[pose] = Vec::Zero(basis_num_);
            for (int n = 0; n < basis_num_; ++n) {
                pose_q_[pose][n] = scale * sqrt(eigen_values[n] / eigen_values[0]) / 4.0 * rand.randNorm();
                q_[n] = pose_q_[pose][n];
            }

            UpdatePosition();
            bool is_inverted = false;
            int t = 0;
            for (; t < tet_number; ++t) {
                int* verts = tet_ + t * 4;
                MapVec3 pos[4] = {
                    MapVec3(X + verts[0] * 3),
                    MapVec3(X + verts[1] * 3),
                    MapVec3(X + verts[2] * 3),
                    MapVec3(X + verts[3] * 3),
                };
                Vec3 v01 = pos[1] - pos[0];
                Vec3 v02 = pos[2] - pos[0];
                Vec3 v03 = pos[3] - pos[0];
                Vec3 cross = v01.cross(v02);
                if (cross.dot(v03) < 0) {
                    is_inverted = true;
                    break;
                }
            }

            if (attempt % 100 == 0)
                P(pose, attempt, q_.norm(), t);
            if (!is_inverted) {
                break;
            }
        }
        for (int n = 0; n < basis_num_; ++n) {
            pose_file << q_[n] << " ";
        }
        pose_file << "\n";
    }
}


int SingleDomainCubature::numTotalPoints() {
    return tet_number;
}

void SingleDomainCubature::evalPointForceDensity(int pointId, VECTOR &q, VECTOR &gOut, int poseIdx) {
    int global_tet_id = pointId;
    int* verts = tet_ + global_tet_id * 4;
    for(int i=0; i<4;i++){
        int v = verts[i];
        MapVec3((double*) & (inv_fem_->u_[v][0])) = vert_basis_[v] * MapVec(&q(0),basis_num_);
    }
    std::vector<int> one_tet(1, global_tet_id);
    inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(one_tet);
    Vec subspace_force = Vec::Zero(basis_num_);
    for (int i = 0; i < 4; ++i) {
      int v = verts[i];
      // vega internal force is in opposite direction
      subspace_force -= vert_basis_transpose_[v] * MapVec3(inv_fem_->element_force_ + global_tet_id * 12 + i * 3);
    }
    memcpy(&gOut(0), &subspace_force[0], basis_num_* sizeof(double));
}

void SingleDomainCubature::handleCubature(std::vector<int> &selectedPoints, VECTOR &weights, Real relErr) {
    std::string error_file = dj::Format("%z/error.txt",output_folder_);
    std::ofstream out(error_file, std::ios::app);
    int numP = selectedPoints.size();
    ASSERT(out.is_open(), P(error_file));
    out << dj::Format("%d/%d\t %z\n",selectedPoints.size(), numTotalPoints(), relErr);
    out.close();

    std::string cubature_file = dj::Format("%z/cubature_pnts.txt",output_folder_);
    std::ofstream cub_out(cubature_file);
    ASSERT(cub_out.is_open(), P(cubature_file));
    cub_out << selectedPoints.size() << std::endl;
    for (int i = 0; i < int(selectedPoints.size()); ++i) {
        cub_out << selectedPoints[i] << " " << weights(i) << std::endl;
    }
    cub_out.close();
    L(dj::Format("Write cubature optimization result to file %z", cubature_file));

    //Copy them to this data
    cubPoints1_.clear();
    cubWeights1_.clear();
    cubPoints1_.resize(numP);
    cubWeights1_.resize(numP);
    std::copy(selectedPoints.begin(),selectedPoints.end(),cubPoints1_.begin());
    std::copy(weights.data(),weights.data()+numP,cubWeights1_.begin());
}

void SingleDomainCubature::GenerateCubature(int max_cubature_point, double relative_error) {
    std::string file_prefix_ = "partition";
    //current_domain_ = p;
    //LoadDomainPartitionInfo(file_prefix_);
    P("Optimizing cubature....");
    std::vector<VECTOR*> training_set;
    VECTOR training_force(pose_num_ * basis_num_);
    // generate training data
    for (int pose = 0; pose < pose_num_; ++pose) {
        // update local offset
        for (int v =0 ; v< vertex_num_; v++) {
            MapVec3((double*) & (inv_fem_->u_[v][0])) = vert_basis_[v] * pose_q_[pose];
        }
        // compute force
       // inv_fem_->ComputePartialInternalForceAndTangentStiffnessMatrix(local_tet2global_tet_);
        inv_fem_->ComputeInternalForceAndTangentStiffnessMatrix(1.0); //Dummy value of dt
        // assenble subspace force
        Vec subspace_force = Vec::Zero(basis_num_);
        for (int t = 0 ; t < tet_number; t++) {
            int* verts = tet_ + t * 4;
            for (int i = 0; i < 4; ++i) {
                // vega internal force is in negative direction
                subspace_force -= vert_basis_transpose_[verts[i]] * MapVec3(inv_fem_->element_force_ + t * 12 + i * 3);
            }
        }
        memcpy(&training_force(pose * basis_num_), &subspace_force[0], sizeof(double) *basis_num_);
        VECTOR* q = new VECTOR(basis_num_);
        memcpy(&(*q)(0), &pose_q_[pose][0], sizeof(double) * basis_num_);
        training_set.push_back(q);
        // verify eval
        /*
        if (0) {
            int basis_size = basis_num_;
            VECTOR gout(basis_size);
            Vec sub = Vec::Zero(basis_size);
            for (int t = 0; t < int(local_tet2global_tet_.size()); ++t) {
                evalPointForceDensity(t, *q, gout, pose);
                sub += MapVec(&gout(0), basis_size);
            }
            Vec diff = subspace_force - sub;
            ASSERT(diff.norm() < 1e-5, P(diff.norm()));
        }
        */
    }
    int validation_sample_num = int(0.05 * pose_num_);
    if (validation_sample_num == 0) validation_sample_num = 3;
    run(training_set, training_force, relative_error, max_cubature_point, 1000, 120, pose_num_ - validation_sample_num);
    {
        std::string file_name = dj::Format("%z/cubature/error.txt", GetDataFolder());
        std::ofstream out(file_name, std::ios::app);
        out << std::endl;
        out.close();
    }
    for (VECTOR * q : training_set) {
        delete q;
    }
}
